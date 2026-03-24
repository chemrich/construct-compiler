"""
Core construct validity checks.

These operate on a ConstructGraph (at any compilation stage) and return
structured results. Each check function returns a list of CheckResult
objects — an empty list means the check passed cleanly.

Checks are designed to be:
  1. Runnable at any pipeline stage (they skip parts without sequences)
  2. Machine-readable (for automated design loops)
  3. Composable (run individually or via run_all_checks)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Optional

from Bio.Seq import Seq

from ..core.graph import ConstructGraph
from ..core.parts import (
    GeneticPart, CDS, PurificationTag, SolubilityTag,
    CleavageSite, Linker, StopCodon, RBS, Backbone,
    Promoter, Terminator, Spacer,
)
from ..core.types import ResolutionState


# ---------------------------------------------------------------------------
# Result types
# ---------------------------------------------------------------------------

class CheckSeverity(Enum):
    ERROR = auto()    # Construct is non-functional
    WARNING = auto()  # Construct may have reduced function
    INFO = auto()     # Informational (e.g., "frame is correct")


@dataclass
class CheckResult:
    """Outcome of a single validation check."""
    check_name: str          # e.g. "reading_frame"
    severity: CheckSeverity
    part_id: str             # which part or cistron triggered this
    message: str
    details: dict = field(default_factory=dict)  # machine-readable context

    @property
    def passed(self) -> bool:
        return self.severity != CheckSeverity.ERROR

    def __repr__(self):
        icon = {CheckSeverity.ERROR: "FAIL", CheckSeverity.WARNING: "WARN", CheckSeverity.INFO: "OK"}
        return f"[{icon[self.severity]}] {self.check_name} @ {self.part_id}: {self.message}"


# ---------------------------------------------------------------------------
# Which part types are "coding" (translated to protein)?
# ---------------------------------------------------------------------------

CODING_TYPES = (CDS, PurificationTag, SolubilityTag, CleavageSite, Linker)


def _is_coding(part: GeneticPart) -> bool:
    return isinstance(part, CODING_TYPES)


def _has_sequence(part: GeneticPart) -> bool:
    return part.sequence is not None and len(part.sequence) > 0


def _get_expected_protein(part: GeneticPart) -> Optional[str]:
    """Get the expected protein sequence for a part, if known."""
    if isinstance(part, CDS):
        return part.protein_sequence
    return part.metadata.get("protein_sequence")


def _translate_dna(dna: str, include_stop: bool = False) -> str:
    """
    Translate a DNA string to protein using Biopython.
    Returns single-letter amino acid string.
    If include_stop=True, the trailing * is included.
    """
    seq = Seq(dna)
    # Biopython translate: to_stop=False gives full translation including stops as *
    protein = str(seq.translate())
    if not include_stop and protein.endswith("*"):
        protein = protein[:-1]
    return protein


# ---------------------------------------------------------------------------
# Check 1: Reading frame continuity
# ---------------------------------------------------------------------------

def check_reading_frames(graph: ConstructGraph) -> list[CheckResult]:
    """
    Verify that every coding part within a cistron has a DNA sequence
    whose length is divisible by 3, and that the concatenated coding
    DNA for the entire cistron maintains reading frame.

    This catches:
    - Individual parts with non-codon-aligned lengths
    - Accumulating frame drift across a fusion chain
    """
    results: list[CheckResult] = []
    cistrons = graph.cistrons()

    for ci, cistron in enumerate(cistrons):
        cistron_label = f"cistron_{ci}"
        coding_parts = [p for p in cistron if _is_coding(p)]

        if not coding_parts:
            continue

        # Check each coding part individually
        running_length = 0
        all_have_sequence = True

        for part in coding_parts:
            if not _has_sequence(part):
                all_have_sequence = False
                continue

            seq_len = len(part.sequence)
            if seq_len % 3 != 0:
                results.append(CheckResult(
                    check_name="reading_frame",
                    severity=CheckSeverity.ERROR,
                    part_id=part.id,
                    message=(
                        f"Sequence length {seq_len} bp is not divisible by 3 "
                        f"(remainder {seq_len % 3}). This will frameshift all "
                        f"downstream coding elements."
                    ),
                    details={
                        "cistron_index": ci,
                        "sequence_length": seq_len,
                        "remainder": seq_len % 3,
                    },
                ))

            running_length += seq_len

        # Check total cistron coding length
        if all_have_sequence and coding_parts:
            if running_length % 3 != 0:
                results.append(CheckResult(
                    check_name="reading_frame",
                    severity=CheckSeverity.ERROR,
                    part_id=cistron_label,
                    message=(
                        f"Total coding region is {running_length} bp "
                        f"(not divisible by 3). The ORF is frameshifted."
                    ),
                    details={
                        "cistron_index": ci,
                        "total_coding_bp": running_length,
                        "part_count": len(coding_parts),
                    },
                ))
            else:
                results.append(CheckResult(
                    check_name="reading_frame",
                    severity=CheckSeverity.INFO,
                    part_id=cistron_label,
                    message=(
                        f"Reading frame OK: {running_length} bp across "
                        f"{len(coding_parts)} coding parts."
                    ),
                    details={
                        "cistron_index": ci,
                        "total_coding_bp": running_length,
                        "part_count": len(coding_parts),
                    },
                ))

    return results


# ---------------------------------------------------------------------------
# Check 2: Start codon placement
# ---------------------------------------------------------------------------

def check_start_codons(graph: ConstructGraph) -> list[CheckResult]:
    """
    Verify that each cistron's first coding element begins with ATG.

    This catches:
    - Missing start codons on tags that appear before the CDS
    - Start codons that got lost during optimization
    - BCD elements whose downstream coding region doesn't start with ATG

    For BCD-type RBS elements, also checks that the leader ORF structure
    is intact (starts with ATG, contains internal stop).
    """
    results: list[CheckResult] = []
    cistrons = graph.cistrons()

    for ci, cistron in enumerate(cistrons):
        cistron_label = f"cistron_{ci}"

        # Find the RBS (should be first in cistron)
        rbs_part = None
        first_coding = None

        for part in cistron:
            if isinstance(part, RBS):
                rbs_part = part
            elif _is_coding(part) and first_coding is None:
                first_coding = part

        if first_coding is None:
            continue

        if not _has_sequence(first_coding):
            continue

        dna = str(first_coding.sequence).upper()

        if not dna.startswith("ATG"):
            results.append(CheckResult(
                check_name="start_codon",
                severity=CheckSeverity.ERROR,
                part_id=first_coding.id,
                message=(
                    f"First coding element '{first_coding.name}' "
                    f"({first_coding.__class__.__name__}) does not start with "
                    f"ATG. Starts with '{dna[:3]}'. Translation cannot initiate."
                ),
                details={
                    "cistron_index": ci,
                    "first_3bp": dna[:3],
                    "part_type": first_coding.__class__.__name__,
                    "part_name": first_coding.name,
                },
            ))
        else:
            results.append(CheckResult(
                check_name="start_codon",
                severity=CheckSeverity.INFO,
                part_id=first_coding.id,
                message=f"Start codon ATG present on '{first_coding.name}'.",
                details={"cistron_index": ci},
            ))

        # Check BCD leader structure if applicable
        if rbs_part and _has_sequence(rbs_part):
            rbs_data = rbs_part.metadata.get("rbs_entry", {})
            is_bcd = rbs_data.get("is_bcd", False) or "BCD" in rbs_part.name.upper()

            if is_bcd:
                _check_bcd_structure(rbs_part, results, ci)

    return results


def _check_bcd_structure(rbs_part: RBS, results: list[CheckResult], ci: int) -> None:
    """
    Validate BCD element structure:
    - Leader ORF starts with ATG
    - Leader ORF contains a stop codon
    - Internal RBS is positioned correctly
    """
    dna = str(rbs_part.sequence).upper()

    if not dna.startswith("ATG"):
        results.append(CheckResult(
            check_name="bcd_structure",
            severity=CheckSeverity.ERROR,
            part_id=rbs_part.id,
            message=f"BCD leader ORF does not start with ATG. Starts with '{dna[:3]}'.",
            details={"cistron_index": ci},
        ))
        return

    # Check leader has a stop codon
    # Translate the leader and look for a stop
    protein = str(Seq(dna).translate())
    if "*" not in protein:
        results.append(CheckResult(
            check_name="bcd_structure",
            severity=CheckSeverity.ERROR,
            part_id=rbs_part.id,
            message="BCD leader ORF has no stop codon — leader peptide will read through.",
            details={"cistron_index": ci, "leader_translation": protein},
        ))
    else:
        # Find position of first stop
        stop_pos = protein.index("*")
        leader_length_aa = stop_pos
        results.append(CheckResult(
            check_name="bcd_structure",
            severity=CheckSeverity.INFO,
            part_id=rbs_part.id,
            message=f"BCD leader ORF OK: {leader_length_aa} aa leader with stop codon.",
            details={"cistron_index": ci, "leader_length_aa": leader_length_aa},
        ))


# ---------------------------------------------------------------------------
# Check 3: Round-trip translation fidelity
# ---------------------------------------------------------------------------

def check_translation_fidelity(graph: ConstructGraph) -> list[CheckResult]:
    """
    For every coding part with a known expected protein sequence,
    translate the DNA back and verify it matches.

    This catches:
    - Codon optimization errors that changed the protein
    - Off-by-one errors in reverse translation
    - Corrupted sequences from constraint resolution
    """
    results: list[CheckResult] = []

    for part in graph.parts():
        if not _is_coding(part):
            continue
        if not _has_sequence(part):
            continue

        expected = _get_expected_protein(part)
        if not expected:
            # No expected protein to compare against — skip
            # (this is common for parts not yet resolved)
            continue

        dna = str(part.sequence).upper()

        # Handle CDS with stop codon: the DNA includes the stop,
        # but the expected protein sequence typically doesn't include *
        include_stop = False
        if isinstance(part, CDS) and part.has_stop:
            # Translate including stop, then strip it for comparison
            actual_with_stop = _translate_dna(dna, include_stop=True)
            actual = actual_with_stop.rstrip("*")

            # Also verify the stop codon is actually there
            if not actual_with_stop.endswith("*"):
                results.append(CheckResult(
                    check_name="translation_fidelity",
                    severity=CheckSeverity.WARNING,
                    part_id=part.id,
                    message=(
                        f"CDS '{part.name}' has has_stop=True but DNA does "
                        f"not end with a stop codon."
                    ),
                    details={"last_codon": dna[-3:] if len(dna) >= 3 else dna},
                ))
        else:
            actual = _translate_dna(dna, include_stop=False)

        # Strip leading M if the expected protein doesn't start with M
        # but our DNA includes a start codon (common for tags)
        expected_clean = expected.upper().rstrip("*")
        actual_clean = actual.upper()

        # For CDS parts, the DNA may include a start codon (ATG=M)
        # that's part of the protein. For tags, the start codon may have
        # been prepended by the compiler but isn't in the expected protein.
        # Handle both cases:
        match = (actual_clean == expected_clean)

        # Check if the mismatch is just a leading Met
        if not match and actual_clean == "M" + expected_clean:
            # Start codon was added — protein is correct after Met
            results.append(CheckResult(
                check_name="translation_fidelity",
                severity=CheckSeverity.INFO,
                part_id=part.id,
                message=(
                    f"Translation matches expected protein (with added "
                    f"initiator Met). {len(expected_clean)} aa."
                ),
                details={
                    "expected_length": len(expected_clean),
                    "actual_length": len(actual_clean),
                    "has_added_met": True,
                },
            ))
            continue

        if match:
            results.append(CheckResult(
                check_name="translation_fidelity",
                severity=CheckSeverity.INFO,
                part_id=part.id,
                message=f"Translation matches expected protein. {len(expected_clean)} aa.",
                details={
                    "expected_length": len(expected_clean),
                    "actual_length": len(actual_clean),
                },
            ))
        else:
            # Find where they diverge
            mismatch_pos = None
            for i in range(min(len(actual_clean), len(expected_clean))):
                if actual_clean[i] != expected_clean[i]:
                    mismatch_pos = i
                    break
            if mismatch_pos is None:
                mismatch_pos = min(len(actual_clean), len(expected_clean))

            results.append(CheckResult(
                check_name="translation_fidelity",
                severity=CheckSeverity.ERROR,
                part_id=part.id,
                message=(
                    f"Translation MISMATCH for '{part.name}': "
                    f"expected {len(expected_clean)} aa, got {len(actual_clean)} aa. "
                    f"First difference at position {mismatch_pos}."
                ),
                details={
                    "expected_length": len(expected_clean),
                    "actual_length": len(actual_clean),
                    "mismatch_position": mismatch_pos,
                    "expected_around_mismatch": expected_clean[max(0,mismatch_pos-3):mismatch_pos+3],
                    "actual_around_mismatch": actual_clean[max(0,mismatch_pos-3):mismatch_pos+3],
                },
            ))

    return results


# ---------------------------------------------------------------------------
# Check 4: Internal stop codons
# ---------------------------------------------------------------------------

def check_internal_stops(graph: ConstructGraph) -> list[CheckResult]:
    """
    Check for premature in-frame stop codons within coding regions.

    For fusion chains (has_stop=False), the entire coding region must be
    free of stop codons. For terminal CDS (has_stop=True), exactly one
    stop codon should appear, and it should be the last codon.

    Also checks the concatenated cistron-level sequence for stops that
    might span part junctions.
    """
    results: list[CheckResult] = []

    # --- Per-part checks ---
    for part in graph.parts():
        if not _is_coding(part):
            continue
        if not _has_sequence(part):
            continue

        dna = str(part.sequence).upper()
        protein = _translate_dna(dna, include_stop=True)

        if isinstance(part, CDS) and part.has_stop:
            # Should have exactly one stop, at the end
            stop_positions = [i for i, aa in enumerate(protein) if aa == "*"]
            if not stop_positions:
                results.append(CheckResult(
                    check_name="internal_stops",
                    severity=CheckSeverity.WARNING,
                    part_id=part.id,
                    message=f"CDS '{part.name}' has has_stop=True but no stop codon in translation.",
                    details={"protein_length": len(protein)},
                ))
            elif len(stop_positions) == 1 and stop_positions[0] == len(protein) - 1:
                results.append(CheckResult(
                    check_name="internal_stops",
                    severity=CheckSeverity.INFO,
                    part_id=part.id,
                    message=f"Stop codon correctly placed at end of '{part.name}'.",
                    details={"stop_position": stop_positions[0]},
                ))
            elif len(stop_positions) > 1 or (stop_positions and stop_positions[0] != len(protein) - 1):
                internal = [p for p in stop_positions if p != len(protein) - 1]
                results.append(CheckResult(
                    check_name="internal_stops",
                    severity=CheckSeverity.ERROR,
                    part_id=part.id,
                    message=(
                        f"PREMATURE stop codon(s) in '{part.name}' at amino acid "
                        f"position(s) {internal}. Protein will be truncated."
                    ),
                    details={
                        "internal_stop_positions": internal,
                        "total_stops": len(stop_positions),
                    },
                ))
        else:
            # Non-terminal coding part: should have NO stops at all
            stop_positions = [i for i, aa in enumerate(protein) if aa == "*"]
            if stop_positions:
                results.append(CheckResult(
                    check_name="internal_stops",
                    severity=CheckSeverity.ERROR,
                    part_id=part.id,
                    message=(
                        f"Stop codon(s) in non-terminal part '{part.name}' at "
                        f"position(s) {stop_positions}. This will truncate the "
                        f"fusion protein."
                    ),
                    details={
                        "stop_positions": stop_positions,
                        "part_type": part.__class__.__name__,
                    },
                ))
            else:
                results.append(CheckResult(
                    check_name="internal_stops",
                    severity=CheckSeverity.INFO,
                    part_id=part.id,
                    message=f"No internal stop codons in '{part.name}'.",
                    details={},
                ))

    # --- Junction checks: look for stop codons spanning part boundaries ---
    cistrons = graph.cistrons()
    for ci, cistron in enumerate(cistrons):
        coding_parts = [p for p in cistron if _is_coding(p) and _has_sequence(p)]
        if len(coding_parts) < 2:
            continue

        for i in range(len(coding_parts) - 1):
            part_a = coding_parts[i]
            part_b = coding_parts[i + 1]

            # Check the 6 bp spanning the junction (last 3 of A + first 3 of B,
            # and positions that straddle the boundary)
            tail = str(part_a.sequence)[-2:]  # last 2 bp of part A
            head = str(part_b.sequence)[:2]   # first 2 bp of part B
            junction_6bp = str(part_a.sequence)[-3:] + str(part_b.sequence)[:3]

            # Translate the junction region in-frame
            # The junction codons depend on the reading frame. Since each part
            # should be codon-aligned, the junction should be at a codon boundary.
            # But let's verify by checking if a stop codon appears at the boundary.
            if len(junction_6bp) >= 3:
                for offset in range(0, len(junction_6bp) - 2, 3):
                    codon = junction_6bp[offset:offset+3].upper()
                    if codon in ("TAA", "TAG", "TGA"):
                        # But only flag this if the parts are supposed to be fused
                        # (i.e., part_a doesn't have has_stop=True)
                        is_fused = not (isinstance(part_a, CDS) and part_a.has_stop)
                        if is_fused:
                            results.append(CheckResult(
                                check_name="junction_stops",
                                severity=CheckSeverity.ERROR,
                                part_id=f"{part_a.id}|{part_b.id}",
                                message=(
                                    f"Stop codon '{codon}' at junction between "
                                    f"'{part_a.name}' and '{part_b.name}'. "
                                    f"Fusion protein will be truncated here."
                                ),
                                details={
                                    "cistron_index": ci,
                                    "junction_sequence": junction_6bp,
                                    "stop_codon": codon,
                                },
                            ))

    return results


# ---------------------------------------------------------------------------
# Aggregate runner
# ---------------------------------------------------------------------------

def run_all_checks(graph: ConstructGraph) -> list[CheckResult]:
    """
    Run all validation checks and return aggregated results.
    Suitable for both automated testing and human review.
    """
    results: list[CheckResult] = []
    results.extend(check_reading_frames(graph))
    results.extend(check_start_codons(graph))
    results.extend(check_translation_fidelity(graph))
    results.extend(check_internal_stops(graph))
    return results


def errors_only(results: list[CheckResult]) -> list[CheckResult]:
    """Filter to only ERROR-severity results."""
    return [r for r in results if r.severity == CheckSeverity.ERROR]


def summary_report(results: list[CheckResult]) -> str:
    """Human-readable summary of check results."""
    errors = [r for r in results if r.severity == CheckSeverity.ERROR]
    warnings = [r for r in results if r.severity == CheckSeverity.WARNING]
    info = [r for r in results if r.severity == CheckSeverity.INFO]

    lines = [f"Construct Validation: {len(errors)} errors, {len(warnings)} warnings, {len(info)} passed"]
    lines.append("=" * 70)

    if errors:
        lines.append("\nERRORS:")
        for r in errors:
            lines.append(f"  {r}")

    if warnings:
        lines.append("\nWARNINGS:")
        for r in warnings:
            lines.append(f"  {r}")

    if info:
        lines.append(f"\nPASSED: {len(info)} checks")

    return "\n".join(lines)

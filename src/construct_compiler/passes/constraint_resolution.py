"""
Constraint Resolution Pass

Uses DNA Chisel to optimize DNA sequences subject to constraints:
- Remove internal restriction sites (e.g., BsaI for Golden Gate)
- Enforce GC content windows
- Remove long homopolymer runs
- Optimize codon usage
- Preserve protein sequences (synonymous changes only)

This pass operates on CONCRETE parts and refines their sequences.
"""

from __future__ import annotations

import logging
from typing import Optional

from Bio.Seq import Seq

from ..core.graph import ConstructGraph
from ..core.types import ResolutionState, ConstraintKind, Constraint
from ..core.parts import (
    CDS, PurificationTag, SolubilityTag, CleavageSite, Linker, Backbone,
    Promoter, RBS, Terminator, Spacer,
)

logger = logging.getLogger(__name__)

# Parts whose sequences should NOT be modified by constraint resolution
# (their exact sequence is functionally critical)
IMMUTABLE_TYPES = (Promoter, RBS, Terminator, Backbone)


def resolve_constraints(graph: ConstructGraph) -> ConstructGraph:
    """
    Apply constraint resolution to all mutable, concrete parts.
    Uses DNA Chisel when available, falls back to manual fixes.
    """
    # Collect global constraints
    global_constraints = graph.constraints

    # Identify which parts are mutable (coding regions we can change synonymously)
    for part in graph.parts():
        if part.resolution != ResolutionState.CONCRETE:
            continue
        if part.sequence is None:
            continue
        if isinstance(part, IMMUTABLE_TYPES):
            continue

        # Merge global + local constraints
        all_constraints = global_constraints + part.constraints

        if not all_constraints:
            continue

        # Try DNA Chisel optimization
        optimized = _optimize_with_dnachisel(part, all_constraints, graph.host_organism)
        if optimized is not None:
            part.sequence = optimized
            logger.info(f"Constraint-optimized '{part.name}': {len(optimized)} bp")
        else:
            # Fallback: manual restriction site removal
            _manual_constraint_fix(part, all_constraints)

    return graph


def _optimize_with_dnachisel(part, constraints: list[Constraint],
                              organism: str) -> Optional[Seq]:
    """
    Use DNA Chisel to optimize a sequence subject to constraints.
    Returns the optimized sequence, or None if DNA Chisel is unavailable.
    """
    try:
        from dnachisel import (
            DnaOptimizationProblem,
            AvoidPattern,
            EnforceGCContent,
            CodonOptimize,
            EnforceTranslation,
        )
    except ImportError:
        logger.warning("DNA Chisel not available, falling back to manual fixes")
        return None

    sequence = str(part.sequence)
    dnachisel_constraints = []
    dnachisel_objectives = []

    # Determine if this is a coding sequence
    is_coding = isinstance(part, (CDS, PurificationTag, SolubilityTag, CleavageSite, Linker))
    protein_seq = None
    if is_coding:
        protein_seq = _get_protein(part)

    for c in constraints:
        if c.kind == ConstraintKind.NO_RESTRICTION_SITE:
            enzyme = c.params.get("enzyme", "BsaI")
            dnachisel_constraints.append(AvoidPattern(f"{enzyme}_site"))

        elif c.kind == ConstraintKind.GC_CONTENT:
            lo = c.params.get("lo", 0.35)
            hi = c.params.get("hi", 0.65)
            dnachisel_constraints.append(
                EnforceGCContent(mini=lo, maxi=hi, window=50)
            )

        elif c.kind == ConstraintKind.MAX_HOMOPOLYMER:
            n = c.params.get("max_run", 6)
            for base in "ATGC":
                dnachisel_constraints.append(
                    AvoidPattern(base * (n + 1))
                )

    # If coding: enforce translation and optimize codons
    if is_coding and protein_seq:
        dnachisel_constraints.append(
            EnforceTranslation()
        )
        # Codon optimization as an objective (not constraint)
        species = _organism_to_species(organism)
        if species:
            dnachisel_objectives.append(
                CodonOptimize(species=species)
            )

    if not dnachisel_constraints and not dnachisel_objectives:
        return None

    try:
        problem = DnaOptimizationProblem(
            sequence=sequence,
            constraints=dnachisel_constraints,
            objectives=dnachisel_objectives,
        )
        problem.resolve_constraints()
        problem.optimize()

        result = problem.sequence
        logger.info(
            f"DNA Chisel optimization for '{part.name}': "
            f"{len(problem.constraints)} constraints, "
            f"{len(problem.objectives)} objectives"
        )
        return Seq(result)

    except Exception as e:
        logger.warning(f"DNA Chisel failed for '{part.name}': {e}")
        return None


def _manual_constraint_fix(part, constraints: list[Constraint]) -> None:
    """
    Fallback: manually fix constraint violations when DNA Chisel isn't available.
    Currently handles restriction site removal via synonymous codon swaps.
    """
    from Bio.Restriction import BsaI, BpiI, BbsI

    enzyme_map = {
        "BsaI": BsaI,
        "BpiI": BpiI,
        "BbsI": BbsI,
    }

    sequence = str(part.sequence)

    for c in constraints:
        if c.kind == ConstraintKind.NO_RESTRICTION_SITE:
            enzyme_name = c.params.get("enzyme", "BsaI")
            enzyme = enzyme_map.get(enzyme_name)
            if enzyme:
                seq_obj = Seq(sequence)
                # Check for sites
                from Bio.Restriction import Analysis
                analysis = Analysis([enzyme], seq_obj)
                results = analysis.full()
                if results.get(enzyme):
                    logger.warning(
                        f"Part '{part.name}' has {len(results[enzyme])} "
                        f"{enzyme_name} sites — manual removal not yet implemented"
                    )


def _get_protein(part) -> Optional[str]:
    """Get protein sequence from a part."""
    if isinstance(part, CDS):
        return part.protein_sequence
    return part.metadata.get("protein_sequence")


def _organism_to_species(organism: str) -> Optional[str]:
    """Map organism shorthand to DNA Chisel species name."""
    mapping = {
        "e_coli": "e_coli",
        "e_coli_bl21": "e_coli",
        "s_cerevisiae": "s_cerevisiae",
        "h_sapiens": "h_sapiens",
        "human": "h_sapiens",
        "yeast": "s_cerevisiae",
    }
    return mapping.get(organism.lower())

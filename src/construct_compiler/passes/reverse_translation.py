"""
Reverse Translation Pass

Converts protein sequences to DNA sequences using E. coli preferred codons.
This is a first-pass "rough draft" — the constraint resolution pass will
refine it further (removing restriction sites, adjusting GC, etc.).

Parts that are already CONCRETE are left untouched.
Parts that are RESOLVED with protein sequences get translated to DNA.
"""

from __future__ import annotations

import logging

from Bio.Seq import Seq

from ..core.graph import ConstructGraph
from ..core.types import ResolutionState
from ..core.parts import (
    CDS, PurificationTag, SolubilityTag, CleavageSite, Linker, StopCodon,
)
from ..data.parts_db import ECOLI_PREFERRED_CODONS

logger = logging.getLogger(__name__)


def reverse_translate(graph: ConstructGraph) -> ConstructGraph:
    """
    Reverse-translate all RESOLVED protein-coding parts to DNA.
    Uses E. coli preferred codons as the initial codon selection.
    """
    # Identify which parts are the first coding element in each cistron.
    # These need a start codon (ATG) regardless of part type — a His-tag
    # that leads a fusion chain needs ATG just as much as a CDS does.
    coding_types = (CDS, PurificationTag, SolubilityTag, CleavageSite, Linker)
    cistron_leaders = _find_cistron_leaders(graph, coding_types)

    for part in graph.parts():
        if part.resolution == ResolutionState.CONCRETE:
            continue
        if part.resolution != ResolutionState.RESOLVED:
            continue

        protein_seq = _get_protein_sequence(part)
        if not protein_seq:
            continue

        dna = _protein_to_dna(protein_seq, graph.host_organism)

        # Add start codon only if this part initiates translation in its
        # cistron — i.e., it's the first coding element after the RBS.
        # Mid-chain parts (tags, CDS in a fusion) must NOT get ATG forced
        # onto them; their first codon encodes part of the fusion protein.
        needs_start = part.id in cistron_leaders
        if needs_start and not dna.startswith("ATG"):
            # PREPEND ATG rather than replacing the first codon.
            # For CDS parts whose protein starts with Met, the DNA already
            # starts with ATG so this branch is skipped. For tags/other parts
            # whose protein doesn't include the initiator Met (e.g., 6xHis =
            # "HHHHHH"), we prepend ATG so the expressed protein is
            # Met + tag_protein. The initiator Met may be removed in vivo
            # by methionine aminopeptidase (MAP) depending on the second
            # residue, but the DNA must include it for translation initiation.
            dna = "ATG" + dna

        # Add stop codon to terminal CDS
        if isinstance(part, CDS) and part.has_stop and not _ends_with_stop(dna):
            dna += "TAA"  # TAA preferred in E. coli

        part.sequence = Seq(dna)
        part.resolution = ResolutionState.CONCRETE
        logger.info(
            f"{part.__class__.__name__} '{part.name}' reverse-translated: "
            f"{len(protein_seq)} aa -> {len(dna)} bp"
        )

    return graph


def _find_cistron_leaders(graph: ConstructGraph, coding_types: tuple) -> set[str]:
    """
    Walk each cistron and return the part IDs of the first coding element
    in each one. These parts need a start codon prepended.

    A "cistron leader" is the first part after the RBS that will actually
    be translated — typically a PurificationTag, SolubilityTag, or CDS.

    Special case: when a catalog backbone provides the RBS (so no RBS node
    exists in the insert graph for the first cistron), cistrons() misses it.
    We detect this and add the first coding element in the cassette as a
    leader, since the backbone's RBS expects the next ORF to start with ATG.
    """
    from ..core.parts import Backbone, RBS as RBSPart

    leaders: set[str] = set()

    # Standard path: RBS-delimited cistrons
    for cistron in graph.cistrons():
        for part in cistron:
            if isinstance(part, coding_types):
                leaders.add(part.id)
                break  # only the first one per cistron

    # Catalog backbone path: backbone provides RBS, so the first cistron
    # in the insert may not have an RBS node. The first coding element
    # before any insert-provided RBS is still a cistron leader.
    backbone_provides_rbs = False
    for part in graph.parts():
        if isinstance(part, Backbone) and part.is_catalog_vector:
            if "rbs" in part.provides:
                backbone_provides_rbs = True
                break

    if backbone_provides_rbs:
        # Find the first coding element that appears before any RBS in the
        # insert graph. This is the leader of the backbone-provided cistron.
        for part in graph.parts():
            if isinstance(part, Backbone):
                continue
            if isinstance(part, RBSPart):
                break  # hit an insert RBS — everything after is already handled
            if isinstance(part, coding_types):
                leaders.add(part.id)
                break  # only the first one

    return leaders


def _get_protein_sequence(part) -> str | None:
    """Extract protein sequence from a part."""
    if isinstance(part, CDS):
        return part.protein_sequence
    return part.metadata.get("protein_sequence")


def _protein_to_dna(protein: str, organism: str = "e_coli") -> str:
    """
    Reverse-translate a protein sequence to DNA using preferred codons.

    Currently only supports E. coli. Future: add codon tables for
    yeast, human, insect cells, etc.
    """
    codons = ECOLI_PREFERRED_CODONS
    dna_codons = []
    for aa in protein.upper():
        if aa == "*":
            continue  # skip stop symbols in protein sequences
        codon = codons.get(aa)
        if codon is None:
            raise ValueError(f"Unknown amino acid: '{aa}'")
        dna_codons.append(codon)
    return "".join(dna_codons)


def _ends_with_stop(dna: str) -> bool:
    """Check if a DNA sequence ends with a stop codon."""
    if len(dna) < 3:
        return False
    last_codon = dna[-3:].upper()
    return last_codon in ("TAA", "TAG", "TGA")

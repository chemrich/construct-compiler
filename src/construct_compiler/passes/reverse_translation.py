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
    for part in graph.parts():
        if part.resolution == ResolutionState.CONCRETE:
            continue
        if part.resolution != ResolutionState.RESOLVED:
            continue

        protein_seq = _get_protein_sequence(part)
        if not protein_seq:
            continue

        dna = _protein_to_dna(protein_seq, graph.host_organism)

        # Add start codon if this is a CDS or the first coding element after an RBS
        if isinstance(part, CDS):
            if not dna.startswith("ATG"):
                dna = "ATG" + dna[3:] if len(dna) >= 3 else "ATG" + dna
            # Add stop codon if needed
            if part.has_stop and not _ends_with_stop(dna):
                dna += "TAA"  # TAA preferred in E. coli

        part.sequence = Seq(dna)
        part.resolution = ResolutionState.CONCRETE
        logger.info(
            f"{part.__class__.__name__} '{part.name}' reverse-translated: "
            f"{len(protein_seq)} aa -> {len(dna)} bp"
        )

    return graph


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

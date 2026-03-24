"""
Shared fixtures for construct validation tests.

Provides compiled graphs at each pipeline stage:
  1. post_parse     — fresh from YAML, all parts ABSTRACT
  2. post_resolve   — parts resolved (sequences fetched/looked up)
  3. post_revtrans  — protein→DNA reverse translation done
  4. post_constrain — DNA Chisel constraint resolution done

Also provides hand-built "known-good" and "known-bad" constructs
for unit-testing individual validators.
"""

from __future__ import annotations

import copy
from pathlib import Path

import pytest
from Bio.Seq import Seq

from construct_compiler.frontend.parser import parse_spec
from construct_compiler.passes.part_resolution import resolve_parts
from construct_compiler.passes.reverse_translation import reverse_translate
from construct_compiler.passes.constraint_resolution import resolve_constraints
from construct_compiler.core.graph import ConstructGraph
from construct_compiler.core.parts import (
    Backbone, Promoter, RBS, CDS, PurificationTag, SolubilityTag,
    CleavageSite, Linker, StopCodon, Terminator, Spacer,
    Origin, RegulationType, RBSDesignMethod, TagPosition,
)
from construct_compiler.core.types import ResolutionState


# ---------------------------------------------------------------------------
# Path to example specs
# ---------------------------------------------------------------------------

EXAMPLES_DIR = Path(__file__).parent.parent / "examples"
EXAMPLE_SPEC = EXAMPLES_DIR / "his_tev_mbp_egfp.yaml"


# ---------------------------------------------------------------------------
# Pipeline stage fixtures (from the example polycistronic spec)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def parsed_graph() -> ConstructGraph:
    """Graph immediately after YAML parsing (ABSTRACT parts)."""
    return parse_spec(EXAMPLE_SPEC)


@pytest.fixture(scope="session")
def resolved_graph() -> ConstructGraph:
    """Graph after part resolution (RESOLVED — sequences looked up)."""
    graph = parse_spec(EXAMPLE_SPEC)
    return resolve_parts(graph)


@pytest.fixture(scope="session")
def revtrans_graph() -> ConstructGraph:
    """Graph after reverse translation (CONCRETE — DNA sequences assigned)."""
    graph = parse_spec(EXAMPLE_SPEC)
    graph = resolve_parts(graph)
    return reverse_translate(graph)


@pytest.fixture(scope="session")
def constrained_graph() -> ConstructGraph:
    """Graph after constraint resolution (CONCRETE — codon-optimized)."""
    graph = parse_spec(EXAMPLE_SPEC)
    graph = resolve_parts(graph)
    graph = reverse_translate(graph)
    return resolve_constraints(graph)


# ---------------------------------------------------------------------------
# Minimal hand-built fixtures for unit testing validators
# ---------------------------------------------------------------------------

@pytest.fixture
def good_single_cistron() -> ConstructGraph:
    """
    A minimal correct single-cistron construct:
    Promoter → RBS → 6xHis(+ATG) → CDS(with stop) → Terminator

    All sequences are codon-aligned and correct.
    """
    graph = ConstructGraph(name="good_single", host_organism="e_coli")

    graph.add_part(Backbone(id="bb", name="pBR322", origin=Origin.PBR322))

    graph.add_part(Promoter(
        id="prom", name="T7", regulation=RegulationType.INDUCIBLE,
        sequence=Seq("TAATACGACTCACTATAGGG"),
        resolution=ResolutionState.CONCRETE,
    ))

    graph.add_part(RBS(
        id="rbs", name="B0034", part_name="BBa_B0034",
        sequence=Seq("AAAGAGGAGAAA"),
        resolution=ResolutionState.CONCRETE,
    ))

    # 6xHis: ATG + HHHHHH = ATG + CAT*6 = 21 bp (divisible by 3? 21/3=7 ✓)
    # Wait: we need initiator Met + His*6 = M + HHHHHH
    # ATG CAT CAT CAT CAT CAT CAT = 21 bp ✓
    graph.add_part(PurificationTag(
        id="his", name="6xHis", tag_type="6xHis",
        sequence=Seq("ATGCATCATCATCATCATCAT"),  # M + 6xHis = 21 bp
        resolution=ResolutionState.CONCRETE,
        metadata={"protein_sequence": "HHHHHH"},
    ))

    # Small test CDS: protein "MVKF" + stop
    # M=ATG, V=GTG, K=AAA, F=TTC + TAA = 15 bp
    graph.add_part(CDS(
        id="cds", name="test_gene",
        protein_sequence="MVKF",
        has_stop=True,
        sequence=Seq("ATGGTGAAATTCTAA"),  # MVKF*
        resolution=ResolutionState.CONCRETE,
    ))

    graph.add_part(Terminator(
        id="term", name="rrnB_T1",
        sequence=Seq("CAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCG"),
        resolution=ResolutionState.CONCRETE,
    ))

    graph.connect_linear()
    return graph


@pytest.fixture
def bad_frameshift() -> ConstructGraph:
    """
    A construct with a frameshifted tag (17 bp instead of 18).
    This should be caught by the reading frame check.
    """
    graph = ConstructGraph(name="bad_frameshift", host_organism="e_coli")

    graph.add_part(Backbone(id="bb", name="pBR322"))

    graph.add_part(Promoter(
        id="prom", name="T7",
        sequence=Seq("TAATACGACTCACTATAGGG"),
        resolution=ResolutionState.CONCRETE,
    ))

    graph.add_part(RBS(
        id="rbs", name="B0034",
        sequence=Seq("AAAGAGGAGAAA"),
        resolution=ResolutionState.CONCRETE,
    ))

    # BAD: 17 bp tag — not divisible by 3!
    graph.add_part(PurificationTag(
        id="his", name="bad_tag", tag_type="bad",
        sequence=Seq("ATGCATCATCATCATCAC"),  # 18 bp — wait, let me make it 17
        resolution=ResolutionState.CONCRETE,
        metadata={"protein_sequence": "HHHH"},
    ))
    # Actually set it to 17 bp
    graph.get_part("his").sequence = Seq("ATGCATCATCATCATCA")  # 17 bp

    graph.add_part(CDS(
        id="cds", name="test_gene",
        protein_sequence="MVKF",
        has_stop=True,
        sequence=Seq("ATGGTGAAATTCTAA"),
        resolution=ResolutionState.CONCRETE,
    ))

    graph.add_part(Terminator(
        id="term", name="rrnB_T1",
        sequence=Seq("CAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCG"),
        resolution=ResolutionState.CONCRETE,
    ))

    graph.connect_linear()
    return graph


@pytest.fixture
def bad_no_start_codon() -> ConstructGraph:
    """
    A construct where the first coding element (tag) lacks a start codon.
    This mirrors the likely bug in the current pipeline.
    """
    graph = ConstructGraph(name="bad_no_start", host_organism="e_coli")

    graph.add_part(Backbone(id="bb", name="pBR322"))

    graph.add_part(Promoter(
        id="prom", name="T7",
        sequence=Seq("TAATACGACTCACTATAGGG"),
        resolution=ResolutionState.CONCRETE,
    ))

    graph.add_part(RBS(
        id="rbs", name="B0034",
        sequence=Seq("AAAGAGGAGAAA"),
        resolution=ResolutionState.CONCRETE,
    ))

    # Tag without start codon — just HHHHHH → CAT*6 = 18 bp
    graph.add_part(PurificationTag(
        id="his", name="6xHis", tag_type="6xHis",
        sequence=Seq("CATCATCATCATCATCAT"),  # No ATG!
        resolution=ResolutionState.CONCRETE,
        metadata={"protein_sequence": "HHHHHH"},
    ))

    graph.add_part(CDS(
        id="cds", name="test_gene",
        protein_sequence="MVKF",
        has_stop=True,
        sequence=Seq("ATGGTGAAATTCTAA"),
        resolution=ResolutionState.CONCRETE,
    ))

    graph.add_part(Terminator(
        id="term", name="rrnB_T1",
        sequence=Seq("CAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCG"),
        resolution=ResolutionState.CONCRETE,
    ))

    graph.connect_linear()
    return graph


@pytest.fixture
def bad_internal_stop() -> ConstructGraph:
    """
    A construct with a premature stop codon inside a non-terminal CDS.
    """
    graph = ConstructGraph(name="bad_internal_stop", host_organism="e_coli")

    graph.add_part(Backbone(id="bb", name="pBR322"))
    graph.add_part(Promoter(
        id="prom", name="T7",
        sequence=Seq("TAATACGACTCACTATAGGG"),
        resolution=ResolutionState.CONCRETE,
    ))
    graph.add_part(RBS(
        id="rbs", name="B0034",
        sequence=Seq("AAAGAGGAGAAA"),
        resolution=ResolutionState.CONCRETE,
    ))

    # Non-terminal CDS (has_stop=False) but DNA contains a stop codon
    # ATG GTG TAA TTC = M V * F — premature stop!
    graph.add_part(CDS(
        id="cds1", name="bad_gene",
        protein_sequence="MVF",
        has_stop=False,
        sequence=Seq("ATGGTGTAATTC"),  # Contains TAA at position 6-8
        resolution=ResolutionState.CONCRETE,
    ))

    graph.add_part(CDS(
        id="cds2", name="downstream_gene",
        protein_sequence="MVKF",
        has_stop=True,
        sequence=Seq("ATGGTGAAATTCTAA"),
        resolution=ResolutionState.CONCRETE,
    ))

    graph.add_part(Terminator(
        id="term", name="rrnB_T1",
        sequence=Seq("CAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCG"),
        resolution=ResolutionState.CONCRETE,
    ))

    graph.connect_linear()
    return graph

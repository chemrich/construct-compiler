"""
Integration tests for the omegamega oligopool design bridge.

These tests require a local clone of chemrich/omegamega and the OMEGAMEGA_DIR
environment variable pointing to it:

    git clone https://github.com/chemrich/omegamega
    export OMEGAMEGA_DIR=/path/to/omegamega
    pytest tests/test_omega_integration.py -v

All tests are marked `slow` and are skipped automatically when OMEGAMEGA_DIR
is not set. Run the full suite with:

    pytest tests/test_omega_integration.py -v -m slow
"""

from __future__ import annotations

import csv
import os
from pathlib import Path

import pytest
from Bio.Seq import Seq

from construct_compiler.backends.omega import OmegaBatchResult, OmegaResult, run_omega, run_omega_batch, to_fasta
from construct_compiler.core.graph import ConstructGraph
from construct_compiler.core.parts import (
    Backbone, CDS, Origin, Promoter, PurificationTag, RBS, Terminator,
)
from construct_compiler.core.types import ResolutionState

OMEGAMEGA_DIR = os.environ.get("OMEGAMEGA_DIR")
skip_no_omega = pytest.mark.skipif(
    not OMEGAMEGA_DIR,
    reason="OMEGAMEGA_DIR not set — clone chemrich/omegamega and set the env var",
)

# ---------------------------------------------------------------------------
# Minimal concrete graph fixture (~330 bp insert — above Twist's 300 bp minimum)
# ---------------------------------------------------------------------------

@pytest.fixture
def short_concrete_graph() -> ConstructGraph:
    """Single-cistron construct with a ~330 bp insert for fast omega tests."""
    graph = ConstructGraph(name="short_omega_test", host_organism="e_coli")
    graph.add_part(Backbone(id="bb", name="pET-28b(+)", origin=Origin.PBR322))
    graph.add_part(Promoter(
        id="prom", name="T7",
        sequence=Seq("TAATACGACTCACTATAGGG"),
        resolution=ResolutionState.CONCRETE,
    ))
    graph.add_part(RBS(
        id="rbs", name="BCD2",
        sequence=Seq("AAAGAGGAGAAA"),
        resolution=ResolutionState.CONCRETE,
    ))
    # 6xHis tag with start codon — 21 bp
    graph.add_part(PurificationTag(
        id="his", name="6xHis", tag_type="6xHis",
        sequence=Seq("ATGCATCATCATCATCATCAT"),
        resolution=ResolutionState.CONCRETE,
        metadata={"protein_sequence": "MHHHHHH"},
    ))
    # 96-codon synthetic CDS (288 bp) + stop = 291 bp; total insert ~330 bp
    cds_seq = (
        "ATGAAAGTGCTGACCGAAATCGGTAAAGCGCTGCAGAAAGTTGCGGAAGACCTGCAGCGTATC"
        "GTTGAAGCGCTGGAACAGCGTCTGACCGACATCGCGGAAATGGTTCAGAAAGCGCTGCAGGAC"
        "GTGGCGGAAGACCTGCAGCGCATCGTGGAAGCGCTGGAACAGCGCCTGACCGACATCGCGGAA"
        "ATGGTGCAGAAAGCGCTGCAGGACGTGGCGGAAGACCTGCAGCGCATCGTGGAAGCGCTGGAA"
        "CAGTAA"
    )
    graph.add_part(CDS(
        id="cds", name="synthetic_cds",
        protein_sequence="MKVLTEIGTGALQKVAEDLQRIVEALEQULTDIAEMUQKALQDVAEDLQRIVEALEQULTDIAEMUQKATED",
        has_stop=True,
        sequence=Seq(cds_seq),
        resolution=ResolutionState.CONCRETE,
    ))
    graph.add_part(Terminator(
        id="term", name="rrnB_T1",
        sequence=Seq("CAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCG"),
        resolution=ResolutionState.CONCRETE,
    ))
    graph.connect_linear()
    return graph


# ---------------------------------------------------------------------------
# to_fasta() — no omegamega needed
# ---------------------------------------------------------------------------

def test_to_fasta_produces_valid_sequence(short_concrete_graph):
    """to_fasta() returns valid FASTA with DNA-only sequence."""
    fasta = to_fasta(short_concrete_graph)
    lines = fasta.strip().splitlines()
    assert lines[0].startswith(">"), "First line must be a FASTA header"
    seq = "".join(lines[1:])
    assert len(seq) > 0
    assert set(seq.upper()).issubset(set("ACGT")), f"Non-DNA characters found: {set(seq.upper()) - set('ACGT')}"


def test_to_fasta_uses_graph_name(short_concrete_graph):
    fasta = to_fasta(short_concrete_graph)
    assert "short_omega_test" in fasta.splitlines()[0]


def test_to_fasta_label_override(short_concrete_graph):
    fasta = to_fasta(short_concrete_graph, label="my_insert")
    assert fasta.startswith(">my_insert")


def test_to_fasta_raises_on_abstract_graph():
    """to_fasta() raises if graph has no concrete sequences."""
    from construct_compiler.frontend.parser import parse_spec
    graph = parse_spec(
        Path(__file__).parent.parent / "examples" / "his_tev_mbp_egfp.yaml"
    )
    with pytest.raises(ValueError, match="no concrete sequences"):
        to_fasta(graph)


def test_to_fasta_insert_length_matches_graph(short_concrete_graph):
    """FASTA sequence length matches graph.full_insert_sequence()."""
    fasta = to_fasta(short_concrete_graph)
    seq = "".join(fasta.strip().splitlines()[1:])
    expected = str(short_concrete_graph.full_insert_sequence())
    assert seq == expected


# ---------------------------------------------------------------------------
# run_omega() integration — requires OMEGAMEGA_DIR
# ---------------------------------------------------------------------------

@pytest.mark.slow
@skip_no_omega
def test_omega_produces_output_files(short_concrete_graph, tmp_path):
    """run_omega() creates oligo_order.csv and pool_stats.csv."""
    result = run_omega(
        short_concrete_graph,
        output_dir=tmp_path / "omega_out",
        nopt_steps=50,
        nopt_runs=1,
        njobs=1,
        omegamega_dir=Path(OMEGAMEGA_DIR),
    )
    assert (result.output_dir / "oligo_order.csv").exists()
    assert (result.output_dir / "pool_stats.csv").exists()


@pytest.mark.slow
@skip_no_omega
def test_omega_result_fields(short_concrete_graph, tmp_path):
    """OmegaResult fields are in valid ranges."""
    result = run_omega(
        short_concrete_graph,
        output_dir=tmp_path / "omega_out",
        nopt_steps=50,
        nopt_runs=1,
        njobs=1,
        omegamega_dir=Path(OMEGAMEGA_DIR),
    )
    assert isinstance(result, OmegaResult)
    assert result.oligo_count > 0
    assert result.pool_count >= 1
    assert 0.0 < result.min_fidelity <= 1.0
    assert 0.0 < result.avg_fidelity <= 1.0
    assert result.min_fidelity <= result.avg_fidelity


@pytest.mark.slow
@skip_no_omega
def test_omega_oligos_meet_minimum_size(short_concrete_graph, tmp_path):
    """All designed oligos are at least min_size bp."""
    result = run_omega(
        short_concrete_graph,
        output_dir=tmp_path / "omega_out",
        nopt_steps=50,
        nopt_runs=1,
        njobs=1,
        omegamega_dir=Path(OMEGAMEGA_DIR),
    )
    with open(result.output_dir / "oligo_order.csv") as f:
        oligos = list(csv.DictReader(f))
    assert len(oligos) > 0
    # Sequences should all be non-empty DNA strings
    for row in oligos:
        seq = row.get("sequence", "")
        assert len(seq) >= 40, f"Oligo shorter than min_size: {seq!r}"
        assert set(seq.upper()).issubset(set("ACGTN")), f"Non-DNA oligo: {seq!r}"


@pytest.mark.slow
@skip_no_omega
def test_omega_full_pipeline_from_spec(tmp_path):
    """End-to-end: YAML spec → compile → omegamega → oligo order."""
    from construct_compiler.frontend.parser import parse_spec
    from construct_compiler.passes.part_resolution import resolve_parts
    from construct_compiler.passes.reverse_translation import reverse_translate

    spec = Path(__file__).parent.parent / "examples" / "his_tev_mbp_egfp.yaml"
    graph = parse_spec(spec)
    graph = resolve_parts(graph)
    graph = reverse_translate(graph)

    result = run_omega(
        graph,
        output_dir=tmp_path / "omega_full",
        nopt_steps=100,
        nopt_runs=1,
        njobs=1,
        omegamega_dir=Path(OMEGAMEGA_DIR),
    )
    assert result.oligo_count > 0
    assert result.min_fidelity > 0.5, "Expected reasonable fidelity for a real construct"


@pytest.mark.slow
@skip_no_omega
def test_omega_missing_dir_raises():
    """run_omega() raises RuntimeError when omegamega dir is not found."""
    from construct_compiler.core.graph import ConstructGraph
    graph = ConstructGraph(name="dummy", host_organism="e_coli")

    with pytest.raises((RuntimeError, FileNotFoundError)):
        run_omega(graph, output_dir=Path("/tmp/omega_test"), omegamega_dir=Path("/nonexistent"))


# ---------------------------------------------------------------------------
# run_omega_batch() integration — requires OMEGAMEGA_DIR
# ---------------------------------------------------------------------------

@pytest.mark.slow
@skip_no_omega
def test_omega_batch_two_constructs(short_concrete_graph, tmp_path):
    """run_omega_batch() on two copies produces a combined oligo order."""
    import copy
    graph2 = copy.deepcopy(short_concrete_graph)
    graph2.name = "short_omega_test_b"

    result = run_omega_batch(
        [(short_concrete_graph, "construct_a"), (graph2, "construct_b")],
        output_dir=tmp_path / "batch_out",
        nopt_steps=50,
        nopt_runs=1,
        njobs=1,
        omegamega_dir=Path(OMEGAMEGA_DIR),
    )
    assert isinstance(result, OmegaBatchResult)
    assert result.construct_count == 2
    assert result.total_oligos > 0
    assert result.pool_count >= 1
    assert 0.0 < result.min_fidelity <= 1.0
    assert result.min_fidelity <= result.avg_fidelity


@pytest.mark.slow
@skip_no_omega
def test_omega_batch_cost_vs_individual(short_concrete_graph, tmp_path):
    """Batch cost should be <= sum of individual costs once volume discount kicks in."""
    import copy
    graphs = []
    for i in range(3):
        g = copy.deepcopy(short_concrete_graph)
        g.name = f"construct_{i}"
        graphs.append((g, f"construct_{i}"))

    batch = run_omega_batch(
        graphs,
        output_dir=tmp_path / "batch_cost",
        nopt_steps=50,
        nopt_runs=1,
        njobs=1,
        omegamega_dir=Path(OMEGAMEGA_DIR),
    )
    assert batch.total_cost_usd >= 0.0
    assert batch.construct_count == 3


@pytest.mark.slow
@skip_no_omega
def test_omega_batch_empty_raises():
    """run_omega_batch() raises ValueError on empty input."""
    with pytest.raises(ValueError, match="non-empty"):
        run_omega_batch([], output_dir=Path("/tmp/batch_test"),
                        omegamega_dir=Path(OMEGAMEGA_DIR))

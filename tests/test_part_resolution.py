"""Tests for CDS resolution paths in part_resolution.py."""

from unittest.mock import patch, MagicMock

import pytest

from construct_compiler.core.parts import CDS
from construct_compiler.core.types import ResolutionState
from construct_compiler.passes.part_resolution import (
    _resolve_cds,
    _lookup_builtin,
    _is_uniprot_accession,
    _BUILTIN_FP_SEQUENCES,
    _FP_ALIASES,
    _uniprot_gene_cache,
)


# ---------------------------------------------------------------------------
# _is_uniprot_accession
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("accession", [
    "P03023", "P0ACT4", "P0A9E0", "P03034", "P42212",
    "Q9Y4K3", "O15116", "A0A000AAA0",
])
def test_is_accession_true(accession):
    assert _is_uniprot_accession(accession)


@pytest.mark.parametrize("name", [
    "lacI", "tetR", "araC", "GFP", "EGFP", "mEGFP", "dCas9", "rtTA3",
    "cI", "luxR", "", "P1234",
])
def test_is_accession_false(name):
    assert not _is_uniprot_accession(name)


# ---------------------------------------------------------------------------
# _lookup_builtin — aliases
# ---------------------------------------------------------------------------

def test_lookup_gfp_alias():
    result = _lookup_builtin("GFP", _BUILTIN_FP_SEQUENCES, _FP_ALIASES)
    assert result == _BUILTIN_FP_SEQUENCES["mEGFP"]


def test_lookup_egfp_alias():
    result = _lookup_builtin("EGFP", _BUILTIN_FP_SEQUENCES, _FP_ALIASES)
    assert result == _BUILTIN_FP_SEQUENCES["mEGFP"]


def test_lookup_rfp_alias():
    result = _lookup_builtin("RFP", _BUILTIN_FP_SEQUENCES, _FP_ALIASES)
    assert result == _BUILTIN_FP_SEQUENCES["mCherry"]


def test_lookup_alias_case_insensitive():
    result = _lookup_builtin("gfp", _BUILTIN_FP_SEQUENCES, _FP_ALIASES)
    assert result == _BUILTIN_FP_SEQUENCES["mEGFP"]


def test_lookup_direct_case_insensitive():
    result = _lookup_builtin("mscarlet", _BUILTIN_FP_SEQUENCES)
    assert result == _BUILTIN_FP_SEQUENCES["mScarlet"]


def test_lookup_missing_returns_none():
    assert _lookup_builtin("unknownProtein", _BUILTIN_FP_SEQUENCES, _FP_ALIASES) is None


# ---------------------------------------------------------------------------
# _resolve_cds — pre-set protein sequence
# ---------------------------------------------------------------------------

def _make_cds(name, source_db="", source_id=None, protein_sequence=None):
    return CDS(
        id=f"cds_{name}",
        name=name,
        source_db=source_db,
        source_id=source_id or name,
        protein_sequence=protein_sequence,
    )


def test_resolve_cds_already_has_sequence():
    part = _make_cds("myGene", protein_sequence="MAAAK")
    _resolve_cds(part)
    assert part.resolution == ResolutionState.RESOLVED
    assert part.protein_sequence == "MAAAK"


# ---------------------------------------------------------------------------
# _resolve_cds — FP aliases (no network)
# ---------------------------------------------------------------------------

def test_resolve_cds_gfp_alias_no_source():
    part = _make_cds("GFP", source_db="")
    _resolve_cds(part)
    assert part.resolution == ResolutionState.RESOLVED
    assert part.protein_sequence == _BUILTIN_FP_SEQUENCES["mEGFP"]


def test_resolve_cds_egfp_alias_fpbase_source():
    """When fpbase fetch fails, falls back to alias lookup."""
    part = _make_cds("EGFP", source_db="fpbase", source_id="EGFP")
    with patch("construct_compiler.passes.part_resolution._fetch_fpbase", return_value=None):
        _resolve_cds(part)
    assert part.resolution == ResolutionState.RESOLVED
    assert part.protein_sequence == _BUILTIN_FP_SEQUENCES["mEGFP"]


def test_resolve_cds_rfp_alias_no_source():
    part = _make_cds("RFP", source_db="")
    _resolve_cds(part)
    assert part.resolution == ResolutionState.RESOLVED
    assert part.protein_sequence == _BUILTIN_FP_SEQUENCES["mCherry"]


def test_resolve_cds_mcherry_case_insensitive():
    part = _make_cds("mcherry", source_db="")
    _resolve_cds(part)
    assert part.resolution == ResolutionState.RESOLVED
    assert part.protein_sequence == _BUILTIN_FP_SEQUENCES["mCherry"]


# ---------------------------------------------------------------------------
# _resolve_cds — UniProt accession path
# ---------------------------------------------------------------------------

def test_resolve_cds_uniprot_accession_fetch():
    part = _make_cds("GFP", source_db="uniprot", source_id="P42212")
    fake_seq = "MSKGEELFT"
    with patch("construct_compiler.passes.part_resolution._fetch_uniprot", return_value=fake_seq):
        _resolve_cds(part)
    assert part.resolution == ResolutionState.RESOLVED
    assert part.protein_sequence == fake_seq


def test_resolve_cds_uniprot_accession_fallback_to_builtin():
    part = _make_cds("GFP", source_db="uniprot", source_id="P42212")
    with patch("construct_compiler.passes.part_resolution._fetch_uniprot", return_value=None):
        _resolve_cds(part)
    assert part.resolution == ResolutionState.RESOLVED  # P42212 in _BUILTIN_UNIPROT_SEQUENCES


# ---------------------------------------------------------------------------
# _resolve_cds — UniProt gene-name search
# ---------------------------------------------------------------------------

def test_resolve_cds_gene_name_via_uniprot_search():
    """source=uniprot + gene name (not accession) triggers gene-name search."""
    part = _make_cds("lacI", source_db="uniprot", source_id="lacI")
    fake_seq = "MKPVTLYDVAEYAG"
    with patch("construct_compiler.passes.part_resolution._fetch_uniprot_by_gene",
               return_value=fake_seq) as mock_search:
        _resolve_cds(part)
    mock_search.assert_called_once_with("lacI")
    assert part.resolution == ResolutionState.RESOLVED
    assert part.protein_sequence == fake_seq


def test_resolve_cds_gene_name_no_source_last_resort():
    """Gene name with no source falls through to gene-name search as last resort."""
    part = _make_cds("tetR", source_db="", source_id="tetR")
    fake_seq = "MSRLDKSKVINA"
    with patch("construct_compiler.passes.part_resolution._fetch_uniprot_by_gene",
               return_value=fake_seq):
        _resolve_cds(part)
    assert part.resolution == ResolutionState.RESOLVED
    assert part.protein_sequence == fake_seq


def test_resolve_cds_gene_name_search_fails():
    """If all resolution paths fail, part stays unresolved."""
    part = _make_cds("unknownXYZ", source_db="", source_id="unknownXYZ")
    with patch("construct_compiler.passes.part_resolution._fetch_uniprot_by_gene",
               return_value=None):
        _resolve_cds(part)
    assert part.resolution != ResolutionState.RESOLVED


def test_resolve_cds_gene_name_search_cached(monkeypatch):
    """Gene-name search result is cached so the HTTP call happens only once."""
    monkeypatch.setitem(_uniprot_gene_cache, "arac", "MCACALSD")

    part = _make_cds("araC", source_db="", source_id="araC")
    with patch("construct_compiler.passes.part_resolution._fetch_uniprot_by_gene",
               wraps=lambda name: _uniprot_gene_cache.get(name.lower())) as mock_search:
        _resolve_cds(part)
    assert part.protein_sequence == "MCACALSD"

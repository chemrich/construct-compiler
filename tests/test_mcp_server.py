"""
Tests for the MCP server tool handlers.

These test the handler functions directly (no JSON-RPC transport),
ensuring that each MCP tool returns well-formed responses.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from construct_compiler.mcp_server import (
    _handle_compile,
    _handle_check,
    _handle_variants,
    MAX_VARIANTS,
)


EXAMPLE_SPEC_PATH = Path(__file__).parent.parent / "examples" / "his_tev_mbp_egfp.yaml"
EXAMPLE_SPEC = yaml.safe_load(EXAMPLE_SPEC_PATH.read_text())


# ---------------------------------------------------------------------------
# compile_spec
# ---------------------------------------------------------------------------

class TestCompileSpec:
    def test_basic_compile(self):
        result = _handle_compile({"spec": EXAMPLE_SPEC, "skip_constraints": True})
        assert len(result) == 1
        data = json.loads(result[0].text)
        assert data["success"] is True
        assert data["construct_name"] != ""
        assert data["insert_length_bp"] > 0
        assert len(data["parts"]) > 0
        assert len(data["strategies"]) > 0

    def test_compile_includes_costs(self):
        result = _handle_compile({"spec": EXAMPLE_SPEC, "skip_constraints": True})
        data = json.loads(result[0].text)
        for strat in data["strategies"]:
            assert "total_cost" in strat
            assert strat["total_cost"] > 0
            assert "wall_clock_days" in strat
            assert len(strat["wall_clock_days"]) == 2

    def test_compile_with_genbank(self):
        result = _handle_compile({
            "spec": EXAMPLE_SPEC,
            "skip_constraints": True,
            "include_genbank": True,
        })
        data = json.loads(result[0].text)
        assert "genbank" in data
        assert "LOCUS" in data["genbank"]

    def test_compile_without_genbank(self):
        result = _handle_compile({
            "spec": EXAMPLE_SPEC,
            "skip_constraints": True,
            "include_genbank": False,
        })
        data = json.loads(result[0].text)
        assert "genbank" not in data

    def test_compile_bad_spec_does_not_crash(self):
        """Even a minimal/invalid spec should return a result, not crash."""
        result = _handle_compile({
            "spec": {"construct": {"cassette": [{"promoter": "NONEXISTENT"}]}},
            "skip_constraints": True,
        })
        data = json.loads(result[0].text)
        # Should return some kind of structured response (success or error)
        assert isinstance(data, dict)
        assert "error" in data or "success" in data


# ---------------------------------------------------------------------------
# check_spec
# ---------------------------------------------------------------------------

class TestCheckSpec:
    def test_example_passes(self):
        result = _handle_check({
            "spec": EXAMPLE_SPEC,
            "skip_constraints": True,
            "auto_save": False,
        })
        data = json.loads(result[0].text)
        assert data["passed"] is True
        assert data["score"] > 0.0
        assert data["error_count"] == 0

    def test_check_includes_summary(self):
        result = _handle_check({
            "spec": EXAMPLE_SPEC,
            "skip_constraints": True,
            "auto_save": False,
        })
        data = json.loads(result[0].text)
        assert "summary" in data
        assert "PASS" in data["summary"] or "FAIL" in data["summary"]

    def test_check_with_intermediate(self):
        result = _handle_check({
            "spec": EXAMPLE_SPEC,
            "skip_constraints": True,
            "validate_intermediate": True,
            "auto_save": False,
        })
        data = json.loads(result[0].text)
        assert "stages" in data
        assert len(data["stages"]) >= 2

    def test_check_auto_save(self, tmp_path, monkeypatch):
        """When auto_save=True and spec passes, it should be saved."""
        import construct_compiler.mcp_server as mcp_mod
        monkeypatch.setattr(mcp_mod, "AGENT_GENERATED_DIR", tmp_path)

        result = _handle_check({
            "spec": EXAMPLE_SPEC,
            "skip_constraints": True,
            "auto_save": True,
        })
        data = json.loads(result[0].text)
        assert data["passed"] is True
        assert "saved_to" in data
        saved = Path(data["saved_to"])
        assert saved.exists()
        assert saved.suffix == ".yaml"


# ---------------------------------------------------------------------------
# evaluate_variants
# ---------------------------------------------------------------------------

class TestEvaluateVariants:
    def test_expression_variants(self):
        result = _handle_variants({
            "spec": EXAMPLE_SPEC,
            "axes": [
                {
                    "name": "expression",
                    "field_path": "cassette.1.cistron.expression",
                    "values": ["high", "medium", "low"],
                }
            ],
            "skip_constraints": True,
        })
        data = json.loads(result[0].text)
        assert data["total_variants"] == 3
        assert data["passed_count"] + data["failed_count"] == 3
        assert len(data["variants"]) == 3

    def test_variant_cap_enforced(self):
        """Requesting more than MAX_VARIANTS should return an error."""
        # Create axes that produce > MAX_VARIANTS combinations
        result = _handle_variants({
            "spec": EXAMPLE_SPEC,
            "axes": [
                {
                    "name": "a",
                    "field_path": "cassette.1.cistron.expression",
                    "values": list(range(10)),
                },
                {
                    "name": "b",
                    "field_path": "cassette.2.spacer",
                    "values": list(range(10)),
                },
            ],
            "skip_constraints": True,
        })
        data = json.loads(result[0].text)
        assert "error" in data
        assert data["total_combinations"] == 100

    def test_variants_include_parameters(self):
        result = _handle_variants({
            "spec": EXAMPLE_SPEC,
            "axes": [
                {
                    "name": "expression",
                    "field_path": "cassette.1.cistron.expression",
                    "values": ["high", "low"],
                }
            ],
            "skip_constraints": True,
        })
        data = json.loads(result[0].text)
        for v in data["variants"]:
            assert "parameters" in v
            assert "expression" in v["parameters"]

    def test_variants_sorted_by_score(self):
        result = _handle_variants({
            "spec": EXAMPLE_SPEC,
            "axes": [
                {
                    "name": "expression",
                    "field_path": "cassette.1.cistron.expression",
                    "values": ["high", "medium", "low"],
                }
            ],
            "skip_constraints": True,
        })
        data = json.loads(result[0].text)
        scores = [v["score"] for v in data["variants"]]
        assert scores == sorted(scores, reverse=True)

    def test_batch_summary_in_response(self):
        result = _handle_variants({
            "spec": EXAMPLE_SPEC,
            "axes": [
                {
                    "name": "expression",
                    "field_path": "cassette.1.cistron.expression",
                    "values": ["high", "low"],
                }
            ],
            "skip_constraints": True,
        })
        data = json.loads(result[0].text)
        assert "batch_summary" in data
        assert "Batch Evaluation" in data["batch_summary"]

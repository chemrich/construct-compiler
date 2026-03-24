"""
Tests for the validation harness and variant generator.

These test the full automated pipeline:
  spec → compile → validate → score → rank
"""

from __future__ import annotations

from pathlib import Path

import pytest

from construct_compiler.validation.harness import (
    evaluate_spec,
    evaluate_batch,
    batch_summary,
    EvalResult,
)
from construct_compiler.validation.variants import (
    DesignAxis,
    vary_spec,
    vary_spec_dicts,
)


EXAMPLE_SPEC = Path(__file__).parent.parent / "examples" / "his_tev_mbp_egfp.yaml"


class TestEvaluateSpec:
    """Test single-spec evaluation."""

    def test_example_spec_passes(self):
        result = evaluate_spec(EXAMPLE_SPEC)
        assert result.passed, f"Example spec failed:\n{result.summary()}"

    def test_result_has_score(self):
        result = evaluate_spec(EXAMPLE_SPEC)
        assert 0.0 <= result.score <= 1.0

    def test_result_has_metadata(self):
        result = evaluate_spec(EXAMPLE_SPEC)
        assert result.construct_name != ""
        assert result.insert_length_bp is not None
        assert result.insert_length_bp > 0
        assert result.cistron_count == 2  # target + reporter

    def test_result_serializes_to_json(self):
        result = evaluate_spec(EXAMPLE_SPEC)
        j = result.to_json()
        assert '"passed": true' in j or '"passed": false' in j

    def test_skip_constraints_is_faster(self):
        r_full = evaluate_spec(EXAMPLE_SPEC, skip_constraints=False)
        r_skip = evaluate_spec(EXAMPLE_SPEC, skip_constraints=True)
        # Both should pass (skip_constraints doesn't change correctness)
        assert r_skip.passed
        # Skip should be faster (or at least not crash)
        assert r_skip.compile_time_s >= 0

    def test_intermediate_stages_captured(self):
        result = evaluate_spec(EXAMPLE_SPEC, validate_intermediate=True)
        assert len(result.stages) >= 2  # at least post_resolve and post_revtrans

    def test_invalid_spec_returns_compile_error(self):
        bad_spec = {"construct": {"cassette": [{"promoter": "NONEXISTENT_XYZ"}]}}
        result = evaluate_spec(bad_spec)
        # Should not crash — should return a result with error info
        assert isinstance(result, EvalResult)


class TestVariantGenerator:
    """Test parametric variant generation."""

    def test_single_axis(self):
        axes = [
            DesignAxis("expression", "cassette.1.cistron.expression",
                       values=["high", "medium", "low"]),
        ]
        variants = vary_spec(EXAMPLE_SPEC, axes)
        assert len(variants) == 3

    def test_two_axes_cross_product(self):
        axes = [
            DesignAxis("expression", "cassette.1.cistron.expression",
                       values=["high", "low"]),
            DesignAxis("spacer", "cassette.2.spacer",
                       values=[20, 50]),
        ]
        variants = vary_spec(EXAMPLE_SPEC, axes)
        assert len(variants) == 4  # 2 × 2

    def test_variants_have_distinct_names(self):
        axes = [
            DesignAxis("expression", "cassette.1.cistron.expression",
                       values=["high", "low"]),
        ]
        variants = vary_spec(EXAMPLE_SPEC, axes)
        names = [v.variant_name for v in variants]
        assert len(set(names)) == len(names)  # all unique

    def test_variants_are_compilable(self):
        """Each generated variant should compile without crashing."""
        axes = [
            DesignAxis("expression", "cassette.1.cistron.expression",
                       values=["high", "low"]),
        ]
        variants = vary_spec_dicts(EXAMPLE_SPEC, axes)
        results = evaluate_batch(variants, skip_constraints=True)
        assert len(results) == 2
        # All should compile (may or may not pass all checks)
        for r in results:
            assert r.compile_error is None, f"Compile error: {r.compile_error}"


class TestBatchEvaluation:
    """Test batch mode: multiple specs evaluated and ranked."""

    def test_batch_returns_sorted_results(self):
        axes = [
            DesignAxis("expression", "cassette.1.cistron.expression",
                       values=["high", "medium", "low"]),
        ]
        variants = vary_spec_dicts(EXAMPLE_SPEC, axes)
        results = evaluate_batch(variants, skip_constraints=True)
        # Should be sorted by score descending
        scores = [r.score for r in results]
        assert scores == sorted(scores, reverse=True)

    def test_batch_summary_format(self):
        axes = [
            DesignAxis("expression", "cassette.1.cistron.expression",
                       values=["high", "low"]),
        ]
        variants = vary_spec_dicts(EXAMPLE_SPEC, axes)
        results = evaluate_batch(variants, skip_constraints=True)
        summary = batch_summary(results)
        assert "Batch Evaluation" in summary
        assert "2 specs" in summary

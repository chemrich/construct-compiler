"""
Variant evaluation smoke tests.

Ensures the variant generator + batch evaluator work together,
and that known-good parameter combinations still pass. These tests
serve as regression guards for the design space exploration feature.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from construct_compiler.validation.harness import evaluate_batch, batch_summary
from construct_compiler.validation.variants import DesignAxis, vary_spec, vary_spec_dicts


EXAMPLE_SPEC = Path(__file__).parent.parent / "examples" / "his_tev_mbp_egfp.yaml"


class TestExpressionVariants:
    """All expression levels should produce valid constructs."""

    def test_all_expression_levels_pass(self):
        axes = [
            DesignAxis("expression", "cassette.1.cistron.expression",
                       ["high", "medium", "low"]),
        ]
        specs = vary_spec_dicts(EXAMPLE_SPEC, axes)
        results = evaluate_batch(specs, skip_constraints=True)
        for r in results:
            assert r.passed, (
                f"Expression variant failed:\n{r.summary()}"
            )

    def test_expression_scores_are_positive(self):
        axes = [
            DesignAxis("expression", "cassette.1.cistron.expression",
                       ["high", "medium", "low"]),
        ]
        specs = vary_spec_dicts(EXAMPLE_SPEC, axes)
        results = evaluate_batch(specs, skip_constraints=True)
        for r in results:
            assert r.score > 0.0, f"Score should be positive, got {r.score}"


class TestSpacerVariants:
    """All reasonable spacer lengths should produce valid constructs."""

    def test_standard_spacer_lengths_pass(self):
        axes = [
            DesignAxis("spacer", "cassette.2.spacer", [20, 30, 50, 100]),
        ]
        specs = vary_spec_dicts(EXAMPLE_SPEC, axes)
        results = evaluate_batch(specs, skip_constraints=True)
        for r in results:
            assert r.passed, (
                f"Spacer variant failed:\n{r.summary()}"
            )

    def test_minimum_spacer_length(self):
        """Spacer of 20 bp (the documented minimum) should work."""
        axes = [
            DesignAxis("spacer", "cassette.2.spacer", [20]),
        ]
        specs = vary_spec_dicts(EXAMPLE_SPEC, axes)
        results = evaluate_batch(specs, skip_constraints=True)
        assert results[0].passed


class TestCrossProductVariants:
    """Multi-axis variant generation and evaluation."""

    def test_expression_x_spacer(self):
        """2x3 = 6 variants should all compile without error."""
        axes = [
            DesignAxis("expression", "cassette.1.cistron.expression",
                       ["high", "low"]),
            DesignAxis("spacer", "cassette.2.spacer", [20, 50, 100]),
        ]
        variants = vary_spec(EXAMPLE_SPEC, axes)
        assert len(variants) == 6

        specs = [v.spec for v in variants]
        results = evaluate_batch(specs, skip_constraints=True)

        # All should at least compile
        for r in results:
            assert r.compile_error is None, (
                f"Compile error in variant: {r.compile_error}"
            )

    def test_results_are_sorted_by_score(self):
        axes = [
            DesignAxis("expression", "cassette.1.cistron.expression",
                       ["high", "medium", "low"]),
        ]
        specs = vary_spec_dicts(EXAMPLE_SPEC, axes)
        results = evaluate_batch(specs, skip_constraints=True)
        scores = [r.score for r in results]
        assert scores == sorted(scores, reverse=True)

    def test_batch_summary_is_readable(self):
        axes = [
            DesignAxis("expression", "cassette.1.cistron.expression",
                       ["high", "low"]),
        ]
        specs = vary_spec_dicts(EXAMPLE_SPEC, axes)
        results = evaluate_batch(specs, skip_constraints=True)
        summary = batch_summary(results)
        assert "Batch Evaluation" in summary
        assert "2 specs" in summary

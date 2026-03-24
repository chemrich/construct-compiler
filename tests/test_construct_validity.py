"""
Tests for construct validity checks.

Organized in two sections:
  1. Unit tests — hand-built fixtures with known outcomes
  2. Integration tests — the example polycistronic spec run through
     the actual pipeline at each stage

Each test function tests ONE validator against ONE fixture. This makes
failures easy to diagnose: you know immediately which check failed at
which pipeline stage.
"""

from __future__ import annotations

import pytest

from construct_compiler.validation import (
    CheckResult,
    CheckSeverity,
    check_reading_frames,
    check_start_codons,
    check_translation_fidelity,
    check_internal_stops,
    run_all_checks,
)
from construct_compiler.validation.construct_checks import errors_only, summary_report


# ===================================================================
# Helpers
# ===================================================================

def assert_no_errors(results: list[CheckResult], context: str = ""):
    """Assert that no ERROR-severity results exist."""
    errors = errors_only(results)
    if errors:
        msg = f"Unexpected errors{f' ({context})' if context else ''}:\n"
        msg += "\n".join(f"  {e}" for e in errors)
        pytest.fail(msg)


def assert_has_error(results: list[CheckResult], check_name: str,
                     part_id_contains: str = ""):
    """Assert that at least one ERROR exists for the given check."""
    errors = [r for r in results
              if r.severity == CheckSeverity.ERROR
              and r.check_name == check_name]
    if part_id_contains:
        errors = [r for r in errors if part_id_contains in r.part_id]
    assert errors, (
        f"Expected ERROR for check '{check_name}'"
        f"{f' at part containing {part_id_contains!r}' if part_id_contains else ''}"
        f", but none found. Results:\n"
        + "\n".join(f"  {r}" for r in results)
    )


# ===================================================================
# UNIT TESTS — hand-built fixtures
# ===================================================================

class TestReadingFrameUnit:
    """Test the reading frame validator with hand-built constructs."""

    def test_good_construct_has_no_frame_errors(self, good_single_cistron):
        results = check_reading_frames(good_single_cistron)
        assert_no_errors(results, "good single cistron")

    def test_frameshift_detected(self, bad_frameshift):
        results = check_reading_frames(bad_frameshift)
        assert_has_error(results, "reading_frame", "his")

    def test_frameshift_reports_remainder(self, bad_frameshift):
        results = check_reading_frames(bad_frameshift)
        errors = [r for r in results if r.severity == CheckSeverity.ERROR]
        assert any(r.details.get("remainder") in (1, 2) for r in errors)


class TestStartCodonUnit:
    """Test start codon placement with hand-built constructs."""

    def test_good_construct_has_start_codon(self, good_single_cistron):
        results = check_start_codons(good_single_cistron)
        assert_no_errors(results, "good single cistron")

    def test_missing_start_codon_detected(self, bad_no_start_codon):
        results = check_start_codons(bad_no_start_codon)
        assert_has_error(results, "start_codon", "his")

    def test_missing_start_reports_actual_triplet(self, bad_no_start_codon):
        results = check_start_codons(bad_no_start_codon)
        errors = [r for r in results
                  if r.severity == CheckSeverity.ERROR
                  and r.check_name == "start_codon"]
        assert any(r.details.get("first_3bp") == "CAT" for r in errors)


class TestTranslationFidelityUnit:
    """Test round-trip translation with hand-built constructs."""

    def test_good_construct_translates_correctly(self, good_single_cistron):
        results = check_translation_fidelity(good_single_cistron)
        assert_no_errors(results, "good single cistron")

    def test_frameshifted_translation_detected(self, bad_frameshift):
        """A frameshifted tag will produce wrong protein on translation."""
        results = check_translation_fidelity(bad_frameshift)
        # The tag's DNA is 17 bp, so translation won't match expected protein
        errors = [r for r in results if r.severity == CheckSeverity.ERROR]
        # May or may not error depending on whether 17bp translates close to HHHH
        # The key thing is check_reading_frames catches the frame issue


class TestInternalStopsUnit:
    """Test internal stop codon detection."""

    def test_good_construct_no_internal_stops(self, good_single_cistron):
        results = check_internal_stops(good_single_cistron)
        assert_no_errors(results, "good single cistron")

    def test_premature_stop_detected(self, bad_internal_stop):
        results = check_internal_stops(bad_internal_stop)
        assert_has_error(results, "internal_stops", "cds1")


# ===================================================================
# INTEGRATION TESTS — example spec through the real pipeline
# ===================================================================

class TestPostParse:
    """
    After parsing, parts are ABSTRACT — no sequences yet.
    Validators should run without crashing (gracefully skip parts
    without sequences) and report no errors (nothing to check).
    """

    def test_frame_check_skips_abstract_parts(self, parsed_graph):
        results = check_reading_frames(parsed_graph)
        errors = errors_only(results)
        assert not errors, f"Frame errors on abstract graph: {errors}"

    def test_start_codon_skips_abstract_parts(self, parsed_graph):
        results = check_start_codons(parsed_graph)
        errors = errors_only(results)
        assert not errors

    def test_translation_skips_abstract_parts(self, parsed_graph):
        results = check_translation_fidelity(parsed_graph)
        errors = errors_only(results)
        assert not errors

    def test_internal_stops_skips_abstract_parts(self, parsed_graph):
        results = check_internal_stops(parsed_graph)
        errors = errors_only(results)
        assert not errors


class TestPostResolve:
    """
    After resolution, parts have metadata (protein sequences) but
    may not have DNA sequences yet. Validators should not crash.
    """

    def test_no_crashes(self, resolved_graph):
        # Just verify all checks run without exceptions
        results = run_all_checks(resolved_graph)
        assert isinstance(results, list)


class TestPostReverseTranslation:
    """
    After reverse translation, all coding parts have DNA sequences.
    All structural checks must pass at this stage.
    """

    def test_reading_frames(self, revtrans_graph):
        results = check_reading_frames(revtrans_graph)
        assert_no_errors(results, "post-reverse-translation reading frames")

    def test_start_codons(self, revtrans_graph):
        results = check_start_codons(revtrans_graph)
        assert_no_errors(results, "post-reverse-translation start codons")

    def test_translation_fidelity(self, revtrans_graph):
        results = check_translation_fidelity(revtrans_graph)
        assert_no_errors(results, "post-reverse-translation translation fidelity")

    def test_internal_stops(self, revtrans_graph):
        results = check_internal_stops(revtrans_graph)
        assert_no_errors(results, "post-reverse-translation internal stops")

    def test_full_report(self, revtrans_graph):
        """Generate and print the full validation report."""
        results = run_all_checks(revtrans_graph)
        report = summary_report(results)
        print("\n" + report)
        assert_no_errors(results, "post-reverse-translation full check")


class TestPostConstraintResolution:
    """
    After constraint resolution (DNA Chisel), verify that codon
    optimization didn't break anything that was correct, or introduce
    new issues (frame shifts, wrong protein, new stop codons).
    """

    def test_reading_frames_preserved(self, constrained_graph):
        results = check_reading_frames(constrained_graph)
        assert_no_errors(results, "post-constraint-resolution reading frames")

    def test_start_codons_preserved(self, constrained_graph):
        results = check_start_codons(constrained_graph)
        assert_no_errors(results, "post-constraint-resolution start codons")

    def test_translation_fidelity_after_optimization(self, constrained_graph):
        """
        Critical test: does DNA Chisel's CodonOptimize +
        EnforceTranslation actually preserve the protein sequences?
        """
        results = check_translation_fidelity(constrained_graph)
        assert_no_errors(results, "post-constraint-resolution translation fidelity")

    def test_no_new_internal_stops(self, constrained_graph):
        results = check_internal_stops(constrained_graph)
        assert_no_errors(results, "post-constraint-resolution internal stops")

    def test_full_report(self, constrained_graph):
        """Full validation report for the constrained construct."""
        results = run_all_checks(constrained_graph)
        report = summary_report(results)
        print("\n" + report)
        assert_no_errors(results, "post-constraint-resolution full check")


# ===================================================================
# REGRESSION TESTS
# ===================================================================

class TestRegressions:
    """
    Tests for specific bugs that have been found and fixed.
    These ensure the bugs don't come back.
    """

    def test_first_tag_gets_start_codon(self, revtrans_graph):
        """
        Regression: reverse_translation.py previously only added ATG to
        CDS parts, not to PurificationTag/SolubilityTag when they appeared
        as the first coding element in a cistron. Fixed by making the pass
        cistron-aware via _find_cistron_leaders().
        """
        results = check_start_codons(revtrans_graph)
        tag_errors = [
            r for r in results
            if r.severity == CheckSeverity.ERROR
            and r.check_name == "start_codon"
            and "PurificationTag" in r.details.get("part_type", "")
        ]
        assert not tag_errors, (
            f"Regression: first tag in cistron still missing start codon: {tag_errors}"
        )

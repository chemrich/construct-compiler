"""
Regression tests: every example spec must compile and pass all checks.

These tests ensure that pipeline changes don't silently break existing
construct designs. Specs in both examples/ and examples/agent_generated/
are discovered automatically — so every design the agent validates in
conversation becomes a permanent regression test case.

Run fast (skip codon optimization):
    pytest tests/test_harness_regression.py -v

Run full (including codon optimization):
    pytest tests/test_harness_regression.py -v --include-slow
"""

from __future__ import annotations

from pathlib import Path

import pytest

from construct_compiler.validation.harness import evaluate_spec


# ---------------------------------------------------------------------------
# Discover all example specs (hand-written + agent-generated)
# ---------------------------------------------------------------------------

EXAMPLES_DIR = Path(__file__).parent.parent / "examples"
AGENT_DIR = EXAMPLES_DIR / "agent_generated"

# Collect all YAML specs from both directories
EXAMPLE_SPECS = sorted(EXAMPLES_DIR.glob("*.yaml"))
if AGENT_DIR.exists():
    EXAMPLE_SPECS += sorted(AGENT_DIR.glob("*.yaml"))


# Skip the module entirely if there are no specs (avoids pytest warning)
if not EXAMPLE_SPECS:
    pytest.skip("No example specs found", allow_module_level=True)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

@pytest.fixture(params=EXAMPLE_SPECS, ids=lambda p: p.stem)
def example_spec(request):
    return request.param


def test_example_compiles_and_passes(example_spec):
    """Every example spec must compile cleanly and pass all 4 checks.

    This test skips the codon optimization pass for speed. The full
    pipeline is tested separately with the @pytest.mark.slow marker.
    """
    result = evaluate_spec(example_spec, skip_constraints=True)
    assert result.passed, (
        f"{example_spec.name} failed with score {result.score}:\n"
        f"{result.summary()}"
    )
    assert result.error_count == 0


@pytest.mark.slow
def test_example_with_constraints(example_spec):
    """Full pipeline including codon optimization."""
    result = evaluate_spec(example_spec, skip_constraints=False)
    assert result.passed, (
        f"{example_spec.name} failed (with constraints):\n"
        f"{result.summary()}"
    )


def test_example_has_valid_metadata(example_spec):
    """Each spec should produce a named construct with measurable properties."""
    result = evaluate_spec(example_spec, skip_constraints=True)
    if result.compile_error:
        pytest.skip(f"Compile error: {result.compile_error}")
    assert result.construct_name, "Construct name should not be empty"
    assert result.insert_length_bp is not None and result.insert_length_bp > 0
    assert result.cistron_count >= 1

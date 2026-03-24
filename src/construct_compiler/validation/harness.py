"""
Test harness for automated construct validation.

Compiles one or many YAML specs through the full pipeline, runs all
validity checks, and returns structured results suitable for:

  1. Agent-in-the-loop design iteration (compile → validate → score → tweak)
  2. Batch evaluation of design variants
  3. CI/CD gating (exit code reflects pass/fail)

Usage as a library:

    from construct_compiler.validation.harness import evaluate_spec, evaluate_batch

    result = evaluate_spec("spec.yaml")
    print(result.passed)       # True/False
    print(result.score)        # 0.0–1.0 composite score
    print(result.error_count)  # number of ERROR-severity issues

    results = evaluate_batch(["v1.yaml", "v2.yaml", "v3.yaml"])
    ranked = sorted(results, key=lambda r: r.score, reverse=True)

Usage from CLI:

    construct-compiler check spec.yaml
    construct-compiler check variants/*.yaml --json
"""

from __future__ import annotations

import json
import logging
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional, Union

from ..core.graph import ConstructGraph
from ..frontend.parser import parse_spec
from ..passes.part_resolution import resolve_parts
from ..passes.reverse_translation import reverse_translate
from ..passes.constraint_resolution import resolve_constraints
from .construct_checks import (
    CheckResult,
    CheckSeverity,
    check_reading_frames,
    check_start_codons,
    check_translation_fidelity,
    check_internal_stops,
    run_all_checks,
    errors_only,
    summary_report,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Result types
# ---------------------------------------------------------------------------

@dataclass
class StageSnapshot:
    """Validation results captured at a specific pipeline stage."""
    stage: str
    errors: int = 0
    warnings: int = 0
    passed: int = 0
    details: list[dict] = field(default_factory=list)


@dataclass
class EvalResult:
    """Complete evaluation of a single construct spec."""
    spec_path: Optional[str] = None
    construct_name: str = ""
    passed: bool = False
    score: float = 0.0          # 0.0 (broken) to 1.0 (all checks clean)
    error_count: int = 0
    warning_count: int = 0
    check_count: int = 0
    insert_length_bp: Optional[int] = None
    cistron_count: int = 0
    compile_time_s: float = 0.0
    stages: list[StageSnapshot] = field(default_factory=list)
    errors: list[dict] = field(default_factory=list)   # ERROR-severity details
    warnings: list[dict] = field(default_factory=list)  # WARNING-severity details
    compile_error: Optional[str] = None  # if compilation itself failed

    def to_dict(self) -> dict:
        return asdict(self)

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent, default=str)

    def summary(self) -> str:
        status = "PASS" if self.passed else "FAIL"
        lines = [
            f"[{status}] {self.construct_name or self.spec_path or 'unnamed'}",
            f"  Score: {self.score:.2f}  |  "
            f"Errors: {self.error_count}  |  "
            f"Warnings: {self.warning_count}  |  "
            f"Checks: {self.check_count}",
        ]
        if self.insert_length_bp:
            lines.append(
                f"  Insert: {self.insert_length_bp} bp  |  "
                f"Cistrons: {self.cistron_count}  |  "
                f"Compile: {self.compile_time_s:.2f}s"
            )
        if self.compile_error:
            lines.append(f"  Compile error: {self.compile_error}")
        if self.errors:
            lines.append("  Errors:")
            for e in self.errors:
                lines.append(f"    - [{e['check_name']}] {e['part_id']}: {e['message']}")
        if self.warnings:
            lines.append("  Warnings:")
            for w in self.warnings:
                lines.append(f"    - [{w['check_name']}] {w['part_id']}: {w['message']}")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Core evaluation
# ---------------------------------------------------------------------------

def evaluate_spec(
    spec: Union[str, Path, dict],
    skip_constraints: bool = False,
    validate_intermediate: bool = True,
) -> EvalResult:
    """
    Compile a spec through the full pipeline and run all validity checks.

    Args:
        spec: YAML file path, YAML string, or parsed dict
        skip_constraints: skip DNA Chisel pass (faster, less thorough)
        validate_intermediate: also check after each pass (not just final)

    Returns:
        EvalResult with pass/fail, score, and detailed check results
    """
    result = EvalResult()

    if isinstance(spec, (str, Path)) and Path(spec).exists():
        result.spec_path = str(spec)

    t0 = time.monotonic()

    # --- Compile through the pipeline, validating at each stage ---
    try:
        # Stage 1: Parse
        graph = parse_spec(spec)
        result.construct_name = graph.name

        # Stage 2: Resolve parts
        graph = resolve_parts(graph)
        if validate_intermediate:
            _snapshot_stage(graph, "post_resolve", result)

        # Stage 3: Reverse translate
        graph = reverse_translate(graph)
        if validate_intermediate:
            _snapshot_stage(graph, "post_revtrans", result)

        # Stage 4: Constraint resolution
        if not skip_constraints:
            graph = resolve_constraints(graph)
        stage_name = "post_constraint" if not skip_constraints else "post_revtrans_final"
        _snapshot_stage(graph, stage_name, result)

        # Capture construct metadata
        result.insert_length_bp = graph.total_insert_length()
        result.cistron_count = len(graph.cistrons())

    except Exception as e:
        result.compile_error = str(e)
        result.passed = False
        result.score = 0.0
        result.compile_time_s = time.monotonic() - t0
        return result

    result.compile_time_s = time.monotonic() - t0

    # --- Run final validation ---
    all_checks = run_all_checks(graph)
    result.check_count = len(all_checks)
    result.error_count = sum(1 for c in all_checks if c.severity == CheckSeverity.ERROR)
    result.warning_count = sum(1 for c in all_checks if c.severity == CheckSeverity.WARNING)
    result.passed = result.error_count == 0

    result.errors = [
        {"check_name": c.check_name, "part_id": c.part_id,
         "message": c.message, "details": c.details}
        for c in all_checks if c.severity == CheckSeverity.ERROR
    ]
    result.warnings = [
        {"check_name": c.check_name, "part_id": c.part_id,
         "message": c.message, "details": c.details}
        for c in all_checks if c.severity == CheckSeverity.WARNING
    ]

    # Compute composite score
    result.score = _compute_score(all_checks)

    return result


def evaluate_batch(
    specs: list[Union[str, Path, dict]],
    skip_constraints: bool = False,
    validate_intermediate: bool = False,
) -> list[EvalResult]:
    """
    Evaluate multiple specs and return results sorted by score (best first).

    Args:
        specs: list of YAML file paths, strings, or dicts
        skip_constraints: skip DNA Chisel (faster batch runs)
        validate_intermediate: check at each pass (slower but more diagnostic)

    Returns:
        List of EvalResults, sorted by score descending
    """
    results = []
    for i, spec in enumerate(specs):
        logger.info(f"Evaluating spec {i+1}/{len(specs)}...")
        r = evaluate_spec(
            spec,
            skip_constraints=skip_constraints,
            validate_intermediate=validate_intermediate,
        )
        results.append(r)

    # Sort by score descending, then by error count ascending
    results.sort(key=lambda r: (-r.score, r.error_count))
    return results


# ---------------------------------------------------------------------------
# Internals
# ---------------------------------------------------------------------------

def _snapshot_stage(graph: ConstructGraph, stage_name: str, result: EvalResult) -> None:
    """Run checks at a pipeline stage and record the snapshot."""
    checks = run_all_checks(graph)
    snap = StageSnapshot(
        stage=stage_name,
        errors=sum(1 for c in checks if c.severity == CheckSeverity.ERROR),
        warnings=sum(1 for c in checks if c.severity == CheckSeverity.WARNING),
        passed=sum(1 for c in checks if c.severity == CheckSeverity.INFO),
        details=[
            {"check": c.check_name, "severity": c.severity.name,
             "part": c.part_id, "message": c.message}
            for c in checks if c.severity != CheckSeverity.INFO
        ],
    )
    result.stages.append(snap)


def _compute_score(checks: list[CheckResult]) -> float:
    """
    Compute a 0.0–1.0 composite score from check results.

    Scoring:
    - Each check contributes equally
    - ERROR = 0 points, WARNING = 0.5 points, INFO/PASS = 1.0 points
    - Score = total points / number of checks

    A score of 1.0 means all checks passed cleanly.
    A score of 0.0 means every check found an error.
    """
    if not checks:
        return 1.0  # no checks to fail = pass

    points = 0.0
    for c in checks:
        if c.severity == CheckSeverity.INFO:
            points += 1.0
        elif c.severity == CheckSeverity.WARNING:
            points += 0.5
        # ERROR = 0

    return round(points / len(checks), 4)


# ---------------------------------------------------------------------------
# Batch report formatting
# ---------------------------------------------------------------------------

def batch_summary(results: list[EvalResult]) -> str:
    """Human-readable summary of a batch evaluation."""
    lines = [
        f"Batch Evaluation: {len(results)} specs",
        "=" * 60,
    ]

    passed = sum(1 for r in results if r.passed)
    failed = len(results) - passed
    lines.append(f"Passed: {passed}  |  Failed: {failed}")
    lines.append("")

    for i, r in enumerate(results):
        rank = i + 1
        status = "PASS" if r.passed else "FAIL"
        name = r.construct_name or r.spec_path or "unnamed"
        lines.append(
            f"  {rank}. [{status}] {name}  "
            f"score={r.score:.2f}  "
            f"errors={r.error_count}  "
            f"warnings={r.warning_count}"
        )
        if r.compile_error:
            lines.append(f"     compile error: {r.compile_error}")

    return "\n".join(lines)

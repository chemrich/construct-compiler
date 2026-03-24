"""
Construct validation — in silico checks that a compiled construct
is biologically viable before sending it to synthesis.

Validators return structured CheckResult objects that can be:
- Consumed by pytest assertions
- Aggregated into a human-readable report
- Fed into an automated design-evaluate loop
"""

from .construct_checks import (
    CheckResult,
    CheckSeverity,
    check_reading_frames,
    check_start_codons,
    check_translation_fidelity,
    check_internal_stops,
    run_all_checks,
)
from .harness import (
    EvalResult,
    evaluate_spec,
    evaluate_batch,
    batch_summary,
)

__all__ = [
    "CheckResult",
    "CheckSeverity",
    "check_reading_frames",
    "check_start_codons",
    "check_translation_fidelity",
    "check_internal_stops",
    "run_all_checks",
    "EvalResult",
    "evaluate_spec",
    "evaluate_batch",
    "batch_summary",
]

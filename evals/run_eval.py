#!/usr/bin/env python3
"""
Eval runner: natural language prompt → LLM-generated spec → harness validation.

Drives the full loop:
  1. Load prompt corpus from prompt_corpus.yaml
  2. For each prompt, call the Anthropic API to generate a YAML spec
  3. Run the spec through the validation harness
  4. Check expected properties (host, cistron count, required parts, pass/fail)
  5. Write structured results to evals/results/

Usage:
    # Run all prompts (requires ANTHROPIC_API_KEY)
    python evals/run_eval.py

    # Run a single prompt by ID
    python evals/run_eval.py --id basic_gfp

    # Run a category
    python evals/run_eval.py --category polycistronic

    # Skip LLM generation and re-evaluate previously generated specs
    python evals/run_eval.py --reeval

    # Use a specific model
    python evals/run_eval.py --model claude-sonnet-4-20250514

    # Dry run — show what would be tested without calling the API
    python evals/run_eval.py --dry-run

    # Specify number of retries per prompt (LLM generates spec again on failure)
    python evals/run_eval.py --retries 2

    # Run with parallel API calls (default: 1 = sequential)
    python evals/run_eval.py --concurrency 10
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import re
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Any, Optional

# ---------------------------------------------------------------------------
# Rate limiter — prevents API throttling when running concurrent workers
# ---------------------------------------------------------------------------

class _RateLimiter:
    """Simple token-bucket rate limiter (thread-safe)."""

    def __init__(self, requests_per_minute: float):
        self._min_interval = 60.0 / requests_per_minute if requests_per_minute > 0 else 0
        self._lock = threading.Lock()
        self._last_request = 0.0

    def wait(self):
        """Block until it's safe to send the next request."""
        if self._min_interval <= 0:
            return
        with self._lock:
            now = time.monotonic()
            wait_time = self._last_request + self._min_interval - now
            if wait_time > 0:
                time.sleep(wait_time)
            self._last_request = time.monotonic()

# Global rate limiter — set from CLI via --rpm flag
_rate_limiter: _RateLimiter | None = None

import yaml

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from construct_compiler.validation.harness import evaluate_spec, EvalResult

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

EVALS_DIR = Path(__file__).parent
CORPUS_PATH = EVALS_DIR / "prompt_corpus.yaml"
SYSTEM_PROMPT_PATH = EVALS_DIR / "spec_generation_prompt.txt"
RESULTS_DIR = EVALS_DIR / "results"
SPECS_DIR = EVALS_DIR / "generated_specs"


# ---------------------------------------------------------------------------
# Data types
# ---------------------------------------------------------------------------

@dataclass
class PromptCase:
    id: str
    category: str
    difficulty: str
    prompt: str
    expect: dict = field(default_factory=dict)
    notes: str = ""


@dataclass
class EvalOutcome:
    prompt_id: str
    category: str
    difficulty: str
    prompt: str

    # LLM generation
    generated_spec: Optional[dict] = None
    generation_error: Optional[str] = None
    raw_llm_output: str = ""
    model: str = ""
    generation_time_s: float = 0.0

    # Harness validation
    harness_passed: Optional[bool] = None
    harness_score: Optional[float] = None
    harness_error_count: int = 0
    harness_warning_count: int = 0
    compile_error: Optional[str] = None
    harness_errors: list[dict] = field(default_factory=list)
    harness_summary: str = ""

    # Expectation checks
    expectation_results: dict = field(default_factory=dict)
    all_expectations_met: bool = False

    # Meta
    attempt: int = 1
    timestamp: str = ""

    def to_dict(self) -> dict:
        return asdict(self)


# ---------------------------------------------------------------------------
# Load corpus
# ---------------------------------------------------------------------------

def load_corpus(
    corpus_path: Path = CORPUS_PATH,
    filter_id: str | None = None,
    filter_category: str | None = None,
) -> list[PromptCase]:
    """Load prompt cases from the corpus YAML."""
    data = yaml.safe_load(corpus_path.read_text())
    cases = []
    for entry in data["prompts"]:
        case = PromptCase(
            id=entry["id"],
            category=entry["category"],
            difficulty=entry.get("difficulty", "medium"),
            prompt=entry["prompt"],
            expect=entry.get("expect", {}),
            notes=entry.get("notes", ""),
        )
        if filter_id and case.id != filter_id:
            continue
        if filter_category and case.category != filter_category:
            continue
        cases.append(case)
    return cases


# ---------------------------------------------------------------------------
# LLM spec generation
# ---------------------------------------------------------------------------

def generate_spec_via_api(
    prompt: str,
    model: str = "claude-sonnet-4-20250514",
    system_prompt: str | None = None,
) -> tuple[dict | None, str, float, str | None]:
    """
    Call the Anthropic API to generate a YAML spec from a natural language prompt.

    Returns: (parsed_spec_dict, raw_output, time_seconds, error_or_none)
    """
    try:
        import anthropic
    except ImportError:
        return None, "", 0.0, "anthropic package not installed. Run: pip install anthropic"

    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        return None, "", 0.0, "ANTHROPIC_API_KEY environment variable not set"

    if system_prompt is None:
        system_prompt = SYSTEM_PROMPT_PATH.read_text()

    client = anthropic.Anthropic(api_key=api_key)

    # Respect rate limit before firing the request
    if _rate_limiter:
        _rate_limiter.wait()

    t0 = time.monotonic()
    max_api_retries = 5
    for api_attempt in range(max_api_retries):
        try:
            response = client.messages.create(
                model=model,
                max_tokens=4096,
                system=system_prompt,
                messages=[{"role": "user", "content": prompt}],
            )
            raw = response.content[0].text
            elapsed = time.monotonic() - t0
            break
        except Exception as e:
            err_str = str(e).lower()
            is_rate_limit = "rate" in err_str or "429" in err_str or "overloaded" in err_str
            if is_rate_limit and api_attempt < max_api_retries - 1:
                backoff = 2 ** api_attempt * 5  # 5s, 10s, 20s, 40s
                logger.warning(f"Rate limited, backing off {backoff}s (attempt {api_attempt + 1})")
                time.sleep(backoff)
                continue
            return None, "", time.monotonic() - t0, f"API error: {e}"

    # Parse the YAML from the response
    spec, parse_error = _extract_yaml(raw)
    return spec, raw, elapsed, parse_error


def _extract_yaml(raw: str) -> tuple[dict | None, str | None]:
    """Extract and parse YAML from LLM output, handling markdown fences."""
    # Strip markdown code fences if present
    cleaned = raw.strip()
    cleaned = re.sub(r'^```ya?ml\s*\n', '', cleaned, flags=re.MULTILINE)
    cleaned = re.sub(r'\n```\s*$', '', cleaned, flags=re.MULTILINE)
    # Also strip bare ``` at start/end
    cleaned = re.sub(r'^```\s*\n', '', cleaned, flags=re.MULTILINE)
    cleaned = re.sub(r'\n```\s*$', '', cleaned, flags=re.MULTILINE)
    cleaned = cleaned.strip()

    # Handle multi-document YAML: take only the first document
    if "\n---\n" in cleaned or cleaned.startswith("---\n"):
        parts = re.split(r'\n---\n', cleaned)
        # Use the first non-empty document that looks like it has 'construct'
        for part in parts:
            part = part.strip().lstrip("---").strip()
            if part and ("construct" in part or "name" in part):
                cleaned = part
                break
        else:
            cleaned = parts[0].strip().lstrip("---").strip()

    try:
        spec = yaml.safe_load(cleaned)
        if not isinstance(spec, dict):
            return None, f"YAML parsed but is not a dict (got {type(spec).__name__})"
        # Normalize: ensure top-level 'construct' key
        if "construct" not in spec:
            # Maybe the LLM output the construct contents directly
            spec = {"construct": spec}
        return spec, None
    except yaml.YAMLError as e:
        return None, f"YAML parse error: {e}"


# ---------------------------------------------------------------------------
# Expectation checking
# ---------------------------------------------------------------------------

def check_expectations(
    case: PromptCase,
    spec: dict | None,
    eval_result: EvalResult | None,
) -> dict[str, dict]:
    """
    Check whether the generated spec meets the expected properties.

    Returns a dict of {check_name: {expected, actual, passed}}.
    """
    results = {}
    expect = case.expect
    if not expect:
        return results

    # Host check
    if "host" in expect:
        actual_host = ""
        if spec:
            root = spec.get("construct", spec)
            actual_host = root.get("host", "")
        expected = expect["host"]
        passed = _hosts_match(expected, actual_host or "")
        results["host"] = {"expected": expected, "actual": actual_host, "passed": passed}

    # Cistron count
    if "min_cistrons" in expect and eval_result:
        actual = eval_result.cistron_count
        expected = expect["min_cistrons"]
        results["min_cistrons"] = {
            "expected": f">= {expected}", "actual": actual,
            "passed": actual >= expected,
        }
    if "max_cistrons" in expect and eval_result:
        actual = eval_result.cistron_count
        expected = expect["max_cistrons"]
        results["max_cistrons"] = {
            "expected": f"<= {expected}", "actual": actual,
            "passed": actual <= expected,
        }

    # Required part types
    if "must_have_parts" in expect and spec:
        # We check the spec structure since the harness result doesn't list part types
        found_types = _extract_part_types_from_spec(spec)
        for required in expect["must_have_parts"]:
            present = required in found_types
            results[f"has_{required}"] = {
                "expected": True, "actual": present,
                "passed": present,
            }

    # Must pass harness
    if expect.get("must_pass") and eval_result:
        results["must_pass"] = {
            "expected": True, "actual": eval_result.passed,
            "passed": eval_result.passed,
        }

    return results


def _hosts_match(expected: str, actual: str) -> bool:
    """
    Flexible host matching. Handles:
    - Substring: "e_coli" matches "e_coli_bl21_de3"
    - Category equivalence: "mammalian" matches "HEK293", "CHO", etc.
    - Cell line specificity: "HEK293" matches "mammalian", "HEK293T", etc.
    """
    if not actual:
        return False

    # Direct substring match (original behavior)
    if expected in actual or actual in expected:
        return True

    # Normalize
    exp_low = expected.lower().replace("-", "").replace("_", "")
    act_low = actual.lower().replace("-", "").replace("_", "")

    if exp_low in act_low or act_low in exp_low:
        return True

    # Category mappings
    mammalian_hosts = {"mammalian", "hek293", "hek293t", "cho", "hela",
                       "cos7", "nih3t3", "vero", "a549", "jurkat",
                       "k562", "u2os", "sf9", "sf21", "hi5"}
    bacterial_hosts = {"ecoli", "ecolibl21", "ecolibl21de3", "ecolik12",
                       "ecolirosetta", "ecolishuffle", "ecoliarctic"}

    exp_norm = exp_low
    act_norm = act_low

    # If both are in the same category, it's a match
    if exp_norm in mammalian_hosts and act_norm in mammalian_hosts:
        return True
    if exp_norm in bacterial_hosts and act_norm in bacterial_hosts:
        return True

    return False


def _extract_part_types_from_spec(spec: dict) -> set[str]:
    """Walk the spec to find what part types are present."""
    types = set()
    root = spec.get("construct", spec)
    cassette = root.get("cassette", [])

    for element in cassette:
        if not isinstance(element, dict):
            continue
        for key, value in element.items():
            if key == "cistron":
                chain = value.get("chain", []) if isinstance(value, dict) else []
                gene = value.get("gene") if isinstance(value, dict) else None
                n_tag = value.get("n_tag") if isinstance(value, dict) else None
                c_tag = value.get("c_tag") if isinstance(value, dict) else None

                if gene:
                    types.add("CDS")
                if n_tag or c_tag:
                    types.add("PurificationTag")  # simplification

                for item in chain:
                    if isinstance(item, dict):
                        for k in item.keys():
                            if k == "tag":
                                types.add("PurificationTag")
                            elif k == "solubility_tag":
                                types.add("SolubilityTag")
                            elif k == "cleavage_site":
                                types.add("CleavageSite")
                            elif k == "linker":
                                types.add("Linker")
                            elif k == "gene":
                                types.add("CDS")
    return types


# ---------------------------------------------------------------------------
# Core eval loop
# ---------------------------------------------------------------------------

def _eval_single_case(
    case: PromptCase,
    index: int,
    total: int,
    model: str,
    retries: int,
    reeval: bool,
    skip_constraints: bool,
    system_prompt: str,
    print_lock: threading.Lock | None = None,
) -> EvalOutcome | None:
    """
    Evaluate a single prompt case (thread-safe).

    Generates a spec via the API, validates it, and checks expectations.
    """
    def _print(*args, **kwargs):
        if print_lock:
            with print_lock:
                print(*args, **kwargs)
        else:
            print(*args, **kwargs)

    _print(f"\n[{index}/{total}] {case.id} ({case.category}/{case.difficulty})")
    _print(f"  Prompt: {case.prompt[:80]}...")

    spec_path = SPECS_DIR / f"{case.id}.yaml"
    best_outcome = None

    for attempt in range(1, retries + 1):
        outcome = EvalOutcome(
            prompt_id=case.id,
            category=case.category,
            difficulty=case.difficulty,
            prompt=case.prompt,
            model=model,
            attempt=attempt,
            timestamp=datetime.now().isoformat(),
        )

        # --- Generate spec ---
        if reeval and spec_path.exists():
            raw = spec_path.read_text()
            spec, parse_error = _extract_yaml(raw)
            outcome.raw_llm_output = raw
            outcome.generated_spec = spec
            outcome.generation_error = parse_error
            _print(f"  [reeval] Loaded saved spec from {spec_path.name}")
        else:
            spec, raw, elapsed, gen_error = generate_spec_via_api(
                case.prompt, model=model, system_prompt=system_prompt,
            )
            outcome.raw_llm_output = raw
            outcome.generated_spec = spec
            outcome.generation_error = gen_error
            outcome.generation_time_s = elapsed

            if gen_error:
                _print(f"  [FAIL] Generation error: {gen_error}")
            else:
                _print(f"  Generated spec in {elapsed:.1f}s")
                spec_path.write_text(raw)

        # --- Validate with harness ---
        if outcome.generated_spec:
            try:
                eval_result = evaluate_spec(
                    outcome.generated_spec,
                    skip_constraints=skip_constraints,
                )
                outcome.harness_passed = eval_result.passed
                outcome.harness_score = eval_result.score
                outcome.harness_error_count = eval_result.error_count
                outcome.harness_warning_count = eval_result.warning_count
                outcome.compile_error = eval_result.compile_error
                outcome.harness_errors = eval_result.errors
                outcome.harness_summary = eval_result.summary()

                status = "PASS" if eval_result.passed else "FAIL"
                _print(f"  [{status}] score={eval_result.score:.2f} "
                       f"errors={eval_result.error_count} "
                       f"warnings={eval_result.warning_count}")
                if eval_result.compile_error:
                    _print(f"  Compile error: {eval_result.compile_error}")
                for err in eval_result.errors:
                    _print(f"    ERROR: [{err['check_name']}] {err['message']}")
            except Exception as e:
                outcome.compile_error = str(e)
                _print(f"  [CRASH] Harness exception: {e}")

        # --- Check expectations ---
        eval_r = None
        if outcome.harness_passed is not None:
            eval_r = EvalResult()
            eval_r.passed = outcome.harness_passed
            eval_r.cistron_count = 0
            try:
                temp = evaluate_spec(outcome.generated_spec, skip_constraints=True)
                eval_r.cistron_count = temp.cistron_count
            except Exception:
                pass

        outcome.expectation_results = check_expectations(
            case, outcome.generated_spec, eval_r,
        )
        outcome.all_expectations_met = all(
            r["passed"] for r in outcome.expectation_results.values()
        ) if outcome.expectation_results else (outcome.harness_passed or False)

        for check_name, check_result in outcome.expectation_results.items():
            status = "OK" if check_result["passed"] else "MISS"
            _print(f"    [{status}] {check_name}: "
                   f"expected={check_result['expected']} "
                   f"actual={check_result['actual']}")

        best_outcome = outcome

        if outcome.harness_passed and outcome.all_expectations_met:
            break

    return best_outcome


def run_eval(
    cases: list[PromptCase],
    model: str = "claude-sonnet-4-20250514",
    retries: int = 1,
    reeval: bool = False,
    dry_run: bool = False,
    skip_constraints: bool = True,
    concurrency: int = 1,
) -> list[EvalOutcome]:
    """
    Run the eval for a list of prompt cases.

    Args:
        cases: prompt cases to evaluate
        model: Anthropic model to use
        retries: number of LLM generation attempts per prompt
        reeval: if True, use previously saved specs instead of calling the API
        dry_run: if True, just print what would be tested
        skip_constraints: skip codon optimization for speed
        concurrency: number of parallel workers (1 = sequential)
    """
    SPECS_DIR.mkdir(parents=True, exist_ok=True)

    system_prompt = SYSTEM_PROMPT_PATH.read_text()

    if dry_run:
        for i, case in enumerate(cases):
            print(f"\n[{i+1}/{len(cases)}] {case.id} ({case.category}/{case.difficulty})")
            print(f"  Prompt: {case.prompt[:80]}...")
            print("  [DRY RUN] Would generate spec and validate")
        return []

    total = len(cases)

    if concurrency <= 1:
        # Sequential path (original behavior, cleaner output)
        outcomes = []
        for i, case in enumerate(cases):
            outcome = _eval_single_case(
                case, i + 1, total, model, retries, reeval,
                skip_constraints, system_prompt,
            )
            if outcome:
                outcomes.append(outcome)
        return outcomes

    # Parallel path
    print(f"Running with {concurrency} parallel workers")
    print_lock = threading.Lock()
    outcomes: list[EvalOutcome] = []
    completed = 0

    def _worker(args):
        nonlocal completed
        i, case = args
        result = _eval_single_case(
            case, i + 1, total, model, retries, reeval,
            skip_constraints, system_prompt, print_lock,
        )
        with print_lock:
            completed += 1
            pct = completed / total * 100
            passed_so_far = sum(1 for o in outcomes if o.harness_passed) + (1 if result and result.harness_passed else 0)
            print(f"\n  >>> Progress: {completed}/{total} ({pct:.0f}%) — {passed_so_far} passed so far")
        return result

    with ThreadPoolExecutor(max_workers=concurrency) as executor:
        futures = {
            executor.submit(_worker, (i, case)): i
            for i, case in enumerate(cases)
        }
        for future in as_completed(futures):
            result = future.result()
            if result:
                outcomes.append(result)

    # Sort outcomes back to original corpus order
    case_order = {case.id: i for i, case in enumerate(cases)}
    outcomes.sort(key=lambda o: case_order.get(o.prompt_id, 0))

    return outcomes


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def write_results(outcomes: list[EvalOutcome], run_name: str | None = None) -> Path:
    """Write eval results to a JSON file."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    if run_name is None:
        run_name = datetime.now().strftime("%Y%m%d_%H%M%S")

    results_path = RESULTS_DIR / f"eval_{run_name}.json"

    data = {
        "run_name": run_name,
        "timestamp": datetime.now().isoformat(),
        "total_prompts": len(outcomes),
        "summary": _compute_summary(outcomes),
        "outcomes": [o.to_dict() for o in outcomes],
    }

    results_path.write_text(json.dumps(data, indent=2, default=str))
    return results_path


def _compute_summary(outcomes: list[EvalOutcome]) -> dict:
    """Compute aggregate stats from eval outcomes."""
    total = len(outcomes)
    if total == 0:
        return {}

    generated = sum(1 for o in outcomes if o.generated_spec is not None)
    gen_failed = sum(1 for o in outcomes if o.generation_error is not None)
    harness_passed = sum(1 for o in outcomes if o.harness_passed)
    harness_failed = sum(1 for o in outcomes if o.harness_passed is False)
    compile_errors = sum(1 for o in outcomes if o.compile_error)
    expectations_met = sum(1 for o in outcomes if o.all_expectations_met)

    # By category
    categories = {}
    for o in outcomes:
        cat = o.category
        if cat not in categories:
            categories[cat] = {"total": 0, "passed": 0, "expectations_met": 0}
        categories[cat]["total"] += 1
        if o.harness_passed:
            categories[cat]["passed"] += 1
        if o.all_expectations_met:
            categories[cat]["expectations_met"] += 1

    # By difficulty
    difficulties = {}
    for o in outcomes:
        diff = o.difficulty
        if diff not in difficulties:
            difficulties[diff] = {"total": 0, "passed": 0}
        difficulties[diff]["total"] += 1
        if o.harness_passed:
            difficulties[diff]["passed"] += 1

    # Error classification
    error_types = {}
    for o in outcomes:
        for err in o.harness_errors:
            check = err.get("check_name", "unknown")
            error_types[check] = error_types.get(check, 0) + 1

    return {
        "total": total,
        "spec_generated": generated,
        "generation_failed": gen_failed,
        "harness_passed": harness_passed,
        "harness_failed": harness_failed,
        "compile_errors": compile_errors,
        "all_expectations_met": expectations_met,
        "pass_rate": round(harness_passed / total * 100, 1) if total else 0,
        "expectation_rate": round(expectations_met / total * 100, 1) if total else 0,
        "by_category": categories,
        "by_difficulty": difficulties,
        "error_types": error_types,
    }


def print_report(outcomes: list[EvalOutcome]) -> None:
    """Print a human-readable eval report to stdout."""
    summary = _compute_summary(outcomes)
    total = summary["total"]

    print("\n" + "=" * 70)
    print("EVAL REPORT")
    print("=" * 70)
    print(f"Total prompts:        {total}")
    print(f"Specs generated:      {summary['spec_generated']}/{total}")
    print(f"Generation failures:  {summary['generation_failed']}")
    print(f"Harness passed:       {summary['harness_passed']}/{total} "
          f"({summary['pass_rate']}%)")
    print(f"Compile errors:       {summary['compile_errors']}")
    print(f"Expectations met:     {summary['all_expectations_met']}/{total} "
          f"({summary['expectation_rate']}%)")

    print(f"\nBy category:")
    for cat, stats in summary.get("by_category", {}).items():
        print(f"  {cat:20s}  {stats['passed']}/{stats['total']} passed  "
              f"  {stats['expectations_met']}/{stats['total']} expectations met")

    print(f"\nBy difficulty:")
    for diff, stats in summary.get("by_difficulty", {}).items():
        print(f"  {diff:10s}  {stats['passed']}/{stats['total']} passed")

    if summary.get("error_types"):
        print(f"\nError types:")
        for check, count in sorted(summary["error_types"].items(),
                                    key=lambda x: -x[1]):
            print(f"  {check:30s}  {count}")

    # List failures
    failures = [o for o in outcomes if not o.all_expectations_met]
    if failures:
        print(f"\nFailed cases ({len(failures)}):")
        for o in failures:
            print(f"  {o.prompt_id} ({o.category}/{o.difficulty})")
            if o.generation_error:
                print(f"    Generation error: {o.generation_error}")
            if o.compile_error:
                print(f"    Compile error: {o.compile_error}")
            if o.harness_passed is False:
                print(f"    Harness: FAIL (score={o.harness_score})")
                for err in o.harness_errors:
                    print(f"      [{err['check_name']}] {err['message']}")
            for check_name, check_result in o.expectation_results.items():
                if not check_result["passed"]:
                    print(f"    Expectation miss: {check_name} "
                          f"(expected={check_result['expected']}, "
                          f"actual={check_result['actual']})")

    print("\n" + "=" * 70)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Run construct compiler evals: prompt → spec → validate",
    )
    parser.add_argument("--id", help="Run a single prompt by ID")
    parser.add_argument("--category", help="Run prompts in a category")
    parser.add_argument("--model", default="claude-sonnet-4-20250514",
                        help="Anthropic model to use")
    parser.add_argument("--retries", type=int, default=1,
                        help="Number of LLM retries per prompt")
    parser.add_argument("--reeval", action="store_true",
                        help="Re-evaluate previously generated specs")
    parser.add_argument("--dry-run", action="store_true",
                        help="Show what would be tested")
    parser.add_argument("--skip-constraints", action="store_true", default=True,
                        help="Skip codon optimization (default: True)")
    parser.add_argument("--full", action="store_true",
                        help="Run with codon optimization (slow)")
    parser.add_argument("--corpus", help="Path to prompt corpus YAML (default: prompt_corpus.yaml)")
    parser.add_argument("--concurrency", "-j", type=int, default=1,
                        help="Number of parallel workers (default: 1 = sequential)")
    parser.add_argument("--rpm", type=int, default=0,
                        help="Max API requests per minute (0 = unlimited). "
                             "Recommended: 40-50 for most tiers.")
    parser.add_argument("--run-name", help="Name for the results file")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO if args.concurrency > 1 else logging.WARNING)

    # Set up rate limiting
    global _rate_limiter
    if args.rpm > 0:
        _rate_limiter = _RateLimiter(args.rpm)
        print(f"Rate limited to {args.rpm} requests/minute")
    elif args.concurrency > 1:
        # Auto-set a safe default when running concurrent
        auto_rpm = max(30, args.concurrency * 8)
        _rate_limiter = _RateLimiter(auto_rpm)
        print(f"Auto rate limit: {auto_rpm} requests/minute (use --rpm to override)")

    corpus_path = Path(args.corpus) if args.corpus else CORPUS_PATH
    cases = load_corpus(corpus_path=corpus_path, filter_id=args.id, filter_category=args.category)
    if not cases:
        print("No matching prompt cases found.")
        sys.exit(1)

    print(f"Loaded {len(cases)} prompt cases")

    skip = not args.full
    outcomes = run_eval(
        cases,
        model=args.model,
        retries=args.retries,
        reeval=args.reeval,
        dry_run=args.dry_run,
        skip_constraints=skip,
        concurrency=args.concurrency,
    )

    if not args.dry_run and outcomes:
        print_report(outcomes)
        results_path = write_results(outcomes, run_name=args.run_name)
        print(f"\nResults saved to: {results_path}")


if __name__ == "__main__":
    main()

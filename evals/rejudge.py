#!/usr/bin/env python3
"""
rejudge.py — Re-run LLM expectation judging over existing eval results.

Loads one or more results JSON files produced by run_eval.py, re-evaluates
the expectation checks using an LLM judge (Claude or Gemini), and writes
new results files with updated expectation_results, expectations_ok,
expectations_total, expectations_pct, and all_expectations_met fields.

No spec generation is performed — only the judge step is re-run.

Usage:
    # Re-judge a results file with Claude Sonnet (streaming, concurrency 5)
    python evals/rejudge.py evals/results/eval_250_baseline.json \
        --judge-model claude-sonnet-4-20250514

    # Re-judge using the Anthropic Batch API (no concurrency needed, 50% cheaper)
    python evals/rejudge.py evals/results/eval_250_baseline.json \
        --judge-model claude-sonnet-4-20250514 --batch

    # Re-judge multiple files in batch mode
    python evals/rejudge.py evals/results/eval_20260410_22*.json \
        --judge-model claude-sonnet-4-20250514 --batch

    # Re-judge with Gemini (streaming only)
    python evals/rejudge.py evals/results/eval_250_baseline.json \
        --judge-model gemini-2.5-flash

    # Dry-run: print what would be judged without calling the API
    python evals/rejudge.py evals/results/eval_250_baseline.json \
        --judge-model claude-sonnet-4-20250514 --dry-run
"""

import argparse
import json
import logging
import os
import re
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import yaml

logging.basicConfig(level=logging.INFO, format="%(levelname)s:%(name)s:%(message)s")
logger = logging.getLogger(__name__)

EVALS_DIR = Path(__file__).parent
RESULTS_DIR = EVALS_DIR / "results"

# ---------------------------------------------------------------------------
# Import shared helpers from run_eval
# ---------------------------------------------------------------------------

sys.path.insert(0, str(EVALS_DIR.parent / "src"))

from run_eval import (
    PromptCase,
    EvalOutcome,
    check_expectations_with_llm,
    write_results,
    _compute_summary,
    _hosts_match,
)


# ---------------------------------------------------------------------------
# Load results JSON
# ---------------------------------------------------------------------------

def load_results_json(path: Path) -> tuple[list[EvalOutcome], list[PromptCase]]:
    """
    Load outcomes from a results JSON file.
    Returns (outcomes, cases) where cases carry the expect field.
    """
    data = json.loads(path.read_text())
    raw_outcomes = data.get("outcomes", [])

    outcomes = []
    cases = []
    for o in raw_outcomes:
        outcome = EvalOutcome(
            prompt_id=o["prompt_id"],
            category=o.get("category", ""),
            difficulty=o.get("difficulty", ""),
            prompt=o.get("prompt", ""),
            generated_spec=o.get("generated_spec"),
            generation_error=o.get("generation_error"),
            raw_llm_output=o.get("raw_llm_output", ""),
            model=o.get("model", ""),
            generation_time_s=o.get("generation_time_s", 0.0),
            harness_passed=o.get("harness_passed"),
            harness_score=o.get("harness_score"),
            harness_error_count=o.get("harness_error_count", 0),
            harness_warning_count=o.get("harness_warning_count", 0),
            compile_error=o.get("compile_error"),
            harness_errors=o.get("harness_errors", []),
            harness_summary=o.get("harness_summary", ""),
            expectation_results=o.get("expectation_results", {}),
            all_expectations_met=o.get("all_expectations_met", False),
            expectations_ok=o.get("expectations_ok", 0),
            expectations_total=o.get("expectations_total", 0),
            expectations_pct=o.get("expectations_pct"),
            attempt=o.get("attempt", 1),
            timestamp=o.get("timestamp", ""),
        )
        outcomes.append(outcome)

        expect = _reconstruct_expect(o)
        case = PromptCase(
            id=o["prompt_id"],
            category=o.get("category", ""),
            difficulty=o.get("difficulty", ""),
            prompt=o.get("prompt", ""),
            expect=expect,
        )
        cases.append(case)

    return outcomes, cases


def _reconstruct_expect(raw_outcome: dict) -> dict:
    """Reconstruct the expect dict from a saved outcome's expectation_results keys."""
    expect = {}
    exp_res = raw_outcome.get("expectation_results", {})

    must_have = []
    must_not_have = []

    for key, val in exp_res.items():
        if key == "host":
            expect["host"] = val.get("expected", "")
        elif key == "min_cistrons":
            try:
                expect["min_cistrons"] = int(str(val.get("expected", "")).replace(">=", "").strip())
            except ValueError:
                pass
        elif key == "max_cistrons":
            try:
                expect["max_cistrons"] = int(str(val.get("expected", "")).replace("<=", "").strip())
            except ValueError:
                pass
        elif key == "must_pass":
            expect["must_pass"] = True
        elif key.startswith("has_"):
            must_have.append(key[4:])
        elif key.startswith("not_have_"):
            must_not_have.append(key[9:])

    if must_have:
        expect["must_have_parts"] = must_have
    if must_not_have:
        expect["must_not_have_parts"] = must_not_have

    return expect


def _infer_cistron_count(outcome: EvalOutcome) -> int:
    """Best-effort cistron count from saved min_cistrons expectation actual."""
    for key, val in outcome.expectation_results.items():
        if key == "min_cistrons":
            try:
                return int(val.get("actual", 0))
            except (ValueError, TypeError):
                pass
    return 0


# ---------------------------------------------------------------------------
# Prompt building and response parsing (shared between streaming and batch)
# ---------------------------------------------------------------------------

JUDGE_SYSTEM_PROMPT = (
    "You are an expert biological engineer. "
    "Respond only with valid JSON — no markdown fences, no commentary."
)


def _build_judge_prompt(case: PromptCase, outcome: EvalOutcome) -> str | None:
    """
    Build the user-facing judge prompt for a single outcome.
    Returns None if there are no parts to check (deterministic checks only).
    """
    must_have = case.expect.get("must_have_parts", [])
    must_not_have = case.expect.get("must_not_have_parts", [])
    if not must_have and not must_not_have:
        return None
    if outcome.generated_spec is None:
        return None

    spec_yaml = yaml.dump(outcome.generated_spec)

    return f"""You are an expert biological engineer evaluating a generated DNA construct specification against requirements.

User Prompt:
{case.prompt}

Generated YAML Specification:
```yaml
{spec_yaml}
```

We need to check if the following required parts are present or reasonably represented in the specification:
{must_have}

We also need to check if the following parts are NOT present or represented in the specification:
{must_not_have}

For each part in BOTH lists above, determine if the specification includes or reasonably represents it.
IMPORTANT CONSIDERATIONS:
1. The schema used by the generator is limited (e.g., it only has bacterial terminators, specific origins, etc.).
2. If the user asked for something (like a PolyA signal or a specific mammalian promoter) and the model used a reasonable placeholder (like a bacterial terminator or a standard strong promoter) and ideally explained it in a comment, consider it PRESENT.
3. We care about intent and reasonable execution within schema limits.

ADDITIONAL EVALUATION:
4. **Host Compatibility**: Verify if the selected parts (promoters, RBSs, terminators) are generally compatible with the specified host organism. If there is a clear mismatch that would prevent function, note it.
5. **Genetic Circuit Logic**: If the user prompt describes a genetic circuit (e.g., a logic gate, toggle switch, oscillator, feedback loop), evaluate if the generated specification logically implements that circuit.

Respond ONLY with a JSON object mapping each part name from BOTH lists above to a boolean 'passed' and a string 'reason'.
Also include "logic_verification" if the prompt describes a circuit, and "host_compatibility" as a general check.

Example:
{{
  "PolyA_Signal": {{"passed": true, "reason": "Used BBa_B0015 as placeholder for mammalian polyA."}},
  "Resistance": {{"passed": true, "reason": "Kanamycin is in backbone."}},
  "host_compatibility": {{"passed": true, "reason": "Parts are appropriate for E. coli host."}}
}}
"""


def _parse_judge_response(raw: str, case: PromptCase, outcome: EvalOutcome) -> dict:
    """
    Parse the LLM's JSON response and merge with deterministic checks
    (host, cistron count, must_pass) to produce a full expectation_results dict.
    """
    # Strip markdown fences
    stripped = raw.strip()
    if stripped.startswith("```"):
        stripped = re.sub(r"^```[a-zA-Z]*\n?", "", stripped)
        stripped = re.sub(r"\n?```$", "", stripped.strip())

    try:
        judge_results = json.loads(stripped)
    except json.JSONDecodeError as e:
        logger.error(f"JSON parse error for {outcome.prompt_id}: {e}")
        judge_results = {}

    results = {}
    expect = case.expect

    # --- Deterministic checks (same as check_expectations) ---
    if "host" in expect:
        actual_host = ""
        if outcome.generated_spec:
            root = outcome.generated_spec.get("construct", outcome.generated_spec)
            actual_host = root.get("host", "") if isinstance(root, dict) else ""
        passed = _hosts_match(expect["host"], actual_host or "")
        results["host"] = {"expected": expect["host"], "actual": actual_host, "passed": passed}

    if "min_cistrons" in expect:
        actual = _infer_cistron_count(outcome)
        expected = expect["min_cistrons"]
        results["min_cistrons"] = {
            "expected": f">= {expected}", "actual": actual, "passed": actual >= expected,
        }
    if "max_cistrons" in expect:
        actual = _infer_cistron_count(outcome)
        expected = expect["max_cistrons"]
        results["max_cistrons"] = {
            "expected": f"<= {expected}", "actual": actual, "passed": actual <= expected,
        }
    if expect.get("must_pass"):
        results["must_pass"] = {
            "expected": True, "actual": bool(outcome.harness_passed), "passed": bool(outcome.harness_passed),
        }

    # --- LLM-judged part checks ---
    for part in expect.get("must_have_parts", []):
        res = judge_results.get(part, {})
        passed = res.get("passed", False)
        results[f"has_{part}"] = {
            "expected": True, "actual": passed, "passed": passed,
            "reason": res.get("reason", ""),
        }

    for part in expect.get("must_not_have_parts", []):
        res = judge_results.get(part, {})
        passed = res.get("passed", False)
        results[f"not_have_{part}"] = {
            "expected": False, "actual": not passed, "passed": passed,
            "reason": res.get("reason", ""),
        }

    for extra in ["logic_verification", "host_compatibility"]:
        if extra in judge_results:
            res = judge_results[extra]
            passed = res.get("passed", False)
            results[extra] = {
                "expected": True, "actual": passed, "passed": passed,
                "reason": res.get("reason", ""),
            }

    return results


def _apply_judge_results(outcome: EvalOutcome, exp_results: dict) -> EvalOutcome:
    """Update an outcome's expectation fields from a parsed judge result dict."""
    outcome.expectation_results = exp_results
    outcome.expectations_total = len(exp_results)
    outcome.expectations_ok = sum(1 for r in exp_results.values() if r["passed"])
    outcome.expectations_pct = (
        round(outcome.expectations_ok / outcome.expectations_total * 100, 1)
        if outcome.expectations_total else None
    )
    outcome.all_expectations_met = (
        all(r["passed"] for r in exp_results.values()) if exp_results else False
    )
    return outcome


def _print_outcome(outcome: EvalOutcome, lock: threading.Lock) -> None:
    with lock:
        pct_str = f"{outcome.expectations_pct}%" if outcome.expectations_pct is not None else "n/a"
        status = "PASS" if outcome.harness_passed else "FAIL"
        print(f"  [{status}] {outcome.prompt_id}  exp {outcome.expectations_ok}/{outcome.expectations_total} ({pct_str})")
        for k, v in outcome.expectation_results.items():
            icon = "[OK]  " if v["passed"] else "[MISS]"
            reason = f" — {v['reason']}" if v.get("reason") else ""
            print(f"    {icon} {k}{reason}")


# ---------------------------------------------------------------------------
# Streaming mode (concurrent, works with Claude and Gemini)
# ---------------------------------------------------------------------------

def _rejudge_outcome_streaming(
    outcome: EvalOutcome,
    case: PromptCase,
    judge_model: str,
    print_lock: threading.Lock,
) -> EvalOutcome:
    """Re-judge a single outcome via the regular (streaming) API."""
    if not case.expect or outcome.generated_spec is None:
        return outcome

    exp_results = check_expectations_with_llm(
        case=case, spec=outcome.generated_spec, eval_result=None, model=judge_model,
    )
    _apply_judge_results(outcome, exp_results)
    _print_outcome(outcome, print_lock)
    return outcome


def rejudge_file_streaming(
    results_path: Path,
    outcomes: list[EvalOutcome],
    cases: list[PromptCase],
    judge_model: str,
    concurrency: int = 5,
) -> list[EvalOutcome]:
    """Re-judge outcomes using the regular API with a thread pool."""
    print_lock = threading.Lock()
    updated = [None] * len(outcomes)

    with ThreadPoolExecutor(max_workers=concurrency) as pool:
        futures = {
            pool.submit(_rejudge_outcome_streaming, outcome, case, judge_model, print_lock): i
            for i, (outcome, case) in enumerate(zip(outcomes, cases))
        }
        for future in as_completed(futures):
            i = futures[future]
            try:
                updated[i] = future.result()
            except Exception as e:
                logger.error(f"Judge failed for {outcomes[i].prompt_id}: {e}")
                updated[i] = outcomes[i]

    return [o for o in updated if o is not None]


# ---------------------------------------------------------------------------
# Batch mode (Anthropic Batch API only — Claude models)
# ---------------------------------------------------------------------------

def _make_batch_request(
    custom_id: str,
    case: PromptCase,
    outcome: EvalOutcome,
    model: str,
) -> dict | None:
    """Build a single Anthropic batch request object, or None if nothing to judge."""
    prompt = _build_judge_prompt(case, outcome)
    if prompt is None:
        return None

    return {
        "custom_id": custom_id,
        "params": {
            "model": model,
            "max_tokens": 2048,
            "system": JUDGE_SYSTEM_PROMPT,
            "messages": [{"role": "user", "content": prompt}],
        },
    }


def rejudge_file_batch(
    results_path: Path,
    outcomes: list[EvalOutcome],
    cases: list[PromptCase],
    judge_model: str,
    poll_interval: int = 30,
) -> list[EvalOutcome]:
    """
    Re-judge outcomes using the Anthropic Message Batches API.

    All requests are submitted in a single batch. The function polls until
    the batch completes, then applies results. No concurrency setting needed —
    Anthropic processes requests in parallel on their end.
    """
    if "claude" not in judge_model.lower():
        raise ValueError(f"Batch mode is only supported for Claude models, got: {judge_model}")

    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        raise RuntimeError("ANTHROPIC_API_KEY environment variable not set")

    import anthropic
    client = anthropic.Anthropic(api_key=api_key)

    # Build all batch requests, keeping a mapping from custom_id → outcome index
    requests = []
    id_to_index: dict[str, int] = {}

    for i, (outcome, case) in enumerate(zip(outcomes, cases)):
        if not case.expect or outcome.generated_spec is None:
            continue
        custom_id = f"{outcome.prompt_id[:50]}_{i}"  # max 64 chars
        req = _make_batch_request(custom_id, case, outcome, judge_model)
        if req is not None:
            requests.append(req)
            id_to_index[custom_id] = i

    if not requests:
        print("  No judgeable outcomes found.")
        return outcomes

    print(f"  Submitting {len(requests)} requests to Anthropic Batch API...")
    batch = client.messages.batches.create(requests=requests)
    print(f"  Batch ID: {batch.id}  |  Status: {batch.processing_status}")
    print(f"  Polling every {poll_interval}s until complete...")

    # Poll until done
    while batch.processing_status != "ended":
        time.sleep(poll_interval)
        batch = client.messages.batches.retrieve(batch.id)
        counts = batch.request_counts
        print(
            f"  [{datetime.now().strftime('%H:%M:%S')}] {batch.processing_status} — "
            f"processing: {counts.processing}  succeeded: {counts.succeeded}  "
            f"errored: {counts.errored}"
        )

    print(f"  Batch complete. Processing results...")

    # Apply results back to outcomes
    for result in client.messages.batches.results(batch.id):
        i = id_to_index.get(result.custom_id)
        if i is None:
            logger.warning(f"Unknown custom_id in batch results: {result.custom_id}")
            continue

        outcome = outcomes[i]
        case = cases[i]

        if result.result.type == "succeeded":
            raw = result.result.message.content[0].text
            exp_results = _parse_judge_response(raw, case, outcome)
            _apply_judge_results(outcome, exp_results)
            pct_str = f"{outcome.expectations_pct}%" if outcome.expectations_pct is not None else "n/a"
            status = "PASS" if outcome.harness_passed else "FAIL"
            print(f"  [{status}] {outcome.prompt_id}  exp {outcome.expectations_ok}/{outcome.expectations_total} ({pct_str})")
        elif result.result.type == "errored":
            err = result.result.error
            logger.error(f"Batch error for {outcome.prompt_id}: {err.type} — {getattr(err, 'message', '')}")
        elif result.result.type == "expired":
            logger.warning(f"Request expired (24h timeout) for {outcome.prompt_id}")
        else:
            logger.warning(f"Unexpected result type for {outcome.prompt_id}: {result.result.type}")

    return outcomes


# ---------------------------------------------------------------------------
# Top-level orchestrator
# ---------------------------------------------------------------------------

def rejudge_file(
    results_path: Path,
    judge_model: str,
    concurrency: int = 5,
    batch: bool = False,
    poll_interval: int = 30,
    limit: int | None = None,
    dry_run: bool = False,
    run_name: str | None = None,
) -> Path:
    """Load a results JSON, re-judge all outcomes, write new results file."""

    logger.info(f"Loading {results_path}")
    outcomes, cases = load_results_json(results_path)

    if limit:
        outcomes = outcomes[:limit]
        cases = cases[:limit]

    judgeable = sum(
        1 for o, c in zip(outcomes, cases)
        if o.generated_spec is not None and c.expect
    )
    mode_str = "Anthropic Batch API" if batch else f"streaming (concurrency={concurrency})"
    print(f"\n{'='*70}")
    print(f"RE-JUDGING: {results_path.name}")
    print(f"  Outcomes: {len(outcomes)}  |  Judgeable: {judgeable}  |  Judge: {judge_model}")
    print(f"  Mode: {mode_str}")
    print(f"{'='*70}")

    if dry_run:
        print(f"[DRY RUN] Would submit {judgeable} requests — skipping API calls.")
        return results_path

    if batch:
        if "claude" not in judge_model.lower():
            print("Error: --batch is only supported for Claude models.", file=sys.stderr)
            sys.exit(1)
        outcomes = rejudge_file_batch(results_path, outcomes, cases, judge_model, poll_interval)
    else:
        outcomes = rejudge_file_streaming(results_path, outcomes, cases, judge_model, concurrency)

    # Write updated results
    if run_name is None:
        stem = results_path.stem
        mode_tag = "batch" if batch else "streaming"
        run_name = (
            f"{stem}_rejudged_{mode_tag}_{judge_model.replace('/', '_')}"
            f"_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        )

    out_path = write_results(outcomes, run_name=run_name)

    summary = _compute_summary(outcomes)
    print(f"\n{'='*70}")
    print(f"SUMMARY — {run_name}")
    print(f"{'='*70}")
    print(f"Harness passed:       {summary['harness_passed']}/{summary['total']} ({summary['pass_rate']}%)")
    print(f"Expectations met:     {summary['all_expectations_met']}/{summary['total']} "
          f"({summary['expectation_rate']}%) [all-or-nothing]")
    print(f"Individual exp. met:  {summary['expectations_ok']}/{summary['expectations_total']} "
          f"({summary['individual_expectation_rate']}%) [fractional]")
    print(f"\nResults written to: {out_path}")

    return out_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Re-judge expectations in existing eval results using an LLM.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "results",
        nargs="+",
        help="Path(s) to results JSON file(s) to re-judge",
    )
    parser.add_argument(
        "--judge-model",
        default="claude-sonnet-4-20250514",
        help="LLM to use as judge. Claude: ANTHROPIC_API_KEY. Gemini: GEMINI_API_KEY. (default: claude-sonnet-4-20250514)",
    )
    parser.add_argument(
        "--batch",
        action="store_true",
        help="Use the Anthropic Message Batches API (Claude only). No concurrency setting needed — "
             "all requests submitted at once, processed async. 50%% cheaper than streaming.",
    )
    parser.add_argument(
        "--poll-interval",
        type=int,
        default=30,
        help="Seconds between batch status polls (default: 30, only used with --batch)",
    )
    parser.add_argument(
        "--concurrency", "-j",
        type=int,
        default=5,
        help="Parallel judge calls for streaming mode (default: 5, ignored with --batch)",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Only re-judge the first N outcomes per file (useful for testing)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would be judged without making any API calls",
    )
    parser.add_argument(
        "--run-name",
        default=None,
        help="Override the output run name (only applies when judging a single file)",
    )
    args = parser.parse_args()

    paths = [Path(p) for p in args.results]
    missing = [p for p in paths if not p.exists()]
    if missing:
        print(f"Error: file(s) not found: {', '.join(str(p) for p in missing)}", file=sys.stderr)
        sys.exit(1)

    run_name = args.run_name if len(paths) == 1 else None

    for path in paths:
        rejudge_file(
            results_path=path,
            judge_model=args.judge_model,
            concurrency=args.concurrency,
            batch=args.batch,
            poll_interval=args.poll_interval,
            limit=args.limit,
            dry_run=args.dry_run,
            run_name=run_name,
        )


if __name__ == "__main__":
    main()

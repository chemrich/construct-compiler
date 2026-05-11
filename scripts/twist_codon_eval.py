#!/usr/bin/env python3
"""
Run the validation harness on a sample of specs, using the Twist API for
codon optimization in place of the local DNA Chisel pass.

For each spec:
  1. parse → resolve_parts → reverse_translate (gives a graph with DNA in
     every coding part)
  2. for each coding part with a known protein sequence, call
     TwistVendor.optimize_codons() and replace the part's sequence with
     the Twist-optimized DNA (cached by (protein, organism))
  3. run the harness validators on the assembled graph
  4. append a JSONL record to the output file

A global rate limiter throttles every HTTP call (5s default), so this is
safe to run against the live API for many specs at once. Per-protein
caching means common parts (His, MBP, eGFP) only cost one Twist call.

The output JSONL is appended incrementally so the run is resumable: on
restart, specs already present in the output file are skipped.

Usage:
    uv run python scripts/twist_codon_eval.py            # 100 specs, seed=0
    uv run python scripts/twist_codon_eval.py -n 5       # quick sanity run
    uv run python scripts/twist_codon_eval.py --interval 10  # 10s between calls
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import random
import sys
import threading
import time
from pathlib import Path
from typing import Optional

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

try:
    from dotenv import load_dotenv
    load_dotenv(ROOT / ".env")
except ImportError:
    pass

from Bio.Seq import Seq

from construct_compiler.frontend.parser import parse_spec
from construct_compiler.passes.part_resolution import resolve_parts
from construct_compiler.passes.reverse_translation import reverse_translate
from construct_compiler.core.parts import (
    CDS, PurificationTag, SolubilityTag, CleavageSite, Linker,
)
from construct_compiler.core.types import ResolutionState
from construct_compiler.validation.construct_checks import (
    run_all_checks, CheckSeverity,
)
from construct_compiler.vendors.twist import TwistVendor

CODING_TYPES = (CDS, PurificationTag, SolubilityTag, CleavageSite, Linker)


# ---------------------------------------------------------------------------
# Rate-limited Twist wrapper
# ---------------------------------------------------------------------------

class RateLimitedTwist(TwistVendor):
    """TwistVendor with a global min-interval gate on every HTTP request."""

    def __init__(self, *args, min_interval: float = 5.0, **kwargs):
        super().__init__(*args, **kwargs)
        self._min_interval = float(min_interval)
        self._last_call = 0.0
        self._lock = threading.Lock()
        self.calls = 0

    def _request(self, method, url, **kwargs):
        with self._lock:
            now = time.monotonic()
            wait = self._min_interval - (now - self._last_call)
            if wait > 0:
                time.sleep(wait)
            self._last_call = time.monotonic()
            self.calls += 1
        return super()._request(method, url, **kwargs)


# ---------------------------------------------------------------------------
# Organism mapping (compiler shorthand → Twist host name)
# ---------------------------------------------------------------------------

ORGANISM_MAP = {
    "e_coli": "Escherichia coli",
    "e_coli_bl21": "Escherichia coli",
    "ecoli": "Escherichia coli",
    "human": "Homo sapiens",
    "h_sapiens": "Homo sapiens",
    "yeast": "Saccharomyces cerevisiae",
    "s_cerevisiae": "Saccharomyces cerevisiae",
}


def twist_organism(host: str) -> str:
    return ORGANISM_MAP.get(host.lower(), "Escherichia coli")


# ---------------------------------------------------------------------------
# Per-spec evaluation
# ---------------------------------------------------------------------------

def _protein_for(part) -> Optional[str]:
    if isinstance(part, CDS):
        return part.protein_sequence
    return part.metadata.get("protein_sequence")


def evaluate_with_twist(
    spec_path: Path,
    vendor: RateLimitedTwist,
    cache: dict[tuple, str],
) -> dict:
    """Compile one spec, swap Twist-optimized DNA into coding parts, run checks."""
    rec: dict = {
        "spec_path": str(spec_path),
        "spec_name": spec_path.name,
        "passed": False,
        "compile_error": None,
        "twist_calls": 0,
        "twist_errors": [],
        "proteins_optimized": 0,
        "proteins_cached_hits": 0,
        "errors": [],
        "warnings": [],
        "error_count": 0,
        "warning_count": 0,
        "elapsed_s": 0.0,
    }
    t0 = time.monotonic()
    calls_before = vendor.calls

    try:
        graph = parse_spec(spec_path)
        rec["construct_name"] = graph.name
        graph = resolve_parts(graph)
        graph = reverse_translate(graph)
    except Exception as exc:
        rec["compile_error"] = f"pre-twist: {exc!r}"
        rec["elapsed_s"] = round(time.monotonic() - t0, 3)
        return rec

    organism = twist_organism(graph.host_organism)

    # Replace coding parts' DNA with Twist-optimized sequences.
    for part in graph.parts():
        if not isinstance(part, CODING_TYPES):
            continue
        if part.resolution != ResolutionState.CONCRETE:
            continue
        protein = _protein_for(part)
        if not protein:
            continue

        cache_key = (protein, organism)
        if cache_key in cache:
            optimized = cache[cache_key]
            rec["proteins_cached_hits"] += 1
        else:
            try:
                result = vendor.optimize_codons(protein, organism=organism)
            except Exception as exc:
                rec["twist_errors"].append(f"{part.name}: {exc!r}")
                continue
            optimized = result.optimized_sequence or ""
            if not optimized:
                rec["twist_errors"].append(
                    f"{part.name}: empty optimized seq (notes={result.notes})"
                )
                continue
            cache[cache_key] = optimized
            rec["proteins_optimized"] += 1

        # Preserve start/stop semantics: if the original DNA was prefixed
        # with ATG (cistron leader) or suffixed with a stop, keep them on
        # the Twist-optimized sequence too.
        original = str(part.sequence)
        if original.startswith("ATG") and not optimized.startswith("ATG"):
            optimized = "ATG" + optimized
        if isinstance(part, CDS) and part.has_stop:
            tail = original[-3:].upper()
            if tail in ("TAA", "TAG", "TGA"):
                opt_tail = optimized[-3:].upper()
                if opt_tail not in ("TAA", "TAG", "TGA"):
                    optimized = optimized + tail

        part.sequence = Seq(optimized)

    # Run validators on the assembled view (matches harness behavior)
    try:
        assembled = graph.assembled_graph()
        checks = run_all_checks(assembled)
    except Exception as exc:
        rec["compile_error"] = f"post-twist checks: {exc!r}"
        rec["twist_calls"] = vendor.calls - calls_before
        rec["elapsed_s"] = round(time.monotonic() - t0, 3)
        return rec

    errors = [c for c in checks if c.severity == CheckSeverity.ERROR]
    warnings = [c for c in checks if c.severity == CheckSeverity.WARNING]
    rec["error_count"] = len(errors)
    rec["warning_count"] = len(warnings)
    rec["passed"] = len(errors) == 0
    rec["errors"] = [
        {"check": c.check_name, "part_id": c.part_id, "message": c.message}
        for c in errors
    ]
    rec["warnings"] = [
        {"check": c.check_name, "part_id": c.part_id, "message": c.message}
        for c in warnings
    ]
    rec["insert_length_bp"] = graph.total_insert_length()
    rec["cistron_count"] = len(assembled.cistrons())
    rec["twist_calls"] = vendor.calls - calls_before
    rec["elapsed_s"] = round(time.monotonic() - t0, 3)
    return rec


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-n", "--num-specs", type=int, default=100)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--interval", type=float, default=5.0,
                   help="Min seconds between Twist HTTP calls (default 5.0)")
    p.add_argument("--specs-dir", type=Path,
                   default=ROOT / "evals" / "generated_specs")
    p.add_argument("--output", type=Path,
                   default=ROOT / "evals" / "results" / "twist_codon_eval.jsonl")
    p.add_argument("--user-email", default=os.environ.get("TWIST_USER_EMAIL"),
                   required=not os.environ.get("TWIST_USER_EMAIL"),
                   help="Twist account email (or set TWIST_USER_EMAIL)")
    p.add_argument("--sandbox", action="store_true")
    p.add_argument("--specs", nargs="+", type=Path, default=None,
                   help="Explicit spec paths (overrides random sampling).")
    args = p.parse_args()

    logging.basicConfig(
        level=logging.WARNING,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    if not args.specs_dir.exists():
        print(f"specs dir not found: {args.specs_dir}", file=sys.stderr)
        return 2

    if args.specs:
        sample = args.specs
    else:
        all_specs = sorted(args.specs_dir.glob("*.yaml"))
        rng = random.Random(args.seed)
        sample = rng.sample(all_specs, k=min(args.num_specs, len(all_specs)))

    args.output.parent.mkdir(parents=True, exist_ok=True)

    # Resume support: skip specs already in the output JSONL
    done: set[str] = set()
    if args.output.exists():
        with args.output.open() as f:
            for line in f:
                try:
                    rec = json.loads(line)
                    done.add(rec.get("spec_path", ""))
                except json.JSONDecodeError:
                    continue
    remaining = [s for s in sample if str(s) not in done]
    print(f"Sampled {len(sample)} specs; {len(done)} already in {args.output.name}; "
          f"{len(remaining)} to process.")

    vendor = RateLimitedTwist(
        user_email=args.user_email,
        sandbox=args.sandbox,
        min_interval=args.interval,
    )
    if not vendor.authenticated:
        print("Twist credentials missing — set TWIST_JWT_TOKEN, "
              "TWIST_END_USER_TOKEN, TWIST_USER_EMAIL.", file=sys.stderr)
        return 1
    print(f"Twist authenticated as {vendor.user_email} "
          f"({'sandbox' if args.sandbox else 'production'}); "
          f"rate limit: 1 req / {args.interval:.1f}s")

    cache: dict[tuple, str] = {}
    t_start = time.monotonic()

    with args.output.open("a") as f:
        for i, spec in enumerate(remaining, 1):
            print(f"[{i}/{len(remaining)}] {spec.name} ...", flush=True)
            rec = evaluate_with_twist(spec, vendor, cache)
            f.write(json.dumps(rec) + "\n")
            f.flush()
            status = "PASS" if rec["passed"] else "FAIL"
            print(f"    {status}  errors={rec['error_count']} "
                  f"warns={rec['warning_count']} "
                  f"twist_calls={rec['twist_calls']} "
                  f"opt={rec['proteins_optimized']} "
                  f"hit={rec['proteins_cached_hits']} "
                  f"elapsed={rec['elapsed_s']}s")
            if rec["compile_error"]:
                print(f"    compile_error: {rec['compile_error']}")
            if rec["twist_errors"]:
                for te in rec["twist_errors"]:
                    print(f"    twist_error: {te}")

    elapsed = time.monotonic() - t_start
    print(f"\nDone. Total wall time: {elapsed/60:.1f} min, "
          f"total Twist HTTP calls: {vendor.calls}, "
          f"unique proteins cached: {len(cache)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

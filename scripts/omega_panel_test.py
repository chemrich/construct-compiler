#!/usr/bin/env python3
"""
omega_panel_test.py — Run the omegamega integration against a panel of specs.

Compiles each spec through reverse_translate (skipping codon optimisation for
speed), runs the omega oligopool pipeline, and prints a summary table.
A second pass combines all passing constructs into a single batch order to
illustrate Twist's non-linear oligo-pool pricing tiers.

Usage:
    OMEGAMEGA_DIR=/path/to/omegamega uv run python scripts/omega_panel_test.py [--n 100]

    --n N   Number of specs to include (default: 100).  Specs are drawn from
            all generated_specs/**/*.yaml that contain ≥2 cistron blocks.
"""

from __future__ import annotations

import argparse
import os
import random
import sys
import time
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from construct_compiler.backends.omega import OmegaBatchResult, OmegaResult, run_omega, run_omega_batch
from construct_compiler.core.graph import ConstructGraph
from construct_compiler.frontend.parser import parse_spec
from construct_compiler.passes.part_resolution import resolve_parts
from construct_compiler.passes.reverse_translation import reverse_translate

OMEGAMEGA_DIR = Path(os.environ.get("OMEGAMEGA_DIR", ""))
SPECS_DIR = Path(__file__).resolve().parent.parent / "evals" / "generated_specs"


def _polycistronic_specs(n: int, seed: int = 42) -> list[Path]:
    """Return up to *n* specs that contain ≥2 cistron blocks, deterministically shuffled."""
    candidates = [
        p for p in sorted(SPECS_DIR.glob("*.yaml"))
        if p.read_text().count("cistron:") >= 2
    ]
    rng = random.Random(seed)
    rng.shuffle(candidates)
    return candidates[:n]


@dataclass
class PanelRow:
    name: str
    status: str          # ok | compile_fail | omega_fail
    insert_bp: int = 0
    oligos: int = 0
    pools: int = 0
    min_fidelity: float = 0.0
    oligo_cost: float = 0.0
    elapsed_s: float = 0.0
    error: str = ""
    graph: object = None  # ConstructGraph if compile succeeded


def run_one(spec: Path, output_root: Path) -> PanelRow:
    name = spec.stem
    t0 = time.time()
    graph: ConstructGraph | None = None
    try:
        graph = parse_spec(spec)
        graph = resolve_parts(graph)
        graph = reverse_translate(graph)
    except Exception as e:
        return PanelRow(name=name, status="compile_fail",
                        elapsed_s=time.time() - t0,
                        error=str(e)[:300])

    insert_seq = graph.full_insert_sequence()
    insert_bp = len(insert_seq) if insert_seq else 0

    try:
        result: OmegaResult = run_omega(
            graph,
            output_dir=output_root / name,
            nopt_steps=100,
            nopt_runs=1,
            njobs=1,
            omegamega_dir=OMEGAMEGA_DIR,
        )
        return PanelRow(
            name=name,
            status="ok",
            insert_bp=insert_bp,
            oligos=result.oligo_count,
            pools=result.pool_count,
            min_fidelity=result.min_fidelity,
            oligo_cost=result.oligo_cost_usd,
            elapsed_s=time.time() - t0,
            graph=graph,
        )
    except Exception as e:
        return PanelRow(name=name, status="omega_fail",
                        insert_bp=insert_bp,
                        elapsed_s=time.time() - t0,
                        error=str(e)[:300])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--n", type=int, default=100, metavar="N",
                        help="Number of specs to run (default: 100)")
    args = parser.parse_args()

    if not OMEGAMEGA_DIR or not (OMEGAMEGA_DIR / "code" / "omega.py").exists():
        print("ERROR: set OMEGAMEGA_DIR to your omegamega repo root.", file=sys.stderr)
        sys.exit(1)

    panel = _polycistronic_specs(args.n)

    output_root = Path("output/omega_panel")
    output_root.mkdir(parents=True, exist_ok=True)

    print(f"Running omega panel: {len(panel)} polycistronic specs\n")

    rows: list[PanelRow] = []
    for i, spec in enumerate(panel, 1):
        print(f"  [{i:3d}/{len(panel)}] {spec.stem[:50]} ...", end=" ", flush=True)
        row = run_one(spec, output_root)
        rows.append(row)
        if row.status == "ok":
            print(f"ok  {row.insert_bp} bp  {row.oligos} oligos  fidelity={row.min_fidelity:.3f}  ${row.oligo_cost:.2f}  ({row.elapsed_s:.1f}s)")
        else:
            print(f"FAIL  [{row.status}]  {row.error[:80]}")

    # --- Summary table -------------------------------------------------------
    ok = [r for r in rows if r.status == "ok"]
    compile_fail = [r for r in rows if r.status == "compile_fail"]
    omega_fail = [r for r in rows if r.status == "omega_fail"]

    print(f"\n{'='*72}")
    print(f"  Panel results: {len(ok)}/{len(rows)} passed")
    print(f"  Compile failures: {len(compile_fail)}")
    print(f"  Omega failures:   {len(omega_fail)}")

    if ok:
        avg_fidelity = sum(r.min_fidelity for r in ok) / len(ok)
        avg_cost = sum(r.oligo_cost for r in ok) / len(ok)
        avg_bp = sum(r.insert_bp for r in ok) / len(ok)
        print(f"\n  Passed constructs:")
        print(f"    Avg insert:   {avg_bp:.0f} bp")
        print(f"    Avg fidelity: {avg_fidelity:.3f}")
        print(f"    Avg cost:     ${avg_cost:.2f}")

        print(f"\n  {'Construct':<45} {'bp':>6} {'oligos':>7} {'pools':>6} {'fidelity':>9} {'cost':>8}")
        print(f"  {'-'*45} {'-'*6} {'-'*7} {'-'*6} {'-'*9} {'-'*8}")
        for r in sorted(ok, key=lambda r: r.name):
            print(f"  {r.name:<45} {r.insert_bp:>6} {r.oligos:>7} {r.pools:>6} {r.min_fidelity:>9.3f} ${r.oligo_cost:>7.2f}")

    if compile_fail or omega_fail:
        print(f"\n  Failures:")
        for r in compile_fail + omega_fail:
            print(f"    [{r.status:>12}] {r.name}: {r.error[:70]}")

    print(f"{'='*72}")

    # --- Batch pricing comparison --------------------------------------------
    if len(ok) < 2:
        print("\nNeed ≥2 passing constructs for batch pricing comparison.")
        return

    print(f"\n{'='*72}")
    print("  BATCH PRICING COMPARISON")
    print(f"  Treating all {len(ok)} passing constructs as a single Twist order")
    print(f"{'='*72}\n")

    batch_graphs = [(r.graph, r.name) for r in ok]
    batch_output = output_root / "_batch"
    print(f"  Running combined omegamega design ({len(batch_graphs)} constructs) ...", flush=True)
    t_batch = time.time()
    try:
        batch: OmegaBatchResult = run_omega_batch(
            batch_graphs,
            output_dir=batch_output,
            nopt_steps=100,
            nopt_runs=1,
            njobs=1,
            omegamega_dir=OMEGAMEGA_DIR,
        )
        elapsed_batch = time.time() - t_batch
        individual_total = sum(r.oligo_cost for r in ok)
        savings = individual_total - batch.total_cost_usd
        pct = 100 * savings / individual_total if individual_total else 0

        print(f"  Done in {elapsed_batch:.1f}s\n")
        print(f"  {'':45} {'oligos':>7} {'cost':>10}")
        print(f"  {'-'*45} {'-'*7} {'-'*10}")
        for r in sorted(ok, key=lambda r: r.name):
            print(f"  {r.name:<45} {r.oligos:>7} ${r.oligo_cost:>9.2f}")
        print(f"  {'-'*45} {'-'*7} {'-'*10}")
        print(f"  {'Individual orders (sum)':<45} {sum(r.oligos for r in ok):>7} ${individual_total:>9.2f}")
        print(f"  {'Combined batch order':<45} {batch.total_oligos:>7} ${batch.total_cost_usd:>9.2f}")
        print(f"\n  Savings: ${savings:.2f}  ({pct:.0f}% reduction)")
        print(f"  Batch min fidelity: {batch.min_fidelity:.3f}  avg: {batch.avg_fidelity:.3f}")
    except Exception as e:
        print(f"  Batch run failed: {e}")

    print(f"{'='*72}")


if __name__ == "__main__":
    main()

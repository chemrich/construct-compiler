#!/usr/bin/env python3
"""
omega_panel_test.py — Run the omegamega integration against a panel of specs.

Compiles each spec through reverse_translate (skipping codon optimisation for
speed), runs the omega oligopool pipeline, and prints a summary table.

Usage:
    OMEGAMEGA_DIR=/path/to/omegamega uv run python scripts/omega_panel_test.py
"""

from __future__ import annotations

import sys
import time
import traceback
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from construct_compiler.backends.omega import OmegaResult, run_omega
from construct_compiler.frontend.parser import parse_spec
from construct_compiler.passes.part_resolution import resolve_parts
from construct_compiler.passes.reverse_translation import reverse_translate

OMEGAMEGA_DIR = Path(
    __import__("os").environ.get("OMEGAMEGA_DIR", "")
)

PANEL: list[Path] = sorted(
    Path(__file__).resolve().parent.parent.glob("evals/generated_specs/poly_*.yaml")
)[:20]


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


def run_one(spec: Path, output_root: Path) -> PanelRow:
    name = spec.stem
    t0 = time.time()
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
        )
    except Exception as e:
        return PanelRow(name=name, status="omega_fail",
                        insert_bp=insert_bp,
                        elapsed_s=time.time() - t0,
                        error=str(e)[:300])


def main() -> None:
    if not OMEGAMEGA_DIR or not (OMEGAMEGA_DIR / "code" / "omega.py").exists():
        print("ERROR: set OMEGAMEGA_DIR to your omegamega repo root.", file=sys.stderr)
        sys.exit(1)

    output_root = Path("output/omega_panel")
    output_root.mkdir(parents=True, exist_ok=True)

    print(f"Running omega panel: {len(PANEL)} polycistronic specs\n")

    rows: list[PanelRow] = []
    for i, spec in enumerate(PANEL, 1):
        print(f"  [{i:2d}/{len(PANEL)}] {spec.stem} ...", end=" ", flush=True)
        row = run_one(spec, output_root)
        rows.append(row)
        if row.status == "ok":
            print(f"✓  {row.insert_bp} bp  {row.oligos} oligos  fidelity={row.min_fidelity:.3f}  ${row.oligo_cost:.2f}  ({row.elapsed_s:.1f}s)")
        else:
            print(f"✗  [{row.status}]  {row.error[:80]}")

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


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
omega_design.py — compile a construct spec and design an oligopool with omegamega.

Runs the full construct-compiler pipeline (parse → resolve → reverse-translate →
constrain), then hands the codon-optimised insert to omegamega for Golden Gate
oligopool fragmentation.

Usage:
    uv run python scripts/omega_design.py examples/his_tev_mbp_egfp.yaml

    # Skip codon optimisation (faster, for quick iteration)
    uv run python scripts/omega_design.py spec.yaml --fast

    # Custom omegamega parameters
    uv run python scripts/omega_design.py spec.yaml \\
        --enzyme BsaI --njunctions 50 --oligo-len 350 \\
        --nopt-steps 1000 --nopt-runs 5

Environment:
    OMEGAMEGA_DIR   Path to a local clone of chemrich/omegamega (required).
                    git clone https://github.com/chemrich/omegamega
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Ensure src/ is on the path when run directly
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from construct_compiler.frontend.parser import parse_spec
from construct_compiler.passes.part_resolution import resolve_parts
from construct_compiler.passes.reverse_translation import reverse_translate
from construct_compiler.passes.constraint_resolution import resolve_constraints
from construct_compiler.backends.omega import run_omega


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compile a construct spec and design an oligopool with omegamega."
    )
    parser.add_argument("spec", type=Path, help="Path to construct YAML spec")
    parser.add_argument(
        "-o", "--output-dir", type=Path, default=None,
        help="Output directory for omegamega results (default: output/<spec_stem>_omega/)",
    )
    parser.add_argument(
        "--fast", action="store_true",
        help="Skip codon optimisation (pass 3); use raw reverse-translated sequence",
    )
    parser.add_argument("--enzyme", default="BsaI", help="Type IIS enzyme (default: BsaI)")
    parser.add_argument("--njunctions", type=int, default=50, help="GG sites per pool")
    parser.add_argument("--oligo-len", type=int, default=350, help="Oligo length in bp")
    parser.add_argument("--upstream-bbsite", default="AATG", help="Upstream backbone junction")
    parser.add_argument("--downstream-bbsite", default="TTAG", help="Downstream backbone junction")
    parser.add_argument("--nopt-steps", type=int, default=500, help="SA optimisation steps")
    parser.add_argument("--nopt-runs", type=int, default=3, help="SA optimisation runs")
    parser.add_argument("--njobs", type=int, default=1, help="Parallel jobs")
    args = parser.parse_args()

    output_dir = args.output_dir or (
        Path("output") / f"{args.spec.stem}_omega"
    )

    # -- Compile construct ----------------------------------------------------
    print(f"Compiling {args.spec} ...", flush=True)
    graph = parse_spec(args.spec)
    graph = resolve_parts(graph)
    graph = reverse_translate(graph)

    if not args.fast:
        print("Running constraint resolution (codon optimisation) ...", flush=True)
        graph = resolve_constraints(graph)

    insert_seq = graph.full_insert_sequence()
    if insert_seq is None:
        print("ERROR: no concrete sequences after compilation.", file=sys.stderr)
        sys.exit(1)

    print(f"Insert length: {len(insert_seq)} bp", flush=True)

    # -- Run omegamega --------------------------------------------------------
    print(f"Running omegamega oligopool design → {output_dir} ...", flush=True)
    result = run_omega(
        graph,
        output_dir=output_dir,
        enzyme=args.enzyme,
        njunctions=args.njunctions,
        oligo_len=args.oligo_len,
        upstream_bbsite=args.upstream_bbsite,
        downstream_bbsite=args.downstream_bbsite,
        nopt_steps=args.nopt_steps,
        nopt_runs=args.nopt_runs,
        njobs=args.njobs,
    )

    print()
    print(result.summary())
    print()
    print(f"Oligo order: {result.output_dir / 'oligo_order.csv'}")
    print(f"Pool stats:  {result.output_dir / 'pool_stats.csv'}")
    if (result.output_dir / "cost_summary.csv").exists():
        print(f"Cost:        {result.output_dir / 'cost_summary.csv'}")


if __name__ == "__main__":
    main()

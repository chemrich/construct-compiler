#!/usr/bin/env python3
"""
Design-evaluate runner for automated construct optimization.

Takes a base YAML spec, optionally generates variants across design axes,
evaluates all of them through the validation harness, and prints a ranked
summary. Designed to be called by an agent in a design loop.

Usage:
    # Evaluate a single spec
    python scripts/design_evaluate.py examples/his_tev_mbp_egfp.yaml

    # Evaluate with variant axes
    python scripts/design_evaluate.py examples/his_tev_mbp_egfp.yaml \
        --axes "expression=high,medium,low" "spacer=20,30,50,100"

    # JSON output for machine consumption
    python scripts/design_evaluate.py examples/his_tev_mbp_egfp.yaml --json

    # Fast mode (skip codon optimization)
    python scripts/design_evaluate.py examples/his_tev_mbp_egfp.yaml --fast

    # With intermediate stage diagnostics
    python scripts/design_evaluate.py examples/his_tev_mbp_egfp.yaml --intermediate

Exit codes:
    0 — all specs passed
    1 — one or more specs failed validation
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

# Add the src directory to the path so this script works without installation
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from construct_compiler.validation.harness import (
    evaluate_spec,
    evaluate_batch,
    batch_summary,
)
from construct_compiler.validation.variants import (
    DesignAxis,
    vary_spec,
    vary_spec_dicts,
)


# ---------------------------------------------------------------------------
# Axis field path mapping — common axes with human-friendly names
# ---------------------------------------------------------------------------

# Maps short axis names to their YAML field paths.
# These assume the standard cassette layout from the example spec:
#   0: promoter, 1: first cistron, 2: spacer, 3: second cistron, 4: terminator
#
# Agents can use these short names or provide full dot-paths.

KNOWN_AXES = {
    "expression":    "cassette.1.cistron.expression",
    "expression1":   "cassette.1.cistron.expression",
    "expression2":   "cassette.3.cistron.expression",
    "spacer":        "cassette.2.spacer",
    "terminator":    "cassette.4.terminator",
    "promoter":      "cassette.0.promoter",
}


def parse_axis(axis_str: str) -> DesignAxis:
    """
    Parse an axis string like "expression=high,medium,low" into a DesignAxis.

    The name part can be a short name from KNOWN_AXES or a full dot-path.
    Values are comma-separated. Numeric values are auto-detected.
    """
    name, _, values_str = axis_str.partition("=")
    if not values_str:
        raise ValueError(f"Axis must be 'name=val1,val2,...', got: {axis_str!r}")

    # Resolve field path
    field_path = KNOWN_AXES.get(name, name)

    # Parse values, auto-detecting ints
    values = []
    for v in values_str.split(","):
        v = v.strip()
        try:
            values.append(int(v))
        except ValueError:
            try:
                values.append(float(v))
            except ValueError:
                values.append(v)

    return DesignAxis(name=name, field_path=field_path, values=values)


def main():
    parser = argparse.ArgumentParser(
        description="Evaluate construct specs through the validation harness.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("spec", type=Path, help="Base YAML spec file.")
    parser.add_argument("--axes", nargs="*", default=[],
                        help='Design axes: "expression=high,medium,low" "spacer=20,50"')
    parser.add_argument("--json", dest="as_json", action="store_true",
                        help="Output results as JSON.")
    parser.add_argument("--fast", action="store_true",
                        help="Skip constraint resolution (faster, uses unoptimized codons).")
    parser.add_argument("--intermediate", action="store_true",
                        help="Also validate intermediate pipeline stages.")
    args = parser.parse_args()

    # Suppress noisy logs
    logging.basicConfig(level=logging.WARNING, format="%(message)s")
    import os
    os.environ["DNACHISEL_QUIET"] = "1"

    if not args.spec.exists():
        print(f"Error: spec file not found: {args.spec}", file=sys.stderr)
        sys.exit(2)

    # Parse axes
    axes = [parse_axis(a) for a in args.axes]

    if not axes:
        # Single spec evaluation
        result = evaluate_spec(
            str(args.spec),
            skip_constraints=args.fast,
            validate_intermediate=args.intermediate,
        )
        if args.as_json:
            print(result.to_json())
        else:
            print(result.summary())
        sys.exit(0 if result.passed else 1)

    else:
        # Variant generation + batch evaluation
        variants = vary_spec(str(args.spec), axes)

        if not args.as_json:
            print(f"Generated {len(variants)} variants from {len(axes)} axes:")
            for a in axes:
                print(f"  {a.name}: {a.values}")
            print()

        specs = [v.spec for v in variants]
        results = evaluate_batch(
            specs,
            skip_constraints=args.fast,
            validate_intermediate=args.intermediate,
        )

        if args.as_json:
            output = []
            for v, r in zip(variants, sorted(
                zip(variants, results), key=lambda x: -x[1].score
            )):
                output.append({
                    "parameters": v.parameters,
                    **r.to_dict(),
                })
            # Re-sort since zip pairs may not align after batch sort
            output = sorted(
                [{"parameters": v.parameters, **r.to_dict()}
                 for v, r in zip(variants, results)],
                key=lambda x: -x["score"],
            )
            print(json.dumps(output, indent=2, default=str))
        else:
            print(batch_summary(results))
            print()

            # Show top 3 with parameters
            print("Top designs:")
            for i, r in enumerate(results[:3]):
                # Find matching variant
                matching = [v for v in variants if v.variant_name == r.construct_name]
                params = matching[0].parameters if matching else {}
                print(f"  {i+1}. score={r.score:.2f}  {params}")

        all_passed = all(r.passed for r in results)
        sys.exit(0 if all_passed else 1)


if __name__ == "__main__":
    main()

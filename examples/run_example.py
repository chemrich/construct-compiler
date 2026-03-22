#!/usr/bin/env python3
"""
End-to-end example: compile a His-TEV-MBP-mEGFP construct
with a polycistronic mScarlet reporter.
"""

import sys
import os
import logging

# Add src/ directory to path so construct_compiler package is importable
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), "src"))

from construct_compiler import compile_construct, export_genbank
from construct_compiler.vendors.twist import TwistVendor

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)


def main():
    # Load the YAML spec
    spec_path = os.path.join(os.path.dirname(__file__), "his_tev_mbp_egfp.yaml")

    print("=" * 70)
    print("GENETIC CONSTRUCT COMPILER - Example Build")
    print("=" * 70)
    print()

    # Compile
    graph, assembly_plan = compile_construct(
        spec_path,
        skip_constraints=False,
        verbose=True,
    )

    # Export GenBank
    output_dir = os.path.dirname(__file__)
    gb_path = os.path.join(output_dir, f"{graph.name}.gb")
    record = export_genbank(graph, gb_path)

    print()
    print("=" * 70)
    print("GENBANK OUTPUT")
    print("=" * 70)
    print(f"File: {gb_path}")
    print(f"Total sequence length: {len(record.seq)} bp")
    print(f"Features: {len(record.features)}")
    for feat in record.features:
        print(f"  {feat.type:15s} {feat.location}  {feat.qualifiers.get('label', [''])[0]}")

    # Run vendor screening (mock mode — no API key)
    print()
    print("=" * 70)
    print("VENDOR SCREENING (Mock)")
    print("=" * 70)
    twist = TwistVendor()
    insert_seq = str(graph.full_insert_sequence() or "")
    if insert_seq:
        result = twist.screen(insert_seq)
        print(f"Vendor: {result.vendor}")
        print(f"Feasible: {result.feasible}")
        print(f"Length: {result.sequence_length} bp")
        print(f"Estimated price: ${result.estimated_price:.2f}")
        print(f"Turnaround: {result.turnaround_days[0]}-{result.turnaround_days[1]} business days")
        if result.warnings:
            for w in result.warnings:
                print(f"  WARNING: {w}")
        if result.errors:
            for e in result.errors:
                print(f"  ERROR: {e}")

    # Assembly plan summary
    print()
    print("=" * 70)
    print("ASSEMBLY PLAN + COST MODEL")
    print("=" * 70)
    print(assembly_plan.summary())

    if assembly_plan.recommended:
        rec = assembly_plan.recommended
        print(f">>> RECOMMENDED: {rec.name}")
        print(f">>> Total cost: ${rec.total_cost:.2f}")
        print(f">>> Wall clock: {rec.total_wall_clock_days[0]}-{rec.total_wall_clock_days[1]} business days")

    return graph, assembly_plan


if __name__ == "__main__":
    main()

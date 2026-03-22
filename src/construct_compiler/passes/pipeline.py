"""
Compilation pipeline — chains passes together.
"""

from __future__ import annotations

import logging
from typing import Optional

from ..core.graph import ConstructGraph
from ..frontend.parser import parse_spec
from .part_resolution import resolve_parts
from .reverse_translation import reverse_translate
from .constraint_resolution import resolve_constraints
from .assembly_planning import plan_assembly, AssemblyPlan, CostParams

logger = logging.getLogger(__name__)


def compile_construct(
    spec: dict | str,
    cost_params: Optional[CostParams] = None,
    skip_constraints: bool = False,
    verbose: bool = True,
) -> tuple[ConstructGraph, AssemblyPlan]:
    """
    Full compilation pipeline: spec → IR → passes → assembly plan.

    Returns (graph, assembly_plan).
    """
    if verbose:
        logging.basicConfig(level=logging.INFO, format="%(message)s")

    # Frontend: parse spec into IR
    logger.info("=== PARSING SPEC ===")
    graph = parse_spec(spec)
    logger.info(f"Parsed: {graph}")
    logger.info(f"Parts: {len(graph.parts())}")

    # Validate initial structure
    errors = graph.validate()
    for e in errors:
        if e.severity == "error":
            logger.error(f"Validation error: {e}")
        else:
            logger.warning(f"Validation warning: {e}")

    # Pass 1: Resolve parts (fetch sequences from DBs)
    logger.info("\n=== PASS 1: PART RESOLUTION ===")
    graph = resolve_parts(graph)

    # Pass 2: Reverse translate protein → DNA
    logger.info("\n=== PASS 2: REVERSE TRANSLATION ===")
    graph = reverse_translate(graph)

    # Pass 3: Constraint resolution (codon optimization, restriction site removal)
    if not skip_constraints:
        logger.info("\n=== PASS 3: CONSTRAINT RESOLUTION ===")
        graph = resolve_constraints(graph)

    # Pass 4: Assembly planning + cost model
    logger.info("\n=== PASS 4: ASSEMBLY PLANNING ===")
    assembly_plan = plan_assembly(graph, cost_params)

    # Final validation
    logger.info("\n=== FINAL VALIDATION ===")
    errors = graph.validate()
    for e in errors:
        logger.info(f"  {e}")

    logger.info(f"\n{graph.summary()}")
    logger.info(f"\n{assembly_plan.summary()}")

    return graph, assembly_plan

"""
MCP (Model Context Protocol) server for the construct compiler.

Exposes the compiler, validation harness, and variant evaluator as tools
that Claude Code (or any MCP client) can call over stdio.

Run directly:
    python -m construct_compiler.mcp_server

Or via entry point (after pip install):
    construct-compiler-mcp

Configuration for Claude Code (.mcp.json):
    {
        "mcpServers": {
            "construct-compiler": {
                "command": "construct-compiler-mcp"
            }
        }
    }
"""

from __future__ import annotations

import asyncio
import json
import logging
import os
import sys
from dataclasses import asdict
from pathlib import Path
from typing import Any

# Suppress noisy loggers before any imports that trigger them
os.environ["TQDM_DISABLE"] = "1"
logging.basicConfig(level=logging.WARNING)
logging.getLogger("dnachisel").setLevel(logging.WARNING)

from mcp.server import Server
from mcp.server.stdio import stdio_server
from mcp.types import TextContent, Tool

# ---------------------------------------------------------------------------
# MCP Server
# ---------------------------------------------------------------------------

server = Server("construct-compiler")

MAX_VARIANTS = 50
AGENT_GENERATED_DIR = Path(__file__).parent.parent.parent / "examples" / "agent_generated"


# ---------------------------------------------------------------------------
# Tool definitions
# ---------------------------------------------------------------------------

COMPILE_SPEC_SCHEMA = {
    "type": "object",
    "properties": {
        "spec": {
            "type": "object",
            "description": (
                "A construct spec dict (the full YAML structure with a top-level "
                "'construct' key containing name, host, backbone, cassette, "
                "constraints, and synthesis sections)."
            ),
        },
        "skip_constraints": {
            "type": "boolean",
            "description": (
                "Skip the DNA Chisel codon optimization pass. Faster but "
                "less thorough. Default false."
            ),
            "default": False,
        },
        "include_genbank": {
            "type": "boolean",
            "description": (
                "Include the GenBank-formatted sequence in the response. "
                "Default false."
            ),
            "default": False,
        },
    },
    "required": ["spec"],
}

CHECK_SPEC_SCHEMA = {
    "type": "object",
    "properties": {
        "spec": {
            "type": "object",
            "description": (
                "A construct spec dict (same format as compile_spec)."
            ),
        },
        "skip_constraints": {
            "type": "boolean",
            "description": "Skip codon optimization pass. Default false.",
            "default": False,
        },
        "validate_intermediate": {
            "type": "boolean",
            "description": (
                "Also run checks after each pipeline stage (not just final). "
                "Useful for diagnosing where errors first appear. Default false."
            ),
            "default": False,
        },
        "auto_save": {
            "type": "boolean",
            "description": (
                "If the spec passes all checks, auto-save it to "
                "examples/agent_generated/ for the regression test corpus. "
                "Default true."
            ),
            "default": True,
        },
    },
    "required": ["spec"],
}

EVALUATE_VARIANTS_SCHEMA = {
    "type": "object",
    "properties": {
        "spec": {
            "type": "object",
            "description": "Base construct spec dict to generate variants from.",
        },
        "axes": {
            "type": "array",
            "description": (
                "Design axes to vary. Each axis defines a parameter name, "
                "a dot-separated field path into the cassette list, and a "
                "list of values to try. Common field paths:\n"
                "  cassette.1.cistron.expression → expression level\n"
                "  cassette.2.spacer → spacer length (bp)\n"
                "  cassette.0.promoter → promoter choice\n"
                "  cassette.4.terminator → terminator choice"
            ),
            "items": {
                "type": "object",
                "properties": {
                    "name": {
                        "type": "string",
                        "description": "Human-readable axis name (e.g. 'expression')",
                    },
                    "field_path": {
                        "type": "string",
                        "description": "Dot-separated path into the spec (e.g. 'cassette.1.cistron.expression')",
                    },
                    "values": {
                        "type": "array",
                        "description": "Values to substitute at this path",
                        "items": {},
                    },
                },
                "required": ["name", "field_path", "values"],
            },
        },
        "skip_constraints": {
            "type": "boolean",
            "description": "Skip codon optimization (recommended for batch runs). Default true.",
            "default": True,
        },
    },
    "required": ["spec", "axes"],
}


@server.list_tools()
async def list_tools() -> list[Tool]:
    return [
        Tool(
            name="compile_spec",
            description=(
                "Compile a genetic construct YAML spec through the full pipeline: "
                "parse → resolve parts → reverse translate → codon optimize → "
                "assembly planning. Returns construct summary, assembly strategies "
                "with cost breakdowns, and optionally GenBank output. "
                "Use this to turn a spec dict into a buildable design."
            ),
            inputSchema=COMPILE_SPEC_SCHEMA,
        ),
        Tool(
            name="check_spec",
            description=(
                "Run all 4 validity checks on a construct spec: "
                "(1) reading frame continuity — every coding part's DNA length is "
                "divisible by 3, (2) start codon placement — first CDS in each "
                "cistron starts with ATG, (3) translation fidelity — DNA translates "
                "back to expected protein, (4) internal stop codons — no premature "
                "stops. Returns pass/fail, a 0.0–1.0 score, and detailed error "
                "messages for any failures. Use this after compile_spec to verify "
                "biological correctness, and iterate on the spec until it passes."
            ),
            inputSchema=CHECK_SPEC_SCHEMA,
        ),
        Tool(
            name="evaluate_variants",
            description=(
                "Generate a combinatorial design space from a base spec and a set "
                "of axes, evaluate all variants, and return results ranked by score. "
                f"Maximum {MAX_VARIANTS} variants per call. Use this after a base "
                "design passes check_spec to explore optimal parameter combinations "
                "(expression levels, spacer lengths, promoter choices, etc.)."
            ),
            inputSchema=EVALUATE_VARIANTS_SCHEMA,
        ),
    ]


# ---------------------------------------------------------------------------
# Tool handlers
# ---------------------------------------------------------------------------

@server.call_tool()
async def call_tool(name: str, arguments: dict) -> list[TextContent]:
    try:
        if name == "compile_spec":
            return _handle_compile(arguments)
        elif name == "check_spec":
            return _handle_check(arguments)
        elif name == "evaluate_variants":
            return _handle_variants(arguments)
        else:
            return [TextContent(type="text", text=f"Unknown tool: {name}")]
    except Exception as e:
        return [TextContent(
            type="text",
            text=json.dumps({"error": str(e), "type": type(e).__name__}, indent=2),
        )]


def _handle_compile(args: dict) -> list[TextContent]:
    """Compile a spec and return construct summary + assembly strategies."""
    from .passes.pipeline import compile_construct
    from .passes.assembly_planning import CostParams
    from .backends.genbank import export_genbank
    from io import StringIO
    from Bio import SeqIO

    spec = args["spec"]
    skip_constraints = args.get("skip_constraints", False)
    include_genbank = args.get("include_genbank", False)

    graph, plan = compile_construct(
        spec,
        cost_params=CostParams(),
        skip_constraints=skip_constraints,
        verbose=False,
    )

    # Build response
    result: dict[str, Any] = {
        "success": True,
        "construct_name": graph.name,
        "host": graph.host_organism,
        "insert_length_bp": graph.total_insert_length(),
        "cistron_count": len(graph.cistrons()),
        "parts": [],
        "strategies": [],
    }

    # Parts list
    from .core.parts import Backbone
    for part in graph.parts():
        if isinstance(part, Backbone):
            result["backbone"] = {
                "name": part.name,
                "resistance": part.resistance,
                "origin": part.origin.value if part.origin else None,
            }
            continue
        result["parts"].append({
            "id": part.id,
            "name": part.name,
            "type": part.__class__.__name__,
            "length_bp": len(part.sequence) if part.sequence else None,
        })

    # Assembly strategies
    for strat in plan.strategies:
        s = {
            "name": strat.name,
            "recommended": strat.recommended,
            "num_fragments": strat.num_fragments,
            "total_cost": round(strat.total_cost, 2),
            "synthesis_cost": round(strat.synthesis_cost, 2),
            "reagent_cost": round(strat.reagent_cost_with_overhead, 2),
            "researcher_cost": round(strat.researcher_cost, 2),
            "researcher_hours": round(strat.researcher_time_hrs, 1),
            "wall_clock_days": list(strat.total_wall_clock_days),
            "overhangs": strat.overhangs,
            "notes": strat.notes,
        }
        result["strategies"].append(s)

    # Validation issues
    errors_list = []
    warnings_list = []
    for e in graph.validate():
        if e.severity == "error":
            errors_list.append(f"{e.location}: {e.message}")
        else:
            warnings_list.append(f"{e.location}: {e.message}")
    result["errors"] = errors_list
    result["warnings"] = warnings_list

    # GenBank output (inline)
    if include_genbank:
        try:
            records = export_genbank(graph)
            buf = StringIO()
            SeqIO.write(records, buf, "genbank")
            result["genbank"] = buf.getvalue()
        except Exception as e:
            result["genbank_error"] = str(e)

    return [TextContent(type="text", text=json.dumps(result, indent=2))]


def _handle_check(args: dict) -> list[TextContent]:
    """Run validity checks and return structured results."""
    from .validation.harness import evaluate_spec

    spec = args["spec"]
    skip_constraints = args.get("skip_constraints", False)
    validate_intermediate = args.get("validate_intermediate", False)
    auto_save = args.get("auto_save", True)

    eval_result = evaluate_spec(
        spec,
        skip_constraints=skip_constraints,
        validate_intermediate=validate_intermediate,
    )

    result = {
        "passed": eval_result.passed,
        "score": eval_result.score,
        "error_count": eval_result.error_count,
        "warning_count": eval_result.warning_count,
        "check_count": eval_result.check_count,
        "construct_name": eval_result.construct_name,
        "insert_length_bp": eval_result.insert_length_bp,
        "cistron_count": eval_result.cistron_count,
        "compile_time_s": round(eval_result.compile_time_s, 3),
        "errors": eval_result.errors,
        "warnings": eval_result.warnings,
        "summary": eval_result.summary(),
    }

    if validate_intermediate:
        result["stages"] = [asdict(s) for s in eval_result.stages]

    if eval_result.compile_error:
        result["compile_error"] = eval_result.compile_error

    # Auto-save passing specs
    if auto_save and eval_result.passed and not eval_result.compile_error:
        save_path = _auto_save_spec(spec, eval_result.construct_name)
        if save_path:
            result["saved_to"] = str(save_path)

    return [TextContent(type="text", text=json.dumps(result, indent=2))]


def _handle_variants(args: dict) -> list[TextContent]:
    """Generate and evaluate design variants."""
    from .validation.harness import evaluate_batch, batch_summary
    from .validation.variants import DesignAxis, vary_spec

    spec = args["spec"]
    axes_raw = args["axes"]
    skip_constraints = args.get("skip_constraints", True)

    # Parse axes
    axes = [
        DesignAxis(
            name=a["name"],
            field_path=a["field_path"],
            values=a["values"],
        )
        for a in axes_raw
    ]

    # Check variant count
    import functools
    import operator
    num_variants = functools.reduce(operator.mul, [len(a.values) for a in axes], 1)
    if num_variants > MAX_VARIANTS:
        return [TextContent(type="text", text=json.dumps({
            "error": (
                f"Too many variants: {num_variants} (max {MAX_VARIANTS}). "
                f"Reduce the number of values per axis or use fewer axes."
            ),
            "axes": [{"name": a.name, "count": len(a.values)} for a in axes],
            "total_combinations": num_variants,
        }, indent=2))]

    # Generate variants
    variants = vary_spec(spec, axes)
    spec_dicts = [v.spec for v in variants]

    # Evaluate
    results = evaluate_batch(spec_dicts, skip_constraints=skip_constraints)

    # Build response — match results back to variant parameters
    variant_params = {v.variant_name: v.parameters for v in variants}
    response_variants = []
    for r in results:
        entry = {
            "variant_name": r.construct_name,
            "parameters": variant_params.get(r.construct_name, {}),
            "passed": r.passed,
            "score": r.score,
            "error_count": r.error_count,
            "warning_count": r.warning_count,
        }
        if r.compile_error:
            entry["compile_error"] = r.compile_error
        if r.errors:
            entry["errors"] = r.errors
        response_variants.append(entry)

    result = {
        "total_variants": len(variants),
        "passed_count": sum(1 for r in results if r.passed),
        "failed_count": sum(1 for r in results if not r.passed),
        "variants": response_variants,
        "batch_summary": batch_summary(results),
    }

    return [TextContent(type="text", text=json.dumps(result, indent=2))]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _auto_save_spec(spec: dict, construct_name: str) -> Path | None:
    """Save a validated spec to the agent_generated examples directory."""
    import yaml

    try:
        AGENT_GENERATED_DIR.mkdir(parents=True, exist_ok=True)

        # Sanitize filename
        safe_name = "".join(
            c if c.isalnum() or c in "-_" else "_"
            for c in (construct_name or "unnamed")
        ).strip("_")[:80]

        out_path = AGENT_GENERATED_DIR / f"{safe_name}.yaml"

        # Don't overwrite existing files — append a counter
        if out_path.exists():
            counter = 1
            while out_path.exists():
                out_path = AGENT_GENERATED_DIR / f"{safe_name}_{counter}.yaml"
                counter += 1

        out_path.write_text(yaml.dump(spec, default_flow_style=False, sort_keys=False))
        return out_path
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

async def _run():
    async with stdio_server() as (read, write):
        await server.run(read, write, server.create_initialization_options())


def main():
    asyncio.run(_run())


if __name__ == "__main__":
    main()

"""
CLI interface for the genetic construct compiler.

Usage:
    construct-compiler compile spec.yaml -o output/
    construct-compiler compile spec.yaml --cost-only
    construct-compiler validate spec.yaml
    construct-compiler parts --list promoters
"""

from __future__ import annotations

import json
import logging
import sys
from pathlib import Path

import click

from .passes.pipeline import compile_construct
from .passes.assembly_planning import CostParams
from .backends.genbank import export_genbank


@click.group()
@click.version_option(version="0.1.0")
def cli():
    """Genetic construct design compiler.

    Compile declarative YAML specs into annotated DNA sequences,
    assembly plans, and cost estimates.
    """
    pass


@cli.command()
@click.argument("spec", type=click.Path(exists=True))
@click.option("-o", "--output", type=click.Path(), default=None,
              help="Output directory for generated files. Defaults to current directory.")
@click.option("--cost-only", is_flag=True, help="Only show cost comparison, don't generate files.")
@click.option("--format", "fmt", type=click.Choice(["genbank", "json", "summary"]),
              default="genbank", help="Output format.")
@click.option("--researcher-rate", type=float, default=150.0,
              help="Researcher hourly rate in $/hr (default: $150).")
@click.option("--overhead", type=float, default=1.5,
              help="Overhead multiplier on consumables (default: 1.5x).")
@click.option("--twist-gene-bp", type=float, default=None, help="Twist gene fragment $/bp.")
@click.option("--twist-clonal-bp", type=float, default=None, help="Twist clonal gene $/bp.")
@click.option("--idt-gblock-bp", type=float, default=None, help="IDT gBlock $/bp.")
@click.option("--bsai-cost", type=float, default=None, help="BsaI $/reaction.")
@click.option("--competent-cells-cost", type=float, default=None, help="Competent cells $/transformation.")
@click.option("--sequencing-cost", type=float, default=None, help="Plasmidsaurus sequencing $/reaction.")
@click.option("--quiet", "-q", is_flag=True, help="Suppress progress output.")
@click.option("--skip-constraints", is_flag=True,
              help="Skip constraint resolution (faster, uses unoptimized codons).")
def compile(spec, output, cost_only, fmt, researcher_rate, overhead,
            twist_gene_bp, twist_clonal_bp, idt_gblock_bp,
            bsai_cost, competent_cells_cost, sequencing_cost,
            quiet, skip_constraints):
    """Compile a construct spec into DNA sequences and an assembly plan."""
    # Configure logging
    if quiet:
        logging.basicConfig(level=logging.WARNING, format="%(message)s")
        # Suppress DNA Chisel progress bars
        logging.getLogger("dnachisel").setLevel(logging.WARNING)
    else:
        logging.basicConfig(level=logging.INFO, format="%(message)s")

    # Suppress tqdm progress bars from DNA Chisel
    import os
    if quiet:
        os.environ["DNACHISEL_QUIET"] = "1"

    cost_kwargs = dict(
        researcher_hourly_rate=researcher_rate,
        overhead_multiplier=overhead,
    )
    # Apply optional overrides from CLI flags
    if twist_gene_bp is not None:
        cost_kwargs["twist_gene_per_bp"] = twist_gene_bp
    if twist_clonal_bp is not None:
        cost_kwargs["twist_clonal_per_bp"] = twist_clonal_bp
    if idt_gblock_bp is not None:
        cost_kwargs["idt_gblock_per_bp"] = idt_gblock_bp
    if bsai_cost is not None:
        cost_kwargs["bsai_per_reaction"] = bsai_cost
    if competent_cells_cost is not None:
        cost_kwargs["competent_cells_per_transformation"] = competent_cells_cost
    if sequencing_cost is not None:
        cost_kwargs["plasmidsaurus_sequencing"] = sequencing_cost

    cost_params = CostParams(**cost_kwargs)

    try:
        graph, plan = compile_construct(
            spec,
            cost_params=cost_params,
            skip_constraints=skip_constraints,
            verbose=not quiet,
        )
    except Exception as e:
        click.secho(f"Compilation failed: {e}", fg="red", err=True)
        sys.exit(1)

    # Cost comparison output
    click.echo()
    click.secho("COST COMPARISON", fg="cyan", bold=True)
    click.secho("=" * 60, fg="cyan")

    for s in plan.strategies:
        if s.recommended:
            click.secho(f"  {s.name:<40s} ${s.total_cost:>8.2f}  "
                        f"({s.total_wall_clock_days[0]}-{s.total_wall_clock_days[1]}d)  "
                        f"<-- recommended", fg="green", bold=True)
        else:
            click.echo(f"  {s.name:<40s} ${s.total_cost:>8.2f}  "
                       f"({s.total_wall_clock_days[0]}-{s.total_wall_clock_days[1]}d)")

    if cost_only:
        return

    # Generate output files
    output_dir = Path(output) if output else Path(".")
    output_dir.mkdir(parents=True, exist_ok=True)

    name_slug = graph.name.replace(" ", "_").replace("/", "_")

    if fmt in ("genbank", "summary"):
        gb_path = output_dir / f"{name_slug}.gb"
        record = export_genbank(graph, gb_path)
        click.echo()
        click.secho(f"GenBank: {gb_path}", fg="green")
        click.echo(f"  Sequence: {len(record.seq)} bp, {len(record.features)} features")

    if fmt == "json":
        result = _build_json_result(graph, plan)
        json_path = output_dir / f"{name_slug}.json"
        json_path.write_text(json.dumps(result, indent=2))
        click.secho(f"JSON: {json_path}", fg="green")

    # Always write the assembly plan as YAML
    plan_path = output_dir / f"{name_slug}_plan.txt"
    plan_path.write_text(plan.summary())
    click.secho(f"Plan: {plan_path}", fg="green")


@cli.command()
@click.argument("spec", type=click.Path(exists=True))
def validate(spec):
    """Validate a construct spec without compiling."""
    from .frontend.parser import parse_spec

    logging.basicConfig(level=logging.WARNING, format="%(message)s")

    try:
        graph = parse_spec(spec)
    except Exception as e:
        click.secho(f"Parse error: {e}", fg="red")
        sys.exit(1)

    errors = graph.validate()

    if not errors:
        click.secho("Valid construct spec.", fg="green")
        click.echo(f"  {len(graph.parts())} parts")
    else:
        for e in errors:
            color = "red" if e.severity == "error" else "yellow"
            click.secho(f"  [{e.severity.upper()}] {e.location}: {e.message}", fg=color)

        error_count = sum(1 for e in errors if e.severity == "error")
        warning_count = sum(1 for e in errors if e.severity == "warning")
        if error_count:
            click.secho(f"\n{error_count} error(s), {warning_count} warning(s)", fg="red")
            sys.exit(1)
        else:
            click.secho(f"\n{warning_count} warning(s), no errors", fg="yellow")


@cli.command()
@click.option("--list", "list_type",
              type=click.Choice(["promoters", "rbs", "terminators", "tags", "cleavage", "linkers",
                                 "overhangs", "all"]),
              default="all", help="Which parts to list.")
def parts(list_type):
    """List available parts in the built-in database."""
    from .data.parts_db import (
        PROMOTERS, RBS_LIBRARY, TERMINATORS,
        PURIFICATION_TAGS, SOLUBILITY_TAGS, CLEAVAGE_SITES, LINKERS,
        HIGH_FIDELITY_OVERHANGS_BSAI,
    )

    sections = {
        "promoters": ("Promoters", PROMOTERS),
        "rbs": ("RBS / Translation Initiation", RBS_LIBRARY),
        "terminators": ("Terminators", TERMINATORS),
        "tags": ("Purification Tags", PURIFICATION_TAGS),
        "cleavage": ("Cleavage Sites", CLEAVAGE_SITES),
        "linkers": ("Linkers", LINKERS),
    }

    to_show = sections if list_type == "all" else {list_type: sections[list_type]}

    for key, (title, data) in to_show.items():
        click.secho(f"\n{title}", fg="cyan", bold=True)
        click.secho("-" * 40, fg="cyan")
        for name, info in data.items():
            notes = info.get("notes", "")
            click.echo(f"  {name:<20s} {notes}")

    if list_type in ("all", "overhangs"):
        click.secho(f"\nHigh-Fidelity BsaI Overhangs ({len(HIGH_FIDELITY_OVERHANGS_BSAI)})",
                     fg="cyan", bold=True)
        click.echo(f"  {', '.join(HIGH_FIDELITY_OVERHANGS_BSAI)}")


def _build_json_result(graph, plan) -> dict:
    """Build a JSON-serializable result dict."""
    parts_list = []
    for part in graph.parts():
        parts_list.append({
            "id": part.id,
            "name": part.name,
            "type": part.__class__.__name__,
            "length_bp": len(part.sequence) if part.sequence else None,
            "sequence": str(part.sequence) if part.sequence else None,
        })

    strategies = []
    for s in plan.strategies:
        strategies.append({
            "name": s.name,
            "recommended": s.recommended,
            "total_cost": round(s.total_cost, 2),
            "synthesis_cost": round(s.synthesis_cost, 2),
            "reagent_cost": round(s.reagent_cost_with_overhead, 2),
            "researcher_cost": round(s.researcher_cost, 2),
            "researcher_hours": s.researcher_time_hrs,
            "risk_surcharge": round(s.failure_risk_surcharge, 2),
            "wall_clock_days": list(s.total_wall_clock_days),
            "num_fragments": s.num_fragments,
            "overhangs": s.overhangs,
            "notes": s.notes,
        })

    return {
        "construct": {
            "name": graph.name,
            "host": graph.host_organism,
            "insert_length_bp": graph.total_insert_length(),
            "parts": parts_list,
        },
        "assembly_plan": {
            "strategies": strategies,
        },
    }


def main():
    cli()


if __name__ == "__main__":
    main()

"""
YAML frontend — parses a declarative construct spec into an IR graph.

The parser translates high-level intent (e.g., "expression: high")
into abstract IR nodes. Concrete sequences are resolved by later passes.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from ..core.graph import ConstructGraph
from ..core.types import (
    Constraint, ConstraintKind, no_restriction_site, gc_content,
    max_homopolymer,
)
from ..core.parts import (
    Backbone, Promoter, RBS, CDS, PurificationTag, SolubilityTag,
    CleavageSite, StopCodon, Linker, Terminator, Spacer,
    RegulationType, RBSDesignMethod, TagPosition, Origin,
)
from ..data.parts_db import PROMOTERS, RBS_TIER_DEFAULTS


def parse_spec(spec: dict | str | Path) -> ConstructGraph:
    """
    Parse a YAML construct spec into a ConstructGraph.

    Accepts a dict (already parsed YAML), a YAML string, or a Path to a file.
    """
    if isinstance(spec, (str, Path)):
        path = Path(spec)
        if path.exists():
            spec = yaml.safe_load(path.read_text())
        else:
            spec = yaml.safe_load(spec)

    root = spec.get("construct", spec)
    graph = ConstructGraph(
        name=root.get("name", "unnamed_construct"),
        host_organism=root.get("host", "e_coli"),
    )
    graph.metadata = {
        k: v for k, v in root.items()
        if k not in ("name", "host", "backbone", "cassette", "constraints", "synthesis")
    }

    # -- Backbone --
    bb_spec = root.get("backbone", {})
    if bb_spec:
        graph.add_part(_parse_backbone(bb_spec))

    # -- Cassette (expression cassette elements in order) --
    cassette = root.get("cassette", [])
    _part_counter = {"n": 0}

    def _next_id(prefix: str) -> str:
        _part_counter["n"] += 1
        return f"{prefix}_{_part_counter['n']}"

    for element in cassette:
        if isinstance(element, str):
            # Shorthand: just a part name, infer type from context
            # (not yet supported — reserved for future)
            continue

        for key, value in element.items():
            if key == "promoter":
                graph.add_part(_parse_promoter(value, _next_id("prom")))

            elif key == "cistron":
                _parse_cistron(value, graph, _next_id)

            elif key == "terminator":
                graph.add_part(_parse_terminator(value, _next_id("term")))

            elif key == "spacer":
                sp = Spacer(
                    id=_next_id("spacer"),
                    name=f"spacer_{value}" if isinstance(value, int) else "spacer",
                    length_bp=value if isinstance(value, int) else value.get("length", 30),
                )
                graph.add_part(sp)

    # -- Constraints --
    constraints_spec = root.get("constraints", {})
    _apply_constraints(constraints_spec, graph)

    # -- Synthesis preferences --
    synth_spec = root.get("synthesis", {})
    graph.metadata["synthesis"] = synth_spec

    # -- Connect linearly --
    graph.connect_linear()

    return graph


# ---------------------------------------------------------------------------
# Sub-parsers
# ---------------------------------------------------------------------------

def _parse_backbone(spec: dict | str) -> Backbone:
    if isinstance(spec, str):
        return Backbone(id="backbone", name=spec, source=spec)

    origin_str = spec.get("ori", "pBR322")
    try:
        origin = Origin(origin_str)
    except ValueError:
        origin = Origin.CUSTOM

    return Backbone(
        id="backbone",
        name=spec.get("name", "backbone"),
        origin=origin,
        resistance=spec.get("resistance", "kanamycin"),
        source=spec.get("source", ""),
        addgene_id=spec.get("addgene_id"),
    )


def _parse_promoter(spec: str | dict, part_id: str) -> Promoter:
    if isinstance(spec, str):
        prom_data = PROMOTERS.get(spec, {})
        reg_str = prom_data.get("regulation", "inducible")
        return Promoter(
            id=part_id,
            name=spec,
            regulation=RegulationType[reg_str.upper()],
            inducer=prom_data.get("inducer", ""),
        )

    return Promoter(
        id=part_id,
        name=spec.get("name", part_id),
        regulation=RegulationType[spec.get("regulation", "inducible").upper()],
        inducer=spec.get("inducer", ""),
        strength=spec.get("strength", "strong"),
    )


def _parse_terminator(spec: str | dict, part_id: str) -> Terminator:
    if isinstance(spec, str):
        return Terminator(id=part_id, name=spec)
    return Terminator(id=part_id, name=spec.get("name", part_id))


def _parse_cistron(spec: dict, graph: ConstructGraph, next_id) -> None:
    """
    Parse a cistron block, which may include:
    - rbs / expression level
    - N-terminal tags, cleavage sites, solubility tags
    - gene (the main CDS)
    - C-terminal tags
    - fused: true/false (whether this is fused to the previous cistron)
    """
    label = spec.get("label", "cistron")
    fused = spec.get("fused", True)  # default: fused within cistron

    # -- RBS --
    expression = spec.get("expression", "high")
    rbs_spec = spec.get("rbs")

    if rbs_spec and isinstance(rbs_spec, dict):
        method = rbs_spec.get("method", "lookup")
        rbs = RBS(
            id=next_id("rbs"),
            name=f"rbs_{label}",
            design_method=RBSDesignMethod[method.upper()],
            target_expression=str(rbs_spec.get("target_rate", expression)),
            part_name=rbs_spec.get("name", ""),
        )
    elif rbs_spec and isinstance(rbs_spec, str):
        rbs = RBS(
            id=next_id("rbs"),
            name=f"rbs_{label}",
            design_method=RBSDesignMethod.LOOKUP,
            part_name=rbs_spec,
            target_expression=expression,
        )
    else:
        # Resolve from expression tier
        default_rbs = RBS_TIER_DEFAULTS.get(expression, "BCD2")
        rbs = RBS(
            id=next_id("rbs"),
            name=f"rbs_{label}",
            design_method=RBSDesignMethod.BCD if "BCD" in default_rbs else RBSDesignMethod.LOOKUP,
            part_name=default_rbs,
            target_expression=expression,
        )

    graph.add_part(rbs)

    # -- Chain: N-terminal elements --
    chain = spec.get("chain", [])
    gene_spec = spec.get("gene")

    # If chain is used, process it in order
    if chain:
        for item in chain:
            for item_key, item_val in (item.items() if isinstance(item, dict) else [(item, {})]):
                _parse_chain_element(item_key, item_val, graph, next_id, label)
    else:
        # Flat format: n_tag, gene, c_tag
        n_tag = spec.get("n_tag")
        if n_tag:
            _parse_tag_shorthand(n_tag, graph, next_id, label, TagPosition.N_TERM)

        if gene_spec:
            _parse_gene(gene_spec, graph, next_id, label, is_last_in_cistron=not spec.get("c_tag"))

        c_tag = spec.get("c_tag")
        if c_tag:
            _parse_tag_shorthand(c_tag, graph, next_id, label, TagPosition.C_TERM)

    # -- Stop codon (if this is a non-fused, independently translated cistron) --
    # The CDS's has_stop flag determines this; parser sets it based on context.


def _parse_chain_element(key: str, value: Any, graph: ConstructGraph,
                         next_id, label: str) -> None:
    """Parse an element in the cistron chain."""
    from ..data.parts_db import PURIFICATION_TAGS, SOLUBILITY_TAGS, CLEAVAGE_SITES, LINKERS

    if key == "tag" or key in PURIFICATION_TAGS:
        tag_name = value if isinstance(value, str) else key
        graph.add_part(PurificationTag(
            id=next_id("tag"),
            name=tag_name,
            tag_type=tag_name,
        ))

    elif key == "solubility_tag" or key in SOLUBILITY_TAGS:
        tag_name = value if isinstance(value, str) else key
        graph.add_part(SolubilityTag(
            id=next_id("sol"),
            name=tag_name,
            tag_type=tag_name,
        ))

    elif key == "cleavage_site" or key in CLEAVAGE_SITES:
        site_name = value if isinstance(value, str) else key
        graph.add_part(CleavageSite(
            id=next_id("cleave"),
            name=site_name,
            protease=site_name,
        ))

    elif key == "linker":
        linker_name = value if isinstance(value, str) else value.get("type", "GS_flexible")
        repeats = value.get("repeats", 3) if isinstance(value, dict) else 3
        graph.add_part(Linker(
            id=next_id("linker"),
            name=linker_name,
            linker_type=linker_name,
            repeats=repeats,
        ))

    elif key == "gene":
        _parse_gene(value, graph, next_id, label, is_last_in_cistron=True)


def _parse_gene(spec: dict | str | tuple, graph: ConstructGraph,
                next_id, label: str, is_last_in_cistron: bool = True) -> None:
    """Parse a gene/CDS entry."""
    if isinstance(spec, str):
        # Simple: just a name or accession
        graph.add_part(CDS(
            id=next_id("cds"),
            name=spec,
            source_id=spec,
            has_stop=is_last_in_cistron,
        ))
        return

    if isinstance(spec, (list, tuple)):
        # Tuple format: (id, source)
        graph.add_part(CDS(
            id=next_id("cds"),
            name=spec[0],
            source_id=spec[0],
            source_db=spec[1] if len(spec) > 1 else "",
            has_stop=is_last_in_cistron,
        ))
        return

    # Dict format
    source = spec.get("source", "")
    gene_id = spec.get("id", spec.get("name", "unknown"))
    codon_opt = spec.get("codon_optimization", "local")

    graph.add_part(CDS(
        id=next_id("cds"),
        name=gene_id,
        source_db=source,
        source_id=gene_id,
        has_stop=is_last_in_cistron,
        codon_optimization=codon_opt,
    ))


def _parse_tag_shorthand(spec, graph, next_id, label, position):
    """Parse the flat n_tag/c_tag shorthand (tuple or string)."""
    from ..data.parts_db import PURIFICATION_TAGS, CLEAVAGE_SITES, SOLUBILITY_TAGS

    if isinstance(spec, str):
        items = [spec]
    elif isinstance(spec, (list, tuple)):
        items = list(spec)
    else:
        return

    for item in items:
        if item in PURIFICATION_TAGS:
            graph.add_part(PurificationTag(
                id=next_id("tag"), name=item, tag_type=item, position=position,
            ))
        elif item in CLEAVAGE_SITES:
            graph.add_part(CleavageSite(
                id=next_id("cleave"), name=item, protease=item,
            ))
        elif item in SOLUBILITY_TAGS:
            graph.add_part(SolubilityTag(
                id=next_id("sol"), name=item, tag_type=item, position=position,
            ))


def _apply_constraints(spec: dict, graph: ConstructGraph) -> None:
    """Apply global constraints from the spec to the graph."""
    assembly = spec.get("assembly", "golden_gate")
    enzyme = spec.get("enzyme", "BsaI")

    if assembly == "golden_gate":
        graph.constraints.append(no_restriction_site(enzyme))

    gc = spec.get("gc_window")
    if gc:
        graph.constraints.append(gc_content(gc[0], gc[1]))
    else:
        graph.constraints.append(gc_content())  # default 35-65%

    homopoly = spec.get("max_homopolymer")
    if homopoly:
        graph.constraints.append(max_homopolymer(homopoly))
    else:
        graph.constraints.append(max_homopolymer())  # default 6

    graph.metadata["assembly_method"] = assembly
    graph.metadata["assembly_enzyme"] = enzyme
    graph.metadata["codon_optimization"] = spec.get("codon_optimization", "local")

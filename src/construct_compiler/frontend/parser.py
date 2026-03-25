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
from ..data.parts_db import PROMOTERS, RBS_TIER_DEFAULTS, TWIST_CATALOG_VECTORS


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
    backbone = None
    if bb_spec:
        backbone = _parse_backbone(bb_spec)
        graph.add_part(backbone)

    # Track what the catalog vector provides (to skip duplicate elements)
    vector_provides = set(backbone.provides) if backbone else set()

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
                if "promoter" in vector_provides:
                    # Catalog vector already has a promoter — skip
                    continue
                graph.add_part(_parse_promoter(value, _next_id("prom")))

            elif key == "cistron":
                _parse_cistron(value, graph, _next_id, vector_provides)

            elif key == "terminator":
                if "terminator" in vector_provides:
                    # Catalog vector already has a terminator — skip
                    continue
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
        # Check if it's a catalog vector name
        if spec in TWIST_CATALOG_VECTORS:
            return _backbone_from_catalog(spec)
        return Backbone(id="backbone", name=spec, source=spec)

    # Check for catalog vector
    catalog_id = spec.get("catalog_vector") or spec.get("catalog_id")
    if catalog_id and catalog_id in TWIST_CATALOG_VECTORS:
        bb = _backbone_from_catalog(catalog_id)
        # Allow explicit cloning pair override from spec
        cp = spec.get("cloning_pair")
        if cp and isinstance(cp, list) and len(cp) == 2:
            bb.cloning_pair = (cp[0], cp[1])
        return bb

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


def _backbone_from_catalog(catalog_id: str) -> Backbone:
    """Create a Backbone from a Twist catalog vector entry."""
    entry = TWIST_CATALOG_VECTORS[catalog_id]
    origin_str = entry.get("ori", "pBR322")
    try:
        origin = Origin(origin_str)
    except ValueError:
        origin = Origin.CUSTOM

    return Backbone(
        id="backbone",
        name=catalog_id,
        origin=origin,
        resistance=entry.get("resistance", ""),
        source="twist_catalog",
        catalog_id=catalog_id,
        catalog_vendor="twist",
        provides=list(entry.get("provides", [])),
        metadata={
            "catalog_entry": entry,
            "size_bp": entry.get("size_bp", 0),
            "cloning_sites": entry.get("cloning_sites", []),
        },
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


def _parse_cistron(spec: dict, graph: ConstructGraph, next_id,
                   vector_provides: set | None = None) -> None:
    """
    Parse a cistron block, which may include:
    - rbs / expression level
    - N-terminal tags, cleavage sites, solubility tags
    - gene (the main CDS)
    - C-terminal tags
    - fused: true/false (whether this is fused to the previous cistron)

    vector_provides: set of element types the catalog vector already has
                     (e.g. {"promoter", "rbs", "n_tag", "cleavage", "terminator"}).
                     These elements are skipped for the first cistron.
    """
    if vector_provides is None:
        vector_provides = set()

    label = spec.get("label", "cistron")
    fused = spec.get("fused", True)  # default: fused within cistron

    # -- RBS --
    # Skip RBS only for the first cistron if the vector provides it
    skip_rbs = "rbs" in vector_provides
    expression = spec.get("expression", "high")
    rbs_spec = spec.get("rbs")

    if not skip_rbs:
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

    # When the user explicitly specifies a chain, always include all elements —
    # the user's insert parts are distinct from the vector's built-in features.
    # The vector's native tags/cleavage sites appear as backbone annotations.
    # Only skip auto-generated elements in the flat format (no explicit chain).
    if chain:
        for item in chain:
            for item_key, item_val in (item.items() if isinstance(item, dict) else [(item, {})]):
                _parse_chain_element(item_key, item_val, graph, next_id, label)
    else:
        # Flat format: n_tag, gene, c_tag
        # Here we DO skip elements the vector provides, since the user didn't
        # explicitly request them — they'd duplicate the vector's built-in features.
        skip_n_tag = "n_tag" in vector_provides
        skip_cleavage = "cleavage" in vector_provides
        skip_solubility = "solubility_tag" in vector_provides

        n_tag = spec.get("n_tag")
        if n_tag and not skip_n_tag:
            _parse_tag_shorthand(n_tag, graph, next_id, label, TagPosition.N_TERM)

        if gene_spec:
            _parse_gene(gene_spec, graph, next_id, label, is_last_in_cistron=not spec.get("c_tag"))

        c_tag = spec.get("c_tag")
        if c_tag:
            _parse_tag_shorthand(c_tag, graph, next_id, label, TagPosition.C_TERM)

    # -- Stop codon (if this is a non-fused, independently translated cistron) --
    # The CDS's has_stop flag determines this; parser sets it based on context.

    # After the first cistron, remove per-ORF elements from vector_provides
    # so subsequent cistrons add their own RBS/tags. Keep promoter/terminator
    # since those are shared across the whole cassette.
    vector_provides.discard("rbs")
    vector_provides.discard("n_tag")
    vector_provides.discard("cleavage")
    vector_provides.discard("solubility_tag")


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

"""
GenBank export backend.

Generates an annotated GenBank file from a compiled ConstructGraph.
Each part becomes a feature with appropriate qualifiers.
"""

from __future__ import annotations

import logging
from datetime import date
from pathlib import Path
from typing import Optional

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

from ..core.graph import ConstructGraph
from ..core.parts import (
    Backbone, Promoter, RBS, CDS, PurificationTag, SolubilityTag,
    CleavageSite, StopCodon, Linker, Terminator, Spacer,
)

logger = logging.getLogger(__name__)


# Map part types to GenBank feature keys
FEATURE_TYPE_MAP = {
    Promoter: "promoter",
    RBS: "RBS",
    CDS: "CDS",
    PurificationTag: "CDS",
    SolubilityTag: "CDS",
    CleavageSite: "CDS",
    Linker: "CDS",
    StopCodon: "CDS",
    Terminator: "terminator",
    Spacer: "misc_feature",
    Backbone: "rep_origin",
}


def export_genbank(graph: ConstructGraph,
                   output_path: Optional[str | Path] = None,
                   molecule_type: str = "DNA",
                   topology: str = "circular") -> SeqRecord:
    """
    Export the compiled construct as an annotated GenBank file.

    Returns the SeqRecord and optionally writes to a file.
    """
    # Build the full sequence
    full_seq = Seq("")
    part_positions: list[tuple] = []  # (part, start, end)

    for part in graph.parts():
        if isinstance(part, Backbone):
            continue  # backbone handled separately
        if part.sequence is None:
            logger.warning(f"Part '{part.name}' has no sequence, skipping")
            continue

        start = len(full_seq)
        full_seq += part.sequence
        end = len(full_seq)
        part_positions.append((part, start, end))

    # Create SeqRecord
    record = SeqRecord(
        full_seq,
        id=graph.name.replace(" ", "_")[:16],
        name=graph.name.replace(" ", "_")[:16],
        description=f"{graph.name} | {graph.host_organism} expression construct",
        annotations={
            "molecule_type": molecule_type,
            "topology": topology,
            "date": date.today().strftime("%d-%b-%Y").upper(),
            "organism": _organism_name(graph.host_organism),
        },
    )

    # Add features for each part
    for part, start, end in part_positions:
        feature = _make_feature(part, start, end)
        if feature:
            record.features.append(feature)

    # Add a source feature
    record.features.insert(0, SeqFeature(
        FeatureLocation(0, len(full_seq)),
        type="source",
        qualifiers={
            "organism": [_organism_name(graph.host_organism)],
            "mol_type": ["other DNA"],
            "note": [f"Designed by construct_compiler for {graph.host_organism}"],
        },
    ))

    # Write to file if path given
    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            SeqIO.write(record, f, "genbank")
        logger.info(f"GenBank file written: {output_path}")

    return record


def _make_feature(part, start: int, end: int) -> Optional[SeqFeature]:
    """Create a GenBank feature for a part."""
    feature_type = FEATURE_TYPE_MAP.get(type(part), "misc_feature")

    qualifiers = {
        "label": [part.name],
    }

    if isinstance(part, Promoter):
        qualifiers["note"] = [f"Promoter: {part.name}"]
        if part.inducer:
            qualifiers["note"].append(f"Inducer: {part.inducer}")

    elif isinstance(part, RBS):
        qualifiers["note"] = [f"RBS: {part.part_name or part.name}"]
        if part.metadata.get("relative_strength"):
            qualifiers["note"].append(
                f"Relative strength: {part.metadata['relative_strength']}"
            )

    elif isinstance(part, CDS):
        qualifiers["gene"] = [part.name]
        qualifiers["codon_start"] = [1]
        if part.protein_sequence:
            qualifiers["translation"] = [part.protein_sequence]
        if part.source_db:
            qualifiers["db_xref"] = [f"{part.source_db}:{part.source_id}"]

    elif isinstance(part, PurificationTag):
        feature_type = "misc_feature"
        qualifiers["note"] = [f"Purification tag: {part.tag_type}"]

    elif isinstance(part, SolubilityTag):
        feature_type = "misc_feature"
        qualifiers["note"] = [f"Solubility tag: {part.tag_type}"]
        if part.metadata.get("protein_sequence"):
            qualifiers["translation"] = [part.metadata["protein_sequence"]]

    elif isinstance(part, CleavageSite):
        feature_type = "misc_feature"
        qualifiers["note"] = [
            f"Protease cleavage site: {part.protease}",
            f"Recognition: {part.recognition_seq}",
        ]

    elif isinstance(part, Linker):
        feature_type = "misc_feature"
        qualifiers["note"] = [f"Linker: {part.linker_type}"]

    elif isinstance(part, Terminator):
        qualifiers["note"] = [f"Terminator: {part.name}"]

    elif isinstance(part, Spacer):
        qualifiers["note"] = [f"Intergenic spacer ({end - start} bp)"]

    return SeqFeature(
        FeatureLocation(start, end),
        type=feature_type,
        qualifiers=qualifiers,
    )


def _organism_name(host: str) -> str:
    """Convert shorthand to full organism name."""
    mapping = {
        "e_coli": "Escherichia coli",
        "e_coli_bl21": "Escherichia coli BL21(DE3)",
        "s_cerevisiae": "Saccharomyces cerevisiae",
        "h_sapiens": "Homo sapiens",
    }
    return mapping.get(host, host)

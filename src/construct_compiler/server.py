"""
FastAPI server for the construct compiler web UI.

Run with:
    uvicorn construct_compiler.server:app --reload

Or:
    python -m construct_compiler.server
"""

from __future__ import annotations

import logging
import os
from pathlib import Path

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, HTMLResponse
from pydantic import BaseModel

# Suppress DNA Chisel progress bars
os.environ["TQDM_DISABLE"] = "1"
logging.basicConfig(level=logging.WARNING)
logging.getLogger("dnachisel").setLevel(logging.WARNING)

app = FastAPI(
    title="Construct Compiler",
    description="Genetic construct design compiler API",
    version="0.1.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# ---------------------------------------------------------------------------
# Request / response models
# ---------------------------------------------------------------------------

class CostParamsRequest(BaseModel):
    """All configurable cost parameters — mirrors CostParams defaults."""
    researcher_hourly_rate: float = 150.0
    overhead_multiplier: float = 1.5

    # Synthesis (per bp)
    twist_gene_per_bp: float = 0.07
    twist_clonal_per_bp: float = 0.09
    idt_gblock_per_bp: float = 0.08

    # Reagents (per reaction, before overhead)
    bsai_per_reaction: float = 3.00
    t4_ligase_per_reaction: float = 0.25
    competent_cells_per_transformation: float = 10.00
    plates_antibiotics: float = 2.00
    colony_pcr_screening: float = 5.00
    miniprep: float = 5.00
    plasmidsaurus_sequencing: float = 15.00

    # Labor (hours)
    golden_gate_setup_hrs: float = 1.5
    colony_screening_hrs: float = 0.0
    miniprep_sequencing_hrs: float = 1.0
    troubleshooting_hrs: float = 3.0

    # Success rates
    two_part_success_rate: float = 0.90
    three_part_success_rate: float = 0.80


class CheckRequest(BaseModel):
    """A construct spec to validate (compile + run all 4 checks)."""
    spec: dict
    skip_constraints: bool = False


class CheckResultItem(BaseModel):
    check_name: str
    severity: str  # "ERROR", "WARNING", "INFO"
    part_id: str
    message: str
    details: dict = {}


class CheckResponse(BaseModel):
    success: bool
    passed: bool = False
    score: float = 0.0
    error_count: int = 0
    warning_count: int = 0
    check_count: int = 0
    construct_name: str = ""
    insert_length_bp: int | None = None
    cistron_count: int = 0
    compile_time_s: float = 0.0
    checks: list[CheckResultItem] = []
    errors: list[str] = []


class CompileRequest(BaseModel):
    """A construct spec submitted from the web UI."""
    spec: dict
    cost_params: CostParamsRequest = CostParamsRequest()
    skip_constraints: bool = False


class PartInfo(BaseModel):
    id: str
    name: str
    type: str
    length_bp: int | None = None
    sequence: str | None = None
    color: str = "#cccccc"


class LineItemInfo(BaseModel):
    category: str
    description: str
    unit_cost: float = 0.0
    quantity: float = 1.0
    multiplier: float = 1.0
    subtotal: float = 0.0
    hours: float = 0.0


class StrategyInfo(BaseModel):
    name: str
    recommended: bool = False
    total_cost: float
    synthesis_cost: float
    reagent_cost: float
    researcher_cost: float
    researcher_hours: float
    risk_surcharge: float
    wall_clock_days: list[int]
    num_fragments: int
    overhangs: list[str] = []
    notes: list[str] = []
    line_items: list[LineItemInfo] = []


class PlasmidAnnotation(BaseModel):
    """An annotation on the full plasmid for seqviz rendering."""
    name: str
    start: int
    end: int
    direction: int = 1   # 1 = forward, -1 = reverse
    color: str = "#cccccc"
    type: str = ""        # "insert_part", "backbone_feature"

class TranslationRange(BaseModel):
    """A CDS region to show amino acid translation in seqviz."""
    start: int
    end: int
    direction: int = 1

class PlasmidMapData(BaseModel):
    """Everything seqviz needs to render the full plasmid."""
    name: str = ""
    sequence: str = ""           # full plasmid sequence (insert real, backbone N-filled)
    total_length_bp: int = 0
    insert_length_bp: int = 0
    backbone_length_bp: int = 0
    backbone_name: str = ""
    backbone_catalog_id: str = ""
    backbone_resistance: str = ""
    backbone_ori: str = ""
    annotations: list[PlasmidAnnotation] = []
    translations: list[TranslationRange] = []
    enzymes: list[str] = []      # restriction enzyme names to highlight

class CompileResponse(BaseModel):
    success: bool
    construct_name: str = ""
    host: str = ""
    insert_length_bp: int = 0
    parts: list[PartInfo] = []
    plasmid_map: PlasmidMapData = PlasmidMapData()
    strategies: list[StrategyInfo] = []
    genbank: str = ""
    errors: list[str] = []
    warnings: list[str] = []


# ---------------------------------------------------------------------------
# Color palette for part types
# ---------------------------------------------------------------------------

PART_COLORS = {
    "Promoter": "#4CAF50",
    "RBS": "#2196F3",
    "CDS": "#FF9800",
    "PurificationTag": "#9C27B0",
    "SolubilityTag": "#E91E63",
    "CleavageSite": "#F44336",
    "Linker": "#607D8B",
    "Terminator": "#795548",
    "Spacer": "#BDBDBD",
    "Backbone": "#37474F",
    "StopCodon": "#F44336",
}


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------

@app.post("/api/compile", response_model=CompileResponse)
async def compile_endpoint(req: CompileRequest):
    """Compile a construct spec and return results."""
    from .passes.pipeline import compile_construct
    from .passes.assembly_planning import CostParams
    from .backends.genbank import export_genbank
    from io import StringIO
    from Bio import SeqIO

    try:
        cp = req.cost_params
        cost_params = CostParams(
            researcher_hourly_rate=cp.researcher_hourly_rate,
            overhead_multiplier=cp.overhead_multiplier,
            twist_gene_per_bp=cp.twist_gene_per_bp,
            twist_clonal_per_bp=cp.twist_clonal_per_bp,
            idt_gblock_per_bp=cp.idt_gblock_per_bp,
            bsai_per_reaction=cp.bsai_per_reaction,
            t4_ligase_per_reaction=cp.t4_ligase_per_reaction,
            competent_cells_per_transformation=cp.competent_cells_per_transformation,
            plates_antibiotics=cp.plates_antibiotics,
            colony_pcr_screening=cp.colony_pcr_screening,
            miniprep=cp.miniprep,
            plasmidsaurus_sequencing=cp.plasmidsaurus_sequencing,
            golden_gate_setup_hrs=cp.golden_gate_setup_hrs,
            colony_screening_hrs=cp.colony_screening_hrs,
            miniprep_sequencing_hrs=cp.miniprep_sequencing_hrs,
            troubleshooting_hrs=cp.troubleshooting_hrs,
            two_part_success_rate=cp.two_part_success_rate,
            three_part_success_rate=cp.three_part_success_rate,
        )

        graph, plan = compile_construct(
            req.spec,
            cost_params=cost_params,
            skip_constraints=req.skip_constraints,
            verbose=False,
        )

        # Collect validation issues
        errors_list = []
        warnings_list = []
        for e in graph.validate():
            if e.severity == "error":
                errors_list.append(f"{e.location}: {e.message}")
            else:
                warnings_list.append(f"{e.location}: {e.message}")

        # Build parts list & plasmid map annotations
        from .core.parts import Backbone, CDS, SolubilityTag, PurificationTag
        from .data.parts_db import TWIST_CATALOG_VECTORS

        parts = []
        backbone_part = None
        insert_annotations = []
        translations = []
        insert_seq_parts = []   # collect sequence fragments for full insert
        cursor = 0  # track bp position along the insert

        # Part types that encode protein (show translation)
        TRANSLATABLE_TYPES = (CDS, SolubilityTag, PurificationTag)

        for part in graph.parts():
            if isinstance(part, Backbone):
                backbone_part = part
                continue
            part_len = len(part.sequence) if part.sequence else 50
            parts.append(PartInfo(
                id=part.id,
                name=part.name,
                type=part.__class__.__name__,
                length_bp=part_len if part.sequence else None,
                sequence=str(part.sequence)[:200] if part.sequence else None,
                color=PART_COLORS.get(part.__class__.__name__, "#cccccc"),
            ))
            insert_annotations.append(PlasmidAnnotation(
                name=part.name,
                start=cursor,
                end=cursor + part_len,
                direction=1,
                color=PART_COLORS.get(part.__class__.__name__, "#cccccc"),
                type="insert_part",
            ))

            # Collect real sequence for the insert
            if part.sequence:
                insert_seq_parts.append(str(part.sequence))
            else:
                insert_seq_parts.append("N" * part_len)

            # Mark CDS / tag regions for amino acid translation display
            if isinstance(part, TRANSLATABLE_TYPES) and part.sequence:
                translations.append(TranslationRange(
                    start=cursor,
                    end=cursor + part_len,
                    direction=1,
                ))

            cursor += part_len

        insert_length = cursor
        insert_sequence = "".join(insert_seq_parts)

        # Build backbone features with approximate sizes
        backbone_name = ""
        backbone_catalog_id = ""
        backbone_resistance = ""
        backbone_ori = ""
        backbone_length = 3000  # default estimate for custom backbone
        backbone_annotations = []

        BACKBONE_FEATURE_COLORS = {
            "resistance": "#e74c3c",
            "ori": "#3498db",
            "promoter": "#2ecc71",
            "rbs": "#58a6ff",
            "terminator": "#8e44ad",
            "tag": "#f39c12",
            "cleavage": "#e67e22",
            "lacI": "#1abc9c",
            "f1_ori": "#9b59b6",
            "enhancer": "#d4ac0d",
            "selection": "#c0392b",
            "lentiviral": "#8e44ad",
        }

        # Approximate sizes for common backbone elements (bp)
        BACKBONE_FEATURE_SIZES = {
            "resistance": {"ampicillin": 861, "kanamycin": 795, "chloramphenicol": 660,
                           "spectinomycin": 789},
            "ori": {"pBR322": 600, "pUC": 600, "p15A": 400, "ColA": 400,
                    "CloDF13": 400, "SC101": 800, "pMB1": 600},
            "lacI": 1083,
            "f1_ori": 440,
            # Vector-provided expression cassette features (upstream of insert)
            "promoter": {"T7": 20, "T7lac": 44, "tac": 34, "lac": 30,
                         "CMV": 600, "EF1a": 1200, "SFFV": 500, "SV40": 350},
            "rbs": 34,  # typical RBS + spacer
            "terminator_cassette": {"T7": 40, "rrnB": 200, "bGH_polyA": 250,
                                    "SV40_polyA": 200, "hBG_polyA": 400},
            "his_tag": 24,       # 6xHis = 18nt + linker
            "thrombin_site": 18,   # LVPR|GS
            "tev_site": 21,        # ENLYFQ|S
            "enterokinase_site": 15,  # DDDDK
            "s_tag": 45,          # S-tag (15 aa)
            "trx_tag": 330,       # thioredoxin
            "kozak": 6,           # GCCACC
            # Mammalian-specific
            "beta_globin_intron": 600,
            "wpre": 600,
            "orip_element": 1800,
            "neo_cassette": 1600,   # SV40-Neo-polyA
            "hygro_cassette": 1700,
            "puro_cassette": 1200,  # PGK/UbC-Puro
            # Lentiviral elements
            "ltr_5": 600,
            "ltr_3_sin": 250,
            "psi_packaging": 400,
            "rre": 750,
            "cppt": 120,
        }

        if backbone_part:
            backbone_name = backbone_part.name
            backbone_catalog_id = backbone_part.catalog_id or ""
            backbone_resistance = backbone_part.resistance
            backbone_ori = backbone_part.origin.value if backbone_part.origin else ""

            catalog_entry = None
            if backbone_part.catalog_id:
                catalog_entry = TWIST_CATALOG_VECTORS.get(backbone_part.catalog_id)

            if catalog_entry and catalog_entry.get("size_bp"):
                backbone_length = catalog_entry["size_bp"] - insert_length
                if backbone_length < 1000:
                    backbone_length = catalog_entry["size_bp"]

            promoter_name = ""
            if catalog_entry:
                promoter_name = catalog_entry.get("promoter", "")

            # ---------------------------------------------------------------
            # Build backbone annotations in two groups:
            #
            # GROUP A — "downstream" features: placed right after the insert,
            #   reading clockwise from the MCS.
            #   pET layout: [INSERT] → T7 term → f1 ori → Kan(←) → lacI(←) → pBR322 ori
            #
            # GROUP B — "upstream MCS" features: placed at the end of the
            #   backbone, wrapping around to sit just before the insert.
            #   pET layout: ... pBR322 ori → T7 promoter → RBS → His → Thrombin → [INSERT]
            # ---------------------------------------------------------------

            downstream_annotations = []   # Group A
            upstream_annotations = []     # Group B
            bb_cursor = insert_length     # cursor for downstream placement

            # === GROUP A: Downstream of insert ===

            # Terminator on vector (immediately after insert)
            provides = set()
            if catalog_entry:
                provides = set(catalog_entry.get("provides", []))
            if "terminator" in provides:
                term_size = 200
                if promoter_name in ("T7", "T7lac"):
                    term_size = BACKBONE_FEATURE_SIZES["terminator_cassette"].get("T7", 40)
                downstream_annotations.append(PlasmidAnnotation(
                    name="T7 terminator" if "T7" in promoter_name else "Terminator",
                    start=bb_cursor,
                    end=bb_cursor + term_size,
                    direction=1,
                    color=BACKBONE_FEATURE_COLORS["terminator"],
                    type="backbone_feature",
                ))
                bb_cursor += term_size

            # Mammalian-specific post-insert elements (WPRE, enhancers downstream)
            if catalog_entry:
                enhancers = catalog_entry.get("enhancers", [])
                placed_enhancers = set()
                # WPRE goes downstream of insert (post-transcriptional enhancer)
                for enh in enhancers:
                    enh_lower = enh.lower()
                    if "wpre" in enh_lower:
                        enh_size = BACKBONE_FEATURE_SIZES["wpre"]
                        downstream_annotations.append(PlasmidAnnotation(
                            name="WPRE", start=bb_cursor, end=bb_cursor + enh_size,
                            direction=1, color=BACKBONE_FEATURE_COLORS["enhancer"],
                            type="backbone_feature",
                        ))
                        bb_cursor += enh_size
                        placed_enhancers.add(enh)

            # f1 origin (pET vectors)
            is_pet = "pET" in (backbone_part.catalog_id or "")
            if catalog_entry and is_pet:
                f1_size = BACKBONE_FEATURE_SIZES["f1_ori"]
                downstream_annotations.append(PlasmidAnnotation(
                    name="f1 ori", start=bb_cursor, end=bb_cursor + f1_size,
                    direction=-1,  # f1 ori typically reverse on pET
                    color=BACKBONE_FEATURE_COLORS["f1_ori"],
                    type="backbone_feature",
                ))
                bb_cursor += f1_size

            # Resistance cassette (reverse orientation on pET/pBR322 vectors)
            res_name = backbone_resistance or "resistance"
            res_size = BACKBONE_FEATURE_SIZES["resistance"].get(res_name, 800)
            # pBR322-based vectors: resistance is counterclockwise
            res_direction = -1 if backbone_ori in ("pBR322",) else 1
            downstream_annotations.append(PlasmidAnnotation(
                name=f"{res_name.capitalize()} resistance",
                start=bb_cursor, end=bb_cursor + res_size,
                direction=res_direction,
                color=BACKBONE_FEATURE_COLORS["resistance"],
                type="backbone_feature",
            ))
            bb_cursor += res_size

            # lacI repressor (for T7lac / tac — reverse orientation)
            if promoter_name in ("T7lac", "tac", "lacUV5"):
                laci_size = BACKBONE_FEATURE_SIZES["lacI"]
                downstream_annotations.append(PlasmidAnnotation(
                    name="lacI repressor",
                    start=bb_cursor, end=bb_cursor + laci_size,
                    direction=-1,
                    color=BACKBONE_FEATURE_COLORS["lacI"],
                    type="backbone_feature",
                ))
                bb_cursor += laci_size

            # Mammalian selection cassette (downstream, own promoter)
            if catalog_entry:
                mam_sel = catalog_entry.get("mammalian_selection", "")
                if mam_sel:
                    sel_map = {
                        "neomycin": ("Neo/G418", BACKBONE_FEATURE_SIZES["neo_cassette"]),
                        "hygromycin": ("Hygromycin", BACKBONE_FEATURE_SIZES["hygro_cassette"]),
                        "puromycin": ("Puromycin", BACKBONE_FEATURE_SIZES["puro_cassette"]),
                    }
                    sel_name, sel_size = sel_map.get(mam_sel, (mam_sel, 1400))
                    downstream_annotations.append(PlasmidAnnotation(
                        name=f"{sel_name} selection",
                        start=bb_cursor, end=bb_cursor + sel_size,
                        direction=1,
                        color=BACKBONE_FEATURE_COLORS["selection"],
                        type="backbone_feature",
                    ))
                    bb_cursor += sel_size

            # Origin of replication
            ori_name = backbone_ori
            ori_size = BACKBONE_FEATURE_SIZES["ori"].get(ori_name, 600)
            downstream_annotations.append(PlasmidAnnotation(
                name=f"{ori_name} ori",
                start=bb_cursor, end=bb_cursor + ori_size,
                direction=1,
                color=BACKBONE_FEATURE_COLORS["ori"],
                type="backbone_feature",
            ))
            bb_cursor += ori_size

            # Lentiviral structural elements (downstream, part of backbone)
            if catalog_entry:
                lenti_elements = catalog_entry.get("lentiviral_elements", [])
                lenti_sizes = {
                    "5' LTR": ("5' LTR", BACKBONE_FEATURE_SIZES["ltr_5"]),
                    "3' SIN-LTR": ("3' SIN-LTR", BACKBONE_FEATURE_SIZES["ltr_3_sin"]),
                    "Psi": ("Ψ packaging", BACKBONE_FEATURE_SIZES["psi_packaging"]),
                    "RRE": ("RRE", BACKBONE_FEATURE_SIZES["rre"]),
                    "cPPT": ("cPPT/CTS", BACKBONE_FEATURE_SIZES["cppt"]),
                    "WPRE": ("WPRE", BACKBONE_FEATURE_SIZES["wpre"]),
                }
                for elem in lenti_elements:
                    # Skip WPRE if already placed downstream
                    if elem == "WPRE" and elem in [e for e in catalog_entry.get("enhancers", [])
                                                    if "wpre" in e.lower()]:
                        continue
                    if elem == "WPRE" and any("wpre" in e.lower()
                                              for e in catalog_entry.get("enhancers", [])):
                        continue
                    if elem in lenti_sizes:
                        elem_name, elem_size = lenti_sizes[elem]
                    else:
                        elem_name, elem_size = elem, 300
                    downstream_annotations.append(PlasmidAnnotation(
                        name=elem_name,
                        start=bb_cursor, end=bb_cursor + elem_size,
                        direction=1 if "3'" not in elem else -1,
                        color=BACKBONE_FEATURE_COLORS["lentiviral"],
                        type="backbone_feature",
                    ))
                    bb_cursor += elem_size

            # === GROUP B: Upstream of insert (MCS flanking region) ===
            # These wrap around from the end of the backbone to sit just
            # before position 0 (the insert start).

            if catalog_entry:
                vec_promoter = catalog_entry.get("promoter", "")

                # Upstream enhancers (beta-globin intron, OriP — before insert)
                enhancers = catalog_entry.get("enhancers", [])
                for enh in enhancers:
                    if enh in placed_enhancers:
                        continue  # already placed downstream (e.g., WPRE)
                    enh_lower = enh.lower()
                    if "beta-globin" in enh_lower or "betaglobin" in enh_lower:
                        enh_size = BACKBONE_FEATURE_SIZES["beta_globin_intron"]
                        enh_name = "β-globin intron"
                    elif "orip" in enh_lower:
                        enh_size = BACKBONE_FEATURE_SIZES["orip_element"]
                        enh_name = "OriP (EBV)"
                    else:
                        enh_size = 400
                        enh_name = enh
                    upstream_annotations.append(PlasmidAnnotation(
                        name=enh_name,
                        start=bb_cursor, end=bb_cursor + enh_size,
                        direction=1,
                        color=BACKBONE_FEATURE_COLORS["enhancer"],
                        type="backbone_feature",
                    ))
                    bb_cursor += enh_size

                # Promoter on vector
                if vec_promoter and "promoter" in provides:
                    p_size = BACKBONE_FEATURE_SIZES["promoter"].get(vec_promoter, 100)
                    upstream_annotations.append(PlasmidAnnotation(
                        name=f"{vec_promoter} promoter",
                        start=bb_cursor, end=bb_cursor + p_size,
                        direction=1,
                        color=BACKBONE_FEATURE_COLORS["promoter"],
                        type="backbone_feature",
                    ))
                    bb_cursor += p_size

                # Kozak sequence (mammalian vectors)
                if "kozak" in provides:
                    kozak_size = BACKBONE_FEATURE_SIZES["kozak"]
                    upstream_annotations.append(PlasmidAnnotation(
                        name="Kozak",
                        start=bb_cursor, end=bb_cursor + kozak_size,
                        direction=1, color="#58a6ff",
                        type="backbone_feature",
                    ))
                    bb_cursor += kozak_size

                # RBS (prokaryotic vectors)
                if "rbs" in provides:
                    rbs_size = BACKBONE_FEATURE_SIZES["rbs"]
                    upstream_annotations.append(PlasmidAnnotation(
                        name="RBS",
                        start=bb_cursor, end=bb_cursor + rbs_size,
                        direction=1,
                        color=BACKBONE_FEATURE_COLORS["rbs"],
                        type="backbone_feature",
                    ))
                    bb_cursor += rbs_size

                # N-terminal tags (His, S-tag — on the vector backbone)
                n_tags = catalog_entry.get("n_tags", [])
                for tag in n_tags:
                    if tag == "6xHis":
                        tag_size = BACKBONE_FEATURE_SIZES["his_tag"]
                        tag_name = "N-His₆ (vector)"
                    else:
                        tag_size = 60
                        tag_name = f"{tag} (vector)"
                    upstream_annotations.append(PlasmidAnnotation(
                        name=tag_name,
                        start=bb_cursor, end=bb_cursor + tag_size,
                        direction=1,
                        color=BACKBONE_FEATURE_COLORS["tag"],
                        type="backbone_feature",
                    ))
                    bb_cursor += tag_size

                # Solubility tag on vector (e.g., Trx on pET-32a)
                if "solubility_tag" in provides:
                    notes = catalog_entry.get("notes", "").lower()
                    if "trx" in notes or "thioredoxin" in notes:
                        sol_name = "Trx (vector)"
                        sol_size = BACKBONE_FEATURE_SIZES["trx_tag"]
                    elif "s-tag" in notes:
                        sol_name = "S-tag (vector)"
                        sol_size = BACKBONE_FEATURE_SIZES["s_tag"]
                    else:
                        sol_name, sol_size = "Solubility tag (vector)", 200
                    upstream_annotations.append(PlasmidAnnotation(
                        name=sol_name,
                        start=bb_cursor, end=bb_cursor + sol_size,
                        direction=1, color="#f778ba",
                        type="backbone_feature",
                    ))
                    bb_cursor += sol_size

                # Cleavage sites on vector (Thrombin, Enterokinase)
                cleavage_sites = catalog_entry.get("cleavage_sites", [])
                for site in cleavage_sites:
                    site_lower = site.lower()
                    if "thrombin" in site_lower:
                        cs_size = BACKBONE_FEATURE_SIZES["thrombin_site"]
                    elif "tev" in site_lower:
                        cs_size = BACKBONE_FEATURE_SIZES["tev_site"]
                    elif "enterokinase" in site_lower:
                        cs_size = BACKBONE_FEATURE_SIZES["enterokinase_site"]
                    else:
                        cs_size = 18
                    upstream_annotations.append(PlasmidAnnotation(
                        name=f"{site} site (vector)",
                        start=bb_cursor, end=bb_cursor + cs_size,
                        direction=1,
                        color=BACKBONE_FEATURE_COLORS["cleavage"],
                        type="backbone_feature",
                    ))
                    bb_cursor += cs_size

            # Combine: downstream first, then upstream (wraps to before insert)
            backbone_annotations = downstream_annotations + upstream_annotations

            # Adjust backbone_length to fit all placed features
            placed_bb = bb_cursor - insert_length
            if placed_bb > backbone_length:
                backbone_length = placed_bb + 200  # small spacer

        total_length = insert_length + backbone_length
        all_annotations = insert_annotations + backbone_annotations

        # Build full plasmid sequence: real insert + N-filled backbone
        full_sequence = insert_sequence + ("N" * backbone_length)

        # Determine which enzymes to highlight (from the assembly constraint)
        assembly_enzymes = []
        spec_constraints = req.spec.get("construct", {}).get("constraints", {})
        enzyme_name = spec_constraints.get("enzyme", "BsaI")
        if enzyme_name:
            assembly_enzymes.append(enzyme_name)

        plasmid_map = PlasmidMapData(
            name=graph.name,
            sequence=full_sequence,
            total_length_bp=total_length,
            insert_length_bp=insert_length,
            backbone_length_bp=backbone_length,
            backbone_name=backbone_name,
            backbone_catalog_id=backbone_catalog_id,
            backbone_resistance=backbone_resistance,
            backbone_ori=backbone_ori,
            annotations=all_annotations,
            translations=translations,
            enzymes=assembly_enzymes,
        )

        # Build strategies list
        strategies = []
        for s in plan.strategies:
            line_items = [
                LineItemInfo(
                    category=li.category,
                    description=li.description,
                    unit_cost=round(li.unit_cost, 4),
                    quantity=round(li.quantity, 2),
                    multiplier=round(li.multiplier, 2),
                    subtotal=round(li.subtotal, 2),
                    hours=round(li.hours, 2),
                )
                for li in s.line_items
            ]
            strategies.append(StrategyInfo(
                name=s.name,
                recommended=s.recommended,
                total_cost=round(s.total_cost, 2),
                synthesis_cost=round(s.synthesis_cost, 2),
                reagent_cost=round(s.reagent_cost_with_overhead, 2),
                researcher_cost=round(s.researcher_cost, 2),
                researcher_hours=s.researcher_time_hrs,
                risk_surcharge=round(s.failure_risk_surcharge, 2),
                wall_clock_days=list(s.total_wall_clock_days),
                num_fragments=s.num_fragments,
                overhangs=s.overhangs,
                notes=s.notes,
                line_items=line_items,
            ))

        # Generate GenBank string
        genbank_str = ""
        record = export_genbank(graph)
        if record:
            buf = StringIO()
            SeqIO.write(record, buf, "genbank")
            genbank_str = buf.getvalue()

        return CompileResponse(
            success=True,
            construct_name=graph.name,
            host=graph.host_organism,
            insert_length_bp=insert_length,
            parts=parts,
            plasmid_map=plasmid_map,
            strategies=strategies,
            genbank=genbank_str,
            errors=errors_list,
            warnings=warnings_list,
        )

    except Exception as e:
        return CompileResponse(
            success=False,
            errors=[str(e)],
        )


@app.post("/api/check", response_model=CheckResponse)
async def check_endpoint(req: CheckRequest):
    """Compile a construct spec and run all 4 validity checks."""
    from .validation.harness import evaluate_spec

    try:
        result = evaluate_spec(
            req.spec,
            skip_constraints=req.skip_constraints,
        )

        checks = []
        for e in result.errors:
            checks.append(CheckResultItem(
                check_name=e["check_name"],
                severity="ERROR",
                part_id=e["part_id"],
                message=e["message"],
                details=e.get("details", {}),
            ))
        for w in result.warnings:
            checks.append(CheckResultItem(
                check_name=w["check_name"],
                severity="WARNING",
                part_id=w["part_id"],
                message=w["message"],
                details=w.get("details", {}),
            ))

        return CheckResponse(
            success=True,
            passed=result.passed,
            score=round(result.score, 4),
            error_count=result.error_count,
            warning_count=result.warning_count,
            check_count=result.check_count,
            construct_name=result.construct_name,
            insert_length_bp=result.insert_length_bp,
            cistron_count=result.cistron_count,
            compile_time_s=round(result.compile_time_s, 3),
            checks=checks,
        )

    except Exception as e:
        return CheckResponse(
            success=False,
            errors=[str(e)],
        )


@app.get("/api/parts/promoters")
async def get_promoters():
    from .data.parts_db import PROMOTERS
    return {k: {"regulation": v["regulation"], "inducer": v.get("inducer", "")}
            for k, v in PROMOTERS.items()}


@app.get("/api/parts/rbs")
async def get_rbs():
    from .data.parts_db import RBS_LIBRARY
    return {k: {"relative_strength": v.get("relative_strength", 0), "tier": v.get("tier", "")}
            for k, v in RBS_LIBRARY.items()}


@app.get("/api/parts/terminators")
async def get_terminators():
    from .data.parts_db import TERMINATORS
    return {k: {"type": v.get("type", "")} for k, v in TERMINATORS.items()}


@app.get("/api/parts/tags")
async def get_tags():
    from .data.parts_db import PURIFICATION_TAGS, SOLUBILITY_TAGS, SPLIT_FP_DETECTORS
    return {
        "purification": list(PURIFICATION_TAGS.keys()),
        "solubility": list(SOLUBILITY_TAGS.keys()),
        "split_fp_detectors": list(SPLIT_FP_DETECTORS.keys()),
    }


@app.get("/api/parts/cleavage")
async def get_cleavage():
    from .data.parts_db import CLEAVAGE_SITES
    return {k: {"recognition": v["protein_sequence"]} for k, v in CLEAVAGE_SITES.items()}


@app.get("/api/parts/linkers")
async def get_linkers():
    from .data.parts_db import LINKERS
    return {k: {"unit": v["unit"], "default_repeats": v["default_repeats"],
                "notes": v.get("notes", "")}
            for k, v in LINKERS.items()}


@app.get("/api/parts/2a_peptides")
async def get_2a_peptides():
    from .data.parts_db import SELF_CLEAVING_2A
    return {k: {"sequence": v["protein_sequence"],
                "length_aa": len(v["protein_sequence"]),
                "notes": v.get("notes", "")}
            for k, v in SELF_CLEAVING_2A.items()}


def _serialize_vector(v: dict) -> dict:
    """Serialize a catalog vector entry for the API."""
    return {
        "description": v.get("description", ""),
        "category": v.get("category", ""),
        "promoter": v.get("promoter", ""),
        "resistance": v.get("resistance", ""),
        "ori": v.get("ori", ""),
        "host": v.get("host", ""),
        "copy_number": v.get("copy_number", ""),
        "n_tags": v.get("n_tags", []),
        "c_tags": v.get("c_tags", []),
        "cleavage_sites": v.get("cleavage_sites", []),
        "provides": v.get("provides", []),
        "mammalian_selection": v.get("mammalian_selection", ""),
        "enhancers": v.get("enhancers", []),
        "lentiviral_elements": v.get("lentiviral_elements", []),
        "notes": v.get("notes", ""),
        "size_bp": v.get("size_bp", 0),
    }


@app.get("/api/vectors/catalog")
async def get_catalog_vectors():
    """Return all Twist catalog vectors."""
    from .data.parts_db import TWIST_CATALOG_VECTORS
    return {k: _serialize_vector(v) for k, v in TWIST_CATALOG_VECTORS.items()}


@app.get("/api/vectors/catalog/ecoli")
async def get_ecoli_catalog_vectors():
    """Return only E. coli expression vectors from Twist catalog."""
    from .data.parts_db import TWIST_ECOLI_EXPRESSION_VECTORS
    return {k: _serialize_vector(v) for k, v in TWIST_ECOLI_EXPRESSION_VECTORS.items()}


@app.get("/api/vectors/catalog/mammalian")
async def get_mammalian_catalog_vectors():
    """Return mammalian expression vectors from Twist catalog."""
    from .data.parts_db import TWIST_MAMMALIAN_VECTORS
    return {k: _serialize_vector(v) for k, v in TWIST_MAMMALIAN_VECTORS.items()}


@app.get("/api/vectors/catalog/lentiviral")
async def get_lentiviral_catalog_vectors():
    """Return lentiviral transfer vectors from Twist catalog."""
    from .data.parts_db import TWIST_LENTIVIRAL_VECTORS
    return {k: _serialize_vector(v) for k, v in TWIST_LENTIVIRAL_VECTORS.items()}


@app.get("/api/vectors/catalog/cloning")
async def get_cloning_catalog_vectors():
    """Return cloning/Gateway vectors from Twist catalog."""
    from .data.parts_db import TWIST_CLONING_VECTORS
    return {k: _serialize_vector(v) for k, v in TWIST_CLONING_VECTORS.items()}


@app.get("/api/vectors/categories")
async def get_vector_categories():
    """Return all vectors organized by category."""
    from .data.parts_db import TWIST_VECTOR_CATEGORIES
    return {
        cat_name: {k: _serialize_vector(v) for k, v in vectors.items()}
        for cat_name, vectors in TWIST_VECTOR_CATEGORIES.items()
    }


# Serve the React frontend
@app.get("/", response_class=HTMLResponse)
async def serve_frontend():
    """Serve the single-page React app."""
    frontend_path = Path(__file__).parent.parent.parent / "frontend" / "index.html"
    if frontend_path.exists():
        return HTMLResponse(content=frontend_path.read_text())
    return HTMLResponse(content="<h1>Frontend not found. Place index.html in frontend/</h1>")


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8421)

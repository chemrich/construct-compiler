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

class CompileRequest(BaseModel):
    """A construct spec submitted from the web UI."""
    spec: dict
    researcher_rate: float = 150.0
    overhead_multiplier: float = 1.5
    skip_constraints: bool = False


class PartInfo(BaseModel):
    id: str
    name: str
    type: str
    length_bp: int | None = None
    sequence: str | None = None
    color: str = "#cccccc"


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


class CompileResponse(BaseModel):
    success: bool
    construct_name: str = ""
    host: str = ""
    insert_length_bp: int = 0
    parts: list[PartInfo] = []
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
        cost_params = CostParams(
            researcher_hourly_rate=req.researcher_rate,
            overhead_multiplier=req.overhead_multiplier,
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

        # Build parts list
        parts = []
        for part in graph.parts():
            from .core.parts import Backbone
            if isinstance(part, Backbone):
                continue
            parts.append(PartInfo(
                id=part.id,
                name=part.name,
                type=part.__class__.__name__,
                length_bp=len(part.sequence) if part.sequence else None,
                sequence=str(part.sequence)[:200] if part.sequence else None,
                color=PART_COLORS.get(part.__class__.__name__, "#cccccc"),
            ))

        # Build strategies list
        strategies = []
        for s in plan.strategies:
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
            insert_length_bp=graph.total_insert_length() or 0,
            parts=parts,
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
    from .data.parts_db import PURIFICATION_TAGS, SOLUBILITY_TAGS
    return {
        "purification": list(PURIFICATION_TAGS.keys()),
        "solubility": list(SOLUBILITY_TAGS.keys()),
    }


@app.get("/api/parts/cleavage")
async def get_cleavage():
    from .data.parts_db import CLEAVAGE_SITES
    return {k: {"recognition": v["protein_sequence"]} for k, v in CLEAVAGE_SITES.items()}


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
    uvicorn.run(app, host="127.0.0.1", port=8421)

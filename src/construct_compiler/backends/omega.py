"""
Omegamega integration backend.

Bridges a compiled ConstructGraph into omegamega's oligopool design pipeline,
producing a vendor-ready oligo order, pool fidelity scores, and cost estimate.

Requires a local clone of chemrich/omegamega:
    git clone https://github.com/chemrich/omegamega

Point to it via the OMEGAMEGA_DIR environment variable, or pass
``omegamega_dir`` explicitly to ``run_omega()``.

omegamega is a GPL-3.0 scripts project (not an installable package), so it
is invoked via subprocess using its own Python interpreter / uv environment.
"""

from __future__ import annotations

import csv
import os
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import yaml

from ..core.graph import ConstructGraph


# ---------------------------------------------------------------------------
# Result type
# ---------------------------------------------------------------------------

@dataclass
class OmegaResult:
    """Summary of a completed omegamega oligopool design run."""
    oligo_count: int
    pool_count: int
    min_fidelity: float
    avg_fidelity: float
    oligo_cost_usd: float
    output_dir: Path

    def summary(self) -> str:
        lines = [
            f"Oligopool design complete",
            f"  Oligos:        {self.oligo_count}",
            f"  Pools:         {self.pool_count}",
            f"  Min fidelity:  {self.min_fidelity:.3f}",
            f"  Avg fidelity:  {self.avg_fidelity:.3f}",
            f"  Oligo cost:   ${self.oligo_cost_usd:.2f}",
            f"  Output dir:    {self.output_dir}",
        ]
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# FASTA export
# ---------------------------------------------------------------------------

def to_fasta(graph: ConstructGraph, label: Optional[str] = None) -> str:
    """Return the codon-optimized insert sequence as a FASTA string.

    Requires the graph to have been compiled through at least pass 2
    (reverse_translate). Pass 3 (constraint_resolution) is recommended so the
    sequence is codon-optimised before fragmentation.
    """
    seq = graph.full_insert_sequence()
    if seq is None:
        raise ValueError(
            "Graph has no concrete sequences — compile through reverse_translate "
            "(pass 2) or resolve_constraints (pass 3) before calling to_fasta()."
        )
    name = label or getattr(graph, "name", None) or "insert"
    # Sanitise the name for FASTA (no spaces)
    name = name.replace(" ", "_")
    return f">{name}\n{seq}\n"


# ---------------------------------------------------------------------------
# Omegamega runner
# ---------------------------------------------------------------------------

def _locate_omegamega(omegamega_dir: Optional[Path]) -> Path:
    """Return a validated path to the omegamega repo root."""
    if omegamega_dir is None:
        env = os.environ.get("OMEGAMEGA_DIR")
        if not env:
            raise RuntimeError(
                "Set the OMEGAMEGA_DIR environment variable to a local clone of "
                "chemrich/omegamega, or pass omegamega_dir= explicitly.\n"
                "  git clone https://github.com/chemrich/omegamega"
            )
        omegamega_dir = Path(env)

    omega_script = omegamega_dir / "code" / "omega.py"
    if not omega_script.exists():
        raise FileNotFoundError(
            f"omega.py not found at {omega_script}. "
            f"Check that OMEGAMEGA_DIR ({omegamega_dir}) points to the repo root."
        )
    return omegamega_dir


def _read_pool_stats(pool_stats_path: Path) -> tuple[float, float]:
    """Return (min_fidelity, avg_fidelity) from pool_stats.csv."""
    min_f = 1.0
    fidelities = []
    with open(pool_stats_path) as f:
        for row in csv.DictReader(f):
            v = float(row.get("min_fidelity", row.get("fidelity", 1.0)))
            fidelities.append(v)
            min_f = min(min_f, v)
    avg_f = sum(fidelities) / len(fidelities) if fidelities else 0.0
    return min_f, avg_f


def _read_oligo_cost(cost_summary_path: Path) -> float:
    """Return the total oligo synthesis cost from cost_summary.csv."""
    with open(cost_summary_path) as f:
        for row in csv.DictReader(f):
            for key in ("oligo_cost_usd", "total_synth_cost_usd", "oligo_cost"):
                if key in row:
                    return float(row[key])
    return 0.0


def run_omega(
    graph: ConstructGraph,
    output_dir: Path,
    *,
    enzyme: str = "BsaI",
    njunctions: int = 50,
    oligo_len: int = 350,
    upstream_bbsite: str = "AATG",
    downstream_bbsite: str = "TTAG",
    nopt_steps: int = 500,
    nopt_runs: int = 3,
    njobs: int = 1,
    label: Optional[str] = None,
    omegamega_dir: Optional[Path] = None,
) -> OmegaResult:
    """Run the omegamega oligopool design pipeline on a compiled construct.

    Writes the insert FASTA to a temp directory, generates an omegamega config,
    invokes ``code/omega.py genes`` via subprocess, then parses the output CSVs.

    Parameters
    ----------
    graph:
        A ConstructGraph compiled through at least pass 2 (reverse_translate).
    output_dir:
        Directory where omegamega writes its output CSVs.
    enzyme:
        Type IIS restriction enzyme for Golden Gate assembly (default: BsaI).
    njunctions:
        Number of Golden Gate sites per pool.
    oligo_len:
        Oligonucleotide length in bp (default: 350, matches Twist OPools).
    upstream_bbsite / downstream_bbsite:
        4-bp backbone junction overhangs flanking the insert.
    nopt_steps / nopt_runs:
        Simulated annealing parameters; reduce for faster tests.
    njobs:
        Parallel jobs for pool optimisation (joblib).
    label:
        FASTA sequence ID; defaults to graph.name.
    omegamega_dir:
        Path to the omegamega repo root. Falls back to OMEGAMEGA_DIR env var.
    """
    omegamega_dir = _locate_omegamega(omegamega_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # -- Write FASTA ----------------------------------------------------------
    fasta_str = to_fasta(graph, label)
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False, dir=output_dir
    ) as fasta_file:
        fasta_file.write(fasta_str)
        fasta_path = Path(fasta_file.name)

    primers_path = omegamega_dir / "data" / "subramanian_orthogonal.csv"

    # -- Write omegamega config -----------------------------------------------
    config = {
        "input_seqs": str(fasta_path),
        "output_dir": str(output_dir),
        "primers": str(primers_path),
        "upstream_bbsite": upstream_bbsite,
        "downstream_bbsite": downstream_bbsite,
        "enzyme": enzyme,
        "njunctions": njunctions,
        "oligo_len": oligo_len,
        "nopt_steps": nopt_steps,
        "nopt_runs": nopt_runs,
        "njobs": njobs,
        "add_primers": True,
        "pad_oligos": True,
        "ligation_data": "T4_18h_37C",
        "pricing_enabled": True,
        "twist_quote": False,
    }
    config_path = output_dir / "omega_config.yml"
    with open(config_path, "w") as f:
        yaml.dump(config, f)

    # -- Invoke omega.py via omegamega's own uv environment -------------------
    # omegamega is a `package = false` scripts project with its own deps
    # (jsonargparse, etc.) not present in construct-compiler's venv.
    omega_script = omegamega_dir / "code" / "omega.py"
    cmd = ["uv", "run", "python", str(omega_script), "genes", "--config", str(config_path)]

    proc = subprocess.run(
        cmd,
        cwd=str(omegamega_dir),
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"omegamega exited with code {proc.returncode}.\n"
            f"stdout:\n{proc.stdout}\n"
            f"stderr:\n{proc.stderr}"
        )

    # -- Parse outputs --------------------------------------------------------
    oligo_order_path = output_dir / "oligo_order.csv"
    pool_stats_path = output_dir / "pool_stats.csv"
    cost_summary_path = output_dir / "cost_summary.csv"

    for path in (oligo_order_path, pool_stats_path):
        if not path.exists():
            raise FileNotFoundError(
                f"Expected omegamega output not found: {path}\n"
                f"omega.py stdout:\n{proc.stdout}"
            )

    oligo_count = sum(1 for _ in open(oligo_order_path)) - 1  # minus header
    pool_count = sum(1 for _ in open(pool_stats_path)) - 1

    min_fidelity, avg_fidelity = _read_pool_stats(pool_stats_path)
    oligo_cost = _read_oligo_cost(cost_summary_path) if cost_summary_path.exists() else 0.0

    fasta_path.unlink(missing_ok=True)

    return OmegaResult(
        oligo_count=oligo_count,
        pool_count=pool_count,
        min_fidelity=min_fidelity,
        avg_fidelity=avg_fidelity,
        oligo_cost_usd=oligo_cost,
        output_dir=output_dir,
    )

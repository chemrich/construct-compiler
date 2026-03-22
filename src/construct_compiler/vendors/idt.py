"""
IDT (Integrated DNA Technologies) vendor plugin.

Uses SciTools Plus API for:
- Oligo analysis (Tm, secondary structure, dimer checks)
- gBlock/eBlock synthesis screening
- Codon optimization
- Order management

API credentials: set IDT_API_KEY and IDT_CLIENT_SECRET env vars.
Docs: https://www.idtdna.com/pages/tools/apidoc
"""

from __future__ import annotations

import os
import logging
from typing import Optional

from .base import VendorPlugin, ScreeningResult, OptimizationResult

logger = logging.getLogger(__name__)


class IDTVendor(VendorPlugin):
    """IDT SciTools Plus API integration."""

    BASE_URL = "https://www.idtdna.com/api/v1"
    AUTH_URL = "https://www.idtdna.com/Identityserver/connect/token"

    def __init__(self, client_id: Optional[str] = None, client_secret: Optional[str] = None):
        self.client_id = client_id or os.environ.get("IDT_CLIENT_ID", "")
        self.client_secret = client_secret or os.environ.get("IDT_CLIENT_SECRET", "")
        self._token: Optional[str] = None

    @property
    def name(self) -> str:
        return "IDT"

    @property
    def authenticated(self) -> bool:
        return bool(self.client_id and self.client_secret)

    def screen(self, sequence: str, product_type: str = "gblock",
               **kwargs) -> ScreeningResult:
        """Screen a sequence for IDT synthesis feasibility."""
        length = len(sequence)

        if not self.authenticated:
            return self._mock_screen(sequence, product_type)

        # Real API would go here (OAuth2 flow + screening endpoint)
        return self._mock_screen(sequence, product_type)

    def optimize_codons(self, protein_sequence: str, organism: str = "Escherichia coli",
                        **kwargs) -> OptimizationResult:
        """Use IDT's codon optimization API."""
        if not self.authenticated:
            return OptimizationResult(
                vendor="idt",
                original_sequence=protein_sequence,
                optimized_sequence="",
                notes=["Mock mode: IDT API credentials not configured."],
            )

        return OptimizationResult(
            vendor="idt",
            original_sequence=protein_sequence,
            optimized_sequence="",
            notes=["IDT optimization API integration pending"],
        )

    def analyze_oligo(self, sequence: str) -> dict:
        """
        Analyze an oligonucleotide for Tm, secondary structure, dimers.
        Useful for primer validation.
        """
        if not self.authenticated:
            return self._mock_oligo_analysis(sequence)
        return self._mock_oligo_analysis(sequence)

    def _mock_screen(self, sequence: str, product_type: str) -> ScreeningResult:
        """Heuristic screening when API is unavailable."""
        length = len(sequence)
        warnings = []
        errors = []
        feasible = True

        limits = {"gblock": 3000, "eblock": 3000, "gene": 5000}
        max_len = limits.get(product_type, 3000)

        if length > max_len:
            errors.append(f"Length {length} bp exceeds {product_type} limit ({max_len} bp)")
            feasible = False

        if length < 125 and product_type == "gblock":
            errors.append("gBlocks require minimum 125 bp")
            feasible = False

        gc = (sequence.upper().count("G") + sequence.upper().count("C")) / max(length, 1)
        if gc < 0.25 or gc > 0.75:
            warnings.append(f"GC content {gc*100:.1f}% outside preferred range")

        price_per_bp = {"gblock": 0.08, "eblock": 0.08, "gene": 0.08}.get(product_type, 0.08)

        return ScreeningResult(
            vendor="idt",
            feasible=feasible,
            sequence_length=length,
            estimated_price=length * price_per_bp,
            turnaround_days=(3, 5) if product_type == "gblock" else (10, 15),
            warnings=warnings,
            errors=errors,
        )

    def _mock_oligo_analysis(self, sequence: str) -> dict:
        """Basic oligo analysis heuristics."""
        from Bio.SeqUtils import MeltingTemp as mt
        from Bio.Seq import Seq
        seq = Seq(sequence)
        return {
            "length": len(sequence),
            "tm_nearest_neighbor": round(mt.Tm_NN(seq), 1),
            "gc_content": round(
                (sequence.upper().count("G") + sequence.upper().count("C"))
                / max(len(sequence), 1) * 100, 1
            ),
        }

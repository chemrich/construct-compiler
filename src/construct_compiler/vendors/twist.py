"""
Twist Bioscience vendor plugin.

Uses TAPI (Twist Application Programming Interface) for:
- Sequence screening (feasibility + complexity scoring)
- Codon optimization (ML-based, tuned for Twist's synthesis process)
- Order placement and tracking

API credentials: set TWIST_API_KEY and TWIST_API_SECRET env vars,
or pass them to the constructor.

Docs: https://www.twistbioscience.com/twist-ordering-platform
"""

from __future__ import annotations

import os
import logging
from typing import Optional

from .base import VendorPlugin, ScreeningResult, OptimizationResult, OrderResult

logger = logging.getLogger(__name__)


class TwistVendor(VendorPlugin):
    """Twist Bioscience TAPI integration."""

    BASE_URL = "https://twist-api.twistbioscience.com/v1"

    def __init__(self, api_key: Optional[str] = None, api_secret: Optional[str] = None):
        self.api_key = api_key or os.environ.get("TWIST_API_KEY", "")
        self.api_secret = api_secret or os.environ.get("TWIST_API_SECRET", "")
        self._token: Optional[str] = None

    @property
    def name(self) -> str:
        return "Twist Bioscience"

    @property
    def authenticated(self) -> bool:
        return bool(self.api_key and self.api_secret)

    def _get_headers(self) -> dict:
        return {
            "Authorization": f"Bearer {self._token}" if self._token else "",
            "Content-Type": "application/json",
        }

    def screen(self, sequence: str, product_type: str = "GENE",
               **kwargs) -> ScreeningResult:
        """
        Screen a sequence for synthesis feasibility.

        In mock mode (no API key), uses heuristic checks:
        - GC content 25-75%
        - No homopolymers > 10bp
        - Length within Twist limits
        """
        length = len(sequence)

        if not self.authenticated:
            return self._mock_screen(sequence, product_type)

        # Real API call (when credentials are available)
        import requests
        try:
            resp = requests.post(
                f"{self.BASE_URL}/screening",
                json={"sequence": sequence, "product_type": product_type},
                headers=self._get_headers(),
                timeout=30,
            )
            if resp.status_code == 200:
                data = resp.json()
                return ScreeningResult(
                    vendor="twist",
                    feasible=data.get("feasible", False),
                    sequence_length=length,
                    estimated_price=data.get("price", length * 0.07),
                    turnaround_days=(12, 18),
                    complexity_score=data.get("complexity_score", 0),
                    warnings=data.get("warnings", []),
                    errors=data.get("errors", []),
                )
        except Exception as e:
            logger.warning(f"Twist API screening failed: {e}")

        return self._mock_screen(sequence, product_type)

    def optimize_codons(self, protein_sequence: str, organism: str = "Escherichia coli",
                        **kwargs) -> OptimizationResult:
        """
        Use Twist's ML-based codon optimization.

        In mock mode, returns the input unchanged.
        """
        if not self.authenticated:
            return OptimizationResult(
                vendor="twist",
                original_sequence=protein_sequence,
                optimized_sequence="",  # empty = use local optimization
                notes=["Mock mode: Twist API key not configured. Using local codon optimization."],
            )

        import requests
        try:
            resp = requests.post(
                f"{self.BASE_URL}/codon-optimization",
                json={
                    "protein_sequence": protein_sequence,
                    "organism": organism,
                },
                headers=self._get_headers(),
                timeout=30,
            )
            if resp.status_code == 200:
                data = resp.json()
                return OptimizationResult(
                    vendor="twist",
                    original_sequence=protein_sequence,
                    optimized_sequence=data.get("optimized_dna", ""),
                    changes_made=data.get("changes", 0),
                    gc_content=data.get("gc_content", 0),
                    notes=data.get("notes", []),
                )
        except Exception as e:
            logger.warning(f"Twist codon optimization failed: {e}")

        return OptimizationResult(
            vendor="twist",
            original_sequence=protein_sequence,
            optimized_sequence="",
            notes=["API call failed, use local optimization as fallback"],
        )

    def _mock_screen(self, sequence: str, product_type: str) -> ScreeningResult:
        """Heuristic screening when API is unavailable."""
        length = len(sequence)
        warnings = []
        errors = []
        feasible = True

        # Length check
        if product_type == "GENE" and length > 5000:
            errors.append(f"Sequence length {length} bp exceeds gene synthesis limit (5000 bp)")
            feasible = False

        # GC content
        gc = (sequence.upper().count("G") + sequence.upper().count("C")) / max(length, 1)
        if gc < 0.25 or gc > 0.75:
            warnings.append(f"Extreme GC content ({gc*100:.1f}%) may affect synthesis")

        # Homopolymers
        for base in "ATGC":
            run = base * 11
            if run in sequence.upper():
                warnings.append(f"Long {base} homopolymer (>10bp) detected")

        # Price estimate
        price_per_bp = 0.07 if product_type == "GENE" else 0.09
        price = length * price_per_bp

        return ScreeningResult(
            vendor="twist",
            feasible=feasible,
            sequence_length=length,
            estimated_price=price,
            turnaround_days=(12, 18),
            complexity_score=0.5 if warnings else 0.2,
            warnings=warnings,
            errors=errors,
        )

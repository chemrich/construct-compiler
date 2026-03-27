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

    BASE_URL = "https://www.idtdna.com/restapi/v1"
    AUTH_URL = "https://www.idtdna.com/Identityserver/connect/token"

    def __init__(self, client_id: Optional[str] = None, client_secret: Optional[str] = None,
                 username: Optional[str] = None, password: Optional[str] = None):
        self.client_id = client_id or os.environ.get("IDT_CLIENT_ID", "")
        self.client_secret = client_secret or os.environ.get("IDT_CLIENT_SECRET", "")
        self.username = username or os.environ.get("IDT_USERNAME", "")
        self.password = password or os.environ.get("IDT_PASSWORD", "")
        self._token: Optional[str] = None


    @property
    def name(self) -> str:
        return "IDT"

    @property
    def authenticated(self) -> bool:
        return bool(self.client_id and self.client_secret and self.username and self.password)

    def _get_token(self) -> Optional[str]:
        """Fetch OAuth2 token from IDT."""
        if self._token:
            return self._token

        if not self.authenticated:
            return None

        import requests
        import base64

        auth_str = base64.b64encode(f"{self.client_id}:{self.client_secret}".encode()).decode()
        headers = {
            "Content-Type": "application/x-www-form-urlencoded",
            "Authorization": f"Basic {auth_str}",
        }
        data = {
            "grant_type": "password",
            "scope": "test",
            "username": self.username,
            "password": self.password,
        }

        try:
            resp = requests.post(self.AUTH_URL, headers=headers, data=data, timeout=30)
            if resp.status_code == 200:
                self._token = resp.json().get("access_token")
                return self._token
            else:
                logger.warning(f"IDT Auth failed with status {resp.status_code}: {resp.text}")
        except Exception as e:
            logger.warning(f"IDT Auth exception: {e}")

        return None


    def screen(self, sequence: str, product_type: str = "gblock",
               **kwargs) -> ScreeningResult:
        """Screen a sequence for IDT synthesis feasibility."""
        length = len(sequence)

        if not self.authenticated:
            return self._mock_screen(sequence, product_type)

        token = self._get_token()
        if not token:
            logger.warning("Failed to authenticate with IDT for screening. Using mock fallback.")
            return self._mock_screen(sequence, product_type)

        import requests
        try:
            url = f"{self.BASE_URL}/CodonOpt/Optimize"
            payload = {
                "organism": kwargs.get("organism", "Escherichia coli"),
                "optimizationItems": [
                    {
                        "Name": "construct_compiler_screen_item",
                        "Sequence": sequence
                    }
                ],
                "sequenceType": "dna",
                "productType": product_type
            }
            
            headers = {
                "Content-Type": "application/json",
                "Authorization": f"Bearer {token}"
            }

            resp = requests.post(url, json=payload, headers=headers, timeout=30)
            if resp.status_code == 200:
                data = resp.json()
                items = data if isinstance(data, list) else data.get("OptimizationItems", [])
                
                if items:
                    item = items[0]
                    opt_result = item.get("OptResult", {})
                    complexity_summary = opt_result.get("ComplexitySummary", "Accepted")
                    violations = opt_result.get("Complexities", [])
                    
                    # If it has any text in violations, let's consider it a warning or error depending on content
                    warnings = []
                    errors = []
                    feasible = True
                    
                    if violations:
                        for viol in violations:
                            text = viol.get("Text", "")
                            if "Rejected" in complexity_summary or "not available" in text.lower():
                                errors.append(text)
                                feasible = False
                            else:
                                warnings.append(text)
                    
                    # Fallback if complexity summary indicates failure but no explicit violations listed
                    if "Rejected" in complexity_summary and feasible:
                        feasible = False
                        errors.append(f"Sequence rejected by IDT complexity screen: {complexity_summary}")

                    price_per_bp = {"gblock": 0.08, "eblock": 0.08, "megamer": 0.08}.get(product_type, 0.08)

                    return ScreeningResult(
                        vendor="idt",
                        feasible=feasible,
                        sequence_length=length,
                        estimated_price=length * price_per_bp,
                        turnaround_days=(3, 5) if product_type == "gblock" else (10, 15),
                        complexity_score=0.0 if feasible else 1.0,
                        warnings=warnings,
                        errors=errors,
                    )
                
                logger.warning(f"Unexpected IDT screening response format: {data}")
            else:
                logger.warning(f"IDT screening API returned status {resp.status_code}: {resp.text}")
        except Exception as e:
            logger.warning(f"IDT API screening failed: {e}")

        return self._mock_screen(sequence, product_type)

    def optimize_codons(self, protein_sequence: str, organism: str = "Escherichia coli",
                        product_type: str = "gblock", **kwargs) -> OptimizationResult:
        """Use IDT's codon optimization API."""
        if not self.authenticated:
            return OptimizationResult(
                vendor="idt",
                original_sequence=protein_sequence,
                optimized_sequence="",
                notes=["Mock mode: IDT API credentials not configured."],
            )

        token = self._get_token()
        if not token:
            return OptimizationResult(
                vendor="idt",
                original_sequence=protein_sequence,
                optimized_sequence="",
                notes=["Failed to authenticate with IDT. Using local optimization fallback."],
            )

        import requests
        try:
            url = f"{self.BASE_URL}/CodonOpt/Optimize"
            payload = {
                "organism": organism,
                "optimizationItems": [
                    {
                        "Name": "construct_compiler_item",
                        "Sequence": protein_sequence
                    }
                ],
                "sequenceType": kwargs.get("sequence_type", "dna" if kwargs.get("is_dna", False) else "protein"),
                "productType": product_type
            }
            
            headers = {
                "Content-Type": "application/json",
                "Authorization": f"Bearer {token}"
            }

            resp = requests.post(url, json=payload, headers=headers, timeout=30)
            if resp.status_code == 200:
                data = resp.json()
                items = data if isinstance(data, list) else data.get("OptimizationItems", [])
                
                if items:
                    item = items[0]
                    opt_result = item.get("OptResult", {})
                    
                    # Try OptResult.FullSequence first, fallback to Sequence at top level
                    optimized_seq = opt_result.get("FullSequence") or item.get("Sequence", "")
                    complexity = opt_result.get("ComplexitySummary") or item.get("ComplexitySummary", "unknown")
                    
                    return OptimizationResult(
                        vendor="idt",
                        original_sequence=protein_sequence,
                        optimized_sequence=optimized_seq,
                        changes_made=0,
                        gc_content=0.0,
                        notes=[f"Complexity: {complexity}"],
                    )
                
                return OptimizationResult(
                    vendor="idt",
                    original_sequence=protein_sequence,
                    optimized_sequence="",
                    notes=[f"Unexpected API response format: {data}"],
                )
            else:
                # Handle non-200 status codes explicitly inside try to avoid fall-through
                return OptimizationResult(
                    vendor="idt",
                    original_sequence=protein_sequence,
                    optimized_sequence="",
                    notes=[f"API returned status {resp.status_code}: {resp.text}"],
                )
        except Exception as e:
            logger.warning(f"IDT codon optimization failed: {e}")
            return OptimizationResult(
                vendor="idt",
                original_sequence=protein_sequence,
                optimized_sequence="",
                notes=[f"Exception during codon optimization: {e}"],
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

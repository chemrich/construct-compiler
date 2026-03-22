"""
Vendor plugin interface — base class for DNA synthesis vendor integrations.

Each vendor implements screen(), optimize(), and order().
API credentials are read from environment variables or a config file.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class ScreeningResult:
    """Result of vendor synthesis feasibility screening."""
    vendor: str
    feasible: bool
    sequence_length: int
    estimated_price: float = 0.0
    turnaround_days: tuple[int, int] = (0, 0)
    complexity_score: float = 0.0   # 0-1, vendor-specific
    warnings: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)


@dataclass
class OptimizationResult:
    """Result of vendor codon optimization."""
    vendor: str
    original_sequence: str
    optimized_sequence: str
    changes_made: int = 0
    gc_content: float = 0.0
    notes: list[str] = field(default_factory=list)


@dataclass
class OrderResult:
    """Result of placing a synthesis order."""
    vendor: str
    order_id: str
    status: str = "submitted"
    estimated_ship_date: str = ""
    total_price: float = 0.0


class VendorPlugin(ABC):
    """Base class for vendor integrations."""

    @property
    @abstractmethod
    def name(self) -> str:
        ...

    @abstractmethod
    def screen(self, sequence: str, **kwargs) -> ScreeningResult:
        """Check synthesis feasibility, return price and turnaround estimate."""
        ...

    @abstractmethod
    def optimize_codons(self, protein_sequence: str, organism: str,
                        **kwargs) -> OptimizationResult:
        """Use vendor's codon optimization engine."""
        ...

    def order(self, sequences: list[dict], **kwargs) -> OrderResult:
        """Place a synthesis order. Override in vendor-specific plugins."""
        raise NotImplementedError(f"{self.name} order API not yet implemented")

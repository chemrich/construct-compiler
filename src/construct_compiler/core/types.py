"""
Core type system for the genetic construct compiler.

Defines PortKind (the biological "types" that flow between parts),
Port (typed connection points on parts), and Constraint (requirements
that must be satisfied during compilation).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Any, Optional


# ---------------------------------------------------------------------------
# Port kinds – the "types" that flow between genetic parts
# ---------------------------------------------------------------------------

class PortKind(Enum):
    """
    Biological signal types that flow through a construct.

    These form a simple type lattice:
        DNA_CONTEXT          -- raw dsDNA (backbone level)
        TRANSCRIPTION        -- region is being transcribed
        TRANSLATION_INIT     -- ribosome is positioned to start translation
        PEPTIDE_CHAIN        -- growing polypeptide (no stop yet)
        TRANSLATED_PRODUCT   -- finished polypeptide (stop codon reached)
        TRANSCRIPT           -- mRNA region downstream of CDS (before terminator)
    """
    DNA_CONTEXT = auto()
    TRANSCRIPTION = auto()
    TRANSLATION_INIT = auto()
    PEPTIDE_CHAIN = auto()
    TRANSLATED_PRODUCT = auto()
    TRANSCRIPT = auto()


# Compatibility matrix: which output port kinds can connect to which inputs.
# An edge (A -> B) is valid iff PORT_COMPATIBILITY[A] contains B.
PORT_COMPATIBILITY: dict[PortKind, set[PortKind]] = {
    PortKind.DNA_CONTEXT: {
        PortKind.DNA_CONTEXT,
        PortKind.TRANSCRIPTION,
    },
    PortKind.TRANSCRIPTION: {
        PortKind.TRANSCRIPTION,       # e.g. spacer within transcribed region
        PortKind.TRANSLATION_INIT,    # RBS consumes transcript, starts translation
        PortKind.TRANSCRIPT,          # terminator can follow promoter (unusual but valid)
    },
    PortKind.TRANSLATION_INIT: {
        PortKind.TRANSLATION_INIT,    # tag/CDS as first element after RBS
        PortKind.PEPTIDE_CHAIN,       # CDS without stop
        PortKind.TRANSLATED_PRODUCT,  # CDS with stop
    },
    PortKind.PEPTIDE_CHAIN: {
        PortKind.PEPTIDE_CHAIN,       # another tag/domain fused in-frame
        PortKind.TRANSLATION_INIT,    # CDS/tag in fusion chain (accepts TI input)
        PortKind.TRANSLATED_PRODUCT,  # stop codon terminates the chain
    },
    PortKind.TRANSLATED_PRODUCT: {
        PortKind.TRANSCRIPTION,       # next cistron's RBS (polycistronic)
        PortKind.TRANSLATION_INIT,    # coupled translation (BCD)
        PortKind.TRANSCRIPT,          # terminator ends the operon
    },
    PortKind.TRANSCRIPT: {
        PortKind.DNA_CONTEXT,         # back to backbone after terminator
        PortKind.TRANSCRIPTION,       # another transcription unit
    },
}


def ports_compatible(output_kind: PortKind, input_kind: PortKind) -> bool:
    """Check whether an output port can connect to an input port."""
    return input_kind in PORT_COMPATIBILITY.get(output_kind, set())


# ---------------------------------------------------------------------------
# Ports
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Port:
    """A typed connection point on a genetic part."""
    name: str
    kind: PortKind
    direction: str  # "in" or "out"

    def __post_init__(self):
        if self.direction not in ("in", "out"):
            raise ValueError(f"Port direction must be 'in' or 'out', got '{self.direction}'")


# ---------------------------------------------------------------------------
# Resolution state
# ---------------------------------------------------------------------------

class ResolutionState(Enum):
    """How concrete is this part's sequence?"""
    ABSTRACT = auto()   # Intent only ("a strong RBS for E. coli")
    RESOLVED = auto()   # Specific part chosen, protein seq known, DNA maybe not
    CONCRETE = auto()   # Final DNA sequence, codon-optimized, constraint-checked


# ---------------------------------------------------------------------------
# Constraints
# ---------------------------------------------------------------------------

class ConstraintKind(Enum):
    """Categories of constraints on parts or junctions."""
    NO_RESTRICTION_SITE = auto()    # e.g. no BsaI sites within this part
    GC_CONTENT = auto()             # GC% must be in [lo, hi]
    MAX_HOMOPOLYMER = auto()        # no runs of >N identical bases
    CODON_FREQUENCY = auto()        # no codons below threshold frequency
    PRESERVE_PROTEIN = auto()       # DNA changes must be synonymous
    SYNTHESIS_COMPATIBLE = auto()   # must pass vendor screening
    CUSTOM = auto()


@dataclass
class Constraint:
    """A requirement that must be satisfied during compilation."""
    kind: ConstraintKind
    params: dict[str, Any] = field(default_factory=dict)
    scope: str = "local"  # "local" (this part only) or "global" (entire construct)
    description: str = ""

    def __repr__(self):
        return f"Constraint({self.kind.name}, {self.params})"


# ---------------------------------------------------------------------------
# Convenience constraint constructors
# ---------------------------------------------------------------------------

def no_restriction_site(enzyme: str) -> Constraint:
    """Construct must not contain recognition site for the given enzyme."""
    return Constraint(
        kind=ConstraintKind.NO_RESTRICTION_SITE,
        params={"enzyme": enzyme},
        scope="global",
        description=f"No internal {enzyme} recognition sites",
    )


def gc_content(lo: float = 0.35, hi: float = 0.65) -> Constraint:
    return Constraint(
        kind=ConstraintKind.GC_CONTENT,
        params={"lo": lo, "hi": hi},
        scope="local",
        description=f"GC content between {lo*100:.0f}% and {hi*100:.0f}%",
    )


def max_homopolymer(n: int = 6) -> Constraint:
    return Constraint(
        kind=ConstraintKind.MAX_HOMOPOLYMER,
        params={"max_run": n},
        scope="local",
        description=f"No homopolymer runs longer than {n} bp",
    )


def preserve_protein() -> Constraint:
    return Constraint(
        kind=ConstraintKind.PRESERVE_PROTEIN,
        params={},
        scope="local",
        description="All nucleotide changes must be synonymous",
    )

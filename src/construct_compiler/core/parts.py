"""
Genetic part definitions — the node types in the IR graph.

Each part is a dataclass with typed ports, resolution state, and
optional constraints. Parts progress from ABSTRACT -> RESOLVED -> CONCRETE
through compiler passes.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional
from enum import Enum, auto

from Bio.Seq import Seq

from .types import (
    Port, PortKind, ResolutionState, Constraint,
    no_restriction_site, preserve_protein,
)


# ---------------------------------------------------------------------------
# Base part
# ---------------------------------------------------------------------------

@dataclass
class GeneticPart:
    """Base class for all IR nodes."""
    id: str
    name: str = ""
    resolution: ResolutionState = ResolutionState.ABSTRACT
    sequence: Optional[Seq] = None          # DNA sequence (None until CONCRETE)
    constraints: list[Constraint] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def ports(self) -> list[Port]:
        raise NotImplementedError

    @property
    def input_port(self) -> Optional[Port]:
        return next((p for p in self.ports if p.direction == "in"), None)

    @property
    def output_port(self) -> Optional[Port]:
        return next((p for p in self.ports if p.direction == "out"), None)

    @property
    def length(self) -> Optional[int]:
        return len(self.sequence) if self.sequence else None

    def __repr__(self):
        state = self.resolution.name[0]  # A/R/C
        seq_len = f" {len(self.sequence)}bp" if self.sequence else ""
        return f"{self.__class__.__name__}({self.id!r}, {state}{seq_len})"


# ---------------------------------------------------------------------------
# Backbone / vector
# ---------------------------------------------------------------------------

class Origin(Enum):
    PBR322 = "pBR322"
    P15A = "p15A"
    COLA = "ColA"
    RSF1030 = "RSF1030"
    COLDF = "CloDF13"
    SC101 = "SC101"
    PUC = "pUC"
    CUSTOM = "custom"


@dataclass
class Backbone(GeneticPart):
    """Plasmid backbone: origin, resistance, linearization site."""
    origin: Origin = Origin.PBR322
    resistance: str = "kanamycin"
    source: str = ""          # e.g. "addgene:26094" or "local"
    addgene_id: Optional[int] = None
    linearization_site: str = ""  # enzyme used to linearize for cloning

    # Catalog vector support
    catalog_id: Optional[str] = None      # e.g. "pET-28a(+)"
    catalog_vendor: str = ""              # e.g. "twist"
    provides: list[str] = field(default_factory=list)  # elements on vector: ["promoter", "rbs", "n_tag", ...]
    cloning_pair: Optional[tuple[str, str]] = None  # e.g. ("NdeI", "XhoI") — restriction sites defining insert window

    @property
    def is_catalog_vector(self) -> bool:
        return bool(self.catalog_id)

    @property
    def ports(self) -> list[Port]:
        # Backbone provides DNA context; insert goes between its ports
        return [
            Port("out", PortKind.DNA_CONTEXT, "out"),
            Port("in", PortKind.DNA_CONTEXT, "in"),
        ]


# ---------------------------------------------------------------------------
# Promoter
# ---------------------------------------------------------------------------

class RegulationType(Enum):
    CONSTITUTIVE = auto()
    INDUCIBLE = auto()
    REPRESSIBLE = auto()


@dataclass
class Promoter(GeneticPart):
    """Promoter element."""
    regulation: RegulationType = RegulationType.INDUCIBLE
    inducer: str = ""         # e.g. "IPTG", "arabinose", "aTc"
    strength: str = "strong"  # "strong", "medium", "weak" or numeric AU

    @property
    def ports(self) -> list[Port]:
        return [
            Port("in", PortKind.DNA_CONTEXT, "in"),
            Port("out", PortKind.TRANSCRIPTION, "out"),
        ]


# ---------------------------------------------------------------------------
# RBS / translation initiation
# ---------------------------------------------------------------------------

class RBSDesignMethod(Enum):
    LOOKUP = auto()          # use a pre-characterized RBS by name
    SALIS_CALCULATOR = auto()  # compute de novo
    BCD = auto()             # bicistronic design element


@dataclass
class RBS(GeneticPart):
    """Ribosome binding site / translation initiation element."""
    design_method: RBSDesignMethod = RBSDesignMethod.LOOKUP
    target_expression: str = "high"  # "high", "medium", "low" or numeric TIR
    part_name: str = ""              # e.g. "BBa_B0034", "BCD2"

    @property
    def ports(self) -> list[Port]:
        return [
            Port("in", PortKind.TRANSCRIPTION, "in"),
            Port("out", PortKind.TRANSLATION_INIT, "out"),
        ]


# ---------------------------------------------------------------------------
# Coding sequence
# ---------------------------------------------------------------------------

@dataclass
class CDS(GeneticPart):
    """Coding sequence — a gene or ORF."""
    protein_sequence: Optional[str] = None
    source_db: str = ""           # "uniprot", "fpbase", "ncbi", "local"
    source_id: str = ""           # accession or slug
    organism: str = "e_coli"      # target expression host
    has_stop: bool = True
    codon_optimization: str = "local"  # "local", "twist", "idt", "none"

    @property
    def ports(self) -> list[Port]:
        out_kind = (
            PortKind.TRANSLATED_PRODUCT if self.has_stop
            else PortKind.PEPTIDE_CHAIN
        )
        # CDS accepts TRANSLATION_INIT (after RBS or tag).
        # The compatibility matrix allows both TRANSLATION_INIT and
        # PEPTIDE_CHAIN to connect to TRANSLATION_INIT.
        return [
            Port("in", PortKind.TRANSLATION_INIT, "in"),
            Port("out", out_kind, "out"),
        ]


# ---------------------------------------------------------------------------
# Fusion tags and cleavage sites
# ---------------------------------------------------------------------------

class TagPosition(Enum):
    N_TERM = auto()
    C_TERM = auto()


@dataclass
class PurificationTag(GeneticPart):
    """Affinity/purification tag (6xHis, GST, Strep-tag II, etc.)."""
    tag_type: str = "6xHis"
    position: TagPosition = TagPosition.N_TERM

    @property
    def ports(self) -> list[Port]:
        # Tags can be first in the chain (right after RBS) or mid-chain.
        # We use TRANSLATION_INIT as input so they're valid after an RBS;
        # the compatibility matrix also allows PEPTIDE_CHAIN -> TRANSLATION_INIT.
        return [
            Port("in", PortKind.TRANSLATION_INIT, "in"),
            Port("out", PortKind.PEPTIDE_CHAIN, "out"),
        ]


@dataclass
class SolubilityTag(GeneticPart):
    """Solubility-enhancing fusion tag (MBP, SUMO, Trx, etc.)."""
    tag_type: str = "MBP"
    position: TagPosition = TagPosition.N_TERM

    @property
    def ports(self) -> list[Port]:
        return [
            Port("in", PortKind.PEPTIDE_CHAIN, "in"),
            Port("out", PortKind.PEPTIDE_CHAIN, "out"),
        ]


@dataclass
class CleavageSite(GeneticPart):
    """Protease cleavage site (TEV, 3C, Factor Xa, Thrombin, etc.)."""
    protease: str = "TEV"
    recognition_seq: str = ""   # amino acid sequence, e.g. "ENLYFQS" for TEV

    @property
    def ports(self) -> list[Port]:
        return [
            Port("in", PortKind.PEPTIDE_CHAIN, "in"),
            Port("out", PortKind.PEPTIDE_CHAIN, "out"),
        ]

    def __post_init__(self):
        # Default recognition sequences
        if not self.recognition_seq:
            defaults = {
                "TEV": "ENLYFQS",
                "3C": "LEVLFQGP",
                "Factor_Xa": "IEGR",
                "Thrombin": "LVPRGS",
                "SUMO": "",  # SUMO protease recognizes the SUMO fold, not a linear seq
            }
            self.recognition_seq = defaults.get(self.protease, "")
        # Cleavage sites MUST preserve their protein sequence
        self.constraints.append(preserve_protein())


@dataclass
class StopCodon(GeneticPart):
    """Explicit stop codon node — terminates a peptide chain."""
    codon: str = "TAA"  # TAA, TAG, or TGA

    @property
    def ports(self) -> list[Port]:
        return [
            Port("in", PortKind.PEPTIDE_CHAIN, "in"),
            Port("out", PortKind.TRANSLATED_PRODUCT, "out"),
        ]


# ---------------------------------------------------------------------------
# Linker peptides
# ---------------------------------------------------------------------------

@dataclass
class Linker(GeneticPart):
    """Peptide linker between fusion domains."""
    linker_type: str = "GS"      # "GS", "rigid", "flexible", "custom"
    repeats: int = 3             # e.g. (GGGGS)x3
    custom_sequence: str = ""    # amino acid sequence if custom

    @property
    def ports(self) -> list[Port]:
        return [
            Port("in", PortKind.PEPTIDE_CHAIN, "in"),
            Port("out", PortKind.PEPTIDE_CHAIN, "out"),
        ]


# ---------------------------------------------------------------------------
# Terminator
# ---------------------------------------------------------------------------

@dataclass
class Terminator(GeneticPart):
    """Transcription terminator."""
    terminator_type: str = "rho_independent"  # or "rho_dependent"

    @property
    def ports(self) -> list[Port]:
        return [
            Port("in", PortKind.TRANSCRIPT, "in"),
            Port("out", PortKind.DNA_CONTEXT, "out"),
        ]


# ---------------------------------------------------------------------------
# Spacer / intergenic region
# ---------------------------------------------------------------------------

@dataclass
class Spacer(GeneticPart):
    """Intergenic spacer between cistrons or other elements."""
    length_bp: int = 30
    avoid_structure: bool = True

    @property
    def ports(self) -> list[Port]:
        # Spacer is transparent — passes through whatever context it's in
        return [
            Port("in", PortKind.TRANSCRIPT, "in"),
            Port("out", PortKind.TRANSCRIPTION, "out"),
        ]

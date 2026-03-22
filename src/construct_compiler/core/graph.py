"""
ConstructGraph — the intermediate representation.

A directed graph of genetic parts connected by typed edges.
Provides validation, traversal, and manipulation methods.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Iterator
from Bio.Seq import Seq

from .types import PortKind, ResolutionState, Constraint, ports_compatible
from .parts import GeneticPart, Backbone


# ---------------------------------------------------------------------------
# Edges
# ---------------------------------------------------------------------------

@dataclass
class Edge:
    """A connection between two parts in the construct."""
    source_id: str
    target_id: str
    source_port: str = "out"
    target_port: str = "in"

    # Assembly-level annotation (filled in by assembly planner)
    is_junction: bool = False       # is this a physical assembly junction?
    overhang: Optional[str] = None  # 4bp overhang sequence (Golden Gate)
    junction_type: str = ""         # "golden_gate", "gibson", "synthesis_join"


# ---------------------------------------------------------------------------
# Validation result
# ---------------------------------------------------------------------------

@dataclass
class ValidationError:
    """A type-checking or structural error in the construct."""
    severity: str  # "error" or "warning"
    location: str  # part id or edge description
    message: str

    def __repr__(self):
        return f"[{self.severity.upper()}] {self.location}: {self.message}"


# ---------------------------------------------------------------------------
# Construct graph
# ---------------------------------------------------------------------------

@dataclass
class ConstructGraph:
    """
    The core IR: an ordered graph of genetic parts.

    For linear expression cassettes (the common case), this is essentially
    an ordered list with typed connections. The graph representation supports
    future extensions (branching, circular constructs, etc.) but day-to-day
    usage is mostly linear traversal.
    """
    name: str = ""
    description: str = ""
    host_organism: str = "e_coli"

    # Storage
    _parts: dict[str, GeneticPart] = field(default_factory=dict)
    _edges: list[Edge] = field(default_factory=list)
    _part_order: list[str] = field(default_factory=list)  # insertion order

    # Top-level construct constraints (applied globally)
    constraints: list[Constraint] = field(default_factory=list)

    # Metadata
    metadata: dict = field(default_factory=dict)

    # ------------------------------------------------------------------
    # Part management
    # ------------------------------------------------------------------

    def add_part(self, part: GeneticPart) -> ConstructGraph:
        """Add a part to the graph. Returns self for chaining."""
        if part.id in self._parts:
            raise ValueError(f"Part with id '{part.id}' already exists")
        self._parts[part.id] = part
        self._part_order.append(part.id)
        return self

    def get_part(self, part_id: str) -> GeneticPart:
        return self._parts[part_id]

    def parts(self) -> list[GeneticPart]:
        """Return parts in insertion order."""
        return [self._parts[pid] for pid in self._part_order]

    def parts_by_type(self, part_type: type) -> list[GeneticPart]:
        """Return all parts of a given type."""
        return [p for p in self.parts() if isinstance(p, part_type)]

    def replace_part(self, part_id: str, new_part: GeneticPart) -> None:
        """Replace a part in-place (same position in order)."""
        if part_id not in self._parts:
            raise KeyError(f"No part with id '{part_id}'")
        new_part.id = part_id
        self._parts[part_id] = new_part

    # ------------------------------------------------------------------
    # Edge management
    # ------------------------------------------------------------------

    def connect(self, source_id: str, target_id: str,
                source_port: str = "out", target_port: str = "in") -> ConstructGraph:
        """Connect two parts. Returns self for chaining."""
        if source_id not in self._parts:
            raise KeyError(f"Source part '{source_id}' not found")
        if target_id not in self._parts:
            raise KeyError(f"Target part '{target_id}' not found")
        self._edges.append(Edge(source_id, target_id, source_port, target_port))
        return self

    def connect_linear(self) -> ConstructGraph:
        """Connect all parts in insertion order (linear construct)."""
        self._edges.clear()
        for i in range(len(self._part_order) - 1):
            self.connect(self._part_order[i], self._part_order[i + 1])
        return self

    def edges(self) -> list[Edge]:
        return list(self._edges)

    def edges_from(self, part_id: str) -> list[Edge]:
        return [e for e in self._edges if e.source_id == part_id]

    def edges_to(self, part_id: str) -> list[Edge]:
        return [e for e in self._edges if e.target_id == part_id]

    # ------------------------------------------------------------------
    # Traversal
    # ------------------------------------------------------------------

    def walk(self) -> Iterator[tuple[GeneticPart, Optional[Edge]]]:
        """Walk the graph in order, yielding (part, incoming_edge) pairs."""
        for i, pid in enumerate(self._part_order):
            part = self._parts[pid]
            incoming = None
            if i > 0:
                # Find edge from previous part
                prev_id = self._part_order[i - 1]
                incoming = next(
                    (e for e in self._edges
                     if e.source_id == prev_id and e.target_id == pid),
                    None
                )
            yield part, incoming

    def cistrons(self) -> list[list[GeneticPart]]:
        """
        Identify cistrons — groups of parts that form a single
        translation unit (RBS through stop codon).
        """
        from .parts import RBS, CDS, StopCodon, PurificationTag, SolubilityTag, CleavageSite, Linker

        cistron_types = (RBS, CDS, StopCodon, PurificationTag, SolubilityTag, CleavageSite, Linker)
        result = []
        current_cistron: list[GeneticPart] = []

        for part in self.parts():
            if isinstance(part, RBS):
                if current_cistron:
                    result.append(current_cistron)
                current_cistron = [part]
            elif isinstance(part, cistron_types) and current_cistron:
                current_cistron.append(part)
            else:
                if current_cistron:
                    result.append(current_cistron)
                    current_cistron = []

        if current_cistron:
            result.append(current_cistron)

        return result

    # ------------------------------------------------------------------
    # Type checking / validation
    # ------------------------------------------------------------------

    def validate(self) -> list[ValidationError]:
        """
        Type-check the construct graph.
        Returns a list of errors and warnings.
        """
        errors: list[ValidationError] = []

        # 1. Check that all edges connect compatible ports
        for edge in self._edges:
            source = self._parts.get(edge.source_id)
            target = self._parts.get(edge.target_id)
            if not source or not target:
                errors.append(ValidationError(
                    "error", f"{edge.source_id}->{edge.target_id}",
                    "Edge references non-existent part"
                ))
                continue

            # Find the actual port objects
            src_port = next(
                (p for p in source.ports if p.name == edge.source_port), None
            )
            tgt_port = next(
                (p for p in target.ports if p.name == edge.target_port), None
            )

            if not src_port:
                errors.append(ValidationError(
                    "error", source.id,
                    f"No port named '{edge.source_port}'"
                ))
                continue
            if not tgt_port:
                errors.append(ValidationError(
                    "error", target.id,
                    f"No port named '{edge.target_port}'"
                ))
                continue

            if not ports_compatible(src_port.kind, tgt_port.kind):
                errors.append(ValidationError(
                    "error",
                    f"{source.id} -> {target.id}",
                    f"Incompatible ports: {source.__class__.__name__}.{src_port.name} "
                    f"({src_port.kind.name}) cannot connect to "
                    f"{target.__class__.__name__}.{tgt_port.name} ({tgt_port.kind.name})"
                ))

        # 2. Check that every CDS in a fusion chain (has_stop=False) is
        #    eventually followed by a stop or a CDS with has_stop=True
        from .parts import CDS
        for i, pid in enumerate(self._part_order):
            part = self._parts[pid]
            if isinstance(part, CDS) and not part.has_stop:
                # Walk forward to find a stop
                found_stop = False
                for j in range(i + 1, len(self._part_order)):
                    downstream = self._parts[self._part_order[j]]
                    if isinstance(downstream, CDS) and downstream.has_stop:
                        found_stop = True
                        break
                    from .parts import StopCodon
                    if isinstance(downstream, StopCodon):
                        found_stop = True
                        break
                if not found_stop:
                    errors.append(ValidationError(
                        "warning", part.id,
                        "CDS without stop codon has no downstream stop — "
                        "peptide chain is never terminated"
                    ))

        # 3. Warn if no backbone
        if not any(isinstance(p, Backbone) for p in self.parts()):
            errors.append(ValidationError(
                "warning", "construct",
                "No backbone defined — construct has no replication origin"
            ))

        return errors

    # ------------------------------------------------------------------
    # Sequence assembly (post-compilation)
    # ------------------------------------------------------------------

    def full_insert_sequence(self) -> Optional[Seq]:
        """
        Concatenate all concrete part sequences in order.
        Returns None if any part is not yet concrete.
        """
        seqs = []
        for part in self.parts():
            if isinstance(part, Backbone):
                continue  # backbone is separate
            if part.sequence is None:
                return None
            seqs.append(part.sequence)
        return sum(seqs, Seq("")) if seqs else None

    def total_insert_length(self) -> Optional[int]:
        seq = self.full_insert_sequence()
        return len(seq) if seq else None

    # ------------------------------------------------------------------
    # Display
    # ------------------------------------------------------------------

    def summary(self) -> str:
        """Human-readable summary of the construct."""
        lines = [f"Construct: {self.name}", f"Host: {self.host_organism}", ""]
        lines.append("Parts:")
        for part, edge in self.walk():
            prefix = "  "
            if edge and edge.is_junction:
                lines.append(f"  --- junction ({edge.overhang or '?'}) ---")
            lines.append(f"{prefix}{part}")
        lines.append("")

        total = self.total_insert_length()
        if total:
            lines.append(f"Total insert: {total} bp")

        errors = self.validate()
        if errors:
            lines.append("")
            lines.append("Validation:")
            for e in errors:
                lines.append(f"  {e}")

        return "\n".join(lines)

    def __repr__(self):
        n = len(self._parts)
        return f"ConstructGraph({self.name!r}, {n} parts)"

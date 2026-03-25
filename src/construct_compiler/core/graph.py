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
    # Assembled view (insert + backbone features)
    # ------------------------------------------------------------------

    def assembled_graph(self) -> ConstructGraph:
        """
        Build a new graph representing the full expression cassette as
        assembled in the backbone plasmid — including backbone-provided
        elements (promoter, RBS, tags, cleavage sites, terminator) with
        their actual DNA sequences from the GenBank annotations.

        This is the view that validators should check against, since it
        represents what the cell actually sees after cloning.

        If no catalog vector is used, returns a copy of the current graph.
        """
        from .parts import (
            Backbone, Promoter, RBS, CDS, PurificationTag, SolubilityTag,
            CleavageSite, Terminator, RBSDesignMethod, RegulationType,
            TagPosition,
        )

        # Find the backbone
        backbone = None
        for part in self.parts():
            if isinstance(part, Backbone) and part.is_catalog_vector:
                backbone = part
                break

        if not backbone:
            # No catalog vector — return a shallow copy of the current graph
            return self._copy_graph()

        # Determine which cloning pair to use:
        # 1. Explicit cloning_pair on the backbone (from spec or user)
        # 2. Auto-infer from insert features vs available pairs
        # 3. Default (no pair → use all upstream/downstream features)
        from ..data.vector_features import get_insert_context
        cloning_pair = backbone.cloning_pair

        if cloning_pair is None:
            # Auto-infer: try to pick the best cloning pair based on
            # what the insert declares vs what each pair provides
            cloning_pair = self._infer_cloning_pair(backbone)

        context = get_insert_context(backbone.catalog_id, cloning_pair=cloning_pair)
        if context is None:
            return self._copy_graph()

        upstream = context["upstream"]
        downstream = context["downstream"]

        # Build a new graph with backbone features + insert parts
        assembled = ConstructGraph(
            name=self.name,
            description=self.description,
            host_organism=self.host_organism,
            constraints=list(self.constraints),
            metadata=dict(self.metadata),
        )

        # Track IDs to avoid collisions
        _counter = [0]

        def _next_id(prefix: str) -> str:
            _counter[0] += 1
            return f"bb_{prefix}_{_counter[0]}"

        # 1. Add the backbone (as metadata, not in the expression cassette)
        assembled.add_part(Backbone(
            id=backbone.id,
            name=backbone.name,
            origin=backbone.origin,
            resistance=backbone.resistance,
            source=backbone.source,
            catalog_id=backbone.catalog_id,
            provides=backbone.provides,
        ))

        # 2. Add upstream backbone features, including inter-feature DNA.
        #    On physical vectors like pET-28b(+), the start codon (ATG)
        #    and initial codons sit in the gap between the annotated RBS
        #    and the first CDS feature (e.g., pET-28b has ATGGGCAGCAGC
        #    encoding M-G-S-S before the His-tag CATCAT...).
        #    We insert this gap as a small CDS "leader" part so the
        #    assembled view has a proper start codon on the first
        #    coding element.
        vector = context["vector"]
        for i, feat in enumerate(upstream):
            part = self._backbone_feature_to_part(feat, _next_id)
            if part is None:
                continue
            assembled.add_part(part)

            # After the RBS, check if there's a coding gap before the
            # next annotated feature. If so, insert a leader CDS part.
            if (feat.feature_type == "rbs"
                    and i + 1 < len(upstream)
                    and upstream[i + 1].feature_type in ("tag", "cleavage_site", "solubility_tag", "leader")
                    and vector.full_sequence):
                next_feat = upstream[i + 1]
                gap_start = feat.end        # 1-based end = 0-based start of gap
                gap_end = next_feat.start - 1  # 1-based start - 1 = 0-based end
                if gap_end > gap_start:
                    gap_seq = vector.full_sequence[gap_start:gap_end].upper()
                    # Only include if it contains a start codon.
                    # Trim to start at the ATG (the 5' UTR before it
                    # is part of the RBS context, not coding).
                    atg_idx = gap_seq.find("ATG")
                    if atg_idx >= 0:
                        gap_seq = gap_seq[atg_idx:]
                        from Bio.Seq import Seq as BioSeq
                        from .parts import CDS as CDSPart
                        leader_part = CDSPart(
                            id=_next_id("initiation"),
                            name="initiation_leader",
                            sequence=BioSeq(gap_seq),
                            has_stop=False,
                            resolution=ResolutionState.CONCRETE,
                            metadata={
                                "backbone_provided": True,
                                "protein_sequence": str(BioSeq(gap_seq).translate()).rstrip("*"),
                            },
                        )
                        leader_part.protein_sequence = leader_part.metadata["protein_sequence"]
                        assembled.add_part(leader_part)

        # 3. Add all insert parts (everything except the backbone)
        for part in self.parts():
            if isinstance(part, Backbone):
                continue
            assembled.add_part(part)

        # 4. Add downstream backbone features (C-terminal tags, terminator)
        # But only add features that aren't already represented by insert parts
        insert_terminators = [p for p in self.parts() if isinstance(p, Terminator)]
        for feat in downstream:
            # Skip terminator if the insert already has one
            if feat.feature_type == "terminator" and insert_terminators:
                continue
            # Skip C-terminal tags if they duplicate N-terminal features
            # already in upstream (common in pET-28 which has N-His and C-His)
            part = self._backbone_feature_to_part(feat, _next_id)
            if part is not None:
                assembled.add_part(part)

        # 5. Connect linearly
        assembled.connect_linear()

        return assembled

    @staticmethod
    def _backbone_feature_to_part(feat, next_id) -> Optional[GeneticPart]:
        """Convert a BackboneFeature to the appropriate GeneticPart subclass."""
        from .parts import (
            Promoter, RBS, PurificationTag, SolubilityTag, CleavageSite,
            Terminator, CDS, RBSDesignMethod, RegulationType, TagPosition,
        )
        from Bio.Seq import Seq as BioSeq

        seq = BioSeq(feat.sequence) if feat.sequence else None
        ftype = feat.feature_type

        if ftype == "promoter":
            return Promoter(
                id=next_id("promoter"),
                name=feat.label,
                sequence=seq,
                resolution=ResolutionState.CONCRETE,
                metadata={"backbone_provided": True},
            )

        elif ftype == "rbs":
            return RBS(
                id=next_id("rbs"),
                name=feat.label,
                design_method=RBSDesignMethod.LOOKUP,
                part_name=feat.label,
                sequence=seq,
                resolution=ResolutionState.CONCRETE,
                metadata={"backbone_provided": True},
            )

        elif ftype == "tag":
            return PurificationTag(
                id=next_id("tag"),
                name=feat.label,
                tag_type=feat.label,
                position=TagPosition.N_TERM,
                sequence=seq,
                resolution=ResolutionState.CONCRETE,
                metadata={
                    "backbone_provided": True,
                    "protein_sequence": feat.translation or "",
                },
            )

        elif ftype == "solubility_tag":
            return SolubilityTag(
                id=next_id("soltag"),
                name=feat.label,
                tag_type=feat.label,
                position=TagPosition.N_TERM,
                sequence=seq,
                resolution=ResolutionState.CONCRETE,
                metadata={
                    "backbone_provided": True,
                    "protein_sequence": feat.translation or "",
                },
            )

        elif ftype == "cleavage_site":
            return CleavageSite(
                id=next_id("cleavage"),
                name=feat.label,
                protease=feat.label,
                sequence=seq,
                resolution=ResolutionState.CONCRETE,
                metadata={"backbone_provided": True},
            )

        elif ftype == "terminator":
            return Terminator(
                id=next_id("terminator"),
                name=feat.label,
                sequence=seq,
                resolution=ResolutionState.CONCRETE,
                metadata={"backbone_provided": True},
            )

        elif ftype == "lac_operator":
            # Include as metadata but not as a coding part
            return None

        elif ftype == "leader":
            # T7 gene 10 leader — include as a CDS part
            return CDS(
                id=next_id("leader"),
                name=feat.label,
                protein_sequence=feat.translation,
                sequence=seq,
                has_stop=False,
                resolution=ResolutionState.CONCRETE,
                metadata={"backbone_provided": True},
            )

        elif ftype == "start_codon":
            # Annotated start codon — skip, the RBS/tag already handles initiation
            return None

        return None

    def _infer_cloning_pair(self, backbone) -> Optional[tuple[str, str]]:
        """
        Infer the best restriction site cloning pair based on what the
        insert declares vs what each pair provides from the backbone.

        Strategy: for each available pair, build the upstream feature list
        and score it based on overlap with the insert's declared features.
        Features the insert already provides should NOT come from the backbone.
        Features the insert relies on the backbone for (RBS, promoter) should.

        Returns the best pair, or None if no pairs are available.
        """
        from ..data.vector_features import get_insert_context, _get_cloning_pairs
        from .parts import (
            PurificationTag, SolubilityTag, CleavageSite, Backbone,
        )

        pairs = _get_cloning_pairs(backbone.catalog_id)
        if not pairs:
            return None

        # Collect what the insert declares
        insert_has_tags = set()
        insert_has_cleavage = set()
        for part in self.parts():
            if isinstance(part, Backbone):
                continue
            if isinstance(part, PurificationTag):
                insert_has_tags.add(part.name.lower())
            elif isinstance(part, SolubilityTag):
                insert_has_tags.add(part.name.lower())
            elif isinstance(part, CleavageSite):
                insert_has_cleavage.add(part.name.lower())

        best_pair = None
        best_score = -999

        for pair in pairs:
            if len(pair) != 2:
                continue
            ctx = get_insert_context(
                backbone.catalog_id,
                cloning_pair=(pair[0], pair[1]),
            )
            if ctx is None:
                continue

            upstream = ctx["upstream"]
            score = 0

            for feat in upstream:
                if feat.feature_type in ("promoter", "rbs", "lac_operator"):
                    # Always want these from backbone
                    score += 1
                elif feat.feature_type == "tag":
                    # Penalize if insert already has this tag (duplication)
                    if feat.label.lower().replace("x", "x") in insert_has_tags:
                        score -= 2
                    else:
                        score += 0  # neutral — backbone provides extra tag
                elif feat.feature_type == "cleavage_site":
                    if feat.label.lower() in insert_has_cleavage:
                        score -= 2
                    else:
                        score += 0
                elif feat.feature_type in ("leader", "solubility_tag"):
                    if feat.label.lower() in insert_has_tags:
                        score -= 2

            best_pair_candidate = (pair[0], pair[1])
            if score > best_score:
                best_score = score
                best_pair = best_pair_candidate

        return best_pair

    def _copy_graph(self) -> ConstructGraph:
        """Create a shallow copy of this graph."""
        new_graph = ConstructGraph(
            name=self.name,
            description=self.description,
            host_organism=self.host_organism,
            constraints=list(self.constraints),
            metadata=dict(self.metadata),
        )
        for part in self.parts():
            new_graph.add_part(part)
        new_graph._edges = list(self._edges)
        return new_graph

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

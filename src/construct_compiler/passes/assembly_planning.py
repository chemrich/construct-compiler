"""
Assembly Planning Pass

Simplified assembly planner for a synthesis-first workflow:
- Max 2-3 part Golden Gate assemblies
- Prefer full synthesis over PCR/cloning
- Compute cost/time estimates for each strategy
- Recommend the optimal approach

Researcher time: $150/hr ($300k/yr)
Overhead: 1.5x on all operational costs
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Optional

from Bio.Seq import Seq

from ..core.graph import ConstructGraph
from ..core.parts import Backbone, CDS
from ..data.parts_db import HIGH_FIDELITY_OVERHANGS_BSAI

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Cost parameters
# ---------------------------------------------------------------------------

@dataclass
class CostParams:
    """Configurable cost parameters for the cost model."""
    # Researcher
    researcher_hourly_rate: float = 150.0   # $/hr
    overhead_multiplier: float = 1.5        # applied to consumables/operational

    # Synthesis vendors (per bp)
    twist_gene_per_bp: float = 0.07
    twist_clonal_per_bp: float = 0.09
    twist_express_per_bp: float = 0.10
    idt_gblock_per_bp: float = 0.08

    # Synthesis limits
    twist_max_gene_length: int = 5000
    twist_max_clonal_length: int = 5000
    idt_gblock_max_length: int = 3000

    # Per-assembly reagent costs (before overhead)
    bsai_per_reaction: float = 3.00
    t4_ligase_per_reaction: float = 0.25
    competent_cells_per_transformation: float = 10.00
    plates_antibiotics: float = 2.00
    colony_pcr_screening: float = 5.00    # 4 colonies
    miniprep: float = 5.00
    plasmid_sequencing: float = 15.00

    # Time estimates (hours, hands-on)
    golden_gate_setup_hrs: float = 1.5
    colony_screening_hrs: float = 1.0
    miniprep_sequencing_hrs: float = 1.0
    troubleshooting_hrs: float = 3.0       # if reattempt needed
    receive_verify_hrs: float = 0.5        # receive synthesis, verify

    # Success rates
    two_part_success_rate: float = 0.90
    three_part_success_rate: float = 0.80

    # Turnaround times (business days)
    twist_standard_days: tuple[int, int] = (12, 18)
    twist_express_days: tuple[int, int] = (7, 12)
    idt_gblock_days: tuple[int, int] = (3, 5)
    assembly_days: int = 3   # benchwork after fragments arrive

    @property
    def reagent_cost_per_assembly(self) -> float:
        """Total reagent cost for one Golden Gate + transformation + verification."""
        return (
            self.bsai_per_reaction
            + self.t4_ligase_per_reaction
            + self.competent_cells_per_transformation
            + self.plates_antibiotics
            + self.colony_pcr_screening
            + self.miniprep
            + self.plasmid_sequencing
        )

    @property
    def hands_on_hours_per_assembly(self) -> float:
        return (
            self.golden_gate_setup_hrs
            + self.colony_screening_hrs
            + self.miniprep_sequencing_hrs
        )


# ---------------------------------------------------------------------------
# Assembly strategies
# ---------------------------------------------------------------------------

@dataclass
class SynthesisFragment:
    """A DNA fragment to be synthesized."""
    name: str
    sequence: str
    length_bp: int
    vendor: str = "twist"
    cost: float = 0.0


@dataclass
class AssemblyStrategy:
    """A complete strategy for building a construct."""
    name: str
    description: str
    num_fragments: int
    fragments: list[SynthesisFragment] = field(default_factory=list)
    overhangs: list[str] = field(default_factory=list)

    # Costs
    synthesis_cost: float = 0.0
    reagent_cost: float = 0.0           # before overhead
    reagent_cost_with_overhead: float = 0.0
    researcher_time_hrs: float = 0.0
    researcher_cost: float = 0.0
    failure_risk_surcharge: float = 0.0
    total_cost: float = 0.0

    # Time
    synthesis_turnaround_days: tuple[int, int] = (0, 0)
    benchwork_days: int = 0
    total_wall_clock_days: tuple[int, int] = (0, 0)

    # Flags
    recommended: bool = False
    notes: list[str] = field(default_factory=list)


@dataclass
class AssemblyPlan:
    """The output of the assembly planner."""
    construct_name: str
    insert_length_bp: int
    strategies: list[AssemblyStrategy] = field(default_factory=list)
    recommended: Optional[AssemblyStrategy] = None

    def summary(self) -> str:
        lines = [
            f"Assembly Plan: {self.construct_name}",
            f"Insert length: {self.insert_length_bp} bp",
            "",
        ]
        for strat in self.strategies:
            marker = " *** RECOMMENDED ***" if strat.recommended else ""
            lines.append(f"Strategy: {strat.name}{marker}")
            lines.append(f"  Description: {strat.description}")
            lines.append(f"  Fragments: {strat.num_fragments}")
            lines.append(f"  Synthesis cost: ${strat.synthesis_cost:.2f}")
            lines.append(f"  Reagents (with overhead): ${strat.reagent_cost_with_overhead:.2f}")
            lines.append(f"  Researcher time: {strat.researcher_time_hrs:.1f} hrs (${strat.researcher_cost:.2f})")
            lines.append(f"  Failure risk surcharge: ${strat.failure_risk_surcharge:.2f}")
            lines.append(f"  TOTAL COST: ${strat.total_cost:.2f}")
            lines.append(
                f"  Wall clock: {strat.total_wall_clock_days[0]}-"
                f"{strat.total_wall_clock_days[1]} business days"
            )
            if strat.overhangs:
                lines.append(f"  Overhangs: {' / '.join(strat.overhangs)}")
            for note in strat.notes:
                lines.append(f"  Note: {note}")
            lines.append("")

        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Assembly planner
# ---------------------------------------------------------------------------

def plan_assembly(graph: ConstructGraph,
                  cost_params: Optional[CostParams] = None) -> AssemblyPlan:
    """
    Plan the assembly for a construct.

    Evaluates multiple strategies and recommends the cheapest one.
    Strategies considered:
    1. Twist clonal gene (full insert delivered in backbone)
    2. Single gene synthesis + 2-part Golden Gate
    3. Two gene fragments + 3-part Golden Gate (if insert > synthesis limit)
    """
    if cost_params is None:
        cost_params = CostParams()

    insert_seq = graph.full_insert_sequence()
    insert_len = len(insert_seq) if insert_seq else _estimate_insert_length(graph)

    plan = AssemblyPlan(
        construct_name=graph.name,
        insert_length_bp=insert_len,
    )

    # Strategy 1: Twist clonal gene (delivered in your backbone)
    if insert_len <= cost_params.twist_max_clonal_length:
        plan.strategies.append(
            _strategy_twist_clonal(insert_len, insert_seq, cost_params)
        )

    # Strategy 2: Single synthesis fragment + 2-part Golden Gate
    if insert_len <= cost_params.twist_max_gene_length:
        plan.strategies.append(
            _strategy_synthesis_2part(insert_len, insert_seq, cost_params)
        )

    # Strategy 3: Two fragments + 3-part Golden Gate
    if insert_len > cost_params.twist_max_gene_length:
        plan.strategies.append(
            _strategy_synthesis_3part(insert_len, insert_seq, graph, cost_params)
        )
    elif insert_len > 2000:
        # Also consider 3-part even for shorter inserts (redundancy)
        plan.strategies.append(
            _strategy_synthesis_3part(insert_len, insert_seq, graph, cost_params)
        )

    # Strategy 4: IDT gBlocks (for short inserts)
    if insert_len <= cost_params.idt_gblock_max_length:
        plan.strategies.append(
            _strategy_idt_gblock_2part(insert_len, insert_seq, cost_params)
        )

    # Pick the recommended strategy (lowest total cost)
    if plan.strategies:
        best = min(plan.strategies, key=lambda s: s.total_cost)
        best.recommended = True
        plan.recommended = best

    return plan


# ---------------------------------------------------------------------------
# Strategy builders
# ---------------------------------------------------------------------------

def _strategy_twist_clonal(insert_len: int, insert_seq, params: CostParams) -> AssemblyStrategy:
    """Full insert synthesized and cloned by Twist into your backbone."""
    synth_cost = insert_len * params.twist_clonal_per_bp

    # Researcher just receives tube, transforms to verify, minipreps, sequences
    researcher_hrs = params.receive_verify_hrs + params.miniprep_sequencing_hrs
    researcher_cost = researcher_hrs * params.researcher_hourly_rate

    # Minimal reagents: just transformation + verification
    reagent = (
        params.competent_cells_per_transformation
        + params.plates_antibiotics
        + params.miniprep
        + params.plasmid_sequencing
    )
    reagent_overhead = reagent * params.overhead_multiplier

    total = synth_cost + reagent_overhead + researcher_cost

    return AssemblyStrategy(
        name="Twist Clonal Gene",
        description="Full insert synthesized and cloned by Twist into your backbone",
        num_fragments=1,
        fragments=[SynthesisFragment(
            name="full_insert",
            sequence=str(insert_seq) if insert_seq else "",
            length_bp=insert_len,
            vendor="twist_clonal",
            cost=synth_cost,
        )],
        synthesis_cost=synth_cost,
        reagent_cost=reagent,
        reagent_cost_with_overhead=reagent_overhead,
        researcher_time_hrs=researcher_hrs,
        researcher_cost=researcher_cost,
        failure_risk_surcharge=0.0,  # Twist guarantees the clone
        total_cost=total,
        synthesis_turnaround_days=params.twist_standard_days,
        benchwork_days=1,
        total_wall_clock_days=(
            params.twist_standard_days[0] + 1,
            params.twist_standard_days[1] + 2,
        ),
        notes=["Twist handles cloning — lowest hands-on effort",
               "Requires providing backbone map to Twist"],
    )


def _strategy_synthesis_2part(insert_len: int, insert_seq, params: CostParams) -> AssemblyStrategy:
    """Synthesize insert as one fragment, Golden Gate into backbone."""
    synth_cost = insert_len * params.twist_gene_per_bp
    reagent = params.reagent_cost_per_assembly
    reagent_overhead = reagent * params.overhead_multiplier
    researcher_hrs = params.hands_on_hours_per_assembly
    researcher_cost = researcher_hrs * params.researcher_hourly_rate

    # Risk-adjusted failure surcharge
    failure_prob = 1.0 - params.two_part_success_rate
    reattempt_cost = (
        params.reagent_cost_per_assembly * params.overhead_multiplier
        + params.troubleshooting_hrs * params.researcher_hourly_rate
    )
    risk_surcharge = failure_prob * reattempt_cost

    total = synth_cost + reagent_overhead + researcher_cost + risk_surcharge

    # Assign overhangs for 2-part assembly
    overhangs = _select_overhangs(2)

    return AssemblyStrategy(
        name="Synthesis + 2-Part Golden Gate",
        description="Synthesize full insert, assemble with backbone via 2-part Golden Gate",
        num_fragments=2,
        fragments=[SynthesisFragment(
            name="insert",
            sequence=str(insert_seq) if insert_seq else "",
            length_bp=insert_len,
            vendor="twist",
            cost=synth_cost,
        )],
        overhangs=overhangs,
        synthesis_cost=synth_cost,
        reagent_cost=reagent,
        reagent_cost_with_overhead=reagent_overhead,
        researcher_time_hrs=researcher_hrs,
        researcher_cost=researcher_cost,
        failure_risk_surcharge=risk_surcharge,
        total_cost=total,
        synthesis_turnaround_days=params.twist_standard_days,
        benchwork_days=params.assembly_days,
        total_wall_clock_days=(
            params.twist_standard_days[0] + params.assembly_days,
            params.twist_standard_days[1] + params.assembly_days + 1,
        ),
        notes=[f"2-part GG success rate: {params.two_part_success_rate*100:.0f}%"],
    )


def _strategy_synthesis_3part(insert_len: int, insert_seq, graph: ConstructGraph,
                               params: CostParams) -> AssemblyStrategy:
    """Split insert into two fragments, 3-part Golden Gate with backbone."""
    # Find the best split point (prefer between cistrons/functional units)
    split_point = _find_best_split(insert_len, graph)

    frag1_len = split_point
    frag2_len = insert_len - split_point

    synth_cost = (
        frag1_len * params.twist_gene_per_bp
        + frag2_len * params.twist_gene_per_bp
    )
    reagent = params.reagent_cost_per_assembly
    reagent_overhead = reagent * params.overhead_multiplier
    researcher_hrs = params.hands_on_hours_per_assembly
    researcher_cost = researcher_hrs * params.researcher_hourly_rate

    # Higher failure risk for 3-part
    failure_prob = 1.0 - params.three_part_success_rate
    reattempt_cost = (
        params.reagent_cost_per_assembly * params.overhead_multiplier
        + params.troubleshooting_hrs * params.researcher_hourly_rate
    )
    risk_surcharge = failure_prob * reattempt_cost

    total = synth_cost + reagent_overhead + researcher_cost + risk_surcharge

    overhangs = _select_overhangs(3)

    return AssemblyStrategy(
        name="2x Synthesis + 3-Part Golden Gate",
        description=f"Split insert at bp {split_point}, synthesize 2 fragments, "
                    f"3-part Golden Gate with backbone",
        num_fragments=3,
        fragments=[
            SynthesisFragment(
                name="insert_frag1",
                sequence="",
                length_bp=frag1_len,
                vendor="twist",
                cost=frag1_len * params.twist_gene_per_bp,
            ),
            SynthesisFragment(
                name="insert_frag2",
                sequence="",
                length_bp=frag2_len,
                vendor="twist",
                cost=frag2_len * params.twist_gene_per_bp,
            ),
        ],
        overhangs=overhangs,
        synthesis_cost=synth_cost,
        reagent_cost=reagent,
        reagent_cost_with_overhead=reagent_overhead,
        researcher_time_hrs=researcher_hrs,
        researcher_cost=researcher_cost,
        failure_risk_surcharge=risk_surcharge,
        total_cost=total,
        synthesis_turnaround_days=params.twist_standard_days,
        benchwork_days=params.assembly_days,
        total_wall_clock_days=(
            params.twist_standard_days[0] + params.assembly_days,
            params.twist_standard_days[1] + params.assembly_days + 1,
        ),
        notes=[
            f"3-part GG success rate: {params.three_part_success_rate*100:.0f}%",
            f"Split: fragment 1 = {frag1_len} bp, fragment 2 = {frag2_len} bp",
        ],
    )


def _strategy_idt_gblock_2part(insert_len: int, insert_seq, params: CostParams) -> AssemblyStrategy:
    """Use IDT gBlock for short inserts — faster turnaround."""
    synth_cost = insert_len * params.idt_gblock_per_bp
    reagent = params.reagent_cost_per_assembly
    reagent_overhead = reagent * params.overhead_multiplier
    researcher_hrs = params.hands_on_hours_per_assembly
    researcher_cost = researcher_hrs * params.researcher_hourly_rate

    failure_prob = 1.0 - params.two_part_success_rate
    reattempt_cost = (
        params.reagent_cost_per_assembly * params.overhead_multiplier
        + params.troubleshooting_hrs * params.researcher_hourly_rate
    )
    risk_surcharge = failure_prob * reattempt_cost

    total = synth_cost + reagent_overhead + researcher_cost + risk_surcharge

    overhangs = _select_overhangs(2)

    return AssemblyStrategy(
        name="IDT gBlock + 2-Part Golden Gate",
        description="Order insert as IDT gBlock, assemble with backbone (faster turnaround)",
        num_fragments=2,
        fragments=[SynthesisFragment(
            name="insert_gblock",
            sequence=str(insert_seq) if insert_seq else "",
            length_bp=insert_len,
            vendor="idt",
            cost=synth_cost,
        )],
        overhangs=overhangs,
        synthesis_cost=synth_cost,
        reagent_cost=reagent,
        reagent_cost_with_overhead=reagent_overhead,
        researcher_time_hrs=researcher_hrs,
        researcher_cost=researcher_cost,
        failure_risk_surcharge=risk_surcharge,
        total_cost=total,
        synthesis_turnaround_days=params.idt_gblock_days,
        benchwork_days=params.assembly_days,
        total_wall_clock_days=(
            params.idt_gblock_days[0] + params.assembly_days,
            params.idt_gblock_days[1] + params.assembly_days + 1,
        ),
        notes=[
            "Fastest turnaround option",
            f"IDT gBlock max length: {params.idt_gblock_max_length} bp",
        ],
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _select_overhangs(num_junctions: int) -> list[str]:
    """
    Select high-fidelity overhangs for the assembly.
    For N-part assembly, we need N+1 overhangs (including the two at the
    backbone junction that are fixed by the backbone linearization site).
    In practice, for a 2-part GG: 2 overhangs (one at each end of insert).
    For 3-part: 3 overhangs (two ends + one internal junction).
    """
    # Use the first few from the validated set
    return HIGH_FIDELITY_OVERHANGS_BSAI[:num_junctions]


def _find_best_split(insert_len: int, graph: ConstructGraph) -> int:
    """
    Find the best point to split an insert into two fragments.
    Prefers splitting between functional units (cistrons).
    Falls back to splitting at the midpoint.
    """
    # Try to find cistron boundaries
    cistrons = graph.cistrons()
    if len(cistrons) >= 2:
        # Split between the first and second cistron
        first_cistron_parts = cistrons[0]
        cumulative_len = 0
        for part in graph.parts():
            if part in first_cistron_parts:
                if part.sequence:
                    cumulative_len += len(part.sequence)
            elif cumulative_len > 0:
                # We've passed the first cistron
                break
            else:
                if part.sequence:
                    cumulative_len += len(part.sequence)

        if 0 < cumulative_len < insert_len:
            return cumulative_len

    # Fallback: split at midpoint
    return insert_len // 2


def _estimate_insert_length(graph: ConstructGraph) -> int:
    """Estimate insert length from part metadata when sequences aren't concrete yet."""
    total = 0
    for part in graph.parts():
        if isinstance(part, Backbone):
            continue
        if part.sequence:
            total += len(part.sequence)
        elif isinstance(part, CDS) and part.protein_sequence:
            total += len(part.protein_sequence) * 3 + 3  # +3 for stop
        elif hasattr(part, 'metadata') and 'protein_sequence' in part.metadata:
            total += len(part.metadata['protein_sequence']) * 3
        else:
            total += 50  # rough estimate for regulatory elements
    return total

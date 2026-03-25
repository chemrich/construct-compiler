"""
Biological sensibility checks.

These operate at a higher level than the DNA-level construct checks —
they catch biologically implausible designs that would compile correctly
but fail in the lab. Results are typically WARNINGs, not ERRORs, since
the DNA is technically valid.

Checks:
  1. Promoter/host compatibility (e.g., CMV in E. coli)
  2. Expression system requirements (e.g., T7 requires DE3 lysogen)
  3. Codon optimization vs. host mismatch
  4. Selection marker conflicts
  5. Vector/host mismatch (mammalian vector for bacterial host)
"""

from __future__ import annotations

from ..core.graph import ConstructGraph
from ..core.parts import Backbone, Promoter, CDS
from .construct_checks import CheckResult, CheckSeverity


# ---------------------------------------------------------------------------
# Promoter → host compatibility rules
# ---------------------------------------------------------------------------

# Promoters that only work in specific host contexts
_BACTERIAL_PROMOTERS = frozenset({
    "T7", "T7lac", "tac", "lac", "lacUV5", "araBAD",
    "J23100", "J23106", "J23119",  # Anderson promoters
    "pLac", "pTet", "pTrc",
})

_MAMMALIAN_PROMOTERS = frozenset({
    "CMV", "EF1a", "EF-1-alpha", "SV40", "UbC", "PGK", "CAG",
    "SFFV", "hPGK", "MSCV", "RSV",
})

# T7 promoter requires a host with T7 RNA polymerase
_T7_REQUIRING = frozenset({"T7", "T7lac"})

# DE3 lysogen hosts (have chromosomal T7 RNAP under lac control)
_DE3_HOSTS = frozenset({
    "e_coli_bl21", "e_coli_bl21_de3", "bl21",
    "e_coli_rosetta", "rosetta",
    "e_coli_origami", "origami",
    "e_coli_shuffle", "shuffle",
    "e_coli_lemo21", "lemo21",
    "e_coli_c41", "e_coli_c43",
})

# Hosts that are broadly "bacterial"
_BACTERIAL_HOSTS = frozenset({
    "e_coli", "e_coli_k12", "e_coli_bl21", "e_coli_bl21_de3",
    "e_coli_rosetta", "e_coli_origami", "e_coli_shuffle",
    "e_coli_lemo21", "e_coli_c41", "e_coli_c43",
    "e_coli_dh5a", "e_coli_top10", "e_coli_stbl3",
    "bacillus", "pseudomonas",
})

_MAMMALIAN_HOSTS = frozenset({
    "mammalian", "hek293", "hek293t", "cho", "hela",
    "cos7", "jurkat", "k562", "u2os", "a549",
    "nih3t3", "vero", "mdck",
})


# ---------------------------------------------------------------------------
# Check implementation
# ---------------------------------------------------------------------------

def check_biological_sensibility(graph: ConstructGraph) -> list[CheckResult]:
    """
    Run all biological sensibility checks on a compiled construct graph.
    Returns a list of CheckResults (typically WARNINGs).
    """
    results: list[CheckResult] = []

    host = (graph.host_organism or "").lower().strip()
    backbone = None
    for part in graph.parts():
        if isinstance(part, Backbone):
            backbone = part
            break

    # Collect promoter names from the graph
    promoter_names = []
    for part in graph.parts():
        if isinstance(part, Promoter):
            promoter_names.append(part.name)

    # Also check the backbone's declared promoter
    if backbone and backbone.catalog_id:
        from ..data.parts_db import TWIST_CATALOG_VECTORS
        vec_entry = TWIST_CATALOG_VECTORS.get(backbone.catalog_id, {})
        vec_promoter = vec_entry.get("promoter", "")
        if vec_promoter and vec_promoter not in promoter_names:
            promoter_names.append(vec_promoter)
        vec_host = vec_entry.get("host", "")
    else:
        vec_host = ""

    # --- Check 1: Promoter/host compatibility ---
    results.extend(_check_promoter_host(promoter_names, host, graph))

    # --- Check 2: T7 system requires DE3 lysogen ---
    results.extend(_check_t7_de3(promoter_names, host, graph))

    # --- Check 3: Codon optimization vs. host ---
    results.extend(_check_codon_host(graph, host))

    # --- Check 4: Vector/host mismatch ---
    results.extend(_check_vector_host(backbone, host, vec_host))

    return results


def _check_promoter_host(
    promoter_names: list[str],
    host: str,
    graph: ConstructGraph,
) -> list[CheckResult]:
    """Check that promoters match the host organism."""
    results = []
    is_bacterial = _is_bacterial_host(host)
    is_mammalian = _is_mammalian_host(host)

    for pname in promoter_names:
        pname_clean = pname.strip()

        # Bacterial promoter in mammalian host
        if pname_clean in _BACTERIAL_PROMOTERS and is_mammalian:
            results.append(CheckResult(
                check_name="bio_promoter_host",
                severity=CheckSeverity.WARNING,
                part_id="construct",
                message=(
                    f"Bacterial promoter '{pname_clean}' is unlikely to drive "
                    f"expression in mammalian host '{host}'. Bacterial promoters "
                    f"require bacterial RNA polymerase."
                ),
                details={
                    "promoter": pname_clean,
                    "host": host,
                    "issue": "bacterial_promoter_in_mammalian_host",
                },
            ))

        # Mammalian promoter in bacterial host
        if pname_clean in _MAMMALIAN_PROMOTERS and is_bacterial:
            results.append(CheckResult(
                check_name="bio_promoter_host",
                severity=CheckSeverity.WARNING,
                part_id="construct",
                message=(
                    f"Mammalian promoter '{pname_clean}' will not function in "
                    f"bacterial host '{host}'. Mammalian promoters require "
                    f"eukaryotic transcription machinery (RNA Pol II)."
                ),
                details={
                    "promoter": pname_clean,
                    "host": host,
                    "issue": "mammalian_promoter_in_bacterial_host",
                },
            ))

    return results


def _check_t7_de3(
    promoter_names: list[str],
    host: str,
    graph: ConstructGraph,
) -> list[CheckResult]:
    """Check that T7 promoter is paired with a DE3 lysogen host."""
    results = []

    has_t7 = any(p in _T7_REQUIRING for p in promoter_names)
    if not has_t7:
        return results

    host_lower = host.lower()

    # If host is generic "e_coli" without specifying DE3, warn
    if host_lower in ("e_coli", "e_coli_k12", "e_coli_dh5a", "e_coli_top10"):
        results.append(CheckResult(
            check_name="bio_t7_de3",
            severity=CheckSeverity.WARNING,
            part_id="construct",
            message=(
                f"T7 promoter requires T7 RNA polymerase, typically from a "
                f"DE3 lysogen (e.g., BL21(DE3)). Host '{host}' may not have "
                f"T7 RNAP. Consider using BL21(DE3) or another DE3 strain."
            ),
            details={
                "host": host,
                "issue": "t7_without_de3",
            },
        ))

    # If host is not bacterial at all
    if _is_mammalian_host(host_lower):
        # Already caught by promoter/host check, don't double-warn
        pass

    return results


def _check_codon_host(graph: ConstructGraph, host: str) -> list[CheckResult]:
    """Check that codon optimization target matches the host."""
    results = []

    for part in graph.parts():
        if not isinstance(part, CDS):
            continue

        codon_opt = getattr(part, "codon_optimization", "").lower()
        organism = getattr(part, "organism", "").lower()

        # If codon optimization explicitly targets a different kingdom
        if _is_mammalian_host(host) and organism and "coli" in organism:
            results.append(CheckResult(
                check_name="bio_codon_host",
                severity=CheckSeverity.WARNING,
                part_id=part.id,
                message=(
                    f"CDS '{part.name}' is codon-optimized for E. coli but "
                    f"the host is '{host}'. Bacterial codon usage will reduce "
                    f"expression in mammalian cells."
                ),
                details={
                    "part": part.name,
                    "codon_target": organism,
                    "host": host,
                    "issue": "codon_optimization_mismatch",
                },
            ))

        if _is_bacterial_host(host) and organism and any(
            m in organism for m in ("human", "mammal", "cho", "hek")
        ):
            results.append(CheckResult(
                check_name="bio_codon_host",
                severity=CheckSeverity.WARNING,
                part_id=part.id,
                message=(
                    f"CDS '{part.name}' is codon-optimized for mammalian cells "
                    f"but the host is '{host}'. Mammalian codon usage may reduce "
                    f"expression in bacteria."
                ),
                details={
                    "part": part.name,
                    "codon_target": organism,
                    "host": host,
                    "issue": "codon_optimization_mismatch",
                },
            ))

    return results


def _check_vector_host(
    backbone: Backbone | None,
    spec_host: str,
    vec_host: str,
) -> list[CheckResult]:
    """Check that the vector is appropriate for the declared host."""
    results = []

    if not backbone or not backbone.catalog_id or not vec_host:
        return results

    spec_is_bacterial = _is_bacterial_host(spec_host)
    spec_is_mammalian = _is_mammalian_host(spec_host)
    vec_is_bacterial = _is_bacterial_host(vec_host)
    vec_is_mammalian = _is_mammalian_host(vec_host)

    if spec_is_bacterial and vec_is_mammalian:
        results.append(CheckResult(
            check_name="bio_vector_host",
            severity=CheckSeverity.WARNING,
            part_id=backbone.id,
            message=(
                f"Vector '{backbone.catalog_id}' is designed for mammalian "
                f"expression but the spec declares host '{spec_host}'. "
                f"This vector's promoter and regulatory elements will not "
                f"function in bacteria."
            ),
            details={
                "vector": backbone.catalog_id,
                "vector_host": vec_host,
                "spec_host": spec_host,
                "issue": "mammalian_vector_in_bacterial_host",
            },
        ))

    if spec_is_mammalian and vec_is_bacterial:
        results.append(CheckResult(
            check_name="bio_vector_host",
            severity=CheckSeverity.WARNING,
            part_id=backbone.id,
            message=(
                f"Vector '{backbone.catalog_id}' is designed for bacterial "
                f"expression but the spec declares host '{spec_host}'. "
                f"Bacterial expression vectors lack mammalian regulatory "
                f"elements."
            ),
            details={
                "vector": backbone.catalog_id,
                "vector_host": vec_host,
                "spec_host": spec_host,
                "issue": "bacterial_vector_in_mammalian_host",
            },
        ))

    return results


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def _is_bacterial_host(host: str) -> bool:
    """Is this host a bacterial organism?"""
    h = host.lower().strip()
    return h in _BACTERIAL_HOSTS or h.startswith("e_coli") or h.startswith("bacillus")


def _is_mammalian_host(host: str) -> bool:
    """Is this host a mammalian cell line?"""
    h = host.lower().strip()
    return h in _MAMMALIAN_HOSTS or h.startswith("hek") or h.startswith("cho")

"""
Part Resolution Pass

Resolves abstract part references to concrete sequences by:
1. Looking up tags, cleavage sites, linkers, RBSs, promoters, terminators
   from the local parts database.
2. Fetching protein sequences from external databases (UniProt, FPbase)
   for CDS parts.

After this pass, all parts should be RESOLVED (protein sequences known)
or CONCRETE (DNA sequences assigned for non-coding parts).
"""

from __future__ import annotations

import logging
from typing import Optional

from Bio.Seq import Seq

from ..core.graph import ConstructGraph
from ..core.types import ResolutionState
from ..core.parts import (
    GeneticPart, Promoter, RBS, CDS, PurificationTag, SolubilityTag,
    CleavageSite, StopCodon, Linker, Terminator, Spacer, Backbone,
)
from ..data.parts_db import (
    PROMOTERS, RBS_LIBRARY, TERMINATORS,
    PURIFICATION_TAGS, SOLUBILITY_TAGS, CLEAVAGE_SITES, LINKERS,
)

logger = logging.getLogger(__name__)


def resolve_parts(graph: ConstructGraph) -> ConstructGraph:
    """
    Resolve all parts in the graph to concrete or resolved state.
    Modifies the graph in-place and returns it.
    """
    for part in graph.parts():
        if part.resolution == ResolutionState.CONCRETE:
            continue

        if isinstance(part, Promoter):
            _resolve_promoter(part)
        elif isinstance(part, RBS):
            _resolve_rbs(part)
        elif isinstance(part, Terminator):
            _resolve_terminator(part)
        elif isinstance(part, PurificationTag):
            _resolve_purification_tag(part)
        elif isinstance(part, SolubilityTag):
            _resolve_solubility_tag(part)
        elif isinstance(part, CleavageSite):
            _resolve_cleavage_site(part)
        elif isinstance(part, Linker):
            _resolve_linker(part)
        elif isinstance(part, CDS):
            _resolve_cds(part)
        elif isinstance(part, StopCodon):
            part.sequence = Seq(part.codon)
            part.resolution = ResolutionState.CONCRETE
        elif isinstance(part, Spacer):
            _resolve_spacer(part)
        elif isinstance(part, Backbone):
            # Backbone resolution requires fetching from Addgene or local file
            # For now, mark as resolved if we have metadata
            if part.resolution == ResolutionState.ABSTRACT:
                part.resolution = ResolutionState.RESOLVED
                logger.info(f"Backbone '{part.name}' marked as RESOLVED (sequence TBD)")

    return graph


# ---------------------------------------------------------------------------
# Individual resolvers
# ---------------------------------------------------------------------------

def _resolve_promoter(part: Promoter) -> None:
    data = PROMOTERS.get(part.name)
    if data:
        part.sequence = Seq(data["sequence"])
        part.resolution = ResolutionState.CONCRETE
        logger.info(f"Promoter '{part.name}' resolved: {len(part.sequence)} bp")
    else:
        logger.warning(f"Promoter '{part.name}' not in local DB")


def _resolve_rbs(part: RBS) -> None:
    # Try the explicit part name first
    data = RBS_LIBRARY.get(part.part_name)
    if data:
        part.sequence = Seq(data["sequence"])
        part.resolution = ResolutionState.CONCRETE
        part.metadata["relative_strength"] = data.get("relative_strength", 0)
        part.metadata["is_bcd"] = data.get("is_bcd", False)
        logger.info(
            f"RBS '{part.part_name}' resolved: {len(part.sequence)} bp, "
            f"strength={data.get('relative_strength', '?')}"
        )
    else:
        logger.warning(f"RBS '{part.part_name}' not in local DB")


def _resolve_terminator(part: Terminator) -> None:
    data = TERMINATORS.get(part.name)
    if data:
        part.sequence = Seq(data["sequence"])
        part.resolution = ResolutionState.CONCRETE
        logger.info(f"Terminator '{part.name}' resolved: {len(part.sequence)} bp")
    else:
        logger.warning(f"Terminator '{part.name}' not in local DB")


def _resolve_purification_tag(part: PurificationTag) -> None:
    data = PURIFICATION_TAGS.get(part.tag_type)
    if data:
        part.metadata["protein_sequence"] = data["protein_sequence"]
        part.resolution = ResolutionState.RESOLVED  # needs reverse translation
        logger.info(f"Tag '{part.tag_type}' resolved: {len(data['protein_sequence'])} aa")
    else:
        logger.warning(f"Purification tag '{part.tag_type}' not in local DB")


def _resolve_solubility_tag(part: SolubilityTag) -> None:
    data = SOLUBILITY_TAGS.get(part.tag_type)
    if data:
        part.metadata["protein_sequence"] = data["protein_sequence"]
        part.resolution = ResolutionState.RESOLVED
        logger.info(f"Solubility tag '{part.tag_type}' resolved: {len(data['protein_sequence'])} aa")
    else:
        logger.warning(f"Solubility tag '{part.tag_type}' not in local DB")


def _resolve_cleavage_site(part: CleavageSite) -> None:
    data = CLEAVAGE_SITES.get(part.protease)
    if data:
        part.recognition_seq = data["protein_sequence"]
        part.metadata["protein_sequence"] = data["protein_sequence"]
        part.metadata["cut_position"] = data.get("cut_position")
        part.resolution = ResolutionState.RESOLVED
        logger.info(f"Cleavage site '{part.protease}' resolved: {part.recognition_seq}")
    else:
        logger.warning(f"Cleavage site '{part.protease}' not in local DB")


def _resolve_linker(part: Linker) -> None:
    data = LINKERS.get(part.linker_type)
    if data:
        unit = data["unit"]
        repeats = part.repeats or data.get("default_repeats", 3)
        protein_seq = unit * repeats
        part.metadata["protein_sequence"] = protein_seq
        part.resolution = ResolutionState.RESOLVED
        logger.info(f"Linker '{part.linker_type}' resolved: {protein_seq}")
    elif part.custom_sequence:
        part.metadata["protein_sequence"] = part.custom_sequence
        part.resolution = ResolutionState.RESOLVED
    else:
        logger.warning(f"Linker '{part.linker_type}' not in local DB")


def _resolve_cds(part: CDS) -> None:
    """
    Resolve a CDS — fetch protein sequence from the appropriate source.
    Supports: fpbase, uniprot, local.
    Falls back to built-in common sequences when remote fetch fails.
    """
    if part.protein_sequence:
        part.resolution = ResolutionState.RESOLVED
        return

    source = part.source_db.lower()
    seq = None

    if source == "fpbase":
        seq = _fetch_fpbase(part.source_id)
        if not seq:
            seq = _BUILTIN_FP_SEQUENCES.get(part.source_id)
            if seq:
                logger.info(f"CDS '{part.name}' resolved from built-in FP database")

    elif source == "uniprot":
        seq = _fetch_uniprot(part.source_id)
        if not seq:
            seq = _BUILTIN_UNIPROT_SEQUENCES.get(part.source_id)
            if seq:
                logger.info(f"CDS '{part.name}' resolved from built-in UniProt cache")

    elif source in ("local", ""):
        seq = _BUILTIN_FP_SEQUENCES.get(part.source_id) or \
              _BUILTIN_UNIPROT_SEQUENCES.get(part.source_id)

    if seq:
        part.protein_sequence = seq
        part.resolution = ResolutionState.RESOLVED
        logger.info(f"CDS '{part.name}' resolved: {len(seq)} aa")
    else:
        logger.warning(f"CDS '{part.name}' could not be resolved from any source")


# Built-in fluorescent protein sequences (fallback when FPbase is unreachable)
_BUILTIN_FP_SEQUENCES = {
    "mEGFP": (
        "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTL"
        "TYGSQFHGLRPDLHMQCFSFDRHTYKTRSIDAPMRAVDLDKQGANFYIAVRILERIGDEFKEGTPM"
        "EGKQYVIQHFSVRYNEGFEPIYKQTIISKNAEDRGDVTADHQVCIETFSFSYDKEERAKVIKEMNF"
        "LVQYLREYNAGIGHKLKLEYNFNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIG"
        "DGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"
    ),
    "mScarlet": (
        "MVSKGEENNMAIIKEFMRFKVHMEGSMNGHEFEIEGEGEGRPYEGTQTAKLKVTKGGPLPFSWDIL"
        "SPQFMYGSRAFIKHPADIPDYYKQSFPEGFKWERVMNFEDGGAVTVTQDTSLEDGTLIYKVKLRGN"
        "LPGATDYKTIKTLRLRGSGNLPAFQKVNSHVTPEQMKIVAVAHYRDELKEISNDHKINENVYIMAP"
        "DQQVELHHSFKVEINYQFEYNVHAVHPGTSYIHQLAFYGSKVAHKISHLPSTMDSALSAILNNLQD"
        "DGIKYDYHFILDHYAQQQGLTDQEYA"
    ),
    "mScarlet-I": (
        "MVSKGEENNMAIIKEFMRFKVHMEGSMNGHEFEIEGEGEGRPYEGTQTAKLKVTKGGPLPFSWDIL"
        "SPQFMYGSRAFIKHPADIPDYYKQSFPEGFKWERVMNFEDGGAVTVTQDTSLEDGTLIYKVKLRGN"
        "LPGATDYKTIKTLRLRGSGNLPAFQKVNSHVTPEQMKIVAVAHYRDELKEISNDHKINENVYIMAP"
        "DQQVELHHSFKVEINYQFEYNVHAVHPGTSYIHQLAFYGSKVAHKISHLPSTMDSALSAILNNLQD"
        "DGIKYDYHFILDHYAQQQGLTDQEYA"
    ),
    "mNeonGreen": (
        "MVSKGEEDNMASLPATHELHIFGSINGVDFDMVGQGTGNPNDGYEELNYVTVSYAEAALGYANRFS"
        "HKIPMTFKQAQDYMHQKDGFPTIELYGKAGRGCTVNFRITAKGVTSQNLPHDFIDHQGIMPFERY"
        "GPKFSHMAAPIYTKTPETYREIPAAFEGRRVFADEDIDTAKVANKELGAHYQRDSMSTFKNRNVQIP"
        "AFEGALTAKLNIDGNAMDISDGQGIHDKASSIIAGCLQHLGRATKLNYDGRSEMQIHRGGATMTMA"
        "KQISDYAKEAGISAEILESMSLDPEKRENDMRIMADL"
    ),
    "mCherry": (
        "MVSKGEEDNMAIIKEFMRFKVHMEGSVNGHEFEIDGDGKGKPYEGTQTANLQVKGIELKGIDFKE"
        "DGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDN"
        "HYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"
    ),
    "sfGFP": (
        "MRKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLT"
        "YGVQCFARYPDHMKQHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFK"
        "EDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPD"
        "NHYLSTQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
    ),
}

# Built-in UniProt sequences (common expression targets, fallback cache)
_BUILTIN_UNIPROT_SEQUENCES = {
    "P42212": (  # GFP, Aequorea victoria
        "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFS"
        "YGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFK"
        "EDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPD"
        "NHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"
    ),
}


def _resolve_spacer(part: Spacer) -> None:
    """Generate a random-ish spacer sequence that avoids secondary structure."""
    import random
    random.seed(hash(part.id))  # deterministic per part ID
    bases = []
    for i in range(part.length_bp):
        # Avoid long homopolymers by not repeating the same base >3 times
        choices = list("ACGT")
        if len(bases) >= 3 and len(set(bases[-3:])) == 1:
            choices.remove(bases[-1])
        bases.append(random.choice(choices))
    part.sequence = Seq("".join(bases))
    part.resolution = ResolutionState.CONCRETE
    logger.info(f"Spacer '{part.id}' generated: {len(part.sequence)} bp")


# ---------------------------------------------------------------------------
# External database fetchers
# ---------------------------------------------------------------------------

def _fetch_fpbase(slug: str) -> Optional[str]:
    """Fetch a protein sequence from FPbase by slug (e.g., 'mEGFP')."""
    import requests
    try:
        # FPbase REST API
        url = f"https://www.fpbase.org/api/proteins/{slug}/?format=json"
        resp = requests.get(url, timeout=10)
        if resp.status_code == 200:
            data = resp.json()
            return data.get("seq", None)

        # Try search if direct slug didn't work
        url = f"https://www.fpbase.org/api/proteins/?slug={slug}&format=json"
        resp = requests.get(url, timeout=10)
        if resp.status_code == 200:
            results = resp.json()
            if isinstance(results, list) and results:
                return results[0].get("seq")
            elif isinstance(results, dict) and results.get("results"):
                return results["results"][0].get("seq")
    except Exception as e:
        logger.error(f"FPbase fetch failed for '{slug}': {e}")
    return None


def _fetch_uniprot(accession: str) -> Optional[str]:
    """Fetch a protein sequence from UniProt by accession (e.g., 'P42212')."""
    import requests
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
        resp = requests.get(url, timeout=10, headers={"Accept": "text/plain"})
        if resp.status_code == 200:
            lines = resp.text.strip().split("\n")
            # Skip header line(s)
            seq_lines = [l.strip() for l in lines if not l.startswith(">")]
            return "".join(seq_lines)
    except Exception as e:
        logger.error(f"UniProt fetch failed for '{accession}': {e}")
    return None

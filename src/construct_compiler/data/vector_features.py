"""
Backbone vector feature extraction from GenBank files.

Parses the Twist catalog vector GenBank files to extract annotated
expression-region features (promoter, RBS, tags, cleavage sites,
terminators) with their actual DNA sequences. These are used to build
an "assembled view" of the full construct — insert + backbone elements.

The assembled view is what the cell actually sees, and is the correct
target for validation checks like reading frame continuity, start codon
placement, and cistron counting.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from Bio import SeqIO
from Bio.Seq import Seq

# ---------------------------------------------------------------------------
# Data directory
# ---------------------------------------------------------------------------

# GenBank files live at project_root/data/twist_vectors/
# Navigate up from src/construct_compiler/data/ to project root
_VECTOR_DIR = Path(__file__).parent.parent.parent.parent / "data" / "twist_vectors"

# ---------------------------------------------------------------------------
# Feature types we care about in the expression region
# ---------------------------------------------------------------------------

# GenBank feature types that are part of the expression cassette
_EXPR_FEATURE_TYPES = {"promoter", "protein_bind", "RBS", "CDS", "terminator"}

# CDS labels that are NOT expression features (resistance, housekeeping, etc.)
_SKIP_CDS_LABELS = frozenset({
    "KanR", "AmpR", "CmR", "NeoR/KanR", "PuroR", "HygR",
    "lacI", "rop", "lacZ-alpha",
    "psi", "Major splice donor", "pMB1 Origin",
    "Ampicillin (AmpR)",
})

# Promoter labels that are NOT the main expression promoter
_SKIP_PROMOTER_LABELS = frozenset({
    "AmpR promoter", "lacI promoter", "tet promoter", "cat promoter",
    "SV40 promoter", "UbC promoter", "lac promoter",
})

# protein_bind labels that are NOT relevant to expression
_SKIP_BIND_LABELS = frozenset({
    "CAP binding site",
})


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass
class BackboneFeature:
    """A single annotated feature from the backbone vector."""
    feature_type: str       # "promoter", "rbs", "tag", "cleavage_site", "terminator", etc.
    label: str              # e.g., "T7 promoter", "6xHis", "thrombin site"
    sequence: str           # DNA sequence (5'→3')
    translation: Optional[str] = None  # protein sequence if CDS
    start: int = 0          # 1-based start position in backbone
    end: int = 0            # 1-based end position in backbone
    genbank_type: str = ""  # original GenBank feature type
    note: str = ""

    @property
    def length_bp(self) -> int:
        return len(self.sequence)

    def __repr__(self):
        return (
            f"BackboneFeature({self.feature_type!r}, {self.label!r}, "
            f"{self.length_bp}bp, pos={self.start}..{self.end})"
        )


@dataclass
class VectorRecord:
    """Parsed backbone vector with expression-region features."""
    name: str                       # e.g., "pET-28b(+)"
    full_sequence: str              # complete backbone DNA sequence
    size_bp: int
    features: list[BackboneFeature] = field(default_factory=list)

    # Expression cassette region (promoter start → terminator end)
    cassette_start: int = 0         # 1-based
    cassette_end: int = 0           # 1-based

    @property
    def promoter(self) -> Optional[BackboneFeature]:
        """The main expression promoter."""
        return next((f for f in self.features if f.feature_type == "promoter"), None)

    @property
    def rbs(self) -> Optional[BackboneFeature]:
        """The RBS / translation initiation element."""
        return next((f for f in self.features if f.feature_type == "rbs"), None)

    @property
    def terminator(self) -> Optional[BackboneFeature]:
        """The primary terminator (first one after the MCS)."""
        return next((f for f in self.features if f.feature_type == "terminator"), None)

    @property
    def n_terminal_features(self) -> list[BackboneFeature]:
        """Tags and cleavage sites before the insert (N-terminal fusions)."""
        return [f for f in self.features
                if f.feature_type in ("tag", "cleavage_site")
                and self.rbs and f.start > (self.rbs.end if self.rbs else 0)
                and self.terminator and f.end < (self.terminator.start if self.terminator else 9999999)]

    @property
    def c_terminal_features(self) -> list[BackboneFeature]:
        """Tags after the MCS region (C-terminal, e.g., C-terminal His)."""
        # These appear after the main MCS/insert but before the terminator
        # We identify them by position relative to other CDS features
        n_feats = self.n_terminal_features
        if not n_feats:
            return []
        last_n_end = max(f.end for f in n_feats)
        return [f for f in self.features
                if f.feature_type == "tag"
                and f.start > last_n_end + 50  # gap suggests MCS in between
                and self.terminator and f.end < self.terminator.start]

    def upstream_of_insert(self) -> list[BackboneFeature]:
        """
        Features that come before the insert cloning site, in order.
        Typically: promoter, lac operator, RBS, N-terminal tag(s), cleavage site(s).
        """
        # Everything before the MCS gap
        result = []
        for f in self.features:
            if f.feature_type in ("promoter", "lac_operator", "rbs", "tag", "cleavage_site"):
                result.append(f)
            elif f.feature_type == "terminator":
                break  # past the insert region
        return result

    def downstream_of_insert(self) -> list[BackboneFeature]:
        """
        Features that come after the insert cloning site, in order.
        Typically: C-terminal tag (optional), terminator.
        """
        result = []
        # Find the terminator and any C-terminal tags
        for f in reversed(self.features):
            if f.feature_type in ("terminator", "tag"):
                result.insert(0, f)
            elif f.feature_type in ("promoter", "rbs"):
                break
        return result

    def summary(self) -> str:
        lines = [f"Vector: {self.name} ({self.size_bp} bp)"]
        for f in self.features:
            lines.append(f"  {f.feature_type:20s} {f.label:30s} {f.start:>5d}..{f.end:>5d}  {f.length_bp:>4d}bp")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Classification: map GenBank annotations → our feature_type taxonomy
# ---------------------------------------------------------------------------

def _classify_cds(label: str, translation: str) -> tuple[str, str]:
    """
    Classify a CDS feature into our taxonomy.
    Returns (feature_type, normalized_label).
    """
    label_lower = label.lower()

    # His tags
    if "his" in label_lower and "x" in label_lower:
        return "tag", label
    if label_lower == "6xhis":
        return "tag", "6xHis"
    if label_lower == "8xhis":
        return "tag", "8xHis"
    if label_lower == "9xhis":
        return "tag", "9xHis"

    # S-Tag
    if "s-tag" in label_lower or label_lower == "s-tag":
        return "tag", "S-Tag"

    # Cleavage sites
    if "thrombin" in label_lower:
        return "cleavage_site", "Thrombin"
    if "factor xa" in label_lower:
        return "cleavage_site", "Factor_Xa"
    if "enterokinase" in label_lower:
        return "cleavage_site", "Enterokinase"
    if "tev" in label_lower:
        return "cleavage_site", "TEV"
    if "3c" in label_lower or "prescission" in label_lower:
        return "cleavage_site", "3C"

    # Solubility/fusion tags
    if label_lower == "trxa" or "thioredoxin" in label_lower:
        return "solubility_tag", "Trx"
    if label_lower == "gst":
        return "solubility_tag", "GST"
    if label_lower == "mbp":
        return "solubility_tag", "MBP"
    if label_lower == "sumo":
        return "solubility_tag", "SUMO"

    # Signal peptides
    if "signal" in label_lower or "pelb" in label_lower:
        return "signal_peptide", label

    # T7 tag / gene 10 leader
    if "t7 tag" in label_lower or "gene 10" in label_lower:
        return "leader", "T7_tag"

    # Start codon annotation (some vectors annotate ATG separately)
    if label_lower == "start codon" or (translation and translation == "M" and len(translation) == 1):
        return "start_codon", "ATG"

    # Puromycin resistance (in mammalian selection cassettes)
    if "puromycin" in label_lower:
        return "selection", label

    return "cds_other", label


def _classify_feature(genbank_feature) -> Optional[BackboneFeature]:
    """
    Convert a BioPython SeqFeature to a BackboneFeature, or None if
    it should be skipped.
    """
    ftype = genbank_feature.type
    quals = genbank_feature.qualifiers
    label = quals.get("label", [""])[0]
    note = quals.get("note", [""])[0]
    translation = quals.get("translation", [None])[0]

    start = int(genbank_feature.location.start) + 1  # 1-based
    end = int(genbank_feature.location.end)
    sequence = str(genbank_feature.extract(genbank_feature._seq) if hasattr(genbank_feature, '_seq') else "")

    if ftype == "source":
        return None

    if ftype not in _EXPR_FEATURE_TYPES:
        return None

    # Skip non-expression features
    if ftype == "CDS" and label in _SKIP_CDS_LABELS:
        return None
    if ftype == "promoter" and label in _SKIP_PROMOTER_LABELS:
        return None
    if ftype == "protein_bind" and label in _SKIP_BIND_LABELS:
        return None

    # Classify
    if ftype == "promoter":
        feature_type = "promoter"
        norm_label = label
    elif ftype == "RBS":
        feature_type = "rbs"
        norm_label = label or "RBS"
    elif ftype == "terminator":
        feature_type = "terminator"
        norm_label = label or "terminator"
    elif ftype == "protein_bind":
        if "lac" in label.lower() and "operator" in label.lower():
            feature_type = "lac_operator"
            norm_label = "lac_operator"
        elif "att" in label.lower():
            feature_type = "recombination_site"
            norm_label = label
        else:
            feature_type = "protein_bind"
            norm_label = label
    elif ftype == "CDS":
        feature_type, norm_label = _classify_cds(label, translation)
    else:
        return None

    return BackboneFeature(
        feature_type=feature_type,
        label=norm_label,
        sequence=sequence,
        translation=translation,
        start=start,
        end=end,
        genbank_type=ftype,
        note=note[:200] if note else "",
    )


# ---------------------------------------------------------------------------
# GenBank parsing
# ---------------------------------------------------------------------------

def parse_vector_genbank(filepath: str | Path) -> VectorRecord:
    """
    Parse a single GenBank file into a VectorRecord with classified features.
    """
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        record = SeqIO.read(str(filepath), "genbank")

    # Try to get vector name from KEYWORDS or filename
    name = ""
    if record.annotations.get("keywords"):
        kw = record.annotations["keywords"]
        if isinstance(kw, list):
            name = kw[0] if kw else ""
        else:
            name = str(kw)
    if not name or name.lower() == "." or name.lower() == "exported":
        # Fall back to filename
        name = Path(filepath).stem
        name = name.replace("+", " ").replace("%28", "(").replace("%29", ")")

    full_seq = str(record.seq).upper()

    # Store the parent sequence on each feature for extraction
    features = []
    for f in record.features:
        f._seq = record.seq  # attach for sequence extraction
        bf = _classify_feature(f)
        if bf is not None:
            # Extract the actual sequence from the record
            bf.sequence = str(f.extract(record.seq)).upper()
            features.append(bf)

    # Sort by position
    features.sort(key=lambda f: f.start)

    # Determine cassette boundaries
    cassette_start = 0
    cassette_end = 0
    promoters = [f for f in features if f.feature_type == "promoter"]
    terminators = [f for f in features if f.feature_type == "terminator"]
    if promoters:
        cassette_start = promoters[0].start
    if terminators:
        cassette_end = terminators[0].end

    return VectorRecord(
        name=name,
        full_sequence=full_seq,
        size_bp=len(full_seq),
        features=features,
        cassette_start=cassette_start,
        cassette_end=cassette_end,
    )


# ---------------------------------------------------------------------------
# Vector catalog: lazy-loaded registry of all Twist vectors
# ---------------------------------------------------------------------------

_VECTOR_CACHE: dict[str, VectorRecord] = {}


def _load_all_vectors() -> None:
    """Parse all GenBank files in the twist_vectors directory."""
    if _VECTOR_CACHE:
        return
    if not _VECTOR_DIR.exists():
        return
    for gb_file in sorted(_VECTOR_DIR.glob("*.gb")):
        try:
            vec = parse_vector_genbank(gb_file)
            _VECTOR_CACHE[vec.name] = vec
        except Exception:
            continue  # skip malformed files


def get_vector(name: str) -> Optional[VectorRecord]:
    """
    Look up a vector by name. Tries exact match first, then fuzzy matching
    (strips parentheses, normalizes whitespace/plus signs).
    """
    _load_all_vectors()

    # Exact match
    if name in _VECTOR_CACHE:
        return _VECTOR_CACHE[name]

    # Normalize and try again
    def normalize(s: str) -> str:
        return (s.lower()
                .replace("+", " ")
                .replace("(", "").replace(")", "")
                .replace("%28", "").replace("%29", "")
                .strip())

    target = normalize(name)
    for cached_name, vec in _VECTOR_CACHE.items():
        if normalize(cached_name) == target:
            return vec

    return None


def list_vectors() -> list[str]:
    """Return names of all available vectors."""
    _load_all_vectors()
    return sorted(_VECTOR_CACHE.keys())


def get_vector_expression_features(name: str) -> Optional[list[BackboneFeature]]:
    """
    Get the expression-region features for a named vector.
    Returns None if the vector is not found.
    """
    vec = get_vector(name)
    return vec.features if vec else None


# ---------------------------------------------------------------------------
# Helper: get upstream/downstream sequences for assembled view
# ---------------------------------------------------------------------------

def get_insert_context(
    vector_name: str,
    cloning_pair: Optional[tuple[str, str]] = None,
) -> Optional[dict]:
    """
    For a given catalog vector, return the features that flank the insert
    cloning site. This is used by assembled_graph() to reconstruct the
    full expression cassette.

    If cloning_pair is given (e.g., ("NcoI", "XhoI")), the feature list
    is trimmed to only include features that fall outside the cloning
    window — i.e., features the restriction digest doesn't remove.

    Returns a dict with:
        upstream: list[BackboneFeature]   — promoter → RBS → N-tags → cleavage
        downstream: list[BackboneFeature] — C-tags → terminator
        insert_site: (start, end)         — approximate cloning site position
        cloning_pairs: list               — available RE pairs for this vector
    """
    vec = get_vector(vector_name)
    if vec is None:
        return None

    # Load available cloning pairs from restriction_sites.json
    available_pairs = _get_cloning_pairs(vector_name)

    if cloning_pair:
        # Use restriction site positions to determine which features
        # are upstream vs downstream of the insert window
        upstream, downstream, insert_site = _features_for_cloning_pair(
            vec, cloning_pair[0], cloning_pair[1]
        )
    else:
        upstream = vec.upstream_of_insert()
        downstream = vec.downstream_of_insert()
        insert_start = max(f.end for f in upstream) + 1 if upstream else 0
        insert_end = min(f.start for f in downstream) - 1 if downstream else 0
        insert_site = (insert_start, insert_end)

    return {
        "upstream": upstream,
        "downstream": downstream,
        "insert_site": insert_site,
        "vector": vec,
        "cloning_pairs": available_pairs,
    }


# ---------------------------------------------------------------------------
# Restriction site / cloning pair support
# ---------------------------------------------------------------------------

_RESTRICTION_SITES: Optional[dict] = None


def _load_restriction_sites() -> dict:
    """Load restriction_sites.json from the twist_vectors directory."""
    global _RESTRICTION_SITES
    if _RESTRICTION_SITES is not None:
        return _RESTRICTION_SITES

    import json
    rs_path = _VECTOR_DIR / "restriction_sites.json"
    if rs_path.exists():
        with open(rs_path) as f:
            _RESTRICTION_SITES = json.load(f)
    else:
        _RESTRICTION_SITES = {}
    return _RESTRICTION_SITES


def _get_cloning_pairs(vector_name: str) -> list[list[str]]:
    """
    Get available restriction site pairs for a vector from the JSON data.
    Returns a list of [site1, site2] pairs.
    """
    data = _load_restriction_sites()

    # Search across all categories
    for category in ("expression_vectors", "cloning_vectors"):
        section = data.get(category, {})
        if isinstance(section, dict):
            # expression_vectors has subcategories (ecoli, mammalian, etc.)
            for sub_key, sub_val in section.items():
                if isinstance(sub_val, dict):
                    # Could be a vector entry or a subcategory
                    if "pairs" in sub_val:
                        # Direct vector entry (cloning_vectors style)
                        if _names_match(sub_key, vector_name):
                            return sub_val.get("pairs", [])
                    else:
                        # Subcategory (ecoli, mammalian, etc.)
                        for vec_name, vec_data in sub_val.items():
                            if isinstance(vec_data, dict) and _names_match(vec_name, vector_name):
                                return vec_data.get("pairs", [])
    return []


def _names_match(json_name: str, query_name: str) -> bool:
    """Fuzzy match vector names between JSON keys and GenBank/catalog names."""
    def norm(s: str) -> str:
        return (s.lower()
                .replace("+", "")
                .replace("(", "").replace(")", "")
                .replace(" ", "").replace("-", "")
                .strip())
    return norm(json_name) == norm(query_name)


def _features_for_cloning_pair(
    vec: VectorRecord,
    site_5prime: str,
    site_3prime: str,
) -> tuple[list[BackboneFeature], list[BackboneFeature], tuple[int, int]]:
    """
    Given a vector and a restriction site pair, determine which expression
    features are upstream and downstream of the insert window.

    The 5' site defines where the insert begins (everything before it
    that's part of the expression cassette is "upstream"). The 3' site
    defines where the insert ends (everything after it up to the
    terminator is "downstream").
    """
    # Find restriction site positions using BioPython
    try:
        import Bio.Restriction as RestrictionModule
        from Bio.Seq import Seq as BioSeq

        seq = BioSeq(vec.full_sequence)

        # Look up enzymes by name from Bio.Restriction module
        enzyme_5 = getattr(RestrictionModule, site_5prime, None)
        enzyme_3 = getattr(RestrictionModule, site_3prime, None)

        if enzyme_5 is None or enzyme_3 is None:
            # Fall back to default behavior
            return (vec.upstream_of_insert(), vec.downstream_of_insert(),
                    (0, 0))

        # Search for cut sites (returns 1-based positions)
        sites_5 = list(enzyme_5.search(seq, linear=False))
        sites_3 = list(enzyme_3.search(seq, linear=False))

        if not sites_5 or not sites_3:
            return (vec.upstream_of_insert(), vec.downstream_of_insert(),
                    (0, 0))

        # Pick the sites within the expression cassette region
        # (near the MCS, not elsewhere on the plasmid)
        cassette_region = range(vec.cassette_start, vec.cassette_end + 200)
        cut_5 = min((s for s in sites_5 if s in cassette_region),
                    default=sites_5[0])
        cut_3 = min((s for s in sites_3 if s in cassette_region),
                    default=sites_3[0])

        # Features upstream: start before the 5' cut site
        upstream = [f for f in vec.features
                    if f.end <= cut_5
                    and f.feature_type in ("promoter", "lac_operator", "rbs",
                                           "tag", "cleavage_site", "solubility_tag",
                                           "leader")]

        # Features downstream: start after the 3' cut site
        downstream = [f for f in vec.features
                      if f.start >= cut_3
                      and f.feature_type in ("tag", "terminator")]

        return (upstream, downstream, (cut_5, cut_3))

    except Exception:
        # Any import or search failure → fall back to default
        return (vec.upstream_of_insert(), vec.downstream_of_insert(),
                (0, 0))

"""
Local parts database — curated DNA and protein sequences for
common genetic elements used in E. coli expression constructs.

All DNA sequences are 5'->3'. Protein sequences are single-letter AA.
Sources cited in comments.
"""

# ---------------------------------------------------------------------------
# Promoters (DNA sequences)
# ---------------------------------------------------------------------------

PROMOTERS = {
    "T7": {
        "sequence": "TAATACGACTCACTATAGGG",
        "regulation": "inducible",
        "inducer": "IPTG",
        "notes": "T7 RNA polymerase promoter (phi10). Requires T7 RNAP (DE3 lysogen).",
    },
    "T7lac": {
        "sequence": "TAATACGACTCACTATAGGGGAATTGTGAGCGGATAACAATTCC",
        "regulation": "inducible",
        "inducer": "IPTG",
        "notes": "T7 promoter with lac operator for tighter repression.",
    },
    "tac": {
        "sequence": "TTGACAATTAATCATCGGCTCGTATAATGTGTGG",
        "regulation": "inducible",
        "inducer": "IPTG",
        "notes": "Hybrid trp/lac promoter. Strong, IPTG-inducible.",
    },
    "araBAD": {
        "sequence": "CTGACGCTTTTTATCGCAACTCTCTACTGTTTCTCCATACCCGTTTTTTTGGATGGAGTGAAACG",
        "regulation": "inducible",
        "inducer": "arabinose",
        "notes": "Arabinose-inducible pBAD promoter.",
    },
    "lacUV5": {
        "sequence": "TTTACACTTTATGCTTCCGGCTCGTATATTGTGTGGAATTGTGAGCGG",
        "regulation": "inducible",
        "inducer": "IPTG",
        "notes": "UV5 mutant of lac promoter. Catabolite repression insensitive.",
    },
    "J23100": {
        "sequence": "TTGACGGCTAGCTCAGTCCTAGGTACAGTGCTAGC",
        "regulation": "constitutive",
        "inducer": "",
        "notes": "Anderson constitutive promoter (strongest). iGEM BBa_J23100.",
    },
    "J23106": {
        "sequence": "TTTACGGCTAGCTCAGTCCTAGGTATAGTGCTAGC",
        "regulation": "constitutive",
        "inducer": "",
        "notes": "Anderson constitutive promoter (medium). iGEM BBa_J23106.",
    },
}


# ---------------------------------------------------------------------------
# RBS / translation initiation elements
# ---------------------------------------------------------------------------

RBS_LIBRARY = {
    # -- Standard RBSs (iGEM community parts) --
    "BBa_B0034": {
        "sequence": "AAAGAGGAGAAA",
        "relative_strength": 1.0,
        "tier": "high",
        "notes": "Community standard strong RBS. ~1.0 relative strength.",
    },
    "BBa_B0032": {
        "sequence": "TCACACAGGAAAG",
        "relative_strength": 0.3,
        "tier": "medium",
        "notes": "Medium-strength RBS. ~0.3 relative to B0034.",
    },
    "BBa_B0031": {
        "sequence": "TCACACAGGAAAC",
        "relative_strength": 0.07,
        "tier": "low",
        "notes": "Weak RBS. ~0.07 relative to B0034.",
    },
    "BBa_B0033": {
        "sequence": "TCACACAGGACTAG",
        "relative_strength": 0.01,
        "tier": "very_low",
        "notes": "Very weak RBS. ~0.01 relative to B0034.",
    },

    # -- Bicistronic designs (BCDs) from Mutalik et al. 2013 --
    # These use a leader peptide for context-insensitive translation initiation.
    # The sequences here include the leader ORF + internal RBS.
    "BCD2": {
        "sequence": (
            "ATGAACAATAATGATTTCGTAGATGATTTAAATATTGAAGGAGATATAAAT"
            "ATGAAAGAGGAGAAATACTAG"
        ),
        "relative_strength": 1.0,
        "tier": "high",
        "is_bcd": True,
        "notes": "BCD2 from Mutalik et al. 2013. Context-insensitive, strong.",
    },
    "BCD12": {
        "sequence": (
            "ATGAACAATAATGATTTCGTAGATGATTTAAATATTGAAGGAGATATAAAT"
            "ATGTCACACAGGAAAG"
        ),
        "relative_strength": 0.2,
        "tier": "medium",
        "is_bcd": True,
        "notes": "BCD12. Context-insensitive, medium strength.",
    },
    "BCD22": {
        "sequence": (
            "ATGAACAATAATGATTTCGTAGATGATTTAAATATTGAAGGAGATATAAAT"
            "ATGTCACACAGGACTAG"
        ),
        "relative_strength": 0.05,
        "tier": "low",
        "is_bcd": True,
        "notes": "BCD22. Context-insensitive, weak.",
    },
}

# Convenience: map tier names to default RBS choices
RBS_TIER_DEFAULTS = {
    "high": "BCD2",
    "medium": "BCD12",
    "low": "BCD22",
    "very_low": "BBa_B0033",
}


# ---------------------------------------------------------------------------
# Terminators
# ---------------------------------------------------------------------------

TERMINATORS = {
    "rrnB_T1": {
        "sequence": (
            "CAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGTCGGTGAACGCTCTC"
            "TGAGTAGTGACAATAAATCATAAATCATAAGTCAAGTAGGA"
        ),
        "type": "rho_independent",
        "notes": "rrnB T1 terminator. Strong, commonly used in pET vectors.",
    },
    "T7term": {
        "sequence": "CTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTG",
        "type": "rho_independent",
        "notes": "T7 transcription terminator.",
    },
    "BBa_B0015": {
        "sequence": (
            "CCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTT"
            "GTTTGTCGGTGAACGCTCTCTACTAGAGTCACACTGGCTCACCTTCGGGTGGGCCTTTCT"
            "GCGTTTATA"
        ),
        "type": "rho_independent",
        "notes": "Double terminator (B0010 + B0012). iGEM standard, very efficient.",
    },
}


# ---------------------------------------------------------------------------
# Purification & affinity tags (amino acid sequences)
# ---------------------------------------------------------------------------

PURIFICATION_TAGS = {
    "6xHis": {
        "protein_sequence": "HHHHHH",
        "notes": "Hexahistidine tag. IMAC purification (Ni-NTA, TALON).",
    },
    "8xHis": {
        "protein_sequence": "HHHHHHHH",
        "notes": "Octa-His. Stronger binding, useful for low-expression proteins.",
    },
    "Strep-II": {
        "protein_sequence": "WSHPQFEK",
        "notes": "Strep-tag II. Binds Strep-Tactin. Mild elution (desthiobiotin).",
    },
    "FLAG": {
        "protein_sequence": "DYKDDDDK",
        "notes": "FLAG tag. Antibody-based detection/purification.",
    },
    "HA": {
        "protein_sequence": "YPYDVPDYA",
        "notes": "HA epitope tag. Detection by anti-HA antibody.",
    },
    "Twin-Strep": {
        "protein_sequence": "WSHPQFEKGGGSGGGSGGSAWSHPQFEK",
        "notes": "Twin Strep-tag. Higher affinity than single Strep-II.",
    },
}


# ---------------------------------------------------------------------------
# Solubility / fusion tags (amino acid sequences)
# ---------------------------------------------------------------------------

SOLUBILITY_TAGS = {
    "MBP": {
        "protein_sequence": (
            "MKIEEGKLVIWINGDKGYNGLAEVGKKFEKDTGIKVTVEHPDKLEEKFPQVAATGDGPDI"
            "IFWAHDRFGGYAQSGLLAEITPDKAFQDKLYPFTWDAVRYNGKLIAYPIAVEALSLIYNK"
            "DLLPNPPKTWEEIPALDKELKAKGKSALMFNLQEPYFTWPLIAADGGYAFKYENGKYDIE"
            "KVCRDLAFPCEGIMVGTDNILAHSSFTREAVLEQFNNAAGAAPQPGELETDVEVVLKNG"
            "DLLTDRIKEIYNKAGQERDAKAAFEELSKIFDDAFKAGVS"
        ),
        "notes": "Maltose binding protein. Strong solubility enhancer. ~42 kDa.",
    },
    "GST": {
        "protein_sequence": (
            "MSPILGYWKIKGLVQPTRLLLEYLEEKYEEHLYERDEGDKWRNKKFELGLEFPNLPYYI"
            "DGDVKLTQSMAIIRYIADKHNMLGGCPKERAEISMLEGAVLDIRYGVSRIAYSKDFETLK"
            "VDFLSKLPEMLKMFEDRLCHKTYLNGDHVTHPDFMLYDALDVVLYMDPMCLDAFPKLVCF"
            "KKRIEAIPQIDKYLKSSKYIAWPLQGWQATFGGGDHPPK"
        ),
        "notes": "Glutathione S-transferase. ~26 kDa. Purification + solubility.",
    },
    "SUMO": {
        "protein_sequence": (
            "MSDSEVNQEAKPEVKPEVKPETHINLKVSDGSSEIFFKIKKTTPLRRLMEAFAKRQGKEM"
            "DSLRFLYDGIRIQADQTPEDLDMEDNDIIEAHREQIGG"
        ),
        "notes": "Small ubiquitin-like modifier. ~11 kDa. SUMO protease cleaves the fold.",
    },
    "Trx": {
        "protein_sequence": (
            "MSDKIIHLTDDSFDTDVLKADGAILVDFWAEWCGPCKMIAPILDEIADEYQGKLTVAKLN"
            "IDQNPGTAPKYGIRGIPTLLLFKNGEVAATKVGALSKGQLKEFLDANLA"
        ),
        "notes": "Thioredoxin. ~12 kDa. Promotes disulfide bond formation in cytoplasm.",
    },
}


# ---------------------------------------------------------------------------
# Cleavage sites (amino acid sequences)
# ---------------------------------------------------------------------------

CLEAVAGE_SITES = {
    "TEV": {
        "protein_sequence": "ENLYFQS",
        "cut_position": 6,  # cuts between Q and S
        "notes": "TEV protease site. Cleaves ENLYFQ/S. Works at 4-30C.",
    },
    "3C": {
        "protein_sequence": "LEVLFQGP",
        "cut_position": 6,  # cuts between Q and G
        "notes": "HRV 3C (PreScission) protease. Cleaves LEVLFQ/GP.",
    },
    "Factor_Xa": {
        "protein_sequence": "IEGR",
        "cut_position": 4,  # cuts after R
        "notes": "Factor Xa. Cleaves after IEGR. Less specific than TEV.",
    },
    "Thrombin": {
        "protein_sequence": "LVPRGS",
        "cut_position": 4,  # cuts between R and G
        "notes": "Thrombin. Cleaves LVPR/GS.",
    },
    "Enterokinase": {
        "protein_sequence": "DDDDK",
        "cut_position": 5,  # cuts after K
        "notes": "Enterokinase. Cleaves after DDDDK.",
    },
}


# ---------------------------------------------------------------------------
# Linkers (amino acid sequences)
# ---------------------------------------------------------------------------

LINKERS = {
    "GS_flexible": {
        "unit": "GGGGS",
        "default_repeats": 3,
        "notes": "(GGGGS)n flexible linker. n=1-6 typical.",
    },
    "rigid_EAAAK": {
        "unit": "EAAAK",
        "default_repeats": 3,
        "notes": "(EAAAK)n rigid alpha-helical linker.",
    },
    "short_GS": {
        "unit": "GS",
        "default_repeats": 1,
        "notes": "Minimal GS linker.",
    },
}


# ---------------------------------------------------------------------------
# Codon table: E. coli optimal codons
# ---------------------------------------------------------------------------

ECOLI_CODON_TABLE = {
    "F": ["TTT", "TTC"],        "L": ["CTG", "TTA", "TTG", "CTT", "CTC", "CTA"],
    "I": ["ATT", "ATC", "ATA"], "M": ["ATG"],
    "V": ["GTG", "GTT", "GTC", "GTA"],
    "S": ["AGC", "TCT", "TCC", "AGT", "TCA", "TCG"],
    "P": ["CCG", "CCA", "CCT", "CCC"],
    "T": ["ACC", "ACG", "ACT", "ACA"],
    "A": ["GCG", "GCC", "GCT", "GCA"],
    "Y": ["TAT", "TAC"],        "*": ["TAA", "TAG", "TGA"],
    "H": ["CAT", "CAC"],        "Q": ["CAG", "CAA"],
    "N": ["AAC", "AAT"],        "K": ["AAA", "AAG"],
    "D": ["GAT", "GAC"],        "E": ["GAA", "GAG"],
    "C": ["TGC", "TGT"],        "W": ["TGG"],
    "R": ["CGT", "CGC", "CGG", "CGA", "AGA", "AGG"],
    "G": ["GGC", "GGT", "GGG", "GGA"],
}

# Preferred (most frequent) codon per amino acid in E. coli K-12
ECOLI_PREFERRED_CODONS = {
    "F": "TTC", "L": "CTG", "I": "ATT", "M": "ATG", "V": "GTG",
    "S": "AGC", "P": "CCG", "T": "ACC", "A": "GCG", "Y": "TAT",
    "H": "CAT", "Q": "CAG", "N": "AAC", "K": "AAA", "D": "GAT",
    "E": "GAA", "C": "TGC", "W": "TGG", "R": "CGT", "G": "GGC",
    "*": "TAA",
}


# ---------------------------------------------------------------------------
# Golden Gate overhang sets (high-fidelity, from Potapov et al. 2018)
# ---------------------------------------------------------------------------

# 20 validated 4-bp overhangs with >95% correct ligation fidelity (BsaI)
# Source: Potapov et al., "Comprehensive Profiling of Four Base Overhang
#         Ligation Fidelity by T4 DNA Ligase", ACS Synth Biol, 2018.
HIGH_FIDELITY_OVERHANGS_BSAI = [
    "GGAG", "TACT", "AATG", "GCTT", "CGCT",
    "TGCC", "ACTA", "AGCG", "TTCG", "GATC",
    "AGTC", "CCAG", "TTAG", "GCTG", "CAGA",
    "ATCC", "TCGC", "CTGA", "GCAA", "TAGA",
]


# ---------------------------------------------------------------------------
# Twist Bioscience catalog (stock) vectors
# ---------------------------------------------------------------------------
# These are pre-made vectors available from Twist for clonal gene synthesis.
# When a catalog vector is selected, the user only needs to specify the insert
# (gene/ORF); the vector already provides promoter, RBS, tags, resistance, etc.
#
# "provides" lists the element types that the vector already contains and
# should NOT be duplicated in the synthesized insert.
# "cloning_method" describes how the insert is cloned into the vector.

TWIST_CATALOG_VECTORS = {
    # --- pET expression vectors (T7 system, for BL21(DE3)) ---
    "pET-21a(+)": {
        "description": "T7lac expression vector, optional C-terminal His-tag",
        "promoter": "T7lac",
        "resistance": "ampicillin",
        "ori": "pBR322",
        "host": "e_coli_bl21",
        "n_tags": [],
        "c_tags": ["6xHis"],   # optional, via reading frame
        "cleavage_sites": [],
        "cloning_sites": ["NdeI", "BamHI", "XhoI"],
        "provides": ["promoter", "rbs", "terminator"],
        "notes": "No N-terminal tag. C-terminal His optional (read through stop codon).",
        "size_bp": 5443,
    },
    "pET-28a(+)": {
        "description": "T7lac expression vector, N-terminal His-Tag + thrombin site",
        "promoter": "T7lac",
        "resistance": "kanamycin",
        "ori": "pBR322",
        "host": "e_coli_bl21",
        "n_tags": ["6xHis"],
        "c_tags": ["6xHis"],    # optional C-terminal His
        "cleavage_sites": ["Thrombin"],
        "cloning_sites": ["NcoI", "NdeI", "BamHI", "XhoI", "NotI"],
        "provides": ["promoter", "rbs", "n_tag", "cleavage", "terminator"],
        "notes": "N-terminal His + thrombin cleavage. Clone at NdeI/XhoI for N-His only.",
        "size_bp": 5369,
    },
    "pET-28b(+)": {
        "description": "T7lac expression vector, N-terminal His-Tag + thrombin site (alternate MCS)",
        "promoter": "T7lac",
        "resistance": "kanamycin",
        "ori": "pBR322",
        "host": "e_coli_bl21",
        "n_tags": ["6xHis"],
        "c_tags": ["6xHis"],
        "cleavage_sites": ["Thrombin"],
        "cloning_sites": ["NcoI", "NdeI", "BamHI", "XhoI"],
        "provides": ["promoter", "rbs", "n_tag", "cleavage", "terminator"],
        "notes": "Same as pET-28a with shifted MCS reading frame.",
        "size_bp": 5368,
    },
    "pET-32a(+)": {
        "description": "T7lac expression vector, Trx-His-S-tag-enterokinase fusion",
        "promoter": "T7lac",
        "resistance": "ampicillin",
        "ori": "pBR322",
        "host": "e_coli_bl21",
        "n_tags": ["6xHis"],
        "c_tags": [],
        "cleavage_sites": ["Enterokinase", "Thrombin"],
        "cloning_sites": ["NcoI", "BamHI", "XhoI", "NotI"],
        "provides": ["promoter", "rbs", "n_tag", "solubility_tag", "cleavage", "terminator"],
        "notes": "N-terminal Trx + His + S-tag + enterokinase site. Best for insoluble proteins.",
        "size_bp": 5900,
    },

    # --- pRSET vectors (T7, high copy) ---
    "pRSET-A": {
        "description": "T7 expression vector, N-terminal His-tag, high copy (pUC ori)",
        "promoter": "T7",
        "resistance": "ampicillin",
        "ori": "pUC",
        "host": "e_coli_bl21",
        "n_tags": ["6xHis"],
        "c_tags": [],
        "cleavage_sites": ["Enterokinase"],
        "cloning_sites": ["BamHI", "EcoRI", "PstI", "HindIII"],
        "provides": ["promoter", "rbs", "n_tag", "cleavage", "terminator"],
        "notes": "High copy number. Good for screening and small-scale expression.",
        "size_bp": 2897,
    },
    "pRSET-B": {
        "description": "T7 expression vector, N-terminal His-tag (alternate reading frame)",
        "promoter": "T7",
        "resistance": "ampicillin",
        "ori": "pUC",
        "host": "e_coli_bl21",
        "n_tags": ["6xHis"],
        "c_tags": [],
        "cleavage_sites": ["Enterokinase"],
        "cloning_sites": ["BamHI", "EcoRI", "PstI", "HindIII"],
        "provides": ["promoter", "rbs", "n_tag", "cleavage", "terminator"],
        "notes": "Same as pRSET-A, different reading frame for insert.",
        "size_bp": 2898,
    },
    "pRSET-C": {
        "description": "T7 expression vector, N-terminal His-tag (third reading frame)",
        "promoter": "T7",
        "resistance": "ampicillin",
        "ori": "pUC",
        "host": "e_coli_bl21",
        "n_tags": ["6xHis"],
        "c_tags": [],
        "cleavage_sites": ["Enterokinase"],
        "cloning_sites": ["BamHI", "EcoRI", "PstI", "HindIII"],
        "provides": ["promoter", "rbs", "n_tag", "cleavage", "terminator"],
        "notes": "Same as pRSET-A, third reading frame.",
        "size_bp": 2899,
    },

    # --- Cloning vectors ---
    "pUC19": {
        "description": "High-copy cloning vector with lacZ MCS",
        "promoter": "lac",
        "resistance": "ampicillin",
        "ori": "pUC",
        "host": "e_coli",
        "n_tags": [],
        "c_tags": [],
        "cleavage_sites": [],
        "cloning_sites": ["EcoRI", "SacI", "KpnI", "XmaI", "SmaI", "BamHI",
                          "XbaI", "AccI", "HincII", "SalI", "SbfI", "PstI",
                          "SphI", "HindIII"],
        "provides": [],
        "notes": "General purpose cloning vector. Blue/white screening. No expression elements.",
        "size_bp": 2686,
    },

    # ===================================================================
    # Twist cloning vectors (no expression elements)
    # ===================================================================

    "pTwist Amp": {
        "description": "Twist high-copy cloning vector, Ampicillin",
        "category": "cloning",
        "promoter": "",
        "resistance": "ampicillin",
        "ori": "pMB1",
        "host": "e_coli",
        "copy_number": "high",
        "n_tags": [],
        "c_tags": [],
        "cleavage_sites": [],
        "cloning_sites": [],
        "provides": [],
        "notes": "Minimal cloning vector. No expression elements. M13 fwd/rev priming sites flank insert.",
        "size_bp": 2221,
    },
    "pTwist Kan": {
        "description": "Twist medium-copy cloning vector, Kanamycin",
        "category": "cloning",
        "promoter": "",
        "resistance": "kanamycin",
        "ori": "pMB1",
        "host": "e_coli",
        "copy_number": "medium",
        "n_tags": [],
        "c_tags": [],
        "cleavage_sites": [],
        "cloning_sites": [],
        "provides": [],
        "notes": "Minimal cloning vector (Kan). M13 fwd/rev priming sites.",
        "size_bp": 2365,
    },

    # ===================================================================
    # Twist Gateway cloning vectors (ENTR)
    # ===================================================================

    "pTwist ENTR": {
        "description": "Gateway entry vector for shuttling into destination vectors",
        "category": "gateway",
        "promoter": "",
        "resistance": "kanamycin",
        "ori": "pMB1",
        "host": "e_coli",
        "copy_number": "high",
        "n_tags": [],
        "c_tags": [],
        "cleavage_sites": [],
        "cloning_sites": ["attL1", "attL2"],
        "provides": [],
        "notes": "Gateway entry clone. attL1/attL2 flanked insert. Use LR Clonase to shuttle into destination vectors. No start/stop codon provided.",
        "size_bp": 2365,
    },
    "pTwist ENTR Kozak": {
        "description": "Gateway entry vector with Kozak consensus for mammalian expression",
        "category": "gateway",
        "promoter": "",
        "resistance": "kanamycin",
        "ori": "pMB1",
        "host": "e_coli",
        "copy_number": "high",
        "n_tags": [],
        "c_tags": [],
        "cleavage_sites": [],
        "cloning_sites": ["attL1", "attL2"],
        "provides": ["kozak"],
        "notes": "Gateway entry clone with Kozak sequence (GCCACC) upstream of insert start. Optimized for mammalian expression destination vectors.",
        "size_bp": 2421,
    },

    # ===================================================================
    # Twist mammalian expression vectors — CMV promoter series
    # ===================================================================

    "pTwist CMV": {
        "description": "CMV-driven mammalian transient expression vector",
        "category": "mammalian_expression",
        "promoter": "CMV",
        "resistance": "ampicillin",
        "ori": "pUC",
        "host": "mammalian",
        "copy_number": "high",
        "n_tags": [],
        "c_tags": [],
        "cleavage_sites": [],
        "cloning_sites": [],
        "mammalian_selection": "",
        "provides": ["promoter", "terminator"],
        "notes": "CMV immediate-early promoter for high-level transient expression in mammalian cells. SV40 polyA signal. No mammalian selection marker.",
        "size_bp": 4831,
    },
    "pTwist CMV BetaGlobin": {
        "description": "CMV + beta-globin intron for enhanced mammalian expression",
        "category": "mammalian_expression",
        "promoter": "CMV",
        "resistance": "ampicillin",
        "ori": "pUC",
        "host": "mammalian",
        "copy_number": "high",
        "n_tags": [],
        "c_tags": [],
        "cleavage_sites": [],
        "cloning_sites": [],
        "mammalian_selection": "",
        "enhancers": ["beta-globin intron"],
        "provides": ["promoter", "terminator"],
        "notes": "CMV promoter + human beta-globin 5'UTR/intron enhances expression. hBG 3'UTR + bGH polyA downstream. No mammalian selection marker.",
        "size_bp": 4893,
    },
    "pTwist CMV BetaGlobin WPRE Neo": {
        "description": "CMV + beta-globin + WPRE, Neomycin/G418 mammalian selection",
        "category": "mammalian_expression",
        "promoter": "CMV",
        "resistance": "ampicillin",
        "ori": "pUC",
        "host": "mammalian",
        "copy_number": "extremely_high",
        "n_tags": [],
        "c_tags": [],
        "cleavage_sites": [],
        "cloning_sites": [],
        "mammalian_selection": "neomycin",
        "enhancers": ["beta-globin intron", "WPRE"],
        "provides": ["promoter", "terminator"],
        "notes": "CMV + beta-globin intron + WPRE for maximal expression. SV40-driven Neomycin resistance for G418 selection of stable mammalian clones.",
        "size_bp": 6737,
    },
    "pTwist CMV Hygro": {
        "description": "CMV-driven mammalian expression, Hygromycin selection",
        "category": "mammalian_expression",
        "promoter": "CMV",
        "resistance": "ampicillin",
        "ori": "pUC",
        "host": "mammalian",
        "copy_number": "extremely_high",
        "n_tags": [],
        "c_tags": [],
        "cleavage_sites": [],
        "cloning_sites": [],
        "mammalian_selection": "hygromycin",
        "provides": ["promoter", "terminator"],
        "notes": "CMV promoter + Hygromycin B resistance for stable mammalian clone selection. bGH polyA signal.",
        "size_bp": 6694,
    },
    "pTwist CMV Puro": {
        "description": "CMV-driven mammalian expression, Puromycin selection",
        "category": "mammalian_expression",
        "promoter": "CMV",
        "resistance": "ampicillin",
        "ori": "pUC",
        "host": "mammalian",
        "copy_number": "high",
        "n_tags": [],
        "c_tags": [],
        "cleavage_sites": [],
        "cloning_sites": [],
        "mammalian_selection": "puromycin",
        "provides": ["promoter", "terminator"],
        "notes": "CMV promoter with UbC-driven Puromycin resistance for rapid stable clone selection in mammalian cells.",
        "size_bp": 6633,
    },
    "pTwist CMV OriP": {
        "description": "CMV + OriP for episomal replication in EBV-positive cells",
        "category": "mammalian_expression",
        "promoter": "CMV",
        "resistance": "ampicillin",
        "ori": "pUC",
        "host": "mammalian",
        "copy_number": "high",
        "n_tags": [],
        "c_tags": [],
        "cleavage_sites": [],
        "cloning_sites": [],
        "mammalian_selection": "",
        "enhancers": ["OriP"],
        "provides": ["promoter", "terminator"],
        "notes": "CMV promoter + EBV OriP for episomal maintenance in cells expressing EBNA1 (e.g., HEK293E). Sustained expression without integration.",
        "size_bp": 4893,
    },

    # ===================================================================
    # Twist mammalian expression vectors — EF1-alpha promoter series
    # ===================================================================

    "pTwist EF1 Alpha": {
        "description": "EF1-alpha driven mammalian expression, medium-level",
        "category": "mammalian_expression",
        "promoter": "EF1a",
        "resistance": "ampicillin",
        "ori": "pUC",
        "host": "mammalian",
        "copy_number": "medium",
        "n_tags": [],
        "c_tags": [],
        "cleavage_sites": [],
        "cloning_sites": [],
        "mammalian_selection": "",
        "provides": ["promoter", "terminator"],
        "notes": "Human EF1-alpha promoter for moderate, sustained expression. Less prone to silencing than CMV in some cell types. bGH polyA signal.",
        "size_bp": 6633,
    },
    "pTwist EF1 Alpha Puro": {
        "description": "EF1-alpha driven mammalian expression, Puromycin selection",
        "category": "mammalian_expression",
        "promoter": "EF1a",
        "resistance": "ampicillin",
        "ori": "pUC",
        "host": "mammalian",
        "copy_number": "medium",
        "n_tags": [],
        "c_tags": [],
        "cleavage_sites": [],
        "cloning_sites": [],
        "mammalian_selection": "puromycin",
        "provides": ["promoter", "terminator"],
        "notes": "EF1-alpha promoter + Puromycin selection for stable mammalian clones. Preferred over CMV for long-term expression in primary cells.",
        "size_bp": 7200,
    },

    # ===================================================================
    # Twist lentiviral transfer vectors
    # ===================================================================

    "pTwist Lenti SFFV": {
        "description": "Lentiviral transfer vector, SFFV promoter",
        "category": "lentiviral",
        "promoter": "SFFV",
        "resistance": "ampicillin",
        "ori": "pUC",
        "host": "mammalian",
        "copy_number": "high",
        "n_tags": [],
        "c_tags": [],
        "cleavage_sites": [],
        "cloning_sites": [],
        "mammalian_selection": "",
        "lentiviral_elements": ["5' LTR", "3' SIN-LTR", "Psi", "RRE", "cPPT", "WPRE"],
        "provides": ["promoter", "terminator"],
        "notes": "3rd-generation lentiviral transfer vector. SFFV promoter for strong expression in hematopoietic cells. Self-inactivating LTR. Requires packaging plasmids for virus production.",
        "size_bp": 5683,
    },
    "pTwist Lenti SFFV Puro": {
        "description": "Lentiviral transfer vector, SFFV promoter, Puromycin selection",
        "category": "lentiviral",
        "promoter": "SFFV",
        "resistance": "ampicillin",
        "ori": "pUC",
        "host": "mammalian",
        "copy_number": "high",
        "n_tags": [],
        "c_tags": [],
        "cleavage_sites": [],
        "cloning_sites": [],
        "mammalian_selection": "puromycin",
        "lentiviral_elements": ["5' LTR", "3' SIN-LTR", "Psi", "RRE", "cPPT", "WPRE"],
        "provides": ["promoter", "terminator"],
        "notes": "3rd-gen lentiviral transfer vector with SFFV + Puromycin selection cassette. SIN-LTR design.",
        "size_bp": 7100,
    },
    "pTwist Lenti EF1 Alpha": {
        "description": "Lentiviral transfer vector, EF1-alpha promoter",
        "category": "lentiviral",
        "promoter": "EF1a",
        "resistance": "ampicillin",
        "ori": "pUC",
        "host": "mammalian",
        "copy_number": "high",
        "n_tags": [],
        "c_tags": [],
        "cleavage_sites": [],
        "cloning_sites": [],
        "mammalian_selection": "",
        "lentiviral_elements": ["5' LTR", "3' SIN-LTR", "Psi", "RRE", "cPPT", "WPRE"],
        "provides": ["promoter", "terminator"],
        "notes": "3rd-gen lentiviral transfer vector. EF1-alpha for broad mammalian expression with reduced silencing risk.",
        "size_bp": 6800,
    },
}

# ---------------------------------------------------------------------------
# Convenience filters
# ---------------------------------------------------------------------------

# E. coli expression vectors (T7, tac, etc.)
TWIST_ECOLI_EXPRESSION_VECTORS = {
    k: v for k, v in TWIST_CATALOG_VECTORS.items()
    if v.get("host", "").startswith("e_coli") and v.get("promoter")
}

# Mammalian expression vectors (CMV, EF1a, etc.)
TWIST_MAMMALIAN_VECTORS = {
    k: v for k, v in TWIST_CATALOG_VECTORS.items()
    if v.get("host") == "mammalian" and v.get("category") == "mammalian_expression"
}

# Lentiviral transfer vectors
TWIST_LENTIVIRAL_VECTORS = {
    k: v for k, v in TWIST_CATALOG_VECTORS.items()
    if v.get("category") == "lentiviral"
}

# Cloning-only vectors (no expression elements)
TWIST_CLONING_VECTORS = {
    k: v for k, v in TWIST_CATALOG_VECTORS.items()
    if v.get("category") in ("cloning", "gateway")
}

# All vectors grouped by category for UI
TWIST_VECTOR_CATEGORIES = {
    "E. coli Expression": TWIST_ECOLI_EXPRESSION_VECTORS,
    "Mammalian Expression": TWIST_MAMMALIAN_VECTORS,
    "Lentiviral": TWIST_LENTIVIRAL_VECTORS,
    "Cloning / Gateway": TWIST_CLONING_VECTORS,
}

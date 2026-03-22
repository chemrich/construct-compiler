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

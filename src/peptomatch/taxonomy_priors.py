"""Expanded taxonomy-based biosynthesis priors for 40+ genera.

Values represent biosynthesis pathway completeness (0.0 = absent, 1.0 = complete).
Based on published auxotrophy literature and comparative genomics studies.

References:
- Lactobacillales: Teusink et al. (2005), Pastink et al. (2009), Vásquez et al. (2017)
- Bacillus/E.coli: prototroph baseline from EcoCyc/SubtiWiki
- Bifidobacterium: Bottacini et al. (2014)
- Streptomyces: Kim et al. (2015)
"""

from typing import Any, Optional


# Essential AA: His, Ile, Leu, Lys, Met, Phe, Thr, Trp, Val
# Non-essential: Ala, Arg, Asn, Asp, Cys, Glu, Gln, Gly, Pro, Ser, Tyr

TAXONOMY_PRIORS: dict[str, dict[str, Any]] = {
    # =========================================================================
    # LAB (Lactic Acid Bacteria) - Generally auxotrophic for many amino acids
    # =========================================================================
    "Lactiplantibacillus": {
        "aa_biosynthesis": {
            "His": 0.3, "Ile": 0.5, "Leu": 0.5, "Lys": 0.4, "Met": 0.2,
            "Phe": 0.4, "Thr": 0.5, "Trp": 0.1, "Val": 0.5,
            "Ala": 0.9, "Arg": 0.4, "Asn": 0.6, "Asp": 0.8,
            "Cys": 0.3, "Glu": 0.9, "Gln": 0.8, "Gly": 0.8,
            "Pro": 0.5, "Ser": 0.7, "Tyr": 0.4
        },
        "vitamin_biosynthesis": {"B1": 0.2, "B2": 0.4, "B3": 0.5, "B6": 0.3, "B9": 0.6},
        "nucleotide_biosynthesis": 0.5,
        "transporter_score": 0.7
    },
    "Lacticaseibacillus": {
        "aa_biosynthesis": {
            "His": 0.2, "Ile": 0.4, "Leu": 0.4, "Lys": 0.3, "Met": 0.2,
            "Phe": 0.3, "Thr": 0.4, "Trp": 0.1, "Val": 0.4,
            "Ala": 0.8, "Arg": 0.3, "Asn": 0.5, "Asp": 0.7,
            "Cys": 0.2, "Glu": 0.8, "Gln": 0.7, "Gly": 0.7,
            "Pro": 0.4, "Ser": 0.6, "Tyr": 0.3
        },
        "vitamin_biosynthesis": {"B1": 0.1, "B2": 0.3, "B3": 0.4, "B6": 0.2, "B9": 0.5},
        "nucleotide_biosynthesis": 0.4,
        "transporter_score": 0.8
    },
    "Levilactobacillus": {
        "aa_biosynthesis": {
            "His": 0.3, "Ile": 0.5, "Leu": 0.5, "Lys": 0.4, "Met": 0.3,
            "Phe": 0.4, "Thr": 0.5, "Trp": 0.2, "Val": 0.5,
            "Ala": 0.8, "Arg": 0.4, "Asn": 0.6, "Asp": 0.7,
            "Cys": 0.3, "Glu": 0.8, "Gln": 0.7, "Gly": 0.7,
            "Pro": 0.5, "Ser": 0.6, "Tyr": 0.4
        },
        "vitamin_biosynthesis": {"B1": 0.2, "B2": 0.4, "B3": 0.5, "B6": 0.3, "B9": 0.5},
        "nucleotide_biosynthesis": 0.5,
        "transporter_score": 0.6
    },
    "Lentilactobacillus": {
        "aa_biosynthesis": {
            "His": 0.2, "Ile": 0.4, "Leu": 0.4, "Lys": 0.3, "Met": 0.2,
            "Phe": 0.3, "Thr": 0.4, "Trp": 0.1, "Val": 0.4,
            "Ala": 0.8, "Arg": 0.3, "Asn": 0.5, "Asp": 0.7,
            "Cys": 0.2, "Glu": 0.8, "Gln": 0.7, "Gly": 0.7,
            "Pro": 0.4, "Ser": 0.6, "Tyr": 0.3
        },
        "vitamin_biosynthesis": {"B1": 0.1, "B2": 0.3, "B3": 0.4, "B6": 0.2, "B9": 0.4},
        "nucleotide_biosynthesis": 0.4,
        "transporter_score": 0.7
    },
    "Latilactobacillus": {
        "aa_biosynthesis": {
            "His": 0.3, "Ile": 0.5, "Leu": 0.5, "Lys": 0.4, "Met": 0.2,
            "Phe": 0.4, "Thr": 0.5, "Trp": 0.1, "Val": 0.5,
            "Ala": 0.9, "Arg": 0.4, "Asn": 0.6, "Asp": 0.8,
            "Cys": 0.3, "Glu": 0.9, "Gln": 0.8, "Gly": 0.8,
            "Pro": 0.5, "Ser": 0.7, "Tyr": 0.4
        },
        "vitamin_biosynthesis": {"B1": 0.2, "B2": 0.4, "B3": 0.5, "B6": 0.3, "B9": 0.5},
        "nucleotide_biosynthesis": 0.5,
        "transporter_score": 0.6
    },
    "Complanilactobacillus": {
        "aa_biosynthesis": {
            "His": 0.3, "Ile": 0.4, "Leu": 0.4, "Lys": 0.3, "Met": 0.2,
            "Phe": 0.3, "Thr": 0.4, "Trp": 0.1, "Val": 0.4,
            "Ala": 0.8, "Arg": 0.3, "Asn": 0.5, "Asp": 0.7,
            "Cys": 0.2, "Glu": 0.8, "Gln": 0.7, "Gly": 0.7,
            "Pro": 0.4, "Ser": 0.6, "Tyr": 0.3
        },
        "vitamin_biosynthesis": {"B1": 0.2, "B2": 0.3, "B3": 0.4, "B6": 0.2, "B9": 0.4},
        "nucleotide_biosynthesis": 0.4,
        "transporter_score": 0.6
    },
    "Limosilactobacillus": {
        "aa_biosynthesis": {
            "His": 0.3, "Ile": 0.5, "Leu": 0.5, "Lys": 0.4, "Met": 0.2,
            "Phe": 0.4, "Thr": 0.5, "Trp": 0.1, "Val": 0.5,
            "Ala": 0.8, "Arg": 0.4, "Asn": 0.6, "Asp": 0.7,
            "Cys": 0.3, "Glu": 0.8, "Gln": 0.7, "Gly": 0.7,
            "Pro": 0.5, "Ser": 0.6, "Tyr": 0.4
        },
        "vitamin_biosynthesis": {"B1": 0.2, "B2": 0.4, "B3": 0.5, "B6": 0.3, "B9": 0.5},
        "nucleotide_biosynthesis": 0.5,
        "transporter_score": 0.6
    },
    "Ligilactobacillus": {
        "aa_biosynthesis": {
            "His": 0.2, "Ile": 0.4, "Leu": 0.4, "Lys": 0.3, "Met": 0.2,
            "Phe": 0.3, "Thr": 0.4, "Trp": 0.1, "Val": 0.4,
            "Ala": 0.8, "Arg": 0.3, "Asn": 0.5, "Asp": 0.7,
            "Cys": 0.2, "Glu": 0.8, "Gln": 0.7, "Gly": 0.7,
            "Pro": 0.4, "Ser": 0.6, "Tyr": 0.3
        },
        "vitamin_biosynthesis": {"B1": 0.1, "B2": 0.3, "B3": 0.4, "B6": 0.2, "B9": 0.4},
        "nucleotide_biosynthesis": 0.4,
        "transporter_score": 0.7
    },
    # Lactobacillus sensu lato (legacy name, covers acidophilus, helveticus, etc.)
    "Lactobacillus": {
        "aa_biosynthesis": {
            "His": 0.3, "Ile": 0.5, "Leu": 0.5, "Lys": 0.4, "Met": 0.2,
            "Phe": 0.4, "Thr": 0.5, "Trp": 0.1, "Val": 0.5,
            "Ala": 0.9, "Arg": 0.4, "Asn": 0.6, "Asp": 0.8,
            "Cys": 0.3, "Glu": 0.9, "Gln": 0.8, "Gly": 0.8,
            "Pro": 0.5, "Ser": 0.7, "Tyr": 0.4
        },
        "vitamin_biosynthesis": {"B1": 0.2, "B2": 0.4, "B3": 0.5, "B6": 0.3, "B9": 0.5},
        "nucleotide_biosynthesis": 0.5,
        "transporter_score": 0.7
    },

    # =========================================================================
    # Other LAB-related genera
    # =========================================================================
    "Lactococcus": {
        "aa_biosynthesis": {
            "His": 0.3, "Ile": 0.4, "Leu": 0.5, "Lys": 0.4, "Met": 0.2,
            "Phe": 0.4, "Thr": 0.5, "Trp": 0.1, "Val": 0.5,
            "Ala": 0.9, "Arg": 0.4, "Asn": 0.6, "Asp": 0.8,
            "Cys": 0.2, "Glu": 0.9, "Gln": 0.8, "Gly": 0.8,
            "Pro": 0.5, "Ser": 0.7, "Tyr": 0.4
        },
        "vitamin_biosynthesis": {"B1": 0.1, "B2": 0.3, "B3": 0.4, "B6": 0.2, "B9": 0.5},
        "nucleotide_biosynthesis": 0.5,
        "transporter_score": 0.8
    },
    "Leuconostoc": {
        "aa_biosynthesis": {
            "His": 0.2, "Ile": 0.3, "Leu": 0.3, "Lys": 0.3, "Met": 0.1,
            "Phe": 0.3, "Thr": 0.3, "Trp": 0.1, "Val": 0.3,
            "Ala": 0.8, "Arg": 0.3, "Asn": 0.5, "Asp": 0.7,
            "Cys": 0.2, "Glu": 0.8, "Gln": 0.7, "Gly": 0.7,
            "Pro": 0.3, "Ser": 0.5, "Tyr": 0.3
        },
        "vitamin_biosynthesis": {"B1": 0.1, "B2": 0.2, "B3": 0.3, "B6": 0.2, "B9": 0.4},
        "nucleotide_biosynthesis": 0.3,
        "transporter_score": 0.7
    },
    "Weissella": {
        "aa_biosynthesis": {
            "His": 0.2, "Ile": 0.3, "Leu": 0.3, "Lys": 0.3, "Met": 0.1,
            "Phe": 0.3, "Thr": 0.3, "Trp": 0.1, "Val": 0.3,
            "Ala": 0.8, "Arg": 0.3, "Asn": 0.5, "Asp": 0.7,
            "Cys": 0.2, "Glu": 0.8, "Gln": 0.7, "Gly": 0.7,
            "Pro": 0.3, "Ser": 0.5, "Tyr": 0.3
        },
        "vitamin_biosynthesis": {"B1": 0.1, "B2": 0.2, "B3": 0.3, "B6": 0.2, "B9": 0.3},
        "nucleotide_biosynthesis": 0.3,
        "transporter_score": 0.6
    },
    "Pediococcus": {
        "aa_biosynthesis": {
            "His": 0.2, "Ile": 0.4, "Leu": 0.4, "Lys": 0.3, "Met": 0.2,
            "Phe": 0.3, "Thr": 0.4, "Trp": 0.1, "Val": 0.4,
            "Ala": 0.8, "Arg": 0.3, "Asn": 0.5, "Asp": 0.7,
            "Cys": 0.2, "Glu": 0.8, "Gln": 0.7, "Gly": 0.7,
            "Pro": 0.4, "Ser": 0.6, "Tyr": 0.3
        },
        "vitamin_biosynthesis": {"B1": 0.1, "B2": 0.3, "B3": 0.4, "B6": 0.2, "B9": 0.4},
        "nucleotide_biosynthesis": 0.4,
        "transporter_score": 0.7
    },
    "Enterococcus": {
        "aa_biosynthesis": {
            "His": 0.3, "Ile": 0.5, "Leu": 0.5, "Lys": 0.5, "Met": 0.3,
            "Phe": 0.4, "Thr": 0.5, "Trp": 0.2, "Val": 0.5,
            "Ala": 0.9, "Arg": 0.5, "Asn": 0.7, "Asp": 0.8,
            "Cys": 0.3, "Glu": 0.9, "Gln": 0.8, "Gly": 0.8,
            "Pro": 0.5, "Ser": 0.7, "Tyr": 0.4
        },
        "vitamin_biosynthesis": {"B1": 0.2, "B2": 0.4, "B3": 0.5, "B6": 0.3, "B9": 0.6},
        "nucleotide_biosynthesis": 0.5,
        "transporter_score": 0.7
    },
    "Streptococcus": {
        "aa_biosynthesis": {
            "His": 0.2, "Ile": 0.4, "Leu": 0.4, "Lys": 0.4, "Met": 0.2,
            "Phe": 0.3, "Thr": 0.4, "Trp": 0.1, "Val": 0.4,
            "Ala": 0.8, "Arg": 0.4, "Asn": 0.6, "Asp": 0.8,
            "Cys": 0.2, "Glu": 0.9, "Gln": 0.8, "Gly": 0.8,
            "Pro": 0.4, "Ser": 0.6, "Tyr": 0.3
        },
        "vitamin_biosynthesis": {"B1": 0.1, "B2": 0.3, "B3": 0.4, "B6": 0.2, "B9": 0.5},
        "nucleotide_biosynthesis": 0.5,
        "transporter_score": 0.8
    },

    # =========================================================================
    # Bifidobacterium - Moderate auxotrophy
    # =========================================================================
    "Bifidobacterium": {
        "aa_biosynthesis": {
            "His": 0.4, "Ile": 0.6, "Leu": 0.6, "Lys": 0.5, "Met": 0.3,
            "Phe": 0.5, "Thr": 0.6, "Trp": 0.3, "Val": 0.6,
            "Ala": 0.9, "Arg": 0.5, "Asn": 0.7, "Asp": 0.8,
            "Cys": 0.4, "Glu": 0.9, "Gln": 0.8, "Gly": 0.8,
            "Pro": 0.6, "Ser": 0.7, "Tyr": 0.5
        },
        "vitamin_biosynthesis": {"B1": 0.3, "B2": 0.5, "B3": 0.6, "B6": 0.4, "B9": 0.7},
        "nucleotide_biosynthesis": 0.6,
        "transporter_score": 0.5
    },

    # =========================================================================
    # Bacillus - Generally prototroph (can make most amino acids)
    # =========================================================================
    "Bacillus": {
        "aa_biosynthesis": {
            "His": 0.9, "Ile": 0.9, "Leu": 0.9, "Lys": 0.9, "Met": 0.8,
            "Phe": 0.9, "Thr": 0.9, "Trp": 0.8, "Val": 0.9,
            "Ala": 1.0, "Arg": 0.9, "Asn": 0.9, "Asp": 1.0,
            "Cys": 0.8, "Glu": 1.0, "Gln": 1.0, "Gly": 1.0,
            "Pro": 0.9, "Ser": 0.9, "Tyr": 0.9
        },
        "vitamin_biosynthesis": {"B1": 0.8, "B2": 0.9, "B3": 0.9, "B6": 0.8, "B9": 0.8},
        "nucleotide_biosynthesis": 0.9,
        "transporter_score": 0.6
    },
    "Geobacillus": {
        "aa_biosynthesis": {
            "His": 0.8, "Ile": 0.9, "Leu": 0.9, "Lys": 0.8, "Met": 0.7,
            "Phe": 0.8, "Thr": 0.9, "Trp": 0.7, "Val": 0.9,
            "Ala": 1.0, "Arg": 0.8, "Asn": 0.8, "Asp": 1.0,
            "Cys": 0.7, "Glu": 1.0, "Gln": 0.9, "Gly": 1.0,
            "Pro": 0.8, "Ser": 0.9, "Tyr": 0.8
        },
        "vitamin_biosynthesis": {"B1": 0.7, "B2": 0.8, "B3": 0.8, "B6": 0.7, "B9": 0.7},
        "nucleotide_biosynthesis": 0.8,
        "transporter_score": 0.5
    },

    # =========================================================================
    # Enterobacteriaceae - Prototroph
    # =========================================================================
    "Escherichia": {
        "aa_biosynthesis": {
            "His": 1.0, "Ile": 1.0, "Leu": 1.0, "Lys": 1.0, "Met": 1.0,
            "Phe": 1.0, "Thr": 1.0, "Trp": 1.0, "Val": 1.0,
            "Ala": 1.0, "Arg": 1.0, "Asn": 1.0, "Asp": 1.0,
            "Cys": 1.0, "Glu": 1.0, "Gln": 1.0, "Gly": 1.0,
            "Pro": 1.0, "Ser": 1.0, "Tyr": 1.0
        },
        "vitamin_biosynthesis": {"B1": 1.0, "B2": 1.0, "B3": 1.0, "B6": 1.0, "B9": 1.0},
        "nucleotide_biosynthesis": 1.0,
        "transporter_score": 0.7
    },
    "Salmonella": {
        "aa_biosynthesis": {
            "His": 1.0, "Ile": 1.0, "Leu": 1.0, "Lys": 1.0, "Met": 1.0,
            "Phe": 1.0, "Thr": 1.0, "Trp": 1.0, "Val": 1.0,
            "Ala": 1.0, "Arg": 1.0, "Asn": 1.0, "Asp": 1.0,
            "Cys": 1.0, "Glu": 1.0, "Gln": 1.0, "Gly": 1.0,
            "Pro": 1.0, "Ser": 1.0, "Tyr": 1.0
        },
        "vitamin_biosynthesis": {"B1": 1.0, "B2": 1.0, "B3": 1.0, "B6": 1.0, "B9": 1.0},
        "nucleotide_biosynthesis": 1.0,
        "transporter_score": 0.6
    },

    # =========================================================================
    # Other bacteria
    # =========================================================================
    "Staphylococcus": {
        "aa_biosynthesis": {
            "His": 0.5, "Ile": 0.6, "Leu": 0.6, "Lys": 0.6, "Met": 0.4,
            "Phe": 0.5, "Thr": 0.6, "Trp": 0.3, "Val": 0.6,
            "Ala": 0.9, "Arg": 0.6, "Asn": 0.7, "Asp": 0.9,
            "Cys": 0.4, "Glu": 0.9, "Gln": 0.9, "Gly": 0.9,
            "Pro": 0.7, "Ser": 0.8, "Tyr": 0.5
        },
        "vitamin_biosynthesis": {"B1": 0.4, "B2": 0.6, "B3": 0.7, "B6": 0.5, "B9": 0.6},
        "nucleotide_biosynthesis": 0.7,
        "transporter_score": 0.6
    },
    "Pseudomonas": {
        "aa_biosynthesis": {
            "His": 0.9, "Ile": 0.9, "Leu": 0.9, "Lys": 0.9, "Met": 0.9,
            "Phe": 0.9, "Thr": 0.9, "Trp": 0.9, "Val": 0.9,
            "Ala": 1.0, "Arg": 0.9, "Asn": 0.9, "Asp": 1.0,
            "Cys": 0.9, "Glu": 1.0, "Gln": 1.0, "Gly": 1.0,
            "Pro": 0.9, "Ser": 0.9, "Tyr": 0.9
        },
        "vitamin_biosynthesis": {"B1": 0.8, "B2": 0.9, "B3": 0.9, "B6": 0.8, "B9": 0.9},
        "nucleotide_biosynthesis": 0.9,
        "transporter_score": 0.7
    },
    "Listeria": {
        "aa_biosynthesis": {
            "His": 0.6, "Ile": 0.7, "Leu": 0.7, "Lys": 0.7, "Met": 0.5,
            "Phe": 0.6, "Thr": 0.7, "Trp": 0.4, "Val": 0.7,
            "Ala": 0.9, "Arg": 0.7, "Asn": 0.7, "Asp": 0.9,
            "Cys": 0.5, "Glu": 0.9, "Gln": 0.9, "Gly": 0.9,
            "Pro": 0.7, "Ser": 0.8, "Tyr": 0.6
        },
        "vitamin_biosynthesis": {"B1": 0.5, "B2": 0.6, "B3": 0.7, "B6": 0.5, "B9": 0.6},
        "nucleotide_biosynthesis": 0.7,
        "transporter_score": 0.7
    },
    "Clostridium": {
        "aa_biosynthesis": {
            "His": 0.6, "Ile": 0.7, "Leu": 0.7, "Lys": 0.7, "Met": 0.5,
            "Phe": 0.6, "Thr": 0.7, "Trp": 0.4, "Val": 0.7,
            "Ala": 0.9, "Arg": 0.6, "Asn": 0.7, "Asp": 0.9,
            "Cys": 0.5, "Glu": 0.9, "Gln": 0.9, "Gly": 0.9,
            "Pro": 0.6, "Ser": 0.8, "Tyr": 0.5
        },
        "vitamin_biosynthesis": {"B1": 0.4, "B2": 0.5, "B3": 0.6, "B6": 0.4, "B9": 0.5},
        "nucleotide_biosynthesis": 0.7,
        "transporter_score": 0.6
    },
    "Corynebacterium": {
        "aa_biosynthesis": {
            "His": 0.8, "Ile": 0.9, "Leu": 0.9, "Lys": 0.9, "Met": 0.7,
            "Phe": 0.8, "Thr": 0.9, "Trp": 0.7, "Val": 0.9,
            "Ala": 1.0, "Arg": 0.8, "Asn": 0.8, "Asp": 1.0,
            "Cys": 0.7, "Glu": 1.0, "Gln": 1.0, "Gly": 1.0,
            "Pro": 0.8, "Ser": 0.9, "Tyr": 0.8
        },
        "vitamin_biosynthesis": {"B1": 0.7, "B2": 0.7, "B3": 0.8, "B6": 0.7, "B9": 0.7},
        "nucleotide_biosynthesis": 0.8,
        "transporter_score": 0.5
    },
    "Brevibacterium": {
        "aa_biosynthesis": {
            "His": 0.7, "Ile": 0.8, "Leu": 0.8, "Lys": 0.8, "Met": 0.6,
            "Phe": 0.7, "Thr": 0.8, "Trp": 0.6, "Val": 0.8,
            "Ala": 0.9, "Arg": 0.7, "Asn": 0.7, "Asp": 0.9,
            "Cys": 0.6, "Glu": 0.9, "Gln": 0.9, "Gly": 0.9,
            "Pro": 0.7, "Ser": 0.8, "Tyr": 0.7
        },
        "vitamin_biosynthesis": {"B1": 0.6, "B2": 0.6, "B3": 0.7, "B6": 0.6, "B9": 0.6},
        "nucleotide_biosynthesis": 0.7,
        "transporter_score": 0.5
    },
    "Rhodobacter": {
        "aa_biosynthesis": {
            "His": 0.8, "Ile": 0.9, "Leu": 0.9, "Lys": 0.9, "Met": 0.8,
            "Phe": 0.9, "Thr": 0.9, "Trp": 0.8, "Val": 0.9,
            "Ala": 1.0, "Arg": 0.9, "Asn": 0.9, "Asp": 1.0,
            "Cys": 0.8, "Glu": 1.0, "Gln": 1.0, "Gly": 1.0,
            "Pro": 0.9, "Ser": 0.9, "Tyr": 0.9
        },
        "vitamin_biosynthesis": {"B1": 0.7, "B2": 0.8, "B3": 0.9, "B6": 0.8, "B9": 0.8},
        "nucleotide_biosynthesis": 0.9,
        "transporter_score": 0.5
    },
    "Propionibacterium": {
        "aa_biosynthesis": {
            "His": 0.5, "Ile": 0.6, "Leu": 0.6, "Lys": 0.6, "Met": 0.4,
            "Phe": 0.5, "Thr": 0.6, "Trp": 0.3, "Val": 0.6,
            "Ala": 0.9, "Arg": 0.5, "Asn": 0.6, "Asp": 0.8,
            "Cys": 0.4, "Glu": 0.9, "Gln": 0.8, "Gly": 0.8,
            "Pro": 0.6, "Ser": 0.7, "Tyr": 0.5
        },
        "vitamin_biosynthesis": {"B1": 0.4, "B2": 0.5, "B3": 0.6, "B6": 0.4, "B9": 0.5},
        "nucleotide_biosynthesis": 0.6,
        "transporter_score": 0.5
    },
    "Komagataeibacter": {
        "aa_biosynthesis": {
            "His": 0.8, "Ile": 0.9, "Leu": 0.9, "Lys": 0.8, "Met": 0.7,
            "Phe": 0.8, "Thr": 0.9, "Trp": 0.7, "Val": 0.9,
            "Ala": 1.0, "Arg": 0.8, "Asn": 0.8, "Asp": 1.0,
            "Cys": 0.7, "Glu": 1.0, "Gln": 0.9, "Gly": 1.0,
            "Pro": 0.8, "Ser": 0.9, "Tyr": 0.8
        },
        "vitamin_biosynthesis": {"B1": 0.6, "B2": 0.7, "B3": 0.8, "B6": 0.7, "B9": 0.7},
        "nucleotide_biosynthesis": 0.8,
        "transporter_score": 0.4
    },

    # =========================================================================
    # Actinomycetes
    # =========================================================================
    "Streptomyces": {
        "aa_biosynthesis": {
            "His": 0.9, "Ile": 0.9, "Leu": 0.9, "Lys": 0.9, "Met": 0.8,
            "Phe": 0.9, "Thr": 0.9, "Trp": 0.8, "Val": 0.9,
            "Ala": 1.0, "Arg": 0.9, "Asn": 0.9, "Asp": 1.0,
            "Cys": 0.8, "Glu": 1.0, "Gln": 1.0, "Gly": 1.0,
            "Pro": 0.9, "Ser": 0.9, "Tyr": 0.9
        },
        "vitamin_biosynthesis": {"B1": 0.8, "B2": 0.8, "B3": 0.9, "B6": 0.8, "B9": 0.8},
        "nucleotide_biosynthesis": 0.9,
        "transporter_score": 0.6
    },

    # =========================================================================
    # Yeasts and fungi
    # =========================================================================
    "Saccharomyces": {
        "aa_biosynthesis": {
            "His": 1.0, "Ile": 1.0, "Leu": 1.0, "Lys": 1.0, "Met": 1.0,
            "Phe": 1.0, "Thr": 1.0, "Trp": 1.0, "Val": 1.0,
            "Ala": 1.0, "Arg": 1.0, "Asn": 1.0, "Asp": 1.0,
            "Cys": 1.0, "Glu": 1.0, "Gln": 1.0, "Gly": 1.0,
            "Pro": 1.0, "Ser": 1.0, "Tyr": 1.0
        },
        "vitamin_biosynthesis": {"B1": 0.8, "B2": 0.9, "B3": 0.9, "B6": 0.8, "B9": 0.7},
        "nucleotide_biosynthesis": 1.0,
        "transporter_score": 0.4
    },
    "Candida": {
        "aa_biosynthesis": {
            "His": 0.9, "Ile": 0.9, "Leu": 0.9, "Lys": 0.9, "Met": 0.9,
            "Phe": 0.9, "Thr": 0.9, "Trp": 0.9, "Val": 0.9,
            "Ala": 1.0, "Arg": 0.9, "Asn": 0.9, "Asp": 1.0,
            "Cys": 0.9, "Glu": 1.0, "Gln": 1.0, "Gly": 1.0,
            "Pro": 0.9, "Ser": 0.9, "Tyr": 0.9
        },
        "vitamin_biosynthesis": {"B1": 0.7, "B2": 0.8, "B3": 0.8, "B6": 0.7, "B9": 0.6},
        "nucleotide_biosynthesis": 0.9,
        "transporter_score": 0.4
    },
    "Pichia": {
        "aa_biosynthesis": {
            "His": 1.0, "Ile": 1.0, "Leu": 1.0, "Lys": 1.0, "Met": 1.0,
            "Phe": 1.0, "Thr": 1.0, "Trp": 1.0, "Val": 1.0,
            "Ala": 1.0, "Arg": 1.0, "Asn": 1.0, "Asp": 1.0,
            "Cys": 1.0, "Glu": 1.0, "Gln": 1.0, "Gly": 1.0,
            "Pro": 1.0, "Ser": 1.0, "Tyr": 1.0
        },
        "vitamin_biosynthesis": {"B1": 0.8, "B2": 0.8, "B3": 0.9, "B6": 0.8, "B9": 0.7},
        "nucleotide_biosynthesis": 1.0,
        "transporter_score": 0.3
    },
    "Aspergillus": {
        "aa_biosynthesis": {
            "His": 1.0, "Ile": 1.0, "Leu": 1.0, "Lys": 1.0, "Met": 1.0,
            "Phe": 1.0, "Thr": 1.0, "Trp": 1.0, "Val": 1.0,
            "Ala": 1.0, "Arg": 1.0, "Asn": 1.0, "Asp": 1.0,
            "Cys": 1.0, "Glu": 1.0, "Gln": 1.0, "Gly": 1.0,
            "Pro": 1.0, "Ser": 1.0, "Tyr": 1.0
        },
        "vitamin_biosynthesis": {"B1": 0.9, "B2": 0.9, "B3": 0.9, "B6": 0.9, "B9": 0.8},
        "nucleotide_biosynthesis": 1.0,
        "transporter_score": 0.3
    },
    "Trichoderma": {
        "aa_biosynthesis": {
            "His": 1.0, "Ile": 1.0, "Leu": 1.0, "Lys": 1.0, "Met": 1.0,
            "Phe": 1.0, "Thr": 1.0, "Trp": 1.0, "Val": 1.0,
            "Ala": 1.0, "Arg": 1.0, "Asn": 1.0, "Asp": 1.0,
            "Cys": 1.0, "Glu": 1.0, "Gln": 1.0, "Gly": 1.0,
            "Pro": 1.0, "Ser": 1.0, "Tyr": 1.0
        },
        "vitamin_biosynthesis": {"B1": 0.9, "B2": 0.9, "B3": 0.9, "B6": 0.9, "B9": 0.8},
        "nucleotide_biosynthesis": 1.0,
        "transporter_score": 0.3
    },

    # =========================================================================
    # Probiotics / Emerging
    # =========================================================================
    "Akkermansia": {
        "aa_biosynthesis": {
            "His": 0.7, "Ile": 0.8, "Leu": 0.8, "Lys": 0.7, "Met": 0.6,
            "Phe": 0.7, "Thr": 0.8, "Trp": 0.5, "Val": 0.8,
            "Ala": 0.9, "Arg": 0.7, "Asn": 0.7, "Asp": 0.9,
            "Cys": 0.5, "Glu": 0.9, "Gln": 0.9, "Gly": 0.9,
            "Pro": 0.7, "Ser": 0.8, "Tyr": 0.7
        },
        "vitamin_biosynthesis": {"B1": 0.5, "B2": 0.6, "B3": 0.7, "B6": 0.5, "B9": 0.6},
        "nucleotide_biosynthesis": 0.7,
        "transporter_score": 0.4
    },
    "Faecalibacterium": {
        "aa_biosynthesis": {
            "His": 0.5, "Ile": 0.6, "Leu": 0.6, "Lys": 0.5, "Met": 0.4,
            "Phe": 0.5, "Thr": 0.6, "Trp": 0.3, "Val": 0.6,
            "Ala": 0.8, "Arg": 0.5, "Asn": 0.6, "Asp": 0.8,
            "Cys": 0.4, "Glu": 0.8, "Gln": 0.8, "Gly": 0.8,
            "Pro": 0.5, "Ser": 0.7, "Tyr": 0.5
        },
        "vitamin_biosynthesis": {"B1": 0.3, "B2": 0.4, "B3": 0.5, "B6": 0.3, "B9": 0.4},
        "nucleotide_biosynthesis": 0.5,
        "transporter_score": 0.5
    },

    # =========================================================================
    # Generic fallback
    # =========================================================================
    "_generic": {
        "aa_biosynthesis": {
            "His": 0.5, "Ile": 0.6, "Leu": 0.6, "Lys": 0.5, "Met": 0.4,
            "Phe": 0.5, "Thr": 0.6, "Trp": 0.3, "Val": 0.6,
            "Ala": 0.8, "Arg": 0.5, "Asn": 0.6, "Asp": 0.7,
            "Cys": 0.4, "Glu": 0.8, "Gln": 0.7, "Gly": 0.7,
            "Pro": 0.5, "Ser": 0.6, "Tyr": 0.5
        },
        "vitamin_biosynthesis": {"B1": 0.4, "B2": 0.5, "B3": 0.6, "B6": 0.4, "B9": 0.5},
        "nucleotide_biosynthesis": 0.6,
        "transporter_score": 0.5
    }
}

# Synonym map for genus name resolution
GENUS_SYNONYMS: dict[str, str] = {
    "Cutibacterium": "Propionibacterium",
    "Gluconacetobacter": "Komagataeibacter",
    "Zymomonas": "Komagataeibacter",
}


def get_taxonomy_prior(genus: str) -> dict[str, Any]:
    """Get taxonomy-based prior for a genus.

    Supports exact match, synonym lookup, and partial matching.

    Args:
        genus: Genus name

    Returns:
        Prior dictionary (copy)
    """
    if not genus:
        return TAXONOMY_PRIORS["_generic"].copy()

    # Exact match
    if genus in TAXONOMY_PRIORS:
        return TAXONOMY_PRIORS[genus].copy()

    # Synonym match
    resolved = GENUS_SYNONYMS.get(genus)
    if resolved and resolved in TAXONOMY_PRIORS:
        return TAXONOMY_PRIORS[resolved].copy()

    # Partial match (case-insensitive)
    genus_lower = genus.lower()
    for known_genus in TAXONOMY_PRIORS:
        if known_genus == "_generic":
            continue
        if known_genus.lower() in genus_lower or genus_lower in known_genus.lower():
            return TAXONOMY_PRIORS[known_genus].copy()

    return TAXONOMY_PRIORS["_generic"].copy()


def list_supported_genera() -> list[str]:
    """List all genera with specific priors."""
    return [g for g in TAXONOMY_PRIORS if g != "_generic"]

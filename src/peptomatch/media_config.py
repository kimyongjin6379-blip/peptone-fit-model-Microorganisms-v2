"""Base media configurations and peptone responsibility calculations.

Peptone is added at 2.0-2.5% in media, replacing nitrogen sources.
The base medium provides sugar, salts, etc. The scoring weight for each
nutrient should reflect what fraction the PEPTONE contributes vs the base medium.
"""

MEDIA_CONFIGS = {
    "MRS_2.0": {
        "display_name": "MRS (2.0%, Peptone+Beef extract 대체)",
        "description": "Yeast extract 유지, Peptone+Beef extract를 Sempio 펩톤으로 대체",
        "peptone_g_per_L": 20.0,
        "default_for": ["LAB"],  # auto-select for these strain types
        "base_provides": {
            "glucose_g": 20.0,
            "yeast_extract_g": 5.0,
            "Na_acetate_g": 5.0,
            "ammonium_citrate_g": 2.0,
            "K2HPO4_g": 2.0,
            "MgSO4_g": 0.2,
            "MnSO4_g": 0.05,
            "Tween80_mL": 1.0,
        },
        "responsibility": {
            "faa_abundance": 1.0,       # peptone is sole complex nitrogen source
            "taa_abundance": 1.0,
            "mw_low": 1.0,
            "mw_medium": 1.0,
            "mw_high": 1.0,
            "aa_biosynthesis_gap": 1.0,
            "transporter_bonus": 1.0,
            "vitamin_b": 0.3,           # YE still provides most vitamins
            "nucleotides": 0.3,         # YE still provides nucleotides
            "sugar": 0.15,              # base has 20g glucose, peptone adds ~3-4g
            "mineral": 0.25,            # base has K2HPO4, MgSO4, MnSO4
            "orgacid": 0.20,            # base has Na acetate 5g
            "nitrogen_quality": 0.8,
        },
    },
    "MRS_2.5": {
        "display_name": "MRS (2.5%, 전체 질소원 대체)",
        "description": "Peptone+Beef extract+Yeast extract 모두 Sempio 펩톤으로 대체",
        "peptone_g_per_L": 25.0,
        "default_for": ["LAB"],
        "base_provides": {
            "glucose_g": 20.0,
            "yeast_extract_g": 0.0,  # removed!
            "Na_acetate_g": 5.0,
            "ammonium_citrate_g": 2.0,
            "K2HPO4_g": 2.0,
            "MgSO4_g": 0.2,
            "MnSO4_g": 0.05,
            "Tween80_mL": 1.0,
        },
        "responsibility": {
            "faa_abundance": 1.0,
            "taa_abundance": 1.0,
            "mw_low": 1.0,
            "mw_medium": 1.0,
            "mw_high": 1.0,
            "aa_biosynthesis_gap": 1.0,
            "transporter_bonus": 1.0,
            "vitamin_b": 1.0,           # YE removed → peptone must supply vitamins!
            "nucleotides": 0.8,         # YE removed → peptone is main source
            "sugar": 0.18,              # slightly more at 2.5%
            "mineral": 0.30,            # slightly more contribution
            "orgacid": 0.15,
            "nitrogen_quality": 0.8,
        },
    },
    "LB_2.0": {
        "display_name": "LB (2.0%, Tryptone+YE 대체)",
        "description": "Tryptone+Yeast extract를 Sempio 펩톤으로 대체",
        "peptone_g_per_L": 20.0,
        "default_for": ["Enterobacteriaceae"],
        "base_provides": {
            "glucose_g": 0.0,   # LB has no glucose!
            "yeast_extract_g": 0.0,
            "NaCl_g": 10.0,
        },
        "responsibility": {
            "faa_abundance": 1.0,
            "taa_abundance": 1.0,
            "mw_low": 1.0,
            "mw_medium": 1.0,
            "mw_high": 1.0,
            "aa_biosynthesis_gap": 1.0,
            "transporter_bonus": 1.0,
            "vitamin_b": 1.0,           # YE replaced
            "nucleotides": 0.8,
            "sugar": 0.8,              # LB has NO glucose → peptone sugar matters!
            "mineral": 0.5,            # only NaCl in base, peptone provides K/Mg/Ca
            "orgacid": 0.3,
            "nitrogen_quality": 0.8,
        },
    },
    "TSB_2.5": {
        "display_name": "TSB (2.5%, Tryptone+Soytone 대체)",
        "description": "Tryptone+Soytone을 Sempio 펩톤으로 대체",
        "peptone_g_per_L": 25.0,
        "default_for": ["Bacillus"],
        "base_provides": {
            "glucose_g": 2.5,
            "yeast_extract_g": 0.0,
            "NaCl_g": 5.0,
            "K2HPO4_g": 2.5,
        },
        "responsibility": {
            "faa_abundance": 1.0,
            "taa_abundance": 1.0,
            "mw_low": 1.0,
            "mw_medium": 1.0,
            "mw_high": 1.0,
            "aa_biosynthesis_gap": 1.0,
            "transporter_bonus": 1.0,
            "vitamin_b": 0.9,
            "nucleotides": 0.8,
            "sugar": 0.50,             # only 2.5g glucose in base
            "mineral": 0.30,
            "orgacid": 0.2,
            "nitrogen_quality": 0.8,
        },
    },
    "NB_2.0": {
        "display_name": "NB (2.0%, Peptone+Beef extract 대체)",
        "description": "Nutrient Broth의 질소원을 Sempio 펩톤으로 대체",
        "peptone_g_per_L": 20.0,
        "default_for": [],
        "base_provides": {
            "glucose_g": 0.0,
            "yeast_extract_g": 0.0,
            "NaCl_g": 5.0,
        },
        "responsibility": {
            "faa_abundance": 1.0,
            "taa_abundance": 1.0,
            "mw_low": 1.0,
            "mw_medium": 1.0,
            "mw_high": 1.0,
            "aa_biosynthesis_gap": 1.0,
            "transporter_bonus": 1.0,
            "vitamin_b": 1.0,
            "nucleotides": 0.8,
            "sugar": 0.8,             # no glucose in base
            "mineral": 0.5,
            "orgacid": 0.3,
            "nitrogen_quality": 0.8,
        },
    },
    "BHI_2.5": {
        "display_name": "BHI (2.5%, Brain+Heart infusion 대체)",
        "description": "Brain Heart Infusion의 질소원을 Sempio 펩톤으로 대체",
        "peptone_g_per_L": 25.0,
        "default_for": ["Bifidobacterium"],
        "base_provides": {
            "glucose_g": 2.0,
            "yeast_extract_g": 0.0,
            "NaCl_g": 5.0,
            "Na2HPO4_g": 2.5,
        },
        "responsibility": {
            "faa_abundance": 1.0,
            "taa_abundance": 1.0,
            "mw_low": 1.0,
            "mw_medium": 1.0,
            "mw_high": 1.0,
            "aa_biosynthesis_gap": 1.0,
            "transporter_bonus": 1.0,
            "vitamin_b": 1.0,
            "nucleotides": 0.8,
            "sugar": 0.55,
            "mineral": 0.35,
            "orgacid": 0.2,
            "nitrogen_quality": 0.8,
        },
    },
    "Custom": {
        "display_name": "사용자 정의 배지",
        "description": "직접 responsibility 값을 조정",
        "peptone_g_per_L": 25.0,
        "default_for": [],
        "base_provides": {},
        "responsibility": {
            "faa_abundance": 1.0,
            "taa_abundance": 1.0,
            "mw_low": 1.0,
            "mw_medium": 1.0,
            "mw_high": 1.0,
            "aa_biosynthesis_gap": 1.0,
            "transporter_bonus": 1.0,
            "vitamin_b": 0.5,
            "nucleotides": 0.5,
            "sugar": 0.5,
            "mineral": 0.5,
            "orgacid": 0.3,
            "nitrogen_quality": 0.7,
        },
    },
}

# Genus → strain type (for default media selection)
GENUS_TO_STRAIN_TYPE = {
    "Lactobacillus": "LAB", "Lactiplantibacillus": "LAB",
    "Lacticaseibacillus": "LAB", "Levilactobacillus": "LAB",
    "Lentilactobacillus": "LAB", "Latilactobacillus": "LAB",
    "Limosilactobacillus": "LAB", "Ligilactobacillus": "LAB",
    "Lactococcus": "LAB", "Streptococcus": "LAB",
    "Leuconostoc": "LAB", "Pediococcus": "LAB",
    "Enterococcus": "LAB", "Weissella": "LAB",
    "Escherichia": "Enterobacteriaceae", "Salmonella": "Enterobacteriaceae",
    "Klebsiella": "Enterobacteriaceae", "Enterobacter": "Enterobacteriaceae",
    "Serratia": "Enterobacteriaceae",
    "Bacillus": "Bacillus", "Priestia": "Bacillus", "Brevibacillus": "Bacillus",
    "Corynebacterium": "Corynebacterium", "Brevibacterium": "Corynebacterium",
    "Bifidobacterium": "Bifidobacterium",
}


def get_default_media(genus: str) -> str:
    """Get default media config key for a genus."""
    strain_type = GENUS_TO_STRAIN_TYPE.get(genus, "")
    for key, cfg in MEDIA_CONFIGS.items():
        if strain_type in cfg.get("default_for", []):
            return key
    return "MRS_2.5"  # safe default


def get_media_responsibility(media_key: str) -> dict[str, float]:
    """Get the peptone responsibility weights for a media type."""
    cfg = MEDIA_CONFIGS.get(media_key, MEDIA_CONFIGS["MRS_2.5"])
    return cfg["responsibility"]


def get_all_media_keys() -> list[str]:
    """Get all available media config keys."""
    return list(MEDIA_CONFIGS.keys())


def get_media_display_name(media_key: str) -> str:
    """Get display name for a media config."""
    return MEDIA_CONFIGS.get(media_key, {}).get("display_name", media_key)

"""ML feature extraction for PeptoMatch.

Builds a (sample × feature) matrix from a (strain, peptone, media) triple.
Used by both training (ml_train.py) and inference (ml_predict.py) so the
feature vocabulary stays identical between the two.

Feature blocks
--------------
1. composition_supply_*   — 45 KEGG-mapped composition values from the peptone
                             (mmol/L contribution at the given peptone conc).
2. genome_demand_*        — strain prior demands (AA biosynthesis gaps,
                             vitamin biosynthesis gaps) from genome_prior.py.
3. fba_predicted_mu       — single scalar from fba_simulator (if a GEM exists).
4. media_onehot_*         — one-hot for the basal medium (MRS_2.0, LB_2.0, …).
5. peptone_blend_*        — blending fractions if the sample is a 2-peptone mix.
6. peptone_pct            — peptone concentration (g/L scaled to 0~1).

All blocks are dense numpy/pandas. Missing genome priors → 0 (model handles).
The FBA block is optional — if FBA is unavailable we still produce the rest
so the ML pipeline degrades gracefully to "ML without FBA".

NOTE: This module is intentionally pure-functional. No model is loaded here.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd

from .fba_simulator import KEGG_COMPOUND_MAP, peptone_to_mmol_L

logger = logging.getLogger("peptomatch.ml_features")

# Stable column order (used by both train and predict to keep the matrix aligned)
COMPOSITION_COLUMNS: list[str] = [f"comp_supply_{cpd}" for cpd in KEGG_COMPOUND_MAP.values()]

# AA prior keys we read from genome_prior_features.json
AA_PRIOR_KEYS = [
    "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
    "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val",
]
VIT_PRIOR_KEYS = ["B1", "B2", "B3", "B5", "B6", "B7", "B9", "B12"]
GENOME_DEMAND_COLUMNS: list[str] = (
    [f"demand_aa_{k}" for k in AA_PRIOR_KEYS]
    + [f"demand_vit_{k}" for k in VIT_PRIOR_KEYS]
)

# Media one-hot vocabulary (extend when media_config gains more presets)
MEDIA_VOCAB: list[str] = [
    "MRS_2.0", "MRS_2.5", "MRS_3.0",
    "LB_1.0", "LB_2.0",
    "TSB_1.7", "TSB_2.0",
    "M9_basic",
]
MEDIA_COLUMNS: list[str] = [f"media_{k}" for k in MEDIA_VOCAB]

# Misc scalar columns
MISC_COLUMNS: list[str] = [
    "peptone_pct",      # 0~1 (g/L ÷ 100 normalised)
    "is_blend",         # 0/1
    "blend_ratio_1",    # 0~1
    "blend_ratio_2",    # 0~1
    "fba_predicted_mu",
    "fba_available",    # 0/1 — was FBA actually run?
    # ── Control baseline (same-experiment MRS positive control) ──
    # These strip plate/season/strain-condition variance so the model
    # learns the *relative* effect of peptone, not absolute OD drift.
    # All zeros + control_available=0 when no control was recorded.
    "control_max_od",
    "control_mu_max",
    "control_auc",
    "control_lag_time_h",
    "control_final_od",
    "control_available",    # 0/1 — was a control recorded on same plate?
]

ALL_FEATURE_COLUMNS: list[str] = (
    COMPOSITION_COLUMNS + GENOME_DEMAND_COLUMNS + MEDIA_COLUMNS + MISC_COLUMNS
)


# ── Lightweight inputs ─────────────────────────────────────────

@dataclass
class FeatureInputs:
    """Minimum information needed to build one feature row.

    Either `peptone_name` (single) or `peptone_1`/`peptone_2`+ratios (blend).

    `control_*` fields are optional baselines from the same-experiment MRS
    positive control. Fill them from :func:`GrowthDB.get_control_for` at
    dataset assembly time; leave as `None` when no control is available and
    the feature vector will fall back to zeros + `control_available=0`.
    """
    strain_id: int
    peptone_name: str
    peptone_pct: float                              # g/L
    media_key: str                                  # e.g. "MRS_2.0"
    peptone_1: Optional[str] = None
    peptone_2: Optional[str] = None
    ratio_1: float = 100.0
    ratio_2: float = 0.0
    fba_predicted_mu: Optional[float] = None        # None if FBA unavailable
    # ── Control baseline (from GrowthDB.get_control_for) ──
    control_max_od: Optional[float] = None
    control_mu_max: Optional[float] = None
    control_auc: Optional[float] = None
    control_lag_time_h: Optional[float] = None
    control_final_od: Optional[float] = None
    extra: dict[str, Any] = field(default_factory=dict)


# ── Per-block builders ─────────────────────────────────────────

def _composition_supply_row(
    composition_df: pd.DataFrame,
    inputs: FeatureInputs,
) -> dict[str, float]:
    """45 KEGG-mapped supply values (mmol/L), accounting for blending."""
    out = {col: 0.0 for col in COMPOSITION_COLUMNS}

    def _add(name: Optional[str], frac: float):
        if not name or frac <= 0:
            return
        sub = composition_df[composition_df["Sample_name"] == name]
        if sub.empty:
            return
        eff_conc = inputs.peptone_pct * frac
        contribution = peptone_to_mmol_L(sub.iloc[0], eff_conc)
        for cpd, mmol in contribution.items():
            col = f"comp_supply_{cpd}"
            if col in out:
                out[col] += float(mmol)

    if inputs.peptone_2 and inputs.ratio_2 > 0:
        _add(inputs.peptone_1, inputs.ratio_1 / 100.0)
        _add(inputs.peptone_2, inputs.ratio_2 / 100.0)
    else:
        _add(inputs.peptone_name or inputs.peptone_1, 1.0)
    return out


def _genome_demand_row(
    strain_id: int,
    genome_priors: dict[int, dict],
) -> dict[str, float]:
    """Pull AA/vitamin biosynthesis values; 'demand' = 1 - synthesis_score.

    A strain that *cannot* synthesise an AA must take it up from the medium →
    high demand. Score in [0, 1]. Missing priors → 0 (neutral).
    """
    out = {col: 0.0 for col in GENOME_DEMAND_COLUMNS}
    prior = genome_priors.get(int(strain_id)) or {}
    aa = prior.get("aa_biosynthesis") or {}
    vit = prior.get("vitamin_biosynthesis") or {}
    for k in AA_PRIOR_KEYS:
        v = aa.get(k)
        if v is not None:
            out[f"demand_aa_{k}"] = max(0.0, 1.0 - float(v))
    for k in VIT_PRIOR_KEYS:
        v = vit.get(k)
        if v is not None:
            out[f"demand_vit_{k}"] = max(0.0, 1.0 - float(v))
    return out


def _media_onehot_row(media_key: str) -> dict[str, float]:
    out = {col: 0.0 for col in MEDIA_COLUMNS}
    col = f"media_{media_key}"
    if col in out:
        out[col] = 1.0
    return out


def _misc_row(inputs: FeatureInputs) -> dict[str, float]:
    is_blend = bool(inputs.peptone_2 and inputs.ratio_2 > 0)
    # A control baseline is "available" only when at least one of the core
    # metrics was recorded. We key availability on max_od (always produced by
    # `_calc_metrics`) so a row where max_od is None means no control.
    control_available = 1.0 if inputs.control_max_od is not None else 0.0
    return {
        "peptone_pct": float(inputs.peptone_pct) / 100.0,
        "is_blend": 1.0 if is_blend else 0.0,
        "blend_ratio_1": (inputs.ratio_1 / 100.0) if is_blend else 1.0,
        "blend_ratio_2": (inputs.ratio_2 / 100.0) if is_blend else 0.0,
        "fba_predicted_mu": float(inputs.fba_predicted_mu or 0.0),
        "fba_available": 1.0 if inputs.fba_predicted_mu is not None else 0.0,
        "control_max_od": float(inputs.control_max_od or 0.0),
        "control_mu_max": float(inputs.control_mu_max or 0.0),
        "control_auc": float(inputs.control_auc or 0.0),
        "control_lag_time_h": float(inputs.control_lag_time_h or 0.0),
        "control_final_od": float(inputs.control_final_od or 0.0),
        "control_available": control_available,
    }


# ── Public API ─────────────────────────────────────────────────

def build_feature_row(
    inputs: FeatureInputs,
    composition_df: pd.DataFrame,
    genome_priors: Optional[dict[int, dict]] = None,
) -> pd.Series:
    """Return a single Series indexed by ALL_FEATURE_COLUMNS."""
    row: dict[str, float] = {}
    row.update(_composition_supply_row(composition_df, inputs))
    row.update(_genome_demand_row(inputs.strain_id, genome_priors or {}))
    row.update(_media_onehot_row(inputs.media_key))
    row.update(_misc_row(inputs))
    return pd.Series(row, index=ALL_FEATURE_COLUMNS).fillna(0.0)


def build_feature_matrix(
    inputs_list: list[FeatureInputs],
    composition_df: pd.DataFrame,
    genome_priors: Optional[dict[int, dict]] = None,
) -> pd.DataFrame:
    """Stack many rows into a (n × n_features) DataFrame."""
    rows = [
        build_feature_row(i, composition_df, genome_priors).to_dict()
        for i in inputs_list
    ]
    return pd.DataFrame(rows, columns=ALL_FEATURE_COLUMNS).fillna(0.0)


# ── Loaders ────────────────────────────────────────────────────

def load_genome_priors(prior_file: Path) -> dict[int, dict]:
    """Read genome_prior_features.json produced by GenomePriorBuilder."""
    import json
    if not prior_file.exists():
        logger.warning(f"Genome prior file not found: {prior_file}")
        return {}
    try:
        with open(prior_file, "r", encoding="utf-8") as f:
            raw = json.load(f)
        return {int(k): v for k, v in raw.items()}
    except Exception as e:
        logger.error(f"Failed to load genome priors: {e}")
        return {}


def feature_dim() -> int:
    """Return total feature dimension (handy for model construction)."""
    return len(ALL_FEATURE_COLUMNS)

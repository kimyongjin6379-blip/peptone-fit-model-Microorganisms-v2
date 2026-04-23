"""ML training pipeline (skeleton).

Pulls labelled samples out of GrowthDB, builds a feature matrix using
ml_features.py, trains a model, and persists it under models/.

Phase progression (per the plan)
--------------------------------
Phase 4b-1  Ridge / BayesianRidge       (≤ ~30 samples)
Phase 4b-2  Log-target + GP + EI        (~30~100 samples)
Phase 4b-3  XGBoost / LightGBM          (≥ 100 samples)
Phase 4b-4  AMN-style NN                (≥ 200 samples + GEM embedding)

This file ships Ridge + GP only — XGBoost and AMN are stubbed so the API
contract stays stable while we accumulate data this week.

Usage
-----
    python -m peptomatch.ml_train --model ridge --target max_od
    python -m peptomatch.ml_train --model gp    --target log_max_od

Saved artifacts
---------------
    models/peptomatch_{model}_{target}.joblib   — sklearn Pipeline
    models/peptomatch_{model}_{target}.meta.json — feature columns + metrics

Inference uses ml_predict.MLPredictor which loads these files.
"""

from __future__ import annotations

import argparse
import json
import logging
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd

from .ml_features import (
    ALL_FEATURE_COLUMNS,
    FeatureInputs,
    build_feature_matrix,
    load_genome_priors,
)

logger = logging.getLogger("peptomatch.ml_train")

DEFAULT_MODELS_DIR = Path("models")
DEFAULT_PRIOR_FILE = Path("outputs/genome_prior_features.json")

# Raw targets pulled straight from growth_metrics.
RAW_TARGETS = {"max_od", "mu_max", "auc", "final_od", "lag_time_h"}
# Fold targets = raw / control (ratio-scale, plate-normalized).
FOLD_TARGETS = {"max_od_fold", "mu_max_fold", "auc_fold", "final_od_fold"}
# Delta targets = raw - control (difference-scale, for lag-like metrics).
DELTA_TARGETS = {"lag_time_delta", "t_max_od_delta"}

VALID_TARGETS = RAW_TARGETS | FOLD_TARGETS | DELTA_TARGETS
VALID_MODELS = {"ridge", "bayesian_ridge", "gp", "xgboost"}


# Map each derived target back to the (raw_metric, control_key, op) triple.
# op is either "fold" (sample / control) or "delta" (sample - control).
# Both sides skipped when control_key is None in the baseline dict.
_DERIVED_TARGET_SPEC: dict[str, tuple[str, str, str]] = {
    "max_od_fold":     ("max_od",      "control_max_od",      "fold"),
    "mu_max_fold":     ("mu_max",      "control_mu_max",      "fold"),
    "auc_fold":        ("auc",         "control_auc",         "fold"),
    "final_od_fold":   ("final_od",    "control_final_od",    "fold"),
    "lag_time_delta":  ("lag_time_h",  "control_lag_time_h",  "delta"),
    "t_max_od_delta":  ("t_max_od_h",  "control_t_max_od_h",  "delta"),
}


# ── Dataset assembly ──────────────────────────────────────────

def _curve_to_inputs(
    curve: dict,
    strain_id_lookup: dict[str, int],
    control_baselines: Optional[dict[tuple[int, str], dict]] = None,
    media_key_default: str = "MRS_2.0",
) -> Optional[FeatureInputs]:
    """Map one growth_curve row → FeatureInputs.

    When `control_baselines` is provided, the matching (experiment_id,
    strain_name) baseline is injected so feature extraction can see the
    plate's own MRS reference — this is how seasonal/instrument variance
    gets normalized out of the learnt signal.
    """
    strain_name = (curve.get("strain_name") or "").strip()
    sid = strain_id_lookup.get(strain_name)
    if sid is None:
        return None
    pep_name = curve.get("peptone_name") or curve.get("peptone_1") or ""
    if not pep_name:
        return None

    # Resolve control baseline (exact strain match → any same-experiment
    # control → None).
    ctrl: dict = {}
    if control_baselines:
        exp_id = curve.get("experiment_id")
        if exp_id is not None:
            key = (int(exp_id), strain_name)
            if key in control_baselines:
                ctrl = control_baselines[key]
            else:
                for (eid, _s), info in control_baselines.items():
                    if eid == int(exp_id):
                        ctrl = info
                        break

    return FeatureInputs(
        strain_id=int(sid),
        peptone_name=pep_name,
        peptone_pct=float(curve.get("peptone_pct") or 0.0),
        media_key=curve.get("media_key") or media_key_default,
        peptone_1=curve.get("peptone_1") or None,
        peptone_2=curve.get("peptone_2") or None,
        ratio_1=float(curve.get("ratio_1") or 100.0),
        ratio_2=float(curve.get("ratio_2") or 0.0),
        fba_predicted_mu=curve.get("fba_predicted_mu"),  # optional column
        control_max_od=ctrl.get("control_max_od"),
        control_mu_max=ctrl.get("control_mu_max"),
        control_auc=ctrl.get("control_auc"),
        control_lag_time_h=ctrl.get("control_lag_time_h"),
        control_final_od=ctrl.get("control_final_od"),
    )


def _derive_target_value(row: dict, target: str, control: dict) -> Optional[float]:
    """Compute one target value from a curve row + its control baseline.

    Returns None when the target is a derived (fold/delta) metric but no
    control baseline is available for that row — the row is then dropped
    from training (you can't learn `max_od_fold` without a reference).
    """
    if target in RAW_TARGETS:
        v = row.get(target)
        try:
            return float(v) if v is not None else None
        except (TypeError, ValueError):
            return None

    spec = _DERIVED_TARGET_SPEC.get(target)
    if spec is None:
        return None
    raw_key, ctrl_key, op = spec
    raw_v = row.get(raw_key)
    ctrl_v = control.get(ctrl_key) if control else None
    if raw_v is None or ctrl_v is None:
        return None
    try:
        raw_f = float(raw_v)
        ctrl_f = float(ctrl_v)
    except (TypeError, ValueError):
        return None
    if op == "fold":
        if ctrl_f <= 1e-9:
            return None
        return raw_f / ctrl_f
    if op == "delta":
        return raw_f - ctrl_f
    return None


@dataclass
class TrainResult:
    n_samples: int
    n_features: int
    target: str
    model_name: str
    r2_train: float
    r2_cv_mean: float
    r2_cv_std: float
    spearman_cv_mean: float
    top_k_hit_rate: float        # fraction of top-3 predictions in top-3 truth
    feature_columns: list[str]
    target_transform: str        # "identity" or "log1p"
    saved_to: str

    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=2, ensure_ascii=False)


def assemble_dataset(
    growth_db,
    composition_df: pd.DataFrame,
    strain_id_lookup: dict[str, int],
    target: str = "max_od",
    genome_priors: Optional[dict[int, dict]] = None,
) -> tuple[pd.DataFrame, np.ndarray]:
    """Pull joined growth_curves+metrics, build (X, y).

    Parameters
    ----------
    growth_db : GrowthDB
    composition_df : peptone composition DataFrame
    strain_id_lookup : map "Lactobacillus rhamnosus" → strain_id (int)
    target : one of VALID_TARGETS (raw / fold / delta).

        Raw targets train directly on the measured metric.
        Fold/delta targets require a same-experiment control baseline; rows
        without one are dropped from the learning set for that target.

    Returns
    -------
    (X, y) where X has columns = ALL_FEATURE_COLUMNS, y is a numpy array.
    """
    if target not in VALID_TARGETS:
        raise ValueError(f"target must be one of {VALID_TARGETS}")

    rows = growth_db.get_ml_training_data() if hasattr(growth_db, "get_ml_training_data") else []
    if not rows:
        return pd.DataFrame(columns=ALL_FEATURE_COLUMNS), np.array([])

    # Pull control baselines once, reuse for every row.
    control_baselines: dict = {}
    if hasattr(growth_db, "get_control_baselines"):
        try:
            control_baselines = growth_db.get_control_baselines() or {}
        except Exception as e:
            logger.warning(f"Could not load control baselines: {e}")

    def _lookup_ctrl(exp_id, strain) -> dict:
        if not control_baselines or exp_id is None:
            return {}
        key = (int(exp_id), (strain or "").strip())
        if key in control_baselines:
            return control_baselines[key]
        for (eid, _s), info in control_baselines.items():
            if eid == int(exp_id):
                return info
        return {}

    inputs: list[FeatureInputs] = []
    targets: list[float] = []
    skipped_no_ctrl = 0
    for row in rows:
        ctrl = _lookup_ctrl(row.get("experiment_id"), row.get("strain_name"))
        y_f = _derive_target_value(row, target, ctrl)
        if y_f is None:
            if target not in RAW_TARGETS and not ctrl:
                skipped_no_ctrl += 1
            continue
        feat_in = _curve_to_inputs(row, strain_id_lookup, control_baselines)
        if feat_in is None:
            continue
        inputs.append(feat_in)
        targets.append(y_f)

    if skipped_no_ctrl:
        logger.info(
            f"Dropped {skipped_no_ctrl} rows from '{target}' training "
            f"(no same-experiment control to normalize against)"
        )

    if not inputs:
        return pd.DataFrame(columns=ALL_FEATURE_COLUMNS), np.array([])

    X = build_feature_matrix(inputs, composition_df, genome_priors)
    y = np.asarray(targets, dtype=float)
    return X, y


# ── Model factories ───────────────────────────────────────────

def _build_model(name: str):
    """Return an sklearn Pipeline (StandardScaler → estimator)."""
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    if name == "ridge":
        from sklearn.linear_model import Ridge
        est = Ridge(alpha=1.0, random_state=0)
    elif name == "bayesian_ridge":
        from sklearn.linear_model import BayesianRidge
        est = BayesianRidge()
    elif name == "gp":
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import Matern, WhiteKernel
        kernel = Matern(length_scale=1.0, nu=2.5) + WhiteKernel(noise_level=1e-2)
        est = GaussianProcessRegressor(kernel=kernel, normalize_y=True, random_state=0)
    elif name == "xgboost":
        try:
            from xgboost import XGBRegressor
        except ImportError as e:
            raise ImportError("xgboost not installed; pip install xgboost") from e
        est = XGBRegressor(
            n_estimators=300, max_depth=4, learning_rate=0.05,
            subsample=0.9, colsample_bytree=0.9, random_state=0,
        )
    else:
        raise ValueError(f"unknown model: {name}")
    return Pipeline([("scaler", StandardScaler(with_mean=True, with_std=True)), ("est", est)])


# ── Training ──────────────────────────────────────────────────

def train(
    X: pd.DataFrame,
    y: np.ndarray,
    model_name: str = "ridge",
    target: str = "max_od",
    log_target: bool = False,
    cv_splits: int = 5,
    models_dir: Path = DEFAULT_MODELS_DIR,
) -> TrainResult:
    """Train one model and persist it. Returns TrainResult."""
    from sklearn.model_selection import KFold, cross_val_predict
    from sklearn.metrics import r2_score
    from scipy.stats import spearmanr
    import joblib

    if len(y) < 5:
        raise ValueError(
            f"Need ≥5 labelled samples to train (have {len(y)}). "
            "Accumulate more experiments first."
        )

    y_train = np.log1p(y) if log_target else y.copy()
    pipe = _build_model(model_name)

    pipe.fit(X.values, y_train)
    yhat_train = pipe.predict(X.values)
    r2_train = r2_score(y_train, yhat_train)

    n_splits = min(cv_splits, len(y))
    if n_splits >= 3:
        kf = KFold(n_splits=n_splits, shuffle=True, random_state=0)
        yhat_cv = cross_val_predict(pipe, X.values, y_train, cv=kf)
        r2_cv = r2_score(y_train, yhat_cv)
        sp_cv, _ = spearmanr(y_train, yhat_cv)
        # Top-K hit rate: top-3 predicted intersect top-3 actual
        k = min(3, len(y))
        top_pred = set(np.argsort(yhat_cv)[-k:])
        top_true = set(np.argsort(y_train)[-k:])
        top_k_hit = len(top_pred & top_true) / k
        # CV stability: bootstrap a small std
        r2_std = float(np.std([
            r2_score(y_train[v], cross_val_predict(pipe, X.values, y_train, cv=kf)[v])
            for v in [list(range(len(y_train)))]
        ]))
    else:
        r2_cv = float("nan")
        sp_cv = float("nan")
        top_k_hit = float("nan")
        r2_std = float("nan")

    models_dir.mkdir(parents=True, exist_ok=True)
    suffix = f"{model_name}_{target}"
    model_path = models_dir / f"peptomatch_{suffix}.joblib"
    meta_path = models_dir / f"peptomatch_{suffix}.meta.json"

    joblib.dump(pipe, model_path)

    result = TrainResult(
        n_samples=int(len(y)),
        n_features=int(X.shape[1]),
        target=target,
        model_name=model_name,
        r2_train=float(r2_train),
        r2_cv_mean=float(r2_cv) if not np.isnan(r2_cv) else 0.0,
        r2_cv_std=float(r2_std) if not np.isnan(r2_std) else 0.0,
        spearman_cv_mean=float(sp_cv) if not np.isnan(sp_cv) else 0.0,
        top_k_hit_rate=float(top_k_hit) if not np.isnan(top_k_hit) else 0.0,
        feature_columns=list(X.columns),
        target_transform="log1p" if log_target else "identity",
        saved_to=str(model_path),
    )

    meta = asdict(result) | {"trained_at": datetime.now().isoformat()}
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)

    logger.info(f"Trained {model_name} on {len(y)} samples. CV R²={result.r2_cv_mean:.3f}")
    return result


# ── CLI ───────────────────────────────────────────────────────

def _cli():
    parser = argparse.ArgumentParser(description="Train PeptoMatch ML model")
    parser.add_argument("--model", choices=sorted(VALID_MODELS), default="ridge")
    parser.add_argument("--target", choices=sorted(VALID_TARGETS), default="max_od")
    parser.add_argument("--log-target", action="store_true",
                        help="Train on log1p(y) — recommended for max_od distribution")
    parser.add_argument("--models-dir", type=Path, default=DEFAULT_MODELS_DIR)
    parser.add_argument("--prior-file", type=Path, default=DEFAULT_PRIOR_FILE)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    # Lazy imports so this module is importable without full backend
    from .growth_db import GrowthDB
    from .strain_db import StrainDB
    from .io_loaders import load_composition_data
    from .utils import load_config

    cfg = load_config()
    gdb = GrowthDB()
    sdb = StrainDB(Path("data/strains.db"))
    sdf = sdb.get_strain_df()
    strain_lookup = {
        str(r.get("full_name", "")).strip(): int(r["strain_id"])
        for _, r in sdf.iterrows()
    }
    comp_df = load_composition_data(
        Path(cfg["data"]["composition_file"]),
        sheet_name=cfg["data"].get("composition_sheet", "data"),
    )
    priors = load_genome_priors(args.prior_file)

    X, y = assemble_dataset(gdb, comp_df, strain_lookup, target=args.target,
                            genome_priors=priors)
    if len(y) == 0:
        logger.error("No labelled samples found. Ingest growth data first.")
        return

    res = train(X, y, model_name=args.model, target=args.target,
                log_target=args.log_target, models_dir=args.models_dir)
    print(res.to_json())


if __name__ == "__main__":
    _cli()

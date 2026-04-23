"""ML inference for PeptoMatch (skeleton).

Loads a trained sklearn pipeline (saved by ml_train.py) and provides:

* `predict_one(inputs)`              — single point prediction (μ or hat-OD).
* `predict_with_uncertainty(inputs)` — (mean, std) for GP / BayesianRidge.
* `rank_candidates(candidates, ...)` — Bayesian-Optimisation-style ranking.

Ranking strategies
------------------
* "mean"          — pure greedy on predicted target (favours exploitation).
* "ei"            — Expected Improvement vs current best, recommended when a
                    GP model is loaded so σ is informative.
* "ucb"           — Upper Confidence Bound (μ + κ·σ).

Returns 503-style errors via `MLPredictor.is_ready` so the API layer can
report "no model trained yet" without crashing.
"""

from __future__ import annotations

import json
import logging
import math
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd

from .ml_features import (
    ALL_FEATURE_COLUMNS,
    FeatureInputs,
    build_feature_matrix,
    build_feature_row,
)

logger = logging.getLogger("peptomatch.ml_predict")

DEFAULT_MODELS_DIR = Path("models")


# ── Loader ────────────────────────────────────────────────────

class MLPredictor:
    """Loads one trained model + meta and serves predictions.

    Lazily loads on first use so the API can construct without a model file.
    """

    def __init__(
        self,
        model_name: str = "ridge",
        target: str = "max_od",
        models_dir: Path = DEFAULT_MODELS_DIR,
    ):
        self.model_name = model_name
        self.target = target
        self.models_dir = Path(models_dir)
        self._pipe = None
        self._meta: Optional[dict] = None

    # ── lifecycle ─────────────────────────────────────
    @property
    def model_path(self) -> Path:
        return self.models_dir / f"peptomatch_{self.model_name}_{self.target}.joblib"

    @property
    def meta_path(self) -> Path:
        return self.models_dir / f"peptomatch_{self.model_name}_{self.target}.meta.json"

    @property
    def is_ready(self) -> bool:
        return self.model_path.exists() and self.meta_path.exists()

    def _load(self):
        if self._pipe is not None:
            return
        if not self.is_ready:
            raise FileNotFoundError(
                f"Model not trained yet: {self.model_path}. "
                f"Run `python -m peptomatch.ml_train --model {self.model_name} "
                f"--target {self.target}` first."
            )
        import joblib
        self._pipe = joblib.load(self.model_path)
        with open(self.meta_path, "r", encoding="utf-8") as f:
            self._meta = json.load(f)
        logger.info(
            f"Loaded model {self.model_path.name} "
            f"(R²cv={self._meta.get('r2_cv_mean', 0):.3f}, "
            f"n={self._meta.get('n_samples', 0)})"
        )

    # ── feature alignment ─────────────────────────────
    def _align(self, X: pd.DataFrame) -> np.ndarray:
        cols = (self._meta or {}).get("feature_columns", ALL_FEATURE_COLUMNS)
        for c in cols:
            if c not in X.columns:
                X[c] = 0.0
        return X[cols].fillna(0.0).values

    def _untransform(self, y_pred: np.ndarray) -> np.ndarray:
        if self._meta and self._meta.get("target_transform") == "log1p":
            return np.expm1(y_pred)
        return y_pred

    # ── public predictions ────────────────────────────
    def predict_one(
        self,
        inputs: FeatureInputs,
        composition_df: pd.DataFrame,
        genome_priors: Optional[dict[int, dict]] = None,
    ) -> float:
        self._load()
        row = build_feature_row(inputs, composition_df, genome_priors).to_frame().T
        Xv = self._align(row)
        yhat = self._pipe.predict(Xv)
        return float(self._untransform(yhat)[0])

    def predict_with_uncertainty(
        self,
        inputs_list: list[FeatureInputs],
        composition_df: pd.DataFrame,
        genome_priors: Optional[dict[int, dict]] = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return (mean, std). std is 0 for models without native uncertainty."""
        self._load()
        X = build_feature_matrix(inputs_list, composition_df, genome_priors)
        Xv = self._align(X)

        # Try estimator native return_std (GP, BayesianRidge)
        est = self._pipe.named_steps["est"]
        scaler = self._pipe.named_steps["scaler"]
        Xs = scaler.transform(Xv)
        try:
            mean, std = est.predict(Xs, return_std=True)
        except (TypeError, AttributeError):
            mean = est.predict(Xs)
            std = np.zeros_like(mean)

        return self._untransform(np.asarray(mean)), np.asarray(std)

    # ── ranking strategies ────────────────────────────
    def rank_candidates(
        self,
        inputs_list: list[FeatureInputs],
        composition_df: pd.DataFrame,
        genome_priors: Optional[dict[int, dict]] = None,
        strategy: str = "ei",
        current_best: Optional[float] = None,
        ucb_kappa: float = 2.0,
        top_k: int = 10,
    ) -> list[dict]:
        """Return ranked candidates with score + uncertainty.

        strategy:
          - "mean": rank by predicted target (exploit only)
          - "ei":   Expected Improvement (balances exploit/explore via σ)
          - "ucb":  μ + κ·σ
        """
        self._load()
        mean, std = self.predict_with_uncertainty(inputs_list, composition_df, genome_priors)

        if current_best is None:
            current_best = float(np.max(mean))

        if strategy == "mean" or np.allclose(std, 0):
            score = mean.copy()
        elif strategy == "ucb":
            score = mean + ucb_kappa * std
        elif strategy == "ei":
            score = _expected_improvement(mean, std, current_best)
        else:
            raise ValueError(f"unknown strategy: {strategy}")

        order = np.argsort(score)[::-1][:top_k]
        out: list[dict] = []
        for rank, idx in enumerate(order, start=1):
            ip = inputs_list[idx]
            out.append({
                "rank": rank,
                "strain_id": ip.strain_id,
                "peptone_name": ip.peptone_name,
                "peptone_pct": ip.peptone_pct,
                "media_key": ip.media_key,
                "predicted_mean": float(mean[idx]),
                "predicted_std": float(std[idx]),
                "score": float(score[idx]),
                "strategy": strategy,
            })
        return out

    def info(self) -> dict:
        try:
            self._load()
        except FileNotFoundError:
            return {"ready": False, "model_name": self.model_name, "target": self.target}
        return {"ready": True, **(self._meta or {})}


# ── helpers ────────────────────────────────────────────────────

def _expected_improvement(mean: np.ndarray, std: np.ndarray, best: float, xi: float = 0.01) -> np.ndarray:
    """Closed-form EI for maximisation under Gaussian posterior."""
    std = np.maximum(std, 1e-9)
    z = (mean - best - xi) / std
    cdf = 0.5 * (1.0 + np.vectorize(math.erf)(z / math.sqrt(2.0)))
    pdf = np.exp(-0.5 * z * z) / math.sqrt(2.0 * math.pi)
    return (mean - best - xi) * cdf + std * pdf


# ── module-level convenience singleton ─────────────────────────

_DEFAULT: Optional[MLPredictor] = None


def get_default_predictor(
    model_name: str = "ridge",
    target: str = "max_od",
) -> MLPredictor:
    global _DEFAULT
    if _DEFAULT is None or _DEFAULT.model_name != model_name or _DEFAULT.target != target:
        _DEFAULT = MLPredictor(model_name=model_name, target=target)
    return _DEFAULT

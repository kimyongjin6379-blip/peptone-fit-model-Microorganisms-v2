"""FBA vs empirical calibration (skeleton).

Compares FBA-predicted growth rates (μ̂) against empirical measurements
collected in growth_db, then fits a simple correction so the FBA scalar
becomes a useful ML feature even when the GEM is biased.

Pipeline
--------
1. Pull (predicted_mu, empirical_mu_max, empirical_max_od) triples.
2. Compute Spearman ρ — is the rank order useful at all?
   - Target: ρ ≥ 0.7 → FBA alone is a usable signal.
   - ρ < 0.4 → GEM likely needs gap-filling on this strain.
3. Fit a per-strain linear correction:
       μ_corrected = a · μ_FBA + b
   so the scaled FBA can be added as a feature to ml_train without
   distorting absolute magnitudes.
4. Persist coefficients to outputs/calibration_{strain_id}.json.

This module deliberately stays small — proper calibration (isotonic,
GP residual model, hierarchical Bayesian per-genus) is Phase 4b-2+.
"""

from __future__ import annotations

import json
import logging
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Optional

import numpy as np

logger = logging.getLogger("peptomatch.calibration")

DEFAULT_OUT_DIR = Path("outputs")
ACCEPTABLE_RHO = 0.7   # plan target


@dataclass
class CalibrationResult:
    strain_id: int
    n_samples: int
    spearman_rho: float
    pearson_r: float
    slope: float
    intercept: float
    rmse: float
    target_metric: str          # "mu_max" | "max_od"
    notes: str = ""

    def predict(self, mu_fba: float) -> float:
        return self.slope * float(mu_fba) + self.intercept


# ── Core ───────────────────────────────────────────────────────

def calibrate_strain(
    samples: list[dict],
    strain_id: int,
    target_metric: str = "mu_max",
) -> Optional[CalibrationResult]:
    """Fit a strain-specific correction.

    `samples` items must contain:
        - "fba_predicted_mu"  : float
        - "<target_metric>"   : float (e.g. "mu_max" or "max_od")

    Returns None if not enough usable samples.
    """
    xs: list[float] = []
    ys: list[float] = []
    for s in samples:
        fb = s.get("fba_predicted_mu")
        emp = s.get(target_metric)
        try:
            xs.append(float(fb))
            ys.append(float(emp))
        except (TypeError, ValueError):
            continue

    if len(xs) < 4:
        logger.info(f"strain {strain_id}: only {len(xs)} samples — skipping")
        return None

    x = np.asarray(xs)
    y = np.asarray(ys)

    try:
        from scipy.stats import pearsonr, spearmanr
        rho, _ = spearmanr(x, y)
        r, _ = pearsonr(x, y)
    except ImportError:
        rho = float(np.corrcoef(np.argsort(x), np.argsort(y))[0, 1])
        r = float(np.corrcoef(x, y)[0, 1])

    # Linear fit: y = a·x + b
    slope, intercept = np.polyfit(x, y, deg=1)
    yhat = slope * x + intercept
    rmse = float(np.sqrt(np.mean((y - yhat) ** 2)))

    note = ""
    if rho < 0.4:
        note = "Low rank correlation — consider gap-filling the GEM."
    elif rho >= ACCEPTABLE_RHO:
        note = f"FBA usable as direct feature (ρ ≥ {ACCEPTABLE_RHO})."

    return CalibrationResult(
        strain_id=int(strain_id),
        n_samples=int(len(xs)),
        spearman_rho=float(rho) if rho is not None else 0.0,
        pearson_r=float(r) if r is not None else 0.0,
        slope=float(slope),
        intercept=float(intercept),
        rmse=rmse,
        target_metric=target_metric,
        notes=note,
    )


def calibrate_all(
    growth_db,
    fba_lookup: dict[tuple[int, str], float],
    target_metric: str = "mu_max",
    out_dir: Path = DEFAULT_OUT_DIR,
) -> list[CalibrationResult]:
    """Iterate every strain with empirical data and fit a calibration.

    `fba_lookup` maps (strain_id, peptone_name) → predicted μ. Build it from
    `peptomatch.fba_simulator.predict_growth_for_recommendation()` per pair.
    """
    rows = growth_db.get_fba_validation_data() if hasattr(growth_db, "get_fba_validation_data") else []
    if not rows:
        logger.info("No FBA validation data available yet.")
        return []

    by_strain: dict[int, list[dict]] = {}
    for r in rows:
        sid = r.get("strain_id")
        pn = r.get("peptone_name", "")
        key = (int(sid), pn) if sid is not None else None
        if key is None or key not in fba_lookup:
            continue
        r["fba_predicted_mu"] = fba_lookup[key]
        by_strain.setdefault(int(sid), []).append(r)

    out_dir.mkdir(parents=True, exist_ok=True)
    results: list[CalibrationResult] = []
    for sid, samples in by_strain.items():
        cal = calibrate_strain(samples, sid, target_metric=target_metric)
        if cal is None:
            continue
        results.append(cal)
        path = out_dir / f"calibration_strain_{sid}.json"
        meta = asdict(cal) | {"computed_at": datetime.now().isoformat()}
        with open(path, "w", encoding="utf-8") as f:
            json.dump(meta, f, indent=2, ensure_ascii=False)
        logger.info(
            f"calibration strain {sid}: n={cal.n_samples}, "
            f"ρ={cal.spearman_rho:.2f}, slope={cal.slope:.2f}"
        )

    return results


def load_calibration(strain_id: int, out_dir: Path = DEFAULT_OUT_DIR) -> Optional[CalibrationResult]:
    """Re-hydrate a previously saved CalibrationResult, or None if absent."""
    path = out_dir / f"calibration_strain_{strain_id}.json"
    if not path.exists():
        return None
    try:
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)
        # drop non-dataclass keys
        data.pop("computed_at", None)
        return CalibrationResult(**data)
    except Exception as e:
        logger.warning(f"failed to load calibration for {strain_id}: {e}")
        return None


def apply_calibration(
    mu_fba: float,
    strain_id: int,
    out_dir: Path = DEFAULT_OUT_DIR,
) -> float:
    """Apply per-strain correction; falls back to identity if no calibration exists."""
    cal = load_calibration(strain_id, out_dir)
    return cal.predict(mu_fba) if cal else float(mu_fba)

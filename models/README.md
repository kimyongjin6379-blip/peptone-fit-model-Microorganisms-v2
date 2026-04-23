# Trained ML models

This directory holds trained PeptoMatch ML pipelines saved by `ml_train.py`.

## Naming convention

```
peptomatch_{model_name}_{target}.joblib       # sklearn Pipeline
peptomatch_{model_name}_{target}.meta.json    # feature columns + CV metrics
```

## Train a model

```bash
# Phase 4b-1: Ridge baseline (works with ≥5 samples)
python -m peptomatch.ml_train --model ridge --target max_od

# Phase 4b-2: GP with log-transformed target + EI ranking
python -m peptomatch.ml_train --model gp --target max_od --log-target

# Phase 4b-3: XGBoost (≥100 samples recommended)
python -m peptomatch.ml_train --model xgboost --target mu_max
```

## Inference

`gateway.py` exposes:
- `GET  /api/ml/status?model=ridge&target=max_od` — check whether trained
- `POST /api/ml/predict` — single-point prediction
- `POST /api/ml/rank` — Bayesian-Optimisation-style ranking (EI/UCB/mean)

When no model file exists, the API returns HTTP 503 with the expected
training command — this is intentional (skeleton for upcoming experimental
data this week).

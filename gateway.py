"""
PeptoMatch FastAPI Application (Railway-ready)

Pure FastAPI + Jinja2 + Tailwind (CDN) stack. Streamlit has been dropped:
earlier Railway deployments struggled with the subprocess/WebSocket proxy,
so this module replaces the UI with server-rendered HTML pages that reuse
the existing peptomatch Python backend (scoring, strain DB, growth DB).

Endpoints
─────────
UI pages
    GET  /                    Dashboard
    GET  /recommend           Peptone recommendation form + results
    GET  /growth              Growth-curve experiment browser

Growth Data API (called by growth-curve-app)
    POST /api/ingest
    GET  /api/growth/summary
    GET  /api/growth/experiments
    GET  /api/growth/curves
    GET  /api/growth/ml-data
    GET  /api/growth/fba-data

Recommendation API
    POST /api/recommend       JSON in/out, used by /recommend page

Health
    GET  /healthz             Always 200 JSON
    GET  /api/health          Always 200 JSON
"""

from __future__ import annotations

import logging
import os
import sys
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Any, Optional

from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

logger = logging.getLogger("peptomatch.gateway")
logging.basicConfig(level=logging.INFO)

# ── Path setup ────────────────────────────────────────────────
_HERE = Path(__file__).parent.resolve()
_SRC = _HERE / "src"
if _SRC.is_dir() and str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

# Backend imports (wrapped in try/except so the app can still start
# and serve /healthz even if any one module fails)
try:
    from peptomatch.growth_db import GrowthDB
except Exception as e:
    logger.error(f"Failed to import GrowthDB: {e}")
    GrowthDB = None  # type: ignore

try:
    from peptomatch.utils import load_config
    from peptomatch.io_loaders import load_composition_data
    from peptomatch.scoring import PeptoneRecommender
    from peptomatch.explain import RecommendationExplainer
    from peptomatch.strain_db import StrainDB
    from peptomatch.media_config import (
        get_all_media_keys, get_media_display_name, get_default_media,
        MEDIA_CONFIGS,
    )
    _BACKEND_OK = True
except Exception as e:
    logger.error(f"Failed to import peptomatch backend: {e}")
    _BACKEND_OK = False

# ── Globals populated by lifespan ─────────────────────────────
growth_db: Optional["GrowthDB"] = None
strain_db: Optional[Any] = None
comp_df: Optional[Any] = None
app_config: Optional[dict] = None


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Startup: all initialization is best-effort so /healthz always returns 200."""
    global growth_db, strain_db, comp_df, app_config

    # Growth DB (required for API ingestion)
    if GrowthDB is not None:
        try:
            growth_db = GrowthDB()
            logger.info(f"GrowthDB ready: {growth_db.db_path}")
        except Exception as e:
            logger.error(f"GrowthDB init failed: {e}")

    # PeptoMatch backend (composition + strain DB + config)
    if _BACKEND_OK:
        try:
            cfg_path = _HERE / "config" / "config.yaml"
            app_config = load_config(cfg_path) if cfg_path.exists() else load_config()
            logger.info("config.yaml loaded")
        except Exception as e:
            logger.error(f"config load failed: {e}")

        try:
            sdb = StrainDB(_HERE / "data" / "strains.db")
            if sdb.count() == 0 and app_config:
                strain_path = Path(app_config["data"]["strain_file"])
                if strain_path.exists():
                    sdb.load_from_excel(strain_path)
            strain_db = sdb
            logger.info(f"StrainDB ready (count={sdb.count()})")
        except Exception as e:
            logger.error(f"StrainDB init failed: {e}")

        try:
            if app_config:
                comp_df = load_composition_data(
                    Path(app_config["data"]["composition_file"]),
                    sheet_name=app_config["data"].get("composition_sheet", "data"),
                )
                logger.info(f"composition loaded: {len(comp_df)} peptones")
        except Exception as e:
            logger.error(f"composition load failed: {e}")

    logger.info("PeptoMatch gateway startup complete")
    yield

    # Shutdown
    if growth_db is not None:
        try:
            growth_db.close()
        except Exception:
            pass
    if strain_db is not None:
        try:
            strain_db.close()
        except Exception:
            pass


app = FastAPI(title="PeptoMatch", lifespan=lifespan)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# ── Global exception handler (logs + visible error in response) ───
import traceback
from fastapi.exceptions import RequestValidationError
from starlette.exceptions import HTTPException as StarletteHTTPException


@app.exception_handler(Exception)
async def _unhandled_exc(request: Request, exc: Exception):
    tb = traceback.format_exc()
    logger.error(f"Unhandled error on {request.url.path}:\n{tb}")
    debug = os.getenv("PEPTOMATCH_DEBUG", "1") == "1"
    if debug:
        body = (
            f"<html><body style='font-family:monospace;background:#0b1020;color:#f87171;"
            f"padding:20px;'><h2>500 on {request.url.path}</h2><pre>{tb}</pre></body></html>"
        )
        return HTMLResponse(body, status_code=500)
    return JSONResponse(status_code=500, content={"error": str(exc)})

# ── Templates + static ────────────────────────────────────────
TEMPLATES_DIR = _HERE / "templates"
STATIC_DIR = _HERE / "static"
TEMPLATES_DIR.mkdir(exist_ok=True)
STATIC_DIR.mkdir(exist_ok=True)

templates = Jinja2Templates(directory=str(TEMPLATES_DIR))
app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")


# ── Health ────────────────────────────────────────────────────

@app.get("/healthz")
async def healthz():
    return JSONResponse({
        "status": "ok",
        "db_ready": growth_db is not None,
        "strain_db_ready": strain_db is not None,
        "comp_loaded": comp_df is not None,
    })


@app.get("/api/health")
async def api_health():
    return JSONResponse({"status": "ok"})


# ── UI pages ──────────────────────────────────────────────────

def _ctx(request: Request, **extra) -> dict:
    """Base template context."""
    return {
        "request": request,
        "backend_ok": _BACKEND_OK and strain_db is not None and comp_df is not None,
        "growth_ok": growth_db is not None,
        **extra,
    }


@app.get("/", response_class=HTMLResponse)
async def home(request: Request):
    summary = {}
    if growth_db is not None:
        try:
            summary = growth_db.get_summary()
        except Exception as e:
            logger.warning(f"summary failed: {e}")

    strain_count = 0
    peptone_count = 0
    if strain_db is not None:
        try:
            strain_count = strain_db.count()
        except Exception:
            pass
    if comp_df is not None:
        try:
            peptone_count = int(len(comp_df))
        except Exception:
            pass

    return templates.TemplateResponse(
        request=request,
        name="home.html",
        context=_ctx(
            request,
            summary=summary,
            strain_count=strain_count,
            peptone_count=peptone_count,
        ),
    )


@app.get("/recommend", response_class=HTMLResponse)
async def recommend_page(request: Request):
    strains = []
    media_options = []
    if strain_db is not None:
        try:
            df = strain_db.get_strain_df()
            strains = [
                {"id": int(r["strain_id"]), "name": r.get("full_name", "")}
                for _, r in df.iterrows()
            ]
        except Exception as e:
            logger.warning(f"strain list failed: {e}")
    if _BACKEND_OK:
        try:
            media_options = [
                {"key": k, "label": get_media_display_name(k)}
                for k in get_all_media_keys()
            ]
        except Exception:
            pass

    return templates.TemplateResponse(
        request=request,
        name="recommend.html",
        context=_ctx(request, strains=strains, media_options=media_options),
    )


@app.get("/growth", response_class=HTMLResponse)
async def growth_page(request: Request):
    experiments: list[dict] = []
    summary: dict = {}
    if growth_db is not None:
        try:
            experiments = growth_db.get_experiments()
            summary = growth_db.get_summary()
        except Exception as e:
            logger.warning(f"growth_db query failed: {e}")
    return templates.TemplateResponse(
        request=request,
        name="growth.html",
        context=_ctx(request, experiments=experiments, summary=summary),
    )


# ── Recommendation API (used by /recommend page JS) ───────────

@app.post("/api/recommend")
async def api_recommend(payload: dict):
    if not (_BACKEND_OK and strain_db is not None and comp_df is not None):
        return JSONResponse(
            status_code=503,
            content={"error": "backend not initialized"},
        )

    try:
        strain_id = int(payload.get("strain_id"))
        top_k = int(payload.get("top_k", 10))
        media_key = payload.get("media_key") or None
        sempio_only = bool(payload.get("sempio_only", True))
        language = payload.get("language", "ko")

        sdf = strain_db.get_strain_df()
        if sdf.empty:
            return JSONResponse(status_code=400, content={"error": "no strains in DB"})

        # Auto default media from genus if not provided
        if not media_key:
            row = sdf[sdf["strain_id"] == strain_id]
            genus = row.iloc[0].get("genus", "") if not row.empty else ""
            media_key = get_default_media(genus)

        recommender = PeptoneRecommender(comp_df, sdf, app_config)
        pf = app_config.get("peptone_filter") if sempio_only else None
        recs = recommender.recommend(
            strain_id, top_k=top_k, peptone_filter=pf, media_key=media_key,
        )

        explainer = RecommendationExplainer(comp_df, sdf, app_config, language=language)
        recs = explainer.explain_batch(strain_id, recs, top_n_reasons=3)
        summary = explainer.get_strain_summary(strain_id)

        # Convert DataFrame → list[dict] with only the columns the UI needs
        keep_cols = [c for c in ("rank", "peptone", "score", "explanation") if c in recs.columns]
        rows = recs[keep_cols].to_dict(orient="records")

        media_cfg = MEDIA_CONFIGS.get(media_key, {})
        return JSONResponse({
            "status": "ok",
            "strain_id": strain_id,
            "media_key": media_key,
            "media_display": media_cfg.get("display_name", media_key),
            "peptone_g_per_L": media_cfg.get("peptone_g_per_L"),
            "summary": summary,
            "recommendations": rows,
        })
    except Exception as e:
        logger.exception(f"recommend failed: {e}")
        return JSONResponse(status_code=500, content={"error": str(e)})


# ── Growth Data API ───────────────────────────────────────────

def _require_db():
    if growth_db is None:
        return JSONResponse(status_code=503, content={"error": "growth_db not initialized"})
    return None


@app.post("/api/ingest")
async def ingest_growth_data(payload: dict):
    err = _require_db()
    if err:
        return err
    try:
        result = growth_db.ingest(payload)
        logger.info(f"Ingest OK: {result}")
        return JSONResponse({"status": "ok", **result})
    except Exception as e:
        logger.error(f"Ingest error: {e}")
        return JSONResponse(status_code=500, content={"status": "error", "detail": str(e)})


@app.get("/api/growth/summary")
async def growth_summary():
    err = _require_db()
    if err:
        return err
    return JSONResponse(content=growth_db.get_summary())


@app.get("/api/growth/experiments")
async def list_experiments(media_type: Optional[str] = None):
    err = _require_db()
    if err:
        return err
    return JSONResponse(content=growth_db.get_experiments(media_type))


@app.get("/api/growth/curves")
async def list_curves(experiment_id: Optional[int] = None):
    err = _require_db()
    if err:
        return err
    return JSONResponse(content=growth_db.get_curves_with_metrics(experiment_id))


@app.get("/api/growth/ml-data")
async def ml_training_data():
    err = _require_db()
    if err:
        return err
    return JSONResponse(content=growth_db.get_ml_training_data())


@app.get("/api/growth/fba-data")
async def fba_validation_data():
    err = _require_db()
    if err:
        return err
    return JSONResponse(content=growth_db.get_fba_validation_data())

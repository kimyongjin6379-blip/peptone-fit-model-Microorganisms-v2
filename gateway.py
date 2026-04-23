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
    from peptomatch.kegg_viz import KEGGVisualizer
    _BACKEND_OK = True
except Exception as e:
    logger.error(f"Failed to import peptomatch backend: {e}")
    _BACKEND_OK = False

# GEM / FBA / ML modules — best-effort (cobra/sklearn/xgboost are optional)
try:
    from peptomatch.gem_manager import get_default_manager as get_gem_manager
    _GEM_OK = True
except Exception as e:
    logger.error(f"Failed to import GEMManager: {e}")
    _GEM_OK = False

try:
    from peptomatch.fba_simulator import (
        FBASimulator,
        load_basal_media,
        basal_medium_to_mmol_L,
        peptone_to_mmol_L,
        predict_growth_for_recommendation,
    )
    _FBA_OK = True
except Exception as e:
    logger.error(f"Failed to import fba_simulator: {e}")
    _FBA_OK = False

try:
    from peptomatch.ml_predict import MLPredictor
    from peptomatch.ml_features import FeatureInputs, load_genome_priors
    _ML_OK = True
except Exception as e:
    logger.error(f"Failed to import ml modules: {e}")
    _ML_OK = False

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


@app.get("/kegg", response_class=HTMLResponse)
async def kegg_page(request: Request):
    strains = []
    if strain_db is not None:
        try:
            df = strain_db.get_strain_df()
            strains = [
                {"id": int(r["strain_id"]), "name": r.get("full_name", "")}
                for _, r in df.iterrows()
            ]
        except Exception as e:
            logger.warning(f"strain list failed: {e}")
    return templates.TemplateResponse(
        request=request,
        name="kegg.html",
        context=_ctx(request, strains=strains),
    )


# ── KEGG API (used by /kegg page JS) ─────────────────────────

@app.post("/api/kegg/analyze")
async def api_kegg_analyze(payload: dict):
    if not (_BACKEND_OK and strain_db is not None):
        return JSONResponse(status_code=503, content={"error": "backend not initialized"})

    try:
        strain_id = int(payload.get("strain_id"))
        sdf = strain_db.get_strain_df()
        viz = KEGGVisualizer(sdf, app_config)

        # Prior info
        prior = viz.prior_builder.get_prior(strain_id)
        source = prior.get("source", "unknown")
        ko_count = prior.get("ko_count", 0)
        org_code = prior.get("kegg_org_code", "")

        # Charts → Plotly JSON
        import json as _json
        overview_json = _json.loads(viz.overview_chart(strain_id).to_json())
        aa_json = _json.loads(viz.aa_pathway_chart(strain_id).to_json())
        vit_json = _json.loads(viz.vitamin_chart(strain_id).to_json())

        # Deficient AAs/vitamins for flowcharts
        aa_synth = prior.get("aa_biosynthesis", {})
        vit_synth = prior.get("vitamin_biosynthesis", {})

        deficient_aa = [aa for aa, v in sorted(aa_synth.items(), key=lambda x: x[1]) if v < 0.8][:3]
        deficient_vit = [v for v, val in sorted(vit_synth.items(), key=lambda x: x[1]) if val < 0.8][:3]

        aa_flowcharts = {}
        for aa in deficient_aa:
            try:
                fig = viz.pathway_detail_chart(strain_id, aa)
                aa_flowcharts[aa] = _json.loads(fig.to_json())
            except Exception:
                pass

        vit_flowcharts = {}
        for vit in deficient_vit:
            try:
                fig = viz.pathway_detail_chart(strain_id, vit, pathway_source="vitamin")
                vit_flowcharts[vit] = _json.loads(fig.to_json())
            except Exception:
                pass

        # Explorable lists (all AA and vitamin names)
        all_aa = sorted(aa_synth.keys())
        all_vit = sorted(vit_synth.keys()) if vit_synth else []

        return JSONResponse({
            "status": "ok",
            "strain_id": strain_id,
            "source": source,
            "ko_count": ko_count,
            "org_code": org_code,
            "overview": overview_json,
            "aa_chart": aa_json,
            "vit_chart": vit_json,
            "deficient_aa_flowcharts": aa_flowcharts,
            "deficient_vit_flowcharts": vit_flowcharts,
            "all_aa": all_aa,
            "all_vit": all_vit,
        })
    except Exception as e:
        logger.exception(f"kegg analyze failed: {e}")
        return JSONResponse(status_code=500, content={"error": str(e)})


@app.post("/api/kegg/pathway-detail")
async def api_kegg_pathway_detail(payload: dict):
    """Fetch a single pathway flowchart for the explorer."""
    if not (_BACKEND_OK and strain_db is not None):
        return JSONResponse(status_code=503, content={"error": "backend not initialized"})

    try:
        strain_id = int(payload.get("strain_id"))
        pathway = payload.get("pathway", "")
        source = payload.get("pathway_source", "aa")  # "aa" or "vitamin"

        sdf = strain_db.get_strain_df()
        viz = KEGGVisualizer(sdf, app_config)
        fig = viz.pathway_detail_chart(
            strain_id, pathway,
            pathway_source="vitamin" if source == "vitamin" else None,
        )
        import json as _json
        return JSONResponse({
            "status": "ok",
            "chart": _json.loads(fig.to_json()),
        })
    except Exception as e:
        logger.exception(f"pathway detail failed: {e}")
        return JSONResponse(status_code=500, content={"error": str(e)})


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


@app.delete("/api/growth/experiments/{experiment_id}")
async def delete_experiment(experiment_id: int):
    """Delete a single experiment and all related curves/metrics."""
    err = _require_db()
    if err:
        return err
    try:
        result = growth_db.delete_experiment(experiment_id)
        if result["experiments"] == 0:
            return JSONResponse(
                status_code=404,
                content={"status": "not_found", "experiment_id": experiment_id},
            )
        return JSONResponse({"status": "ok", **result})
    except Exception as e:
        logger.error(f"delete_experiment failed: {e}")
        return JSONResponse(status_code=500, content={"status": "error", "detail": str(e)})


@app.delete("/api/growth/reset")
async def reset_all_growth(confirm: str = ""):
    """⚠️ Delete ALL experiments, curves, and metrics. Irreversible.

    Requires query param confirm=yes to actually execute.
    """
    err = _require_db()
    if err:
        return err
    if confirm != "yes":
        return JSONResponse(
            status_code=400,
            content={"status": "error", "detail": "confirm=yes required"},
        )
    return JSONResponse({"status": "ok", **growth_db.reset_all()})


# ── GEM management API ────────────────────────────────────────

def _basal_media_path() -> Path:
    """Where the basal_media_composition.xlsx lives. Returns Path even if absent."""
    return _HERE / "data" / "basal_media_composition.xlsx"


@app.get("/api/gem/list")
async def gem_list():
    """List cached GEM SBML files."""
    if not _GEM_OK:
        return JSONResponse(status_code=503, content={"error": "gem_manager unavailable"})
    try:
        mgr = get_gem_manager()
        return JSONResponse({
            "status": "ok",
            "cache_dir": str(mgr.cache_dir),
            "gems": [g.to_dict() for g in mgr.list_gems()],
        })
    except Exception as e:
        logger.exception("gem list failed")
        return JSONResponse(status_code=500, content={"error": str(e)})


@app.get("/api/gem/validate/{gcf}")
async def gem_validate(gcf: str):
    """Run cobra-based light validation on a cached GEM."""
    if not _GEM_OK:
        return JSONResponse(status_code=503, content={"error": "gem_manager unavailable"})
    try:
        mgr = get_gem_manager()
        return JSONResponse({"status": "ok", "result": mgr.validate_gem(gcf)})
    except Exception as e:
        return JSONResponse(status_code=500, content={"error": str(e)})


@app.post("/api/gem/generate")
async def gem_generate(payload: dict):
    """Trigger gapseq find→draft→fill via WSL/Linux. Long-running.

    Body: {"gcf": "GCF_...", "genome_fasta": "/path/to.fna",
           "bitscore_cutoff": 200, "timeout_seconds": 7200}
    Returns immediately if a cached GEM exists; otherwise runs synchronously
    (use a background task in future when we want non-blocking behaviour).
    """
    if not _GEM_OK:
        return JSONResponse(status_code=503, content={"error": "gem_manager unavailable"})
    try:
        gcf = payload.get("gcf", "").strip()
        fasta = payload.get("genome_fasta", "").strip()
        if not gcf or not fasta:
            return JSONResponse(status_code=400, content={"error": "gcf and genome_fasta required"})
        mgr = get_gem_manager()
        if mgr.has_gem(gcf):
            return JSONResponse({"status": "ok", "cached": True, "path": str(mgr.gem_path_for(gcf))})
        path = mgr.generate_gem(
            gcf=gcf,
            genome_fasta=Path(fasta),
            bitscore_cutoff=int(payload.get("bitscore_cutoff", 200)),
            timeout_seconds=int(payload.get("timeout_seconds", 7200)),
        )
        return JSONResponse({"status": "ok", "cached": False, "path": str(path)})
    except Exception as e:
        logger.exception("gem generate failed")
        return JSONResponse(status_code=500, content={"error": str(e)})


# ── FBA API ───────────────────────────────────────────────────

@app.post("/api/fba/predict")
async def fba_predict(payload: dict):
    """Predict growth rate (μ) for a (gcf, media, peptone) triple.

    Body: {"gcf": "...", "media_id": "MRS", "peptone_name": "SOY-N+",
           "peptone_conc_g_per_L": 25.0}
    """
    if not (_FBA_OK and _GEM_OK):
        return JSONResponse(status_code=503, content={"error": "FBA stack unavailable"})
    if comp_df is None:
        return JSONResponse(status_code=503, content={"error": "composition not loaded"})

    try:
        gcf = payload.get("gcf", "").strip()
        media_id = payload.get("media_id", "MRS")
        peptone_name = payload.get("peptone_name", "")
        conc = float(payload.get("peptone_conc_g_per_L", 25.0))

        mgr = get_gem_manager()
        gem_path = mgr.get_gem_path(gcf)
        if gem_path is None:
            return JSONResponse(status_code=404, content={
                "error": f"No GEM cached for {gcf}. POST /api/gem/generate first."
            })

        basal_path = _basal_media_path()
        if not basal_path.exists():
            return JSONResponse(status_code=503, content={
                "error": f"basal_media_composition.xlsx not found at {basal_path}"
            })
        basal_df = load_basal_media(basal_path)

        result = predict_growth_for_recommendation(
            gem_path=gem_path,
            basal_df=basal_df,
            media_id=media_id,
            composition_df=comp_df,
            peptone_name=peptone_name,
            peptone_conc_g_per_L=conc,
        )
        return JSONResponse({"status": "ok", **result})
    except Exception as e:
        logger.exception("fba predict failed")
        return JSONResponse(status_code=500, content={"error": str(e)})


@app.post("/api/fba/shadow-prices")
async def fba_shadow_prices(payload: dict):
    """Return top-N bottleneck nutrients for a (gcf, media, peptone) triple."""
    if not (_FBA_OK and _GEM_OK):
        return JSONResponse(status_code=503, content={"error": "FBA stack unavailable"})
    if comp_df is None:
        return JSONResponse(status_code=503, content={"error": "composition not loaded"})

    try:
        gcf = payload.get("gcf", "").strip()
        media_id = payload.get("media_id", "MRS")
        peptone_name = payload.get("peptone_name", "")
        conc = float(payload.get("peptone_conc_g_per_L", 25.0))
        top_n = int(payload.get("top_n", 10))

        mgr = get_gem_manager()
        gem_path = mgr.get_gem_path(gcf)
        if gem_path is None:
            return JSONResponse(status_code=404, content={"error": f"No GEM cached for {gcf}"})

        basal_path = _basal_media_path()
        if not basal_path.exists():
            return JSONResponse(status_code=503, content={
                "error": f"basal_media_composition.xlsx not found at {basal_path}"
            })
        basal_df = load_basal_media(basal_path)

        sim = FBASimulator(gem_path)
        basal_mmol = basal_medium_to_mmol_L(basal_df, media_id)
        pep_row = comp_df[comp_df["Sample_name"] == peptone_name]
        peptone_mmol = peptone_to_mmol_L(pep_row.iloc[0], conc) if not pep_row.empty else {}
        sim.set_medium(basal_mmol, peptone_mmol)
        return JSONResponse({
            "status": "ok",
            "shadow_prices": sim.get_shadow_prices(top_n=top_n),
        })
    except Exception as e:
        logger.exception("fba shadow prices failed")
        return JSONResponse(status_code=500, content={"error": str(e)})


@app.post("/api/fba/optimize-medium")
async def fba_optimize_medium(payload: dict):
    """Try each candidate supplement at its max bound and report Δμ."""
    if not (_FBA_OK and _GEM_OK):
        return JSONResponse(status_code=503, content={"error": "FBA stack unavailable"})
    try:
        gcf = payload.get("gcf", "").strip()
        media_id = payload.get("media_id", "MRS")
        peptone_name = payload.get("peptone_name", "")
        conc = float(payload.get("peptone_conc_g_per_L", 25.0))
        # candidates: {"cpd00027": [0, 30], ...}
        cand_raw = payload.get("candidate_supplements", {}) or {}
        candidates = {k: (float(v[0]), float(v[1])) for k, v in cand_raw.items()}
        target_growth = payload.get("target_growth")

        mgr = get_gem_manager()
        gem_path = mgr.get_gem_path(gcf)
        if gem_path is None:
            return JSONResponse(status_code=404, content={"error": f"No GEM cached for {gcf}"})

        basal_path = _basal_media_path()
        basal_df = load_basal_media(basal_path) if basal_path.exists() else None
        sim = FBASimulator(gem_path)
        if basal_df is not None and comp_df is not None:
            basal_mmol = basal_medium_to_mmol_L(basal_df, media_id)
            pep_row = comp_df[comp_df["Sample_name"] == peptone_name]
            peptone_mmol = peptone_to_mmol_L(pep_row.iloc[0], conc) if not pep_row.empty else {}
            sim.set_medium(basal_mmol, peptone_mmol)

        result = sim.optimize_medium(
            candidate_supplements=candidates,
            target_growth=float(target_growth) if target_growth is not None else None,
        )
        return JSONResponse({"status": "ok", **result})
    except Exception as e:
        logger.exception("fba optimize medium failed")
        return JSONResponse(status_code=500, content={"error": str(e)})


# ── ML API (returns 503 until a model is trained) ─────────────

def _resolve_ml_predictor(model_name: str, target: str):
    if not _ML_OK:
        return None, JSONResponse(status_code=503, content={"error": "ml stack unavailable"})
    pred = MLPredictor(model_name=model_name, target=target)
    if not pred.is_ready:
        return None, JSONResponse(
            status_code=503,
            content={
                "error": "model not trained yet",
                "expected_path": str(pred.model_path),
                "hint": f"python -m peptomatch.ml_train --model {model_name} --target {target}",
            },
        )
    return pred, None


@app.get("/api/ml/status")
async def ml_status(model: str = "ridge", target: str = "max_od"):
    """Report whether the requested model is trained + meta."""
    if not _ML_OK:
        return JSONResponse(status_code=503, content={"error": "ml stack unavailable"})
    pred = MLPredictor(model_name=model, target=target)
    return JSONResponse({"status": "ok", **pred.info()})


@app.post("/api/ml/predict")
async def ml_predict(payload: dict):
    """Single-point ML prediction.

    Body: {"strain_id": int, "peptone_name": str, "peptone_pct": float,
           "media_key": "MRS_2.0",
           "model": "ridge", "target": "max_od",
           "fba_predicted_mu": float|None}
    """
    if comp_df is None:
        return JSONResponse(status_code=503, content={"error": "composition not loaded"})
    pred, err = _resolve_ml_predictor(payload.get("model", "ridge"), payload.get("target", "max_od"))
    if err:
        return err
    try:
        priors_path = _HERE / "outputs" / "genome_prior_features.json"
        priors = load_genome_priors(priors_path) if _ML_OK else {}
        inputs = FeatureInputs(
            strain_id=int(payload["strain_id"]),
            peptone_name=payload.get("peptone_name", ""),
            peptone_pct=float(payload.get("peptone_pct", 25.0)),
            media_key=payload.get("media_key", "MRS_2.0"),
            peptone_1=payload.get("peptone_1"),
            peptone_2=payload.get("peptone_2"),
            ratio_1=float(payload.get("ratio_1", 100.0)),
            ratio_2=float(payload.get("ratio_2", 0.0)),
            fba_predicted_mu=payload.get("fba_predicted_mu"),
        )
        yhat = pred.predict_one(inputs, comp_df, priors)
        return JSONResponse({"status": "ok", "prediction": yhat, "target": pred.target})
    except Exception as e:
        logger.exception("ml predict failed")
        return JSONResponse(status_code=500, content={"error": str(e)})


@app.post("/api/ml/rank")
async def ml_rank(payload: dict):
    """Rank candidate (strain, peptone, pct, media) tuples using EI/UCB/mean.

    Body: {"candidates": [{"strain_id":..., "peptone_name":..., ...}, ...],
           "model": "ridge", "target": "max_od",
           "strategy": "ei"|"ucb"|"mean", "top_k": 10}
    """
    if comp_df is None:
        return JSONResponse(status_code=503, content={"error": "composition not loaded"})
    pred, err = _resolve_ml_predictor(payload.get("model", "ridge"), payload.get("target", "max_od"))
    if err:
        return err
    try:
        priors_path = _HERE / "outputs" / "genome_prior_features.json"
        priors = load_genome_priors(priors_path) if _ML_OK else {}
        candidates_raw = payload.get("candidates", []) or []
        if not candidates_raw:
            return JSONResponse(status_code=400, content={"error": "candidates required"})
        inputs_list = [
            FeatureInputs(
                strain_id=int(c["strain_id"]),
                peptone_name=c.get("peptone_name", ""),
                peptone_pct=float(c.get("peptone_pct", 25.0)),
                media_key=c.get("media_key", "MRS_2.0"),
                peptone_1=c.get("peptone_1"),
                peptone_2=c.get("peptone_2"),
                ratio_1=float(c.get("ratio_1", 100.0)),
                ratio_2=float(c.get("ratio_2", 0.0)),
                fba_predicted_mu=c.get("fba_predicted_mu"),
            )
            for c in candidates_raw
        ]
        ranked = pred.rank_candidates(
            inputs_list, comp_df, priors,
            strategy=payload.get("strategy", "ei"),
            top_k=int(payload.get("top_k", 10)),
        )
        return JSONResponse({"status": "ok", "ranked": ranked})
    except Exception as e:
        logger.exception("ml rank failed")
        return JSONResponse(status_code=500, content={"error": str(e)})

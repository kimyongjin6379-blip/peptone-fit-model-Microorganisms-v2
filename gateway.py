"""
PeptoMatch API Gateway (Railway-compatible)

FastAPI gateway that:
1. Handles /api/* endpoints for growth data ingestion
2. Reverse-proxies all other requests to Streamlit (internal subprocess)
3. Proxies WebSocket for Streamlit's real-time updates

Key: Returns 200 loading page (not 503) while Streamlit boots,
so Railway health checks pass immediately.
"""

import asyncio
import logging
import os
import subprocess
import sys
from contextlib import asynccontextmanager

import httpx
from fastapi import FastAPI, Request, Response
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse, JSONResponse

logger = logging.getLogger("peptomatch.gateway")
logging.basicConfig(level=logging.INFO)

# Import GrowthDB (handles both installed package and source paths)
# Add src/ to sys.path as a fallback in case `pip install -e .` didn't run
_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = os.path.join(_HERE, "src")
if os.path.isdir(_SRC) and _SRC not in sys.path:
    sys.path.insert(0, _SRC)

try:
    from peptomatch.growth_db import GrowthDB
except ImportError as e:
    logging.error(f"Failed to import GrowthDB: {e}")
    raise

STREAMLIT_PORT = int(os.getenv("STREAMLIT_INTERNAL_PORT", "8501"))
STREAMLIT_URL = f"http://127.0.0.1:{STREAMLIT_PORT}"

growth_db: GrowthDB = None
streamlit_proc: subprocess.Popen = None
streamlit_ready = False

LOADING_HTML = """<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>PeptoMatch</title>
<meta http-equiv="refresh" content="3">
<style>body{display:flex;justify-content:center;align-items:center;height:100vh;
font-family:sans-serif;background:#f0f2f6;color:#333;}
.box{text-align:center;}.spinner{width:40px;height:40px;border:4px solid #ddd;
border-top:4px solid #D32F2F;border-radius:50%;animation:spin 1s linear infinite;
margin:0 auto 16px;}@keyframes spin{to{transform:rotate(360deg);}}</style>
</head><body><div class="box"><div class="spinner"></div>
<h2>PeptoMatch Loading...</h2><p>Streamlit is starting up. This page will auto-refresh.</p>
</div></body></html>"""


async def _wait_for_streamlit():
    """Background task: poll Streamlit health until ready."""
    global streamlit_ready
    for _ in range(120):
        try:
            async with httpx.AsyncClient() as client:
                r = await client.get(f"{STREAMLIT_URL}/_stcore/health", timeout=3)
                if r.status_code == 200:
                    streamlit_ready = True
                    logger.info("Streamlit is ready!")
                    return
        except Exception:
            pass
        await asyncio.sleep(2)
    logger.error("Streamlit did not become ready within timeout")


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Startup: everything is wrapped in try/except so gateway NEVER fails to start.
    Railway health check on /healthz must succeed immediately.
    """
    global growth_db, streamlit_proc

    # 1. Init DB (fast, catch errors but don't block)
    try:
        growth_db = GrowthDB()
        logger.info(f"Growth DB initialized: {growth_db.db_path}")
    except Exception as e:
        logger.error(f"GrowthDB init failed: {e}")
        growth_db = None

    # 2. Start Streamlit subprocess (non-critical for health check)
    try:
        streamlit_cmd = [
            sys.executable, "-m", "streamlit", "run",
            "app/streamlit_app.py",
            "--server.port", str(STREAMLIT_PORT),
            "--server.address", "127.0.0.1",
            "--server.headless", "true",
            "--browser.gatherUsageStats", "false",
        ]
        logger.info(f"Starting Streamlit: {' '.join(streamlit_cmd)}")
        streamlit_proc = subprocess.Popen(
            streamlit_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        logger.info(f"Streamlit subprocess started (PID {streamlit_proc.pid})")
    except Exception as e:
        logger.error(f"Streamlit subprocess failed to start: {e}")
        streamlit_proc = None

    # 3. Background health-check (non-blocking)
    try:
        asyncio.create_task(_wait_for_streamlit())
    except Exception as e:
        logger.error(f"Failed to start health check task: {e}")

    # Gateway is ready NOW for Railway health check
    logger.info("Gateway lifespan startup complete")
    yield

    # Cleanup
    if streamlit_proc:
        try:
            streamlit_proc.terminate()
            streamlit_proc.wait(timeout=5)
        except Exception:
            try:
                streamlit_proc.kill()
            except Exception:
                pass
    if growth_db:
        try:
            growth_db.close()
        except Exception:
            pass


app = FastAPI(title="PeptoMatch Gateway", lifespan=lifespan)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)


# ── Health check (always 200) ──────────────────────────────────

@app.get("/healthz")
async def healthz():
    """Always return 200 for Railway health check."""
    return JSONResponse({
        "status": "ok",
        "streamlit_ready": streamlit_ready,
        "db_ready": growth_db is not None,
    })


@app.get("/api/health")
async def api_health():
    return JSONResponse({"status": "ok"})


# ── Growth Data API ────────────────────────────────────────────

def _require_db():
    if growth_db is None:
        return JSONResponse(status_code=503, content={"error": "database not initialized"})
    return None


@app.post("/api/ingest")
async def ingest_growth_data(payload: dict):
    err = _require_db()
    if err: return err
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
    if err: return err
    return JSONResponse(content=growth_db.get_summary())


@app.get("/api/growth/experiments")
async def list_experiments(media_type: str = None):
    err = _require_db()
    if err: return err
    return JSONResponse(content=growth_db.get_experiments(media_type))


@app.get("/api/growth/curves")
async def list_curves(experiment_id: int = None):
    err = _require_db()
    if err: return err
    return JSONResponse(content=growth_db.get_curves_with_metrics(experiment_id))


@app.get("/api/growth/ml-data")
async def ml_training_data():
    err = _require_db()
    if err: return err
    return JSONResponse(content=growth_db.get_ml_training_data())


@app.get("/api/growth/fba-data")
async def fba_validation_data():
    err = _require_db()
    if err: return err
    return JSONResponse(content=growth_db.get_fba_validation_data())


# ── Streamlit Reverse Proxy ────────────────────────────────────

@app.get("/")
async def root(request: Request):
    return await _proxy_http(request, "")


@app.api_route(
    "/{path:path}",
    methods=["GET", "POST", "PUT", "DELETE", "OPTIONS", "HEAD", "PATCH"],
)
async def proxy_catchall(request: Request, path: str = ""):
    return await _proxy_http(request, path)


async def _proxy_http(request: Request, path: str) -> Response:
    """Forward HTTP request to Streamlit, or show loading page."""
    if not streamlit_ready:
        return HTMLResponse(LOADING_HTML, status_code=200)

    url = f"{STREAMLIT_URL}/{path}"
    if request.url.query:
        url += f"?{request.url.query}"

    try:
        body = await request.body()
        headers = {k: v for k, v in request.headers.items() if k.lower() != "host"}

        async with httpx.AsyncClient(timeout=30) as client:
            resp = await client.request(
                method=request.method, url=url,
                headers=headers, content=body,
            )

        skip = {"content-encoding", "transfer-encoding", "content-length"}
        resp_headers = {k: v for k, v in resp.headers.items() if k.lower() not in skip}
        return Response(content=resp.content, status_code=resp.status_code, headers=resp_headers)

    except httpx.ConnectError:
        return HTMLResponse(LOADING_HTML, status_code=200)
    except Exception as e:
        logger.error(f"Proxy error [{path}]: {e}")
        return HTMLResponse(f"<h3>Proxy Error</h3><p>{e}</p>", status_code=502)


# ── WebSocket Proxy ────────────────────────────────────────────

from starlette.websockets import WebSocket, WebSocketDisconnect

@app.websocket("/_stcore/stream")
async def ws_proxy(client_ws: WebSocket):
    """Proxy Streamlit WebSocket connection."""
    await client_ws.accept()

    import websockets
    target = f"ws://127.0.0.1:{STREAMLIT_PORT}/_stcore/stream"

    try:
        async with websockets.connect(target) as server_ws:
            async def forward_client():
                try:
                    while True:
                        msg = await client_ws.receive_text()
                        await server_ws.send(msg)
                except (WebSocketDisconnect, Exception):
                    pass

            async def forward_server():
                try:
                    async for msg in server_ws:
                        if isinstance(msg, str):
                            await client_ws.send_text(msg)
                        else:
                            await client_ws.send_bytes(msg)
                except Exception:
                    pass

            done, pending = await asyncio.wait(
                [asyncio.create_task(forward_client()),
                 asyncio.create_task(forward_server())],
                return_when=asyncio.FIRST_COMPLETED,
            )
            for t in pending:
                t.cancel()

    except Exception as e:
        logger.warning(f"WS proxy error: {e}")
    finally:
        try:
            await client_ws.close()
        except Exception:
            pass

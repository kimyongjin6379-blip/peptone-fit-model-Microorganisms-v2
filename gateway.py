"""
PeptoMatch API Gateway

FastAPI gateway that:
1. Handles /api/ingest for growth data ingestion from growth-curve-app
2. Reverse-proxies all other requests to Streamlit (running internally)

Railway exposes one port -> this gateway serves on $PORT,
Streamlit runs on an internal port (8501).
"""

import asyncio
import logging
import os
import subprocess
import sys
import time
from contextlib import asynccontextmanager

import httpx
from fastapi import FastAPI, Request, Response
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, StreamingResponse

from src.peptomatch.growth_db import GrowthDB

logger = logging.getLogger("peptomatch.gateway")
logging.basicConfig(level=logging.INFO)

STREAMLIT_INTERNAL_PORT = 8501
STREAMLIT_BASE = f"http://127.0.0.1:{STREAMLIT_INTERNAL_PORT}"

# Global DB instance
growth_db: GrowthDB = None
streamlit_proc: subprocess.Popen = None


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Start Streamlit as subprocess on app startup, stop on shutdown."""
    global growth_db, streamlit_proc

    # Initialize Growth DB
    growth_db = GrowthDB()
    logger.info(f"Growth DB initialized: {growth_db.db_path}")

    # Start Streamlit as subprocess
    streamlit_cmd = [
        sys.executable, "-m", "streamlit", "run",
        "app/streamlit_app.py",
        "--server.port", str(STREAMLIT_INTERNAL_PORT),
        "--server.address", "127.0.0.1",
        "--server.headless", "true",
        "--browser.gatherUsageStats", "false",
    ]
    streamlit_proc = subprocess.Popen(streamlit_cmd)
    logger.info(f"Streamlit started (PID: {streamlit_proc.pid}) on port {STREAMLIT_INTERNAL_PORT}")

    # Wait for Streamlit to be ready
    for _ in range(30):
        try:
            async with httpx.AsyncClient() as client:
                resp = await client.get(f"{STREAMLIT_BASE}/_stcore/health")
                if resp.status_code == 200:
                    logger.info("Streamlit is ready")
                    break
        except Exception:
            pass
        await asyncio.sleep(1)

    yield

    # Cleanup
    if streamlit_proc:
        streamlit_proc.terminate()
        streamlit_proc.wait(timeout=5)
    if growth_db:
        growth_db.close()


app = FastAPI(title="PeptoMatch Gateway", lifespan=lifespan)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)


# ── API Endpoints ───────────────────────────────────────────────

@app.post("/api/ingest")
async def ingest_growth_data(payload: dict):
    """Receive growth data from growth-curve-app and store in SQLite."""
    try:
        result = growth_db.ingest(payload)
        logger.info(f"Ingest success: {result}")
        return JSONResponse(content={"status": "ok", **result})
    except Exception as e:
        logger.error(f"Ingest failed: {e}")
        return JSONResponse(
            status_code=500,
            content={"status": "error", "detail": str(e)},
        )


@app.get("/api/growth/summary")
async def growth_summary():
    """Return DB summary stats."""
    return JSONResponse(content=growth_db.get_summary())


@app.get("/api/growth/experiments")
async def list_experiments(media_type: str = None):
    """List experiments."""
    return JSONResponse(content=growth_db.get_experiments(media_type))


@app.get("/api/growth/curves")
async def list_curves(experiment_id: int = None):
    """List growth curves with metrics."""
    return JSONResponse(content=growth_db.get_curves_with_metrics(experiment_id))


@app.get("/api/growth/ml-data")
async def ml_training_data():
    """Get ML training data (peptone screening experiments)."""
    return JSONResponse(content=growth_db.get_ml_training_data())


@app.get("/api/growth/fba-data")
async def fba_validation_data():
    """Get FBA validation data (media optimization experiments)."""
    return JSONResponse(content=growth_db.get_fba_validation_data())


# ── Reverse Proxy to Streamlit ──────────────────────────────────

@app.api_route(
    "/{path:path}",
    methods=["GET", "POST", "PUT", "DELETE", "OPTIONS", "HEAD", "PATCH"],
)
async def proxy_to_streamlit(request: Request, path: str = ""):
    """Proxy all non-API requests to Streamlit."""
    target_url = f"{STREAMLIT_BASE}/{path}"

    # Forward query params
    if request.url.query:
        target_url += f"?{request.url.query}"

    try:
        body = await request.body()
        headers = dict(request.headers)
        # Remove host header to avoid confusion
        headers.pop("host", None)

        async with httpx.AsyncClient(timeout=30.0) as client:
            resp = await client.request(
                method=request.method,
                url=target_url,
                headers=headers,
                content=body,
            )

        # Forward response
        excluded_headers = {"content-encoding", "transfer-encoding", "content-length"}
        response_headers = {
            k: v for k, v in resp.headers.items()
            if k.lower() not in excluded_headers
        }

        return Response(
            content=resp.content,
            status_code=resp.status_code,
            headers=response_headers,
        )
    except httpx.ConnectError:
        return JSONResponse(
            status_code=503,
            content={"detail": "Streamlit is starting up, please wait..."},
        )
    except Exception as e:
        logger.error(f"Proxy error: {e}")
        return JSONResponse(
            status_code=502,
            content={"detail": f"Proxy error: {str(e)}"},
        )


# Root path (Streamlit index)
@app.get("/")
async def proxy_root(request: Request):
    """Proxy root to Streamlit."""
    return await proxy_to_streamlit(request, path="")


# ── WebSocket Proxy for Streamlit ───────────────────────────────

from starlette.websockets import WebSocket, WebSocketDisconnect
import websockets


@app.websocket("/_stcore/stream")
async def websocket_proxy(ws: WebSocket):
    """Proxy WebSocket connections to Streamlit's _stcore/stream endpoint."""
    await ws.accept()
    target_ws_url = f"ws://127.0.0.1:{STREAMLIT_INTERNAL_PORT}/_stcore/stream"

    try:
        async with websockets.connect(target_ws_url) as target_ws:

            async def client_to_server():
                try:
                    while True:
                        data = await ws.receive_text()
                        await target_ws.send(data)
                except WebSocketDisconnect:
                    pass

            async def server_to_client():
                try:
                    async for message in target_ws:
                        if isinstance(message, str):
                            await ws.send_text(message)
                        else:
                            await ws.send_bytes(message)
                except Exception:
                    pass

            await asyncio.gather(client_to_server(), server_to_client())

    except Exception as e:
        logger.error(f"WebSocket proxy error: {e}")
    finally:
        try:
            await ws.close()
        except Exception:
            pass

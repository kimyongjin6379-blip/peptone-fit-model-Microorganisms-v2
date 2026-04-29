#!/bin/bash
# ============================================================
# PeptoMatch 시작 스크립트 (Railway)
# Pure FastAPI app — Streamlit 제거됨
# ============================================================

echo "=== PeptoMatch Starting ==="
echo "PORT=${PORT:-8000}"
echo "PWD=$(pwd)"

# venv PATH 주입 (Nixpacks가 /opt/venv에 설치)
if [ -d "/opt/venv" ]; then
    export PATH="/opt/venv/bin:$PATH"
    echo "venv active: $(which python)"
fi

# ── Volume sync: build 시점의 xlsx를 /app/data Volume에 동기화 ──
# Volume mount(/app/data)가 build의 새 xlsx를 가리는 문제 해결.
# DB 파일들(*.db)은 건드리지 않음 — 영속화 유지.
BUILD_DATA="/app/_build_data/data"
if [ -d "$BUILD_DATA" ]; then
    echo "[volume-sync] copying build xlsx files into /app/data ..."
    mkdir -p /app/data
    shopt -s nullglob
    for f in "$BUILD_DATA"/*.xlsx; do
        cp -v "$f" /app/data/
    done
    shopt -u nullglob
    echo "[volume-sync] done"
else
    echo "[volume-sync] $BUILD_DATA missing — skipping"
fi

# KofamScan 설치 (선택 기능, 실패해도 무시)
export KOFAM_DIR="/app/kofam"
if [ -f "/app/scripts/setup_kofamscan.sh" ]; then
    timeout 60 bash /app/scripts/setup_kofamscan.sh || echo "KofamScan setup skipped/timed out (non-critical)"
fi
export KOFAMSCAN_PATH="${KOFAM_DIR}/kofam_scan/exec_annotation"
export KOFAMSCAN_PROFILES="${KOFAM_DIR}/profiles"

echo "Starting FastAPI on 0.0.0.0:${PORT:-8000}..."
exec uvicorn gateway:app \
    --host 0.0.0.0 \
    --port "${PORT:-8000}" \
    --log-level info \
    --access-log

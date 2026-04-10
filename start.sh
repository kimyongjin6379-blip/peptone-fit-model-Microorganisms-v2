#!/bin/bash
# ============================================================
# PeptoMatch 시작 스크립트 (Railway용)
# Gateway (FastAPI) → Streamlit (내부 subprocess)
# ============================================================

echo "=== PeptoMatch Starting ==="
echo "PORT=${PORT:-8000}"
echo "PWD=$(pwd)"
echo "PATH=$PATH"

# venv 확인
if [ -d "/opt/venv" ]; then
    echo "venv found at /opt/venv"
    export PATH="/opt/venv/bin:$PATH"
    which python
    which uvicorn || echo "WARN: uvicorn not found in PATH"
fi

# KofamScan 설치 (첫 실행 시에만, 실패해도 무시)
export KOFAM_DIR="/app/kofam"
if [ -f "/app/scripts/setup_kofamscan.sh" ]; then
    timeout 60 bash /app/scripts/setup_kofamscan.sh || echo "KofamScan setup skipped/timed out (non-critical)"
else
    echo "KofamScan setup script not found, skipping"
fi

export KOFAMSCAN_PATH="${KOFAM_DIR}/kofam_scan/exec_annotation"
export KOFAMSCAN_PROFILES="${KOFAM_DIR}/profiles"

echo "Starting gateway on 0.0.0.0:${PORT:-8000}..."
exec uvicorn gateway:app \
    --host 0.0.0.0 \
    --port "${PORT:-8000}" \
    --log-level info \
    --access-log

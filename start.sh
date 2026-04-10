#!/bin/bash
# ============================================================
# PeptoMatch 시작 스크립트 (Railway용)
# Gateway (FastAPI) → Streamlit (내부 subprocess)
# ============================================================

set -e

echo "=== PeptoMatch Starting ==="

# KofamScan 설치 (첫 실행 시에만 다운로드)
export KOFAM_DIR="/app/kofam"
bash /app/scripts/setup_kofamscan.sh || echo "KofamScan setup skipped (non-critical)"

export KOFAMSCAN_PATH="${KOFAM_DIR}/kofam_scan/exec_annotation"
export KOFAMSCAN_PROFILES="${KOFAM_DIR}/profiles"

echo "Starting gateway on port ${PORT:-8000}..."
exec uvicorn gateway:app \
    --host 0.0.0.0 \
    --port "${PORT:-8000}" \
    --log-level info

#!/bin/bash
# ============================================================
# PeptoMatch 시작 스크립트 (Railway용)
# 1. KofamScan 설치 확인 (없으면 자동 다운로드)
# 2. FastAPI 게이트웨이 실행 (Streamlit을 내부 프록시)
# ============================================================

set -e

echo "=== PeptoMatch Starting ==="

# KofamScan 설치 (첫 실행 시에만 다운로드)
export KOFAM_DIR="/app/kofam"
bash /app/scripts/setup_kofamscan.sh || echo "KofamScan setup skipped (non-critical)"

# 환경변수로 KofamScan 경로 전달
export KOFAMSCAN_PATH="${KOFAM_DIR}/kofam_scan/exec_annotation"
export KOFAMSCAN_PROFILES="${KOFAM_DIR}/profiles"

# FastAPI 게이트웨이 실행 (Streamlit은 게이트웨이 내부에서 subprocess로 시작)
exec uvicorn gateway:app \
    --host 0.0.0.0 \
    --port "${PORT:-8000}" \
    --log-level info

#!/bin/bash
# ============================================================
# PeptoMatch 시작 스크립트 (Railway용)
# 1. KofamScan 설치 확인 (없으면 자동 다운로드)
# 2. Streamlit 앱 실행
# ============================================================

set -e

echo "=== PeptoMatch Starting ==="

# KofamScan 설치 (첫 실행 시에만 다운로드)
export KOFAM_DIR="/app/kofam"
bash /app/scripts/setup_kofamscan.sh || echo "KofamScan setup skipped (non-critical)"

# 환경변수로 KofamScan 경로 전달
export KOFAMSCAN_PATH="${KOFAM_DIR}/kofam_scan/exec_annotation"
export KOFAMSCAN_PROFILES="${KOFAM_DIR}/profiles"

# Streamlit 실행
exec streamlit run app/streamlit_app.py \
    --server.port="${PORT:-8501}" \
    --server.address=0.0.0.0 \
    --server.headless=true

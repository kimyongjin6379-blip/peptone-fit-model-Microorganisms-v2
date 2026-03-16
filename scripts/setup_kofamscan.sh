#!/bin/bash
# ============================================================
# KofamScan 자동 설치 스크립트
# Railway 첫 실행 시 KOfam HMM 프로파일 + KofamScan 다운로드
# ============================================================

set -e

KOFAM_DIR="${KOFAM_DIR:-/app/kofam}"
PROFILE_DIR="${KOFAM_DIR}/profiles"
KO_LIST="${KOFAM_DIR}/ko_list"
KOFAMSCAN_BIN="${KOFAM_DIR}/kofam_scan/exec_annotation"

echo "=== KofamScan Setup ==="
echo "KOFAM_DIR: ${KOFAM_DIR}"

# 이미 설치되어 있으면 스킵
if [ -f "${KOFAMSCAN_BIN}" ] && [ -d "${PROFILE_DIR}" ] && [ -f "${KO_LIST}" ]; then
    echo "KofamScan already installed. Skipping setup."
    exit 0
fi

mkdir -p "${KOFAM_DIR}"
cd "${KOFAM_DIR}"

# 1. KofamScan 다운로드 (Ruby 스크립트, ~100KB)
if [ ! -f "${KOFAMSCAN_BIN}" ]; then
    echo "Downloading KofamScan..."
    wget -q "https://www.genome.jp/ftp/tools/kofam_scan/kofam_scan-1.3.0.tar.gz" -O kofam_scan.tar.gz
    tar xzf kofam_scan.tar.gz
    rm kofam_scan.tar.gz
    # Rename to consistent directory
    mv kofam_scan-* kofam_scan 2>/dev/null || true
    chmod +x "${KOFAMSCAN_BIN}"
    echo "KofamScan installed."
fi

# 2. ko_list 다운로드 (~3MB)
if [ ! -f "${KO_LIST}" ]; then
    echo "Downloading ko_list..."
    wget -q "https://www.genome.jp/ftp/db/kofam/ko_list.gz" -O ko_list.gz
    gunzip ko_list.gz
    echo "ko_list downloaded."
fi

# 3. HMM 프로파일 다운로드 (~4GB, 시간 소요)
if [ ! -d "${PROFILE_DIR}" ]; then
    echo "Downloading KOfam HMM profiles (~4GB)... This will take a while."
    wget -q --show-progress "https://www.genome.jp/ftp/db/kofam/profiles.tar.gz" -O profiles.tar.gz
    echo "Extracting profiles..."
    tar xzf profiles.tar.gz
    rm profiles.tar.gz
    echo "Profiles installed."
fi

# 4. KofamScan config 파일 생성
cat > "${KOFAM_DIR}/kofam_scan/config.yml" << EOF
profile: ${PROFILE_DIR}
ko_list: ${KO_LIST}
cpu: 2
hmmer: $(which hmmsearch)
EOF

echo "=== KofamScan setup complete ==="
echo "Binary: ${KOFAMSCAN_BIN}"
echo "Profiles: ${PROFILE_DIR}"
echo "ko_list: ${KO_LIST}"

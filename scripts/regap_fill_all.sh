#!/bin/bash
# ────────────────────────────────────────────────────────────────
# 10균주 gapseq gap-fill 재실행 (균주별 적정 배지)
#
# find/draft 단계는 건너뛰고 fill만 재실행합니다.
# (각 균주 폴더의 {strain}-draft.RDS 를 재사용)
#
# 전제:
#   - 기존 작업 디렉토리: /mnt/d/folder1/gapseq_linux/gapseq_work/{STRAIN}/
#     안에 {strain}-draft.RDS 있음
#   - 우리 custom media CSV:
#     /mnt/d/folder1/peptomatch/data/gapseq_media/MRS_rich.csv, BHI_rich.csv
#   - gapseq 내장: ~/gapseq/dat/media/TSBmed.csv
#
# 실행:
#   bash regap_fill_all.sh               # rich 버전 (권장)
#   bash regap_fill_all.sh strict        # strict 버전 (gap-fill 실패 디버깅용)
#
# 소요 시간:
#   균주당 10~20분, 순차 10개 → 2~3시간
#   (병렬 2로 돌리려면 `sudo apt install parallel` 후 parallel 버전 따로)
# ────────────────────────────────────────────────────────────────

set -u

MODE="${1:-rich}"   # rich 또는 strict

# ── 경로 설정 ─────────────────────────────────────────────
# WSL에서 본 Windows D드라이브 경로
WORK=/mnt/d/folder1/gapseq_linux/gapseq_work
CUSTOM_MEDIA=/mnt/d/folder1/peptomatch/data/gapseq_media
GAPSEQ_MEDIA=$HOME/gapseq/dat/media
PEPTOMATCH_GEM=/mnt/d/folder1/peptomatch/outputs/gem_cache

LOG_DIR="$WORK/logs_v2"
mkdir -p "$LOG_DIR"

# ── 균주별 배지 매핑 ─────────────────────────────────────
# key  : 균주 약어
# value: 배지 CSV 전체 경로
declare -A STRAIN_MEDIA=(
  [LP]="$CUSTOM_MEDIA/MRS_${MODE}.csv"
  [LR]="$CUSTOM_MEDIA/MRS_${MODE}.csv"
  [LA]="$CUSTOM_MEDIA/MRS_${MODE}.csv"
  [LC]="$CUSTOM_MEDIA/MRS_${MODE}.csv"
  [LPC]="$CUSTOM_MEDIA/MRS_${MODE}.csv"
  [LS]="$CUSTOM_MEDIA/MRS_${MODE}.csv"
  [STT]="$CUSTOM_MEDIA/MRS_${MODE}.csv"
  [EF]="$CUSTOM_MEDIA/BHI_${MODE}.csv"
  [BS]="$GAPSEQ_MEDIA/TSBmed.csv"
  [BC]="$GAPSEQ_MEDIA/TSBmed.csv"
)

# ── 사전 점검 ─────────────────────────────────────────────
echo "======================================================================"
echo "  gapseq re-fill (fill stage only)   mode=$MODE"
echo "======================================================================"
echo "  WORK          : $WORK"
echo "  custom media  : $CUSTOM_MEDIA"
echo "  gapseq media  : $GAPSEQ_MEDIA"
echo "  logs          : $LOG_DIR"
echo

missing=0
for strain in LP LR LA LC LPC LS STT EF BS BC; do
  draft="$WORK/$strain/${strain}-draft.RDS"
  media="${STRAIN_MEDIA[$strain]}"
  printf "  %-4s  draft=%-6s  media=%-50s  " "$strain" \
         "$([ -f "$draft" ] && echo OK || echo MISSING)" \
         "$(basename $media) $([ -f "$media" ] && echo OK || echo MISSING)"
  if [ ! -f "$draft" ] || [ ! -f "$media" ]; then
    missing=$((missing+1))
    echo "[SKIP]"
  else
    echo "[READY]"
  fi
done

if [ $missing -gt 0 ]; then
  echo
  echo "WARNING: $missing strain(s) will be skipped. Continue? (Ctrl+C to abort, Enter to go)"
  read -r _
fi

# ── 실행 ──────────────────────────────────────────────────
echo
echo "======================================================================"
echo "  Starting re-fill..."
echo "======================================================================"

START_ALL=$(date +%s)

for strain in LP LR LA LC LPC LS STT EF BS BC; do
  draft="$WORK/$strain/${strain}-draft.RDS"
  media="${STRAIN_MEDIA[$strain]}"
  log="$LOG_DIR/${strain}.log"

  if [ ! -f "$draft" ] || [ ! -f "$media" ]; then
    echo "[$strain] SKIP (missing draft or media)"
    continue
  fi

  # 이미 이번 세션에서 완료된 균주 자동 스킵
  # (새 xml이 존재하고, allmed 백업보다 최신이면 건너뜀)
  new_xml="$WORK/$strain/${strain}.xml"
  bak_xml="$WORK/$strain/${strain}.allmed.xml"
  if [ -s "$new_xml" ] && [ -f "$bak_xml" ] && [ "$new_xml" -nt "$bak_xml" ]; then
    mu=$(grep -oP 'Final growth rate:\s*\K[0-9.]+' "$LOG_DIR/${strain}.log" 2>/dev/null | tail -1)
    echo "[$strain] SKIP (already done, μ=${mu:-?})"
    continue
  fi

  echo
  echo "[$(date +%H:%M:%S)] ▶ $strain  media=$(basename $media)"

  cd "$WORK/$strain"

  # 기존 .xml / .RDS 백업 (ALLmed 버전 보존)
  for ext in xml RDS; do
    if [ -f "${strain}.${ext}" ] && [ ! -f "${strain}.allmed.${ext}" ]; then
      cp "${strain}.${ext}" "${strain}.allmed.${ext}"
      echo "    backup: ${strain}.allmed.${ext}"
    fi
  done

  # fill only
  start_ts=$(date +%s)
  gapseq fill -m "${strain}-draft.RDS" -n "$media" > "$log" 2>&1
  ec=$?
  dur=$(( $(date +%s) - start_ts ))

  if [ $ec -eq 0 ] && [ -s "${strain}.xml" ]; then
    # 생장률 추출
    mu=$(grep -oP 'Final growth rate:\s*\K[0-9.]+' "$log" | tail -1)
    echo "    [OK]  ${strain}.xml  μ=${mu:-?}  (${dur}s)"
  else
    echo "    [FAIL]  (exit=$ec, ${dur}s)  see: $log"
    tail -10 "$log" | sed 's/^/        /'
  fi
done

DUR=$(( $(date +%s) - START_ALL ))
echo
echo "======================================================================"
echo "  TOTAL: ${DUR}s ($((DUR/60))m $((DUR%60))s)"
echo "======================================================================"
echo
echo "Summary:"
for strain in LP LR LA LC LPC LS STT EF BS BC; do
  xml="$WORK/$strain/${strain}.xml"
  if [ -s "$xml" ]; then
    mu=$(grep -oP 'Final growth rate:\s*\K[0-9.]+' "$LOG_DIR/${strain}.log" 2>/dev/null | tail -1)
    echo "  ✓ $strain  μ=${mu:-?}  $(du -h "$xml" | cut -f1)"
  else
    echo "  ✗ $strain  (missing output)"
  fi
done

echo
echo "======================================================================"
echo "  Next: copy new GEMs to peptomatch"
echo "======================================================================"
echo "  # (backup current ALLmed versions first)"
echo "  mkdir -p $PEPTOMATCH_GEM/archive_allmed"
echo "  mv $PEPTOMATCH_GEM/*.xml $PEPTOMATCH_GEM/archive_allmed/"
echo "  cp $WORK/{LP,LR,LA,LC,LPC,LS,STT,EF,BS,BC}/*.xml $PEPTOMATCH_GEM/"
echo "  (LR.xml, LP.xml 등이 복사됨 — .allmed.xml은 제외)"

#!/bin/bash
# ────────────────────────────────────────────────────────────────
# Tier 1 균주 6종 전체 gapseq 파이프라인 실행
#   (LRE, LB, LF, LH, LG, EC)
#
# 전제:
#   - FASTA: /mnt/d/folder1/gapseq_linux/gapseq_work/{STRAIN}/{STRAIN}.fna
#     (download_tier1_genomes.sh 로 다운받음)
#   - Media CSV: /mnt/d/folder1/peptomatch/data/gapseq_media/{MRS,LB}_rich*.csv
#   - gapseq 설치되어 있고 PATH 잡혀 있음
#
# 소요 시간:
#   균주당 find ~60~120분 + find-transport ~15분 + draft ~10분 + fill ~15분
#   ≈ 2시간/균주 × 6 = 약 12시간 (야간 실행 권장)
#
# 실행:
#   bash run_tier1_pipeline.sh
#
# 재개:
#   이미 완료된 단계는 자동 스킵 (output 파일 존재 여부로 판단)
# ────────────────────────────────────────────────────────────────

set -u

WORK=/mnt/d/folder1/gapseq_linux/gapseq_work
CUSTOM_MEDIA=/mnt/d/folder1/peptomatch/data/gapseq_media
LOG_DIR="$WORK/logs_tier1"
mkdir -p "$LOG_DIR"

# strain → media CSV 매핑
declare -A STRAIN_MEDIA=(
  [LRE]="$CUSTOM_MEDIA/MRS_rich.csv"
  [LB]="$CUSTOM_MEDIA/MRS_rich.csv"
  [LF]="$CUSTOM_MEDIA/MRS_rich.csv"
  [LH]="$CUSTOM_MEDIA/MRS_rich.csv"
  [LG]="$CUSTOM_MEDIA/MRS_rich_cys.csv"     # cysteine 보충
  [EC]="$CUSTOM_MEDIA/LB_rich.csv"
)

echo "======================================================================"
echo "  Tier 1 full gapseq pipeline"
echo "======================================================================"
echo "  Work dir : $WORK"
echo "  Logs     : $LOG_DIR"
echo "  Start    : $(date)"
echo

# ── 사전 점검 ─────────────────────────────────────────────
for strain in LRE LB LF LH LG EC; do
  fna="$WORK/$strain/${strain}.fna"
  media="${STRAIN_MEDIA[$strain]}"
  printf "  %-4s  fna=%-8s  media=%s  %s\n" "$strain" \
    "$([ -s "$fna" ] && echo OK || echo MISSING)" \
    "$(basename $media)" \
    "$([ -s "$media" ] && echo OK || echo MEDIA-MISSING)"
done
echo

START_ALL=$(date +%s)

# ── 파이프라인 실행 ───────────────────────────────────────
for strain in LRE LB LF LH LG EC; do
  fna="$WORK/$strain/${strain}.fna"
  media="${STRAIN_MEDIA[$strain]}"

  if [ ! -s "$fna" ] || [ ! -s "$media" ]; then
    echo "[$strain] SKIP (missing fna or media)"
    continue
  fi

  echo
  echo "======================================================================"
  echo "  [$(date +%H:%M:%S)] ▶ $strain   (genome: $(basename $fna))"
  echo "======================================================================"

  cd "$WORK/$strain"

  # Step 1: find (pathway analysis)
  if [ ! -s "${strain}-all-Pathways.tbl" ]; then
    echo "[$strain] Step 1/4  find (pathway analysis) ..."
    t0=$(date +%s)
    gapseq find -p all -b 200 "$fna" \
      > "$LOG_DIR/${strain}_find.log" 2>&1
    ec=$?
    dur=$(( $(date +%s) - t0 ))
    if [ $ec -eq 0 ]; then
      echo "    [OK]  find  (${dur}s)"
    else
      echo "    [FAIL]  find  (exit=$ec, ${dur}s)  → see $LOG_DIR/${strain}_find.log"
      continue
    fi
  else
    echo "[$strain] Step 1/4  find — SKIP (already done)"
  fi

  # Step 2: find-transport
  if [ ! -s "${strain}-Transporter.tbl" ]; then
    echo "[$strain] Step 2/4  find-transport ..."
    t0=$(date +%s)
    gapseq find-transport -b 200 "$fna" \
      > "$LOG_DIR/${strain}_transport.log" 2>&1
    ec=$?
    dur=$(( $(date +%s) - t0 ))
    if [ $ec -eq 0 ]; then
      echo "    [OK]  transport  (${dur}s)"
    else
      echo "    [FAIL]  transport  (exit=$ec, ${dur}s)"
      continue
    fi
  else
    echo "[$strain] Step 2/4  find-transport — SKIP"
  fi

  # Step 3: draft
  # NOTE: gapseq draft 는 genome FASTA 를 받지 않는다 (-c 플래그 없음).
  #       -r/-t/-p 만으로 충분. biomass 는 기본 "auto" (Gram 자동 판별).
  #       실제 RDS 생성 여부를 exit code 대신 파일로 검증.
  if [ ! -s "${strain}-draft.RDS" ]; then
    echo "[$strain] Step 3/4  draft ..."
    t0=$(date +%s)
    gapseq draft \
      -r "${strain}-all-Reactions.tbl" \
      -t "${strain}-Transporter.tbl" \
      -p "${strain}-all-Pathways.tbl" \
      -n "${strain}" \
      > "$LOG_DIR/${strain}_draft.log" 2>&1
    ec=$?
    dur=$(( $(date +%s) - t0 ))
    if [ $ec -eq 0 ] && [ -s "${strain}-draft.RDS" ]; then
      echo "    [OK]  draft  (${dur}s)"
    else
      echo "    [FAIL]  draft  (exit=$ec, ${dur}s)  RDS=$([ -s "${strain}-draft.RDS" ] && echo OK || echo MISSING)"
      tail -10 "$LOG_DIR/${strain}_draft.log" | sed 's/^/        /'
      continue
    fi
  else
    echo "[$strain] Step 3/4  draft - SKIP"
  fi

  # Step 4: fill
  if [ ! -s "${strain}.xml" ] || [ "${strain}.xml" -ot "${strain}-draft.RDS" ]; then
    echo "[$strain] Step 4/4  fill  (media=$(basename $media)) ..."
    t0=$(date +%s)
    gapseq fill -m "${strain}-draft.RDS" -n "$media" \
      > "$LOG_DIR/${strain}_fill.log" 2>&1
    ec=$?
    dur=$(( $(date +%s) - t0 ))
    if [ $ec -eq 0 ] && [ -s "${strain}.xml" ]; then
      mu=$(grep -oP 'Final growth rate:\s*\K[0-9.]+' "$LOG_DIR/${strain}_fill.log" | tail -1)
      echo "    [OK]  fill  μ=${mu:-?}  (${dur}s)  → ${strain}.xml"
    else
      echo "    [FAIL]  fill  (exit=$ec, ${dur}s)"
      tail -10 "$LOG_DIR/${strain}_fill.log" | sed 's/^/        /'
    fi
  else
    echo "[$strain] Step 4/4  fill — SKIP (xml newer than draft)"
  fi
done

DUR=$(( $(date +%s) - START_ALL ))
echo
echo "======================================================================"
echo "  TOTAL: ${DUR}s ($((DUR/60))m $((DUR%60))s)"
echo "  End  : $(date)"
echo "======================================================================"
echo
echo "Summary:"
for strain in LRE LB LF LH LG EC; do
  xml="$WORK/$strain/${strain}.xml"
  if [ -s "$xml" ]; then
    mu=$(grep -oP 'Final growth rate:\s*\K[0-9.]+' "$LOG_DIR/${strain}_fill.log" 2>/dev/null | tail -1)
    echo "  ✓ $strain  μ=${mu:-?}  $(du -h "$xml" | cut -f1)"
  else
    echo "  ✗ $strain  (missing output)"
  fi
done

echo
echo "Next:"
echo "  cp $WORK/{LRE,LB,LF,LH,LG,EC}/*.xml /mnt/d/folder1/peptomatch/outputs/gem_cache/"

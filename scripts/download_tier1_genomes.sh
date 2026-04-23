#!/bin/bash
# ────────────────────────────────────────────────────────────────
# Tier 1 균주 6종 genome FASTA 다운로드 (NCBI datasets CLI 사용)
#
# 전제:
#   - WSL에 NCBI datasets CLI 설치되어 있음
#     설치: conda install -c conda-forge ncbi-datasets-cli
#     또는: curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
#   - 출력: /mnt/d/folder1/gapseq_linux/gapseq_work/{STRAIN}/{STRAIN}.fna
#
# 실행: bash download_tier1_genomes.sh
# ────────────────────────────────────────────────────────────────

set -u

WORK=/mnt/d/folder1/gapseq_linux/gapseq_work
mkdir -p "$WORK"

# strain → GCF accession 매핑
declare -A GENOMES=(
  [LRE]="GCF_000016825.1"   # Limosilactobacillus reuteri KCCM 40717
  [LB]="GCF_000056065.1"    # Lactobacillus delbrueckii subsp. bulgaricus KCCM 35463
  [LF]="GCF_013394085.1"    # Limosilactobacillus fermentum KCCM 35469
  [LH]="GCF_001434945.1"    # Lactobacillus helveticus KCCM 40989
  [LG]="GCF_040050875.1"    # Lactobacillus gasseri KCTC 3163
  [EC]="GCF_000022665.1"    # Escherichia coli BL21
)

echo "======================================================================"
echo "  Tier 1 genome downloader"
echo "======================================================================"
echo "  Work dir: $WORK"
echo

missing=0
for strain in LRE LB LF LH LG EC; do
  gcf="${GENOMES[$strain]}"
  out_dir="$WORK/$strain"
  fna="$out_dir/${strain}.fna"

  mkdir -p "$out_dir"

  if [ -s "$fna" ]; then
    echo "[$strain] SKIP (already downloaded: $fna)"
    continue
  fi

  echo "[$strain] downloading $gcf ..."

  # datasets CLI로 assembly 다운로드
  tmp_dir="$out_dir/_tmp_dl"
  rm -rf "$tmp_dir"
  mkdir -p "$tmp_dir"

  if datasets download genome accession "$gcf" \
       --include genome \
       --filename "$tmp_dir/${gcf}.zip" 2>&1 | tail -3; then
    unzip -q -o "$tmp_dir/${gcf}.zip" -d "$tmp_dir"
    # find the .fna file
    src_fna=$(find "$tmp_dir" -name "*.fna" | head -1)
    if [ -n "$src_fna" ]; then
      cp "$src_fna" "$fna"
      echo "    [OK]  $fna  ($(du -h "$fna" | cut -f1))"
    else
      echo "    [FAIL] .fna not found in downloaded archive"
      missing=$((missing+1))
    fi
  else
    echo "    [FAIL] datasets download failed"
    missing=$((missing+1))
  fi

  rm -rf "$tmp_dir"
done

echo
echo "======================================================================"
if [ $missing -eq 0 ]; then
  echo "  All 6 genomes downloaded successfully."
else
  echo "  WARNING: $missing genome(s) failed. Check errors above."
fi
echo "======================================================================"
echo
echo "Next: bash run_tier1_pipeline.sh"

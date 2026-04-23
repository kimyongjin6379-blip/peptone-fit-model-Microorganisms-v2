#!/bin/bash
# Download NCBI genomes for all strains in strain_registry.csv
# Usage: bash download_genomes.sh [output_dir]
#   default output_dir: ~/genomes

set -e

OUT_DIR="${1:-$HOME/genomes}"
REGISTRY="/mnt/d/folder1/peptomatch/data/strain_registry.csv"

if [[ ! -f "$REGISTRY" ]]; then
    echo "ERROR: registry not found at $REGISTRY" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"
echo "Downloading genomes into: $OUT_DIR"
echo ""

# NCBI FTP URL builder: splits GCF_000011045.1 → .../GCF/000/011/045/GCF_000011045.1_*
build_url() {
    local gcf="$1"   # e.g. GCF_000011045.1
    # split digits into 3/3/3 groups
    local digits="${gcf#GCF_}"    # 000011045.1
    local num="${digits%.*}"      # 000011045
    local d1="${num:0:3}"
    local d2="${num:3:3}"
    local d3="${num:6:3}"
    # Base dir
    local base="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${d1}/${d2}/${d3}"
    # We need the full directory name GCF_000011045.1_<ASM_NAME>
    # Use curl to list directory and find the match
    local dirname
    dirname=$(curl -s "${base}/" | grep -oE "${gcf}_[A-Za-z0-9._]+" | head -1)
    if [[ -z "$dirname" ]]; then
        return 1
    fi
    echo "${base}/${dirname}/${dirname}_genomic.fna.gz"
}

SKIP_COUNT=0
OK_COUNT=0
FAIL_COUNT=0
FAILED=()

# Skip header, iterate
tail -n +2 "$REGISTRY" | while IFS=',' read -r abbrev genus species strain_id gcf gram medium source has_emp priority notes; do
    [[ -z "$abbrev" || "$abbrev" =~ ^# ]] && continue

    out_file="$OUT_DIR/${abbrev}.fna"
    if [[ -f "$out_file" && -s "$out_file" ]]; then
        echo "  [$abbrev] already exists ($(du -h "$out_file" | cut -f1))"
        SKIP_COUNT=$((SKIP_COUNT+1))
        continue
    fi

    echo -n "  [$abbrev] ($gcf) ... "
    url=$(build_url "$gcf" 2>/dev/null) || {
        echo "FAIL (cannot resolve FTP path)"
        FAILED+=("$abbrev ($gcf)")
        FAIL_COUNT=$((FAIL_COUNT+1))
        continue
    }

    if wget -q -O "${out_file}.gz" "$url" && gunzip -f "${out_file}.gz"; then
        size=$(du -h "$out_file" | cut -f1)
        echo "OK ($size)"
        OK_COUNT=$((OK_COUNT+1))
    else
        echo "FAIL (download error)"
        rm -f "${out_file}.gz" "$out_file"
        FAILED+=("$abbrev ($gcf)")
        FAIL_COUNT=$((FAIL_COUNT+1))
    fi
done

echo ""
echo "=== Summary ==="
echo "Output dir: $OUT_DIR"
ls -lh "$OUT_DIR"/*.fna 2>/dev/null | wc -l | xargs echo "Total FASTA files:"

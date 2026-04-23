#!/bin/bash
# Batch gapseq pipeline for multiple strains using strain_registry.csv
#
# Usage:
#   bash gapseq_batch.sh                    # sequential, all priority=1 strains
#   bash gapseq_batch.sh --parallel 2       # run 2 strains in parallel
#   bash gapseq_batch.sh --priority 1       # only priority=1 (default)
#   bash gapseq_batch.sh --only LR,LA       # only these abbreviations
#   bash gapseq_batch.sh --skip LP          # skip these
#   bash gapseq_batch.sh --force            # re-run even if xml exists
#
# Assumptions:
#   - genomes already downloaded to ~/genomes/{abbrev}.fna
#   - gapseq ships with ALLmed.csv (fallback medium for fill)
#   - registry at /mnt/d/folder1/peptomatch/data/strain_registry.csv

set -u

REGISTRY="/mnt/d/folder1/peptomatch/data/strain_registry.csv"
GENOMES_DIR="$HOME/genomes"
WORK_BASE="$HOME/gapseq_work"
LOG_DIR="$WORK_BASE/logs"
GAPSEQ_MEDIA_DIR="$HOME/gapseq/dat/media"

# Medium fallback order (space-separated; arrays don't survive export to subshells)
# Order: MRS > ALLmed > TSB > LB > M9 > gut (extend as needed)
FALLBACK_MEDIA_STR="MRSmed.csv ALLmed.csv TSBmed.csv LBmed.csv M9.csv gut_low_fiber.csv"

# ── parse args ─────────────────────────────────────────────
PARALLEL=1
PRIORITY=1
ONLY=""
SKIP=""
FORCE=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --parallel) PARALLEL="$2"; shift 2 ;;
        --priority) PRIORITY="$2"; shift 2 ;;
        --only)     ONLY="$2";     shift 2 ;;
        --skip)     SKIP="$2";     shift 2 ;;
        --force)    FORCE=1;       shift ;;
        -h|--help)
            sed -n '2,16p' "$0" | sed 's/^#//'
            exit 0 ;;
        *) echo "Unknown arg: $1" >&2; exit 1 ;;
    esac
done

mkdir -p "$WORK_BASE" "$LOG_DIR"

# ── medium picker ──────────────────────────────────────────
pick_medium() {
    local preferred="$1"  # e.g. "MRS", "BHI", "LB"
    # Try preferred first (case insensitive + "med" suffix)
    for name in "${preferred}med.csv" "${preferred}.csv"; do
        if [[ -f "$GAPSEQ_MEDIA_DIR/$name" ]]; then
            echo "$GAPSEQ_MEDIA_DIR/$name"
            return 0
        fi
    done
    # Fallback chain (string → array split, works across subshells)
    local fb
    for fb in $FALLBACK_MEDIA_STR; do
        if [[ -f "$GAPSEQ_MEDIA_DIR/$fb" ]]; then
            echo "$GAPSEQ_MEDIA_DIR/$fb"
            return 0
        fi
    done
    return 1
}

# ── per-strain worker ──────────────────────────────────────
run_strain() {
    local abbrev="$1" gram="$2" medium="$3"
    local genome="$GENOMES_DIR/${abbrev}.fna"
    local work="$WORK_BASE/$abbrev"
    local log="$LOG_DIR/${abbrev}.log"
    local start_ts=$(date +%s)

    # Skip if already done (unless --force)
    if [[ $FORCE -eq 0 && -f "$work/${abbrev}.xml" && -s "$work/${abbrev}.xml" ]]; then
        echo "[$abbrev] SKIP (already done: $work/${abbrev}.xml)"
        return 0
    fi

    if [[ ! -f "$genome" ]]; then
        echo "[$abbrev] FAIL: genome not found at $genome"
        return 1
    fi

    # Map biomass flag: gram=pos/neg
    local biomass="$gram"

    # Pick medium file
    local med_file
    med_file=$(pick_medium "$medium") || {
        echo "[$abbrev] FAIL: no medium file found (tried $medium + fallbacks)"
        return 1
    }

    mkdir -p "$work" && cd "$work"
    cp -u "$genome" genome.fna

    echo "[$abbrev] START gram=$biomass medium=$(basename $med_file) $(date +%T)" | tee "$log"

    {
        echo "=== [$abbrev] find ==="
        time gapseq find -p all -b 200 genome.fna || exit 10

        echo "=== [$abbrev] find-transport ==="
        time gapseq find-transport -b 200 genome.fna || exit 11

        echo "=== [$abbrev] draft ==="
        gapseq draft \
            -r genome-all-Reactions.tbl \
            -t genome-Transporter.tbl \
            -p genome-all-Pathways.tbl \
            -b "$biomass" \
            -u 200 -l 100 \
            -n "$abbrev" || exit 12

        echo "=== [$abbrev] fill ==="
        time gapseq fill -m "${abbrev}-draft.RDS" -n "$med_file" || exit 13
    } >> "$log" 2>&1

    local ec=$?
    local dur=$(( $(date +%s) - start_ts ))

    if [[ $ec -eq 0 && -f "$work/${abbrev}.xml" && -s "$work/${abbrev}.xml" ]]; then
        echo "[$abbrev] DONE ($dur s) → $work/${abbrev}.xml"
        return 0
    else
        echo "[$abbrev] FAIL (exit=$ec, ${dur}s) — see $log"
        return 1
    fi
}

export -f run_strain pick_medium
export GENOMES_DIR WORK_BASE LOG_DIR GAPSEQ_MEDIA_DIR FORCE
export FALLBACK_MEDIA_STR

# ── dispatch ───────────────────────────────────────────────
echo "===================================================="
echo "gapseq batch runner"
echo "  registry : $REGISTRY"
echo "  parallel : $PARALLEL"
echo "  priority : $PRIORITY"
[[ -n "$ONLY" ]] && echo "  only     : $ONLY"
[[ -n "$SKIP" ]] && echo "  skip     : $SKIP"
[[ $FORCE -eq 1 ]] && echo "  FORCE    : re-run even if xml exists"
echo "===================================================="
echo ""

ONLY_ARR=()
SKIP_ARR=()
[[ -n "$ONLY" ]] && IFS=',' read -ra ONLY_ARR <<< "$ONLY"
[[ -n "$SKIP" ]] && IFS=',' read -ra SKIP_ARR <<< "$SKIP"

contains() {
    local needle="$1"; shift
    for x in "$@"; do [[ "$x" == "$needle" ]] && return 0; done
    return 1
}

# Build job list
JOBS=()
while IFS=',' read -r abbrev genus species strain_id gcf gram medium source has_emp prio notes; do
    [[ -z "$abbrev" || "$abbrev" == "abbrev" ]] && continue
    [[ "$prio" != "$PRIORITY" ]] && continue
    if [[ ${#ONLY_ARR[@]} -gt 0 ]] && ! contains "$abbrev" "${ONLY_ARR[@]}"; then continue; fi
    if [[ ${#SKIP_ARR[@]} -gt 0 ]] && contains "$abbrev" "${SKIP_ARR[@]}"; then continue; fi
    JOBS+=("$abbrev|$gram|$medium")
done < "$REGISTRY"

if [[ ${#JOBS[@]} -eq 0 ]]; then
    echo "No jobs matched filters. Exiting."
    exit 0
fi

echo "Jobs queued: ${#JOBS[@]}"
printf '  - %s\n' "${JOBS[@]}"
echo ""

START_ALL=$(date +%s)

# Run
if [[ $PARALLEL -gt 1 ]]; then
    if ! command -v parallel >/dev/null 2>&1; then
        echo "WARN: 'parallel' not installed — falling back to sequential"
        echo "      (install with: sudo apt install parallel)"
        PARALLEL=1
    fi
fi

if [[ $PARALLEL -gt 1 ]]; then
    printf '%s\n' "${JOBS[@]}" | parallel -j "$PARALLEL" --colsep '\|' \
        'run_strain {1} {2} {3}'
else
    for job in "${JOBS[@]}"; do
        IFS='|' read -r a g m <<< "$job"
        run_strain "$a" "$g" "$m"
    done
fi

DUR_ALL=$(( $(date +%s) - START_ALL ))

# ── report ─────────────────────────────────────────────────
echo ""
echo "===================================================="
echo "BATCH SUMMARY (total ${DUR_ALL}s = $((DUR_ALL/60))m $((DUR_ALL%60))s)"
echo "===================================================="
for job in "${JOBS[@]}"; do
    IFS='|' read -r a g m <<< "$job"
    xml="$WORK_BASE/$a/$a.xml"
    if [[ -f "$xml" && -s "$xml" ]]; then
        sz=$(du -h "$xml" | cut -f1)
        echo "  ✓ $a  ($sz)"
    else
        echo "  ✗ $a  (see $LOG_DIR/$a.log)"
    fi
done
echo "===================================================="
echo ""
echo "Next step: bash scripts/copy_gems.sh"

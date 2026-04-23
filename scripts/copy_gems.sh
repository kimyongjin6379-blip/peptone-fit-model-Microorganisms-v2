#!/bin/bash
# Copy completed gapseq xml files from WSL workdir to Windows PeptoMatch gem_cache.
# Usage: bash copy_gems.sh

set -e

SRC_BASE="$HOME/gapseq_work"
DST_DIR="/mnt/d/folder1/peptomatch/outputs/gem_cache"
REGISTRY="/mnt/d/folder1/peptomatch/data/strain_registry.csv"

mkdir -p "$DST_DIR"

echo "Copying GEM xml files..."
echo "  src: $SRC_BASE"
echo "  dst: $DST_DIR"
echo ""

COPIED=0
SKIPPED=0
MISSING=()

tail -n +2 "$REGISTRY" | while IFS=',' read -r abbrev rest; do
    [[ -z "$abbrev" || "$abbrev" == "abbrev" ]] && continue
    src="$SRC_BASE/$abbrev/$abbrev.xml"
    dst="$DST_DIR/$abbrev.xml"
    if [[ -f "$src" && -s "$src" ]]; then
        # Only copy if different or missing
        if [[ ! -f "$dst" ]] || ! cmp -s "$src" "$dst"; then
            cp "$src" "$dst"
            sz=$(du -h "$dst" | cut -f1)
            echo "  ✓ $abbrev ($sz)"
            COPIED=$((COPIED+1))
        else
            echo "  = $abbrev (unchanged)"
            SKIPPED=$((SKIPPED+1))
        fi
    else
        echo "  ✗ $abbrev (no xml at $src)"
        MISSING+=("$abbrev")
    fi
done

echo ""
echo "=== Summary ==="
echo "Destination contents:"
ls -lh "$DST_DIR"/*.xml 2>/dev/null | awk '{print "  " $9 "  " $5}'

echo ""
echo "Next step (run from Windows):"
echo "  cd D:\\folder1\\peptomatch"
echo "  python -c \"from peptomatch.gem_manager import verify_all; verify_all()\""

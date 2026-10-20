#!/usr/bin/env bash
# Check and extract one AlphaGenome Atlas download zip (downloaded by hand from
# https://deepmind.google.com/science/alphagenome/downloads; the site needs a Google login).
# Usage: prepare_alphagenome_downloads.sh check|extract <zip>
#   check    zip complete, readable central directory, contents listed, sha256 recorded
#   extract  unzip next to the zip into <zip_stem>/ (never overwrites; CRCs are verified while
#            extracting), then print each table's first lines and its index files
# Sidecars written next to the zip: <zip>.contents.txt, <zip>.sha256. The zip itself is kept.
set -euo pipefail

MODE="${1:?usage: prepare_alphagenome_downloads.sh check|extract <zip>}"
ZIP="${2:?usage: prepare_alphagenome_downloads.sh check|extract <zip>}"
DIR="$(cd "$(dirname "$ZIP")" && pwd)"
NAME="$(basename "$ZIP")"
STEM="${NAME%.zip}"

# A browser download still in progress leaves a hidden .part file and a 0-byte target
if compgen -G "$DIR/.*${STEM}*.part*" > /dev/null || [[ ! -s "$ZIP" ]]; then
    echo "$ZIP: download not finished (partial file present or zip empty)" >&2; exit 1
fi

case "$MODE" in
check)
    # Fails on a truncated zip: the central directory sits at the end of the file
    unzip -l "$ZIP" > "$DIR/$NAME.contents.txt"
    cat "$DIR/$NAME.contents.txt"
    (cd "$DIR" && sha256sum "$NAME") > "$DIR/$NAME.sha256"
    echo "size_bytes $(stat -c %s "$ZIP")"
    cat "$DIR/$NAME.sha256"
    ;;
extract)
    [[ -s "$DIR/$NAME.sha256" ]] || { echo "run 'check' first" >&2; exit 1; }
    OUT="$DIR/$STEM"
    mkdir -p "$OUT"
    unzip -n -q "$ZIP" -d "$OUT"
    # Show each extracted table's header and first records, and which index files exist
    find "$OUT" -type f | sort | while read -r f; do
        echo "== $f ($(stat -c %s "$f") bytes)"
        case "$f" in
            *.tbi|*.csi) ;;
            *.gz|*.bgz) zcat "$f" 2>/dev/null | head -5 || true ;;
            *) head -5 "$f" || true ;;
        esac
    done
    ;;
*) echo "unknown mode: $MODE" >&2; exit 1 ;;
esac

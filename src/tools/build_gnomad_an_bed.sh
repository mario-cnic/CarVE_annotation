#!/usr/bin/env bash
#$ -P BIGN
#$ -N build_gnomad_an_bed
#$ -A PGP
#$ -pe smp 4
#$ -l h_vmem=4G
#$ -o _log/build_gnomad_an_bed.stdout
#$ -e _log/build_gnomad_an_bed.stderr
#
# Download gnomAD v4.1 "allele number at all callable sites" (exomes, genomes) and convert each
# to a bgzipped, tabix-indexed BED: chrom, start (0-based), end, AN. Consecutive positions with
# the same AN are merged into one interval; positions missing from the source get no interval.
# Re-running resumes the download and rebuilds a BED only if its index is missing.
#
# Usage (from the repo root): bash src/tools/build_gnomad_an_bed.sh [--out-dir DIR] [--only exomes|genomes]
# Output: $OUT_DIR/gnomad.{exomes,genomes}.v4.1.allele_number_all_sites.bed.gz(.tbi), $OUT_DIR/MANIFEST.txt

set -euo pipefail

BASE_URL="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/tsv"
OUT_DIR="/data_lab_PGP/resources/annotation/gnomad_an/v4.1"
ENV_TABIX="/data_lab_PGP/shared/utils/conda_envs/spliceai_env/bin"
THREADS="${NSLOTS:-4}"
SETS=(exomes genomes)

while [[ $# -gt 0 ]]; do
    case "$1" in
        --out-dir) OUT_DIR="$2"; shift 2 ;;
        --only) SETS=("$2"); shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done

mkdir -p "$OUT_DIR/raw"
echo "Start $(date -Is) host=$(hostname) out=$OUT_DIR sets=${SETS[*]}"

for s in "${SETS[@]}"; do
    name="gnomad.$s.v4.1.allele_number_all_sites"
    raw="$OUT_DIR/raw/$name.tsv.bgz"
    bed="$OUT_DIR/$name.bed.gz"
    url="$BASE_URL/$s/$name.tsv.bgz"

    remote_size=$(curl -sI "$url" | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}')
    local_size=$(stat -c %s "$raw" 2>/dev/null || echo 0)
    if [[ "$local_size" != "$remote_size" ]]; then
        curl -sS --fail --retry 5 -C - -o "$raw" "$url"
        local_size=$(stat -c %s "$raw")
    fi
    [[ "$local_size" == "$remote_size" ]] || { echo "ERROR: $raw is $local_size bytes, remote $remote_size" >&2; exit 1; }
    "$ENV_TABIX/bgzip" -t "$raw"
    echo "$s downloaded: $local_size bytes $(date -Is)"

    if [[ -s "$bed.tbi" ]]; then
        echo "$s BED already indexed, skipped"
        continue
    fi
    # Fails if a chromosome reappears or positions are not strictly increasing.
    zcat "$raw" | awk -F'\t' -v OFS='\t' '
        NR == 1 { if ($1 != "locus" || $2 != "AN") { print "bad header: " $0 > "/dev/stderr"; bad = 1; exit 1 } next }
        {
            i = index($1, ":"); c = substr($1, 1, i - 1); p = substr($1, i + 1) + 0; a = $2
            if (c != cc) {
                if (cc != "") print cc, s - 1, e, an
                if (c in seen) { print "chromosome not contiguous: " c > "/dev/stderr"; bad = 1; exit 1 }
                seen[c] = 1; cc = c; s = p; e = p; an = a; next
            }
            if (p <= e) { print "position not increasing: " $1 > "/dev/stderr"; bad = 1; exit 1 }
            if (p == e + 1 && a == an) { e = p; next }
            print cc, s - 1, e, an; s = p; e = p; an = a
        }
        END { if (bad) exit 1; if (cc != "") print cc, s - 1, e, an }' \
        | "$ENV_TABIX/bgzip" -@ "$THREADS" > "$bed.tmp"
    mv "$bed.tmp" "$bed"
    "$ENV_TABIX/tabix" -p bed "$bed"
    echo "$s BED built $(date -Is)"
done

{
    echo "source	$BASE_URL/{exomes,genomes}/gnomad.*.v4.1.allele_number_all_sites.tsv.bgz"
    echo "release	gnomAD v4.1 (GRCh38)"
    echo "format	BED chrom, start (0-based), end, AN; runs of equal AN merged; no interval = position not in source"
    echo "finished	$(date -Is)"
    echo "host	$(hostname)"
    echo "tabix	$("$ENV_TABIX/tabix" --version | head -1)"
    for f in "$OUT_DIR"/raw/*.tsv.bgz "$OUT_DIR"/*.bed.gz; do
        [[ -e "$f" ]] && echo "file	$(basename "$f")	$(stat -c %s "$f")	sha256:$(sha256sum "$f" | cut -d' ' -f1)"
    done
} > "$OUT_DIR/MANIFEST.txt"
cat "$OUT_DIR/MANIFEST.txt"
echo "Done $(date -Is)"

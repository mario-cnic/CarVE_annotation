#!/usr/bin/env bash
# Download the current ClinVar GRCh38 VCF from NCBI into a local directory, verify it against
# NCBI's md5, and name it clinvar_<fileDate>.vcf.gz (plus .tbi and a .md5 sidecar).
# Usage: build_clinvar_vcf.sh <output_dir>
# Intended output: /data_lab_PGP/resources/annotation/clinvar/GRCh38 (readable by the cluster nodes).
# Never writes to /references.
set -euo pipefail

OUT_DIR="${1:?usage: build_clinvar_vcf.sh <output_dir>}"
BCFTOOLS="${BCFTOOLS:-bcftools}"
BASE_URL="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38"

mkdir -p "$OUT_DIR"
WORK="$(mktemp -d "$OUT_DIR/.clinvar_dl.XXXXXX")"
trap 'rm -rf "$WORK"' EXIT

curl -sSfo "$WORK/clinvar.vcf.gz" "$BASE_URL/clinvar.vcf.gz"
curl -sSfo "$WORK/clinvar.vcf.gz.md5" "$BASE_URL/clinvar.vcf.gz.md5"

# Compare the downloaded file with NCBI's published checksum
expected="$(cut -d' ' -f1 "$WORK/clinvar.vcf.gz.md5")"
actual="$(md5sum "$WORK/clinvar.vcf.gz" | cut -d' ' -f1)"
[[ "$expected" == "$actual" ]] || { echo "md5 mismatch: expected $expected, got $actual" >&2; exit 1; }

# Release date and assembly come from the file's own header
file_date="$("$BCFTOOLS" view -h "$WORK/clinvar.vcf.gz" | sed -n 's/^##fileDate=//p' | tr -d '-')"
reference="$("$BCFTOOLS" view -h "$WORK/clinvar.vcf.gz" | sed -n 's/^##reference=//p')"
[[ -n "$file_date" && "$reference" == "GRCh38" ]] || { echo "unexpected header: fileDate='$file_date' reference='$reference'" >&2; exit 1; }

NAME="clinvar_${file_date}.vcf.gz"
[[ ! -e "$OUT_DIR/$NAME" ]] || { echo "$OUT_DIR/$NAME already exists; not overwriting" >&2; exit 1; }
mv "$WORK/clinvar.vcf.gz" "$OUT_DIR/$NAME"
"$BCFTOOLS" index -t "$OUT_DIR/$NAME"
printf '%s  %s\n' "$actual" "$NAME" > "$OUT_DIR/$NAME.md5"

cat <<EOF
Built $OUT_DIR/$NAME (ClinVar release $file_date, md5 $actual).
Point clinvar_vcf in nextflow.config at that file. A newer release is a new file name; old files stay.
EOF

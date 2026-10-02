#!/bin/bash
# Builds one multi-sample VCF from a WGS VCF and a WES VCF of the same family, to be annotated once
# with annotate_vcf.nf. Runs outside the pipeline; HPC-sized inputs, run it on the cluster.
#
# Steps:
#   1. normalise the WES VCF (bcftools norm -m -any -f FASTA) into a temporary BCF
#   2. keep the WES records that carry an alt allele in any WES sample, plus the WES records with the
#      same CHROM/POS/REF/ALT as a WGS record (their hom-ref / no-call genotypes are evidence)
#   3. bcftools merge -m none WGS + selected WES
#   4. store each source's original FILTER as INFO/FILTER_WGS and INFO/FILTER_WES (FILTER "." is
#      written as UNFILTERED; a record absent from a source has no tag for it), then clear FILTER
#
# A record that one source lacks gets ./. for that source's samples: "not called", not hom-ref.
# Records of the normalised WES file with identical CHROM/POS/REF/ALT (different original records
# that normalise to the same allele) are reduced to the first one; only the first is kept (chr22 test: 806 of 89,733 selected records, none losing an alt allele).
# The WGS VCF is assumed normalised (check it with the pipeline's input check). Nothing is
# overwritten: the script stops if an output already exists.
#
# Usage:
#   prepare_family_union.sh --wgs WGS.vcf.gz --wes WES.vcf.gz --fasta REF.fa --outdir DIR \
#       [--prefix NAME] [--region chr22] [--tmpdir DIR] [--threads 4] [--bcftools PATH]
set -euo pipefail
ORIG_ARGS="$*"

WGS="" WES="" FASTA="" OUTDIR="" PREFIX="family_union" REGION="" TMPBASE="${TMPDIR:-/tmp}" THREADS=4 BCFTOOLS="bcftools"
while [[ $# -gt 0 ]]; do
	case "$1" in
		--wgs) WGS="$2"; shift 2 ;;
		--wes) WES="$2"; shift 2 ;;
		--fasta) FASTA="$2"; shift 2 ;;
		--outdir) OUTDIR="$2"; shift 2 ;;
		--prefix) PREFIX="$2"; shift 2 ;;
		--region) REGION="$2"; shift 2 ;;
		--tmpdir) TMPBASE="$2"; shift 2 ;;
		--threads) THREADS="$2"; shift 2 ;;
		--bcftools) BCFTOOLS="$2"; shift 2 ;;
		-h|--help) sed -n '2,22p' "$0"; exit 0 ;;
		*) echo "unknown argument: $1" >&2; exit 2 ;;
	esac
done
for v in WGS WES FASTA OUTDIR; do
	[[ -n "${!v}" ]] || { echo "missing --$(echo "$v" | tr 'A-Z' 'a-z')" >&2; exit 2; }
done
for f in "$WGS" "$WES" "$FASTA" "$WGS.tbi" "$WES.tbi" "$FASTA.fai"; do
	[[ -s "$f" ]] || { echo "not found or empty: $f" >&2; exit 2; }
done

OUT="$OUTDIR/$PREFIX.vcf.gz"
MANIFEST="$OUTDIR/$PREFIX.PREPARE_MANIFEST.tsv"
mkdir -p "$OUTDIR"
for f in "$OUT" "$OUT.tbi" "$MANIFEST"; do
	[[ ! -e "$f" ]] || { echo "refusing to overwrite existing $f" >&2; exit 2; }
done

TMP="$(mktemp -d "$TMPBASE/prepare_family_union.XXXXXX")"
trap 'rm -rf "$TMP"' EXIT
# bgzip and tabix are expected next to bcftools when it is given by path
[[ "$BCFTOOLS" != */* ]] || PATH="$(dirname "$BCFTOOLS"):$PATH"
bt() { local sub="$1"; shift; "$BCFTOOLS" "$sub" --threads "$THREADS" "$@"; }
RARG=()
[[ -z "$REGION" ]] || RARG=(-r "$REGION")
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
log() { echo "[prepare_family_union] $(date -Iseconds) $*" >&2; }
count() { $BCFTOOLS index -n "$1"; }

log "normalising WES"
bt norm -f "$FASTA" -m -any -c w -Ob ${RARG[@]+"${RARG[@]}"} "$WES" -o "$TMP/wes_norm.bcf"
$BCFTOOLS index -f "$TMP/wes_norm.bcf"
WES_NORM_N=$(count "$TMP/wes_norm.bcf")

log "selecting WES records with an alt allele or identical to a WGS record"
bt view -c1:nref -Ob "$TMP/wes_norm.bcf" -o "$TMP/wes_alt.bcf"
WGS_FOR_ISEC="$WGS"
if [[ -n "$REGION" ]]; then
	bt view -r "$REGION" -Ob "$WGS" -o "$TMP/wgs_region.bcf"
	$BCFTOOLS index -f "$TMP/wgs_region.bcf"
	WGS_FOR_ISEC="$TMP/wgs_region.bcf"
fi
# isec keeps WES records with the same CHROM, POS, REF and ALT as a WGS record (their hom-ref and no-call genotypes)
$BCFTOOLS isec -n=2 -w2 -Ob -o "$TMP/wes_at_wgs.bcf" "$WGS_FOR_ISEC" "$TMP/wes_norm.bcf"
$BCFTOOLS index -f "$TMP/wes_alt.bcf"
$BCFTOOLS index -f "$TMP/wes_at_wgs.bcf"
bt concat -a -D -Ob "$TMP/wes_alt.bcf" "$TMP/wes_at_wgs.bcf" -o "$TMP/wes_sel.bcf"
$BCFTOOLS index -f "$TMP/wes_sel.bcf"
WES_SEL_N=$(count "$TMP/wes_sel.bcf")
WES_ALT_N=$(count "$TMP/wes_alt.bcf")
WES_AT_WGS_N=$(count "$TMP/wes_at_wgs.bcf")

log "original FILTER tables"
for src in WGS WES; do
	if [[ $src == WGS ]]; then in="$WGS"; ra=(${RARG[@]+"${RARG[@]}"}); else in="$TMP/wes_sel.bcf"; ra=(); fi
	$BCFTOOLS query ${ra[@]+"${ra[@]}"} -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\n' "$in" \
		| awk 'BEGIN{OFS="\t"} {if ($5 == ".") $5 = "UNFILTERED"; print}' | bgzip -c > "$TMP/filter_$src.tsv.gz"
	tabix -f -s1 -b2 -e2 "$TMP/filter_$src.tsv.gz"
	printf '##INFO=<ID=FILTER_%s,Number=1,Type=String,Description="FILTER value of this record in the %s VCF (UNFILTERED = FILTER was .)">\n' "$src" "$src" > "$TMP/hdr_$src.txt"
done

log "merging"
bt merge -m none -Ob "$WGS_FOR_ISEC" "$TMP/wes_sel.bcf" -o "$TMP/merged.bcf"
$BCFTOOLS index -f "$TMP/merged.bcf"

log "writing FILTER_WGS / FILTER_WES and clearing FILTER"
$BCFTOOLS annotate -a "$TMP/filter_WGS.tsv.gz" -h "$TMP/hdr_WGS.txt" -c CHROM,POS,REF,ALT,INFO/FILTER_WGS -Ou "$TMP/merged.bcf" \
	| bt annotate -a "$TMP/filter_WES.tsv.gz" -h "$TMP/hdr_WES.txt" -c CHROM,POS,REF,ALT,INFO/FILTER_WES -x FILTER -Oz -o "$OUT" -
tabix -p vcf "$OUT"

OUT_N=$(count "$OUT")
SAMPLES=$($BCFTOOLS query -l "$OUT" | paste -sd,)
log "done: $OUT_N records, samples $SAMPLES"

{
	printf 'key\tvalue\n'
	printf 'generated_at\t%s\n' "$(date -Iseconds)"
	printf 'host\t%s\n' "$(hostname)"
	printf 'command\t%s\n' "$0 $ORIG_ARGS"
	printf 'script_commit\t%s\n' "$(git -C "$REPO_ROOT" rev-parse HEAD 2>/dev/null || echo unknown)"
	printf 'script_dirty\t%s\n' "$(git -C "$REPO_ROOT" status --porcelain -- src/tools/prepare_family_union.sh 2>/dev/null | wc -l)"
	printf 'bcftools\t%s\n' "$($BCFTOOLS --version | head -1)"
	printf 'region\t%s\n' "${REGION:-all}"
	printf 'fasta\t%s\n' "$FASTA"
	for src in WGS WES; do
		f="${!src}"
		printf '%s_path\t%s\n%s_size_bytes\t%s\n%s_sha256\t%s\n%s_records\t%s\n' "$src" "$f" "$src" "$(stat -c %s "$f")" "$src" "$(sha256sum "$f" | cut -d' ' -f1)" "$src" "$(count "$f")"
	done
	printf 'wes_records_after_norm\t%s\n' "$WES_NORM_N"
	printf 'wes_records_selected\t%s\n' "$WES_SEL_N"
	printf 'wes_records_with_alt\t%s\nwes_records_identical_to_wgs\t%s\n' "$WES_ALT_N" "$WES_AT_WGS_N"
	printf 'output\t%s\n' "$OUT"
	printf 'output_records\t%s\n' "$OUT_N"
	printf 'output_sha256\t%s\n' "$(sha256sum "$OUT" | cut -d' ' -f1)"
	printf 'output_samples\t%s\n' "$SAMPLES"
} > "$MANIFEST"
log "manifest: $MANIFEST"

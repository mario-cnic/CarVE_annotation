#!/bin/bash
# Launches annotate_vcf.nf and writes the run-provenance manifest before and after it. Every argument
# is passed to `nextflow` unchanged, except that --run_id is always resolved here and forwarded, so
# the manifest location and Nextflow's output folder use the same value.
#
# Usage (same arguments you'd already pass to `nextflow run annotate_vcf.nf` directly):
#   ./run_annotate_vcf.sh -profile standard --input_vcf sample.vcf.gz [--run_id label] \
#     [--spip_tier restricted] [--output_format pq] [-resume]
set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST_PY="$REPO_ROOT/src/python/write_run_manifest.py"
MAIN_NF="$REPO_ROOT/annotate_vcf.nf"
CONFIG="$REPO_ROOT/nextflow.config"

# `nextflow` is not guaranteed to be on PATH on this cluster (confirmed by hitting
# "nextflow: command not found" on a real run) — it lives at /opt/nextflow/nextflow. Resolve once
# here and use this everywhere a bare `nextflow` call would otherwise silently fail or (worse,
# since this script has no `set -e`) silently produce empty provenance values instead of erroring.
NEXTFLOW_BIN="$(command -v nextflow || true)"
if [[ -z "$NEXTFLOW_BIN" && -x /opt/nextflow/nextflow ]]; then
	NEXTFLOW_BIN=/opt/nextflow/nextflow
fi
if [[ -z "$NEXTFLOW_BIN" ]]; then
	echo "error: nextflow not found on PATH or at /opt/nextflow/nextflow — set NEXTFLOW_BIN or add it to PATH" >&2
	exit 2
fi

ARGS=("$@")
INPUT="" RUN_ID="" PROFILE="standard" SPIP_TIER="" OUTPUT_FORMAT="" RESUME=0
FORWARD_ARGS=()
i=0
while [[ $i -lt ${#ARGS[@]} ]]; do
	case "${ARGS[$i]}" in
	--input_vcf) INPUT="${ARGS[$((i + 1))]}"; FORWARD_ARGS+=("${ARGS[$i]}" "${ARGS[$((i + 1))]}"); i=$((i + 1)) ;;
	# --run_id is deliberately dropped from FORWARD_ARGS here, not forwarded as-typed: it's
	# re-added once, below, from $RUN_ID — the one place this script computes it — so Nextflow
	# never derives its own value in parallel (see the note below this loop for why that matters).
	--run_id) RUN_ID="${ARGS[$((i + 1))]}"; i=$((i + 1)) ;;
	-profile) PROFILE="${ARGS[$((i + 1))]}"; FORWARD_ARGS+=("${ARGS[$i]}" "${ARGS[$((i + 1))]}"); i=$((i + 1)) ;;
	--spip_tier) SPIP_TIER="${ARGS[$((i + 1))]}"; FORWARD_ARGS+=("${ARGS[$i]}" "${ARGS[$((i + 1))]}"); i=$((i + 1)) ;;
	--output_format) OUTPUT_FORMAT="${ARGS[$((i + 1))]}"; FORWARD_ARGS+=("${ARGS[$i]}" "${ARGS[$((i + 1))]}"); i=$((i + 1)) ;;
	-resume) RESUME=1; FORWARD_ARGS+=("${ARGS[$i]}") ;;
	*) FORWARD_ARGS+=("${ARGS[$i]}") ;;
	esac
	i=$((i + 1))
done

if [[ -z "$INPUT" ]]; then
	echo "error: --input_vcf is required (same as annotate_vcf.nf itself needs it)" >&2
	exit 2
fi

# run_id defaults to the input file name without .gz/.vcf; passed to Nextflow explicitly below.
if [[ -z "$RUN_ID" ]]; then
	RUN_ID="$(basename "$INPUT")"
	RUN_ID="${RUN_ID%.gz}"
	RUN_ID="${RUN_ID%.vcf}"
fi
FORWARD_ARGS+=(--run_id "$RUN_ID")

# Resolve the same binary paths/config values this run will actually use, from the SAME
# nextflow.config the real run reads — not re-derived or hardcoded here.
# The config is dumped once (each `nextflow config` call is a JVM start) and queried from memory.
CONFIG_DUMP="$("$NEXTFLOW_BIN" config -profile "$PROFILE" "$REPO_ROOT" 2>/dev/null)"
resolve_param() {
	printf '%s\n' "$CONFIG_DUMP" | \
		grep -E "^\s*$1\s*=" | head -1 | sed -E "s/^\s*$1\s*=\s*'?([^']*)'?\s*\$/\1/"
}
BCFTOOLS="$(resolve_param bcftools)"
VEP_SIF="$(resolve_param vep_sif)"
VEP_CACHE_VERSION="$(resolve_param vep_cache_version)"
SPLICEAI_DISTANCE="$(resolve_param spliceai_distance)"
PANGOLIN_DISTANCE="$(resolve_param pangolin_distance)"

# Identity of everything outside this repo that shapes the output (reference FASTA, Pangolin db,
# SPiP script, LaBranchoR BED, VEP cache dir; the shared/utils code repo; each predictor's conda
# env). Names come from nextflow.config params, so the manifest records what THIS run resolved.
# A param that resolves to nothing is omitted, not recorded as an empty path.
EXTERNAL_ARGS=()
add_external() { # kind (--resource|--code-repo|--env) label param
	local v; v="$(resolve_param "$3")"
	[[ -n "$v" ]] && EXTERNAL_ARGS+=("$1" "$2=$v")
}
for p in fasta pangolin_db spip_script labranchor_bed vep_dir clinvar_vcf gnomad_vcf splicevardb_vcf splicevault_file vep_cadd vep_utrannotator vep_spliceai_snv vep_spliceai_indel vep_revel vep_mastermind vep_loftool vep_pli vep_dbnsfp vep_dbnsfp_fields vep_alphamissense vep_maxentscan_dir vep_remm gene_restriction_bed hgnc_table spliceai_symbol_map enst_spip_map dbnsfp_aligned_columns; do
	add_external --resource "$p" "$p"
done
add_external --code-repo shared_utils shared_utils
for p in python_spliceai pangolin_python spip_rscript bcftools vcf_parser_python datasci_python; do
	add_external --env "$p" "$p"
done

OUTDIR="$REPO_ROOT/nf_work/annotation_out/$RUN_ID"
INFO_DIR="$OUTDIR/pipeline_info"
FORWARD_ARGS+=(--provenance_stamp "$INFO_DIR/TABLE_STAMP.json")
NF_LOG="$INFO_DIR/nextflow_$(date +%Y-%m-%d_%H-%M-%S).log"

RESUME_FLAG=()
[[ $RESUME -eq 1 ]] && RESUME_FLAG=(-resume)

echo "[run_annotate_vcf.sh] writing launch manifest to $INFO_DIR/RUN_MANIFEST.json"
python3 "$MANIFEST_PY" launch \
	--outdir "$OUTDIR" \
	--command-line "$NEXTFLOW_BIN -log $NF_LOG run annotate_vcf.nf ${FORWARD_ARGS[*]}" \
	--nextflow-main "$MAIN_NF" \
	--nextflow-binary "$NEXTFLOW_BIN" \
	--config "$CONFIG" \
	--input "$INPUT" \
	--run-id "$RUN_ID" \
	--profile "$PROFILE" \
	--spip-tier "$SPIP_TIER" \
	--output-format "$OUTPUT_FORMAT" \
	--bcftools "$BCFTOOLS" \
	--vep-sif "$VEP_SIF" \
	--vep-cache-version "$VEP_CACHE_VERSION" \
	--spliceai-distance "$SPLICEAI_DISTANCE" \
	--pangolin-distance "$PANGOLIN_DISTANCE" \
	${EXTERNAL_ARGS[@]+"${EXTERNAL_ARGS[@]}"} \
	"${RESUME_FLAG[@]:+--resume}" || {
	echo "error: failed to write launch manifest — aborting before launching nextflow" >&2
	exit 1
}

echo "[run_annotate_vcf.sh] launching: $NEXTFLOW_BIN -log $NF_LOG run annotate_vcf.nf ${FORWARD_ARGS[*]}"
cd "$REPO_ROOT" && "$NEXTFLOW_BIN" -log "$NF_LOG" run annotate_vcf.nf "${FORWARD_ARGS[@]}"
EXIT_CODE=$?

echo "[run_annotate_vcf.sh] nextflow exited $EXIT_CODE — writing completion manifest"
python3 "$MANIFEST_PY" completed \
	--outdir "$OUTDIR" \
	--exit-code "$EXIT_CODE" \
	--nextflow-log "$NF_LOG"

# Run QC report (record counts, time and memory per process, failed attempts) into <run>/reports/.
# A failure here is reported but does not change the run's exit code.
DATASCI_PY="$(resolve_param datasci_python)"
"${DATASCI_PY:-python3}" "$REPO_ROOT/src/python/build_run_qc_report.py" \
	--run-dir "$OUTDIR" --nextflow-log "$NF_LOG" --bcftools "$BCFTOOLS" ||
	echo "warning: run QC report failed; the run itself is unaffected" >&2

# The manifest, any RUN_DIRTY patch and the Nextflow log/trace/timeline/report stay in the run's
# pipeline_info folder; run records are not added to this repository.
echo "[run_annotate_vcf.sh] done. Manifest: $INFO_DIR/RUN_MANIFEST.json"

exit $EXIT_CODE

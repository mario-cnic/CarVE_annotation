#!/bin/bash
#$ -P NEW
#$ -N vcf2parsed
#$ -A PGP
#$ -pe smp 1
#$ -l h_vmem=120G
#$ -o _log/master/parse.stdout
#$ -e _log/master/parse.stderr

# Unified VCF-to-Parsed Clean Table Stage
# Converts annotated VCF directly to final clean Parquet/TSV/Excel table
# using local scratch storage without retaining massive intermediate TSVs on disk.

set -euo pipefail

INPUT=$1
OUTPUT=$2
GENE_SUBSET=${3:-""}

VCF_PARSER=/data_lab_PGP/shared/utils/src/vcf_parser_pysam.py
FILTER_VARIANTS=/data_lab_PGP/shared/utils/src/filter_variants.py

if [ ! -f "$VCF_PARSER" ]; then
    VCF_PARSER=/home/mruizp/data_lab_PGP/shared/utils/src/vcf_parser_pysam.py
fi
if [ ! -f "$FILTER_VARIANTS" ]; then
    FILTER_VARIANTS=/home/mruizp/data_lab_PGP/shared/utils/src/filter_variants.py
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

if [ -f "$PIPELINE_ROOT/resources/all_but_old_gnomad_vep_cols.txt" ]; then
    VEP_COLS_FILE="$PIPELINE_ROOT/resources/all_but_old_gnomad_vep_cols.txt"
elif [ -f "/data_lab_PGP/shared/utils/data/vcf_all_columns.txt" ]; then
    VEP_COLS_FILE="/data_lab_PGP/shared/utils/data/vcf_all_columns.txt"
else
    VEP_COLS_FILE=/home/mruizp/data_lab_PGP/shared/utils/data/vcf_all_columns.txt
fi

LOG_LEVEL=INFO
FILTER_BY=Feature

# Determine local high-speed scratch space for ephemeral stream
SCRATCH_DIR="/data_tmp"
if [ ! -d "$SCRATCH_DIR" ] || [ ! -w "$SCRATCH_DIR" ]; then
    SCRATCH_DIR="/tmp"
fi

TMP_SCRATCH_TSV="${SCRATCH_DIR}/vcf_parse_${$}_$(date +%s%N).tsv"

# Ensure cleanup of scratch file on exit or error
cleanup() {
    if [ -f "$TMP_SCRATCH_TSV" ]; then
        rm -f "$TMP_SCRATCH_TSV"
    fi
}
trap cleanup EXIT INT TERM

# Ensure parent directory of output exists
mkdir -p "$(dirname "$OUTPUT")"

# Limit threads to prevent resource/process limits exhaustion on concurrent cluster slots
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export ARROW_IO_THREADS=1

VCF_PARSER_PYTHONPATH="/data_lab_PGP/shared/utils/conda_envs/vcf_parser/lib/python3.13/site-packages:/home/mruizp/data_lab_PGP/shared/utils/conda_envs/vcf_parser/lib/python3.13/site-packages:${PYTHONPATH:-}"

# ----------------- Logger Setup -----------------
LOGGER_SCRIPT="src/hpc/logger.sh"
if [ -f "$LOGGER_SCRIPT" ] && [ -n "${LOG_DIR:-}" ]; then
    source "$LOGGER_SCRIPT"
    log_step_start "vcf2parsed" "$INPUT"
    TIME_LOG="$LOG_DIR/vcf2parsed.time"
    TIME_CMD=$(get_time_cmd "$TIME_LOG")
else
    TIME_CMD=""
fi
if ! command -v run_command_timed &>/dev/null; then
    run_command_timed() { "$@"; }
fi
# ------------------------------------------------

VCF_PARSER_PYTHON=/data_lab_PGP/shared/utils/conda_envs/vcf_parser/bin/python
if [ ! -x "$VCF_PARSER_PYTHON" ]; then
    VCF_PARSER_PYTHON=/home/mruizp/data_lab_PGP/shared/utils/conda_envs/vcf_parser/bin/python
fi
if [ ! -x "$VCF_PARSER_PYTHON" ]; then
    VCF_PARSER_PYTHON=python3
fi

DATASCI_PYTHON=/data_lab_PGP/shared/utils/conda_envs/datasci/bin/python
if [ ! -x "$DATASCI_PYTHON" ]; then
    DATASCI_PYTHON=/home/mruizp/data_lab_PGP/shared/utils/conda_envs/datasci/bin/python
fi
if [ ! -x "$DATASCI_PYTHON" ]; then
    DATASCI_PYTHON=python3
fi

echo "============================================================================"
echo " ⚡ Streamlined Direct VCF-to-Parsed Table Converter"
echo "============================================================================"
echo " Input VCF     : $INPUT"
echo " Target Output : $OUTPUT"
echo " Feature Target: ${GENE_SUBSET:-'All transcripts'}"
echo " Scratch Buffer: $TMP_SCRATCH_TSV"
echo "============================================================================"

# Step A: Parse Annotated VCF into fast local scratch buffer
echo "Stage 1/2: Parsing annotated VCF fields..."
PYTHONPATH="$VCF_PARSER_PYTHONPATH" run_command_timed $VCF_PARSER_PYTHON $VCF_PARSER \
    --input "$INPUT" \
    --output "$TMP_SCRATCH_TSV" \
    --vep_columns "$VEP_COLS_FILE" \
    --add_info --add_vep --overwrite \
    --logging_level "$LOG_LEVEL" \
    ${GENE_SUBSET:+--gene_set "$GENE_SUBSET"} \
    --filter_by "$FILTER_BY"

# Step B: Filter and transform directly into final clean Parquet/TSV/Excel table
echo "Stage 2/2: Applying status contracts, splicing metrics, and saving $OUTPUT..."
run_command_timed $DATASCI_PYTHON $FILTER_VARIANTS \
    --input "$TMP_SCRATCH_TSV" \
    --output "$OUTPUT" \
    --logging_level "$LOG_LEVEL" \
    --skip_quality_filter

CMD_EXIT_CODE=$?

# ----------------- Logger End -----------------
if [ -f "$LOGGER_SCRIPT" ] && [ -n "${LOG_DIR:-}" ]; then
    log_step_end "vcf2parsed" "$OUTPUT" "$CMD_EXIT_CODE" "$TIME_LOG"
    log_final
fi
# ----------------------------------------------

echo "Direct VCF-to-Parsed conversion completed successfully. Output saved to $OUTPUT"

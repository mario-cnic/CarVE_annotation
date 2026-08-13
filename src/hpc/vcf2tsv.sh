#!/bin/bash
#$ -P NEW
#$ -N parse
#$ -A PGP
#$ -pe smp 1
#$ -l h_vmem=120G
#$ -o _log/master/parse.stdout
#$ -e _log/master/parse.stderr

INPUT=$1
OUTPUT=$2
GENE_SUBSET=$3

# [AUTO-UPDATE] Path updated to shared/utils
VCF_PARSER=/data_lab_PGP/shared/utils/src/vcf_parser_pysam.py
# [AUTO-UPDATE] Path updated to shared/utils
FILTER_VARIANTS=/data_lab_PGP/shared/utils/src/filter_variants.py
if [ -f "resources/all_but_old_gnomad_vep_cols.txt" ]; then
    COLUMNS="resources/all_but_old_gnomad_vep_cols.txt"
else
    COLUMNS=/data_lab_PGP/shared/utils/data/vcf_all_columns.txt
fi
LOG_LEVEL=INFO
FILTER_BY=Feature

# Filter TSV into final CSV
mkdir -p ./results/master
mkdir -p ./_log/master
# Do not source ~/.bashrc in headless SGE cluster jobs

# Limit threads to prevent resource/process limits exhaustion on concurrent cluster slots
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export ARROW_IO_THREADS=1
# ----------------- Logger Setup -----------------
LOGGER_SCRIPT="src/hpc/logger.sh"
if [ -f "$LOGGER_SCRIPT" ] && [ -n "$LOG_DIR" ]; then
    source "$LOGGER_SCRIPT"
    log_step_start "tsv" "$INPUT"
    TIME_LOG="$LOG_DIR/tsv.time"
    TIME_CMD=$(get_time_cmd "$TIME_LOG")
else
    TIME_CMD=""
fi
if ! command -v run_command_timed &>/dev/null; then
    run_command_timed() { "$@"; }
fi
# ------------------------------------------------

echo "Filtering variants from $INPUT and outputting to $OUTPUT"
# [AUTO-UPDATE] Path updated to shared/utils
run_command_timed /data_lab_PGP/shared/utils/conda_envs/vcf_parser/bin/python $VCF_PARSER \
	--input $INPUT \
 	--output $OUTPUT \
 	--vep_columns $COLUMNS \
 	--add_info --add_vep --overwrite \
 	--logging_level $LOG_LEVEL \
 	--gene_set $GENE_SUBSET \
	--filter_by $FILTER_BY
CMD_EXIT_CODE=$?

echo "Filtering complete. Output saved to $OUTPUT"

# ----------------- Logger End -----------------
if [ -f "$LOGGER_SCRIPT" ] && [ -n "$LOG_DIR" ]; then
    log_step_end "tsv" "$OUTPUT" "$CMD_EXIT_CODE" "$TIME_LOG"
fi
# ----------------------------------------------
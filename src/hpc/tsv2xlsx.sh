#!/bin/bash
#$ -P NEW
#$ -N parse
#$ -A PGP
#$ -pe smp 1
#$ -l h_vmem=120G
#$ -o _log/parse_excel.stdout
#$ -e _log/parse_excel.stderr

INPUT=$1
OUT_CSV=$2
# [AUTO-UPDATE] Path updated to shared/utils
FILTER_VARIANTS=/data_lab_PGP/shared/utils/src/filter_variants.py
LOG_LEVEL=INFO

#  ------------------------------------------------------------------------------------------

CLASSIFICATION=dominant
GENE_CATEGORY=all
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
    log_step_start "tsv2xlsx" "$INPUT"
    TIME_LOG="$LOG_DIR/tsv2xlsx.time"
    TIME_CMD=$(get_time_cmd "$TIME_LOG")
else
    TIME_CMD=""
fi
if ! command -v run_command_timed &>/dev/null; then
    run_command_timed() { "$@"; }
fi
# ------------------------------------------------

# [AUTO-UPDATE] Path updated to shared/utils
run_command_timed /data_lab_PGP/shared/utils/conda_envs/datasci/bin/python $FILTER_VARIANTS \
	--input $INPUT \
	--output $OUT_CSV \
	--logging_level $LOG_LEVEL \
	--skip_quality_filter
CMD_EXIT_CODE=$?

# ----------------- Logger End -----------------
if [ -f "$LOGGER_SCRIPT" ] && [ -n "$LOG_DIR" ]; then
    log_step_end "tsv2xlsx" "$OUT_CSV" "$CMD_EXIT_CODE" "$TIME_LOG"
    # This is the final step for this gene, compile the execution report
    log_final
fi
# ----------------------------------------------

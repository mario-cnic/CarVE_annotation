#!/bin/bash
#$ -P PGP
#$ -N BRANCHPOINT
#$ -A PGP
#$ -l thread=2
#$ -l h_vmem=8G
#$ -o _log/master/annotation_branchpointer.stdout
#$ -e _log/master/annotation_branchpointer.stderr

INPUT=$1
OUTPUT=$2

SPLICEAI_ENV=/data_lab_PGP/shared/utils/conda_envs/spliceai_env
if [ ! -d "$SPLICEAI_ENV" ]; then
    SPLICEAI_ENV=/home/mruizp/conda_envs/spliceai_env
fi

PYTHON_BIN=$SPLICEAI_ENV/bin/python3
if [ ! -x "$PYTHON_BIN" ]; then
    PYTHON_BIN="/data_lab_PGP/shared/utils/conda_envs/datasci/bin/python3"
fi
if [ ! -x "$PYTHON_BIN" ]; then
    PYTHON_BIN="python3"
fi

BGZIP_BIN=$SPLICEAI_ENV/bin/bgzip
if [ ! -x "$BGZIP_BIN" ]; then
    BGZIP_BIN="/data_lab_PGP/shared/utils/conda_envs/genomics/bin/bgzip"
fi
if [ ! -x "$BGZIP_BIN" ]; then
    BGZIP_BIN="bgzip"
fi

TABIX_BIN=$SPLICEAI_ENV/bin/tabix
if [ ! -x "$TABIX_BIN" ]; then
    TABIX_BIN="/data_lab_PGP/shared/utils/conda_envs/genomics/bin/tabix"
fi
if [ ! -x "$TABIX_BIN" ]; then
    TABIX_BIN="tabix"
fi

SCRIPT=src/python/annotate_branchpointer.py
BED_REF=/data_lab_PGP/resources/annotation/labranchor/labranchor_grch38_top.bed.gz
if [ ! -f "$BED_REF" ]; then
    BED_REF=/home/mruizp/data_lab_PGP/resources/annotation/labranchor/labranchor_grch38_top.bed.gz
fi

RAW_OUT="${OUTPUT%.gz}"

# ----------------- Logger Setup -----------------
LOGGER_SCRIPT="src/hpc/logger.sh"
if [ -f "$LOGGER_SCRIPT" ] && [ -n "$LOG_DIR" ]; then
    source "$LOGGER_SCRIPT"
    log_step_start "branchpoint" "$INPUT"
    TIME_LOG="$LOG_DIR/branchpoint.time"
    TIME_CMD=$(get_time_cmd "$TIME_LOG")
else
    TIME_CMD=""
fi
# ------------------------------------------------

echo "Running Branch Point Predictor Engine (Branchpointer + LaBranchoR) on $INPUT..."
if [ -n "$TIME_CMD" ]; then
    $TIME_CMD $PYTHON_BIN $SCRIPT \
        "$INPUT" \
        "$RAW_OUT" \
        --labranchor-bed "$BED_REF"
    CMD_EXIT_CODE=$?
else
    $PYTHON_BIN $SCRIPT \
        "$INPUT" \
        "$RAW_OUT" \
        --labranchor-bed "$BED_REF"
    CMD_EXIT_CODE=$?
fi

echo "Compressing and indexing output VCF..."
if [[ "$OUTPUT" == *.gz ]]; then
    if [ -f "$RAW_OUT" ]; then
        $BGZIP_BIN -f "$RAW_OUT"
    fi
    FINAL_OUT="$OUTPUT"
else
    if [ -f "$OUTPUT" ]; then
        $BGZIP_BIN -f "$OUTPUT"
    fi
    FINAL_OUT="${OUTPUT}.gz"
fi
$TABIX_BIN -f "$FINAL_OUT"

echo "Branchpoint annotation completed. Final output: $FINAL_OUT"

# ----------------- Logger End -----------------
if [ -f "$LOGGER_SCRIPT" ] && [ -n "$LOG_DIR" ]; then
    log_step_end "branchpoint" "$FINAL_OUT" "$CMD_EXIT_CODE" "$TIME_LOG"
fi
# ----------------------------------------------

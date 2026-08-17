#!/bin/bash
#$ -P PGP
#$ -N SPLICEAI
#$ -A PGP
#$ -l thread=4
#$ -l h_vmem=16G
#$ -o _log/master/annotation_spliceai.stdout
#$ -e _log/master/annotation_spliceai.stderr

INPUT=$1
OUTPUT=$2

FASTA=/home/mruizp/data_references/genomes/Homo_sapiens/GATK_bundle/v0/Homo_sapiens_assembly38.fasta
if [ ! -f "$FASTA" ]; then
    FASTA=/references/genomes/Homo_sapiens/GATK_bundle/v0/Homo_sapiens_assembly38.fasta
fi

SPLICEAI_ENV=/data_lab_PGP/shared/utils/conda_envs/spliceai_env
if [ ! -d "$SPLICEAI_ENV" ]; then
    SPLICEAI_ENV=/home/mruizp/conda_envs/spliceai_env
fi

export TF_NUM_INTEROP_THREADS=2
export TF_NUM_INTRAOP_THREADS=2
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2
export TF_ENABLE_ONEDNN_OPTS=0

# ----------------- Logger Setup -----------------
LOGGER_SCRIPT="src/hpc/logger.sh"
if [ -f "$LOGGER_SCRIPT" ] && [ -n "$LOG_DIR" ]; then
    source "$LOGGER_SCRIPT"
    log_step_start "spliceai" "$INPUT"
    TIME_LOG="$LOG_DIR/spliceai.time"
    TIME_CMD=$(get_time_cmd "$TIME_LOG")
else
    TIME_CMD=""
fi
# ------------------------------------------------

RAW_OUT="${OUTPUT%.gz}"
PYTHON_BIN=$SPLICEAI_ENV/bin/python3
SPLICEAI_SCRIPT=src/python/annotate_spliceai.py

echo "Running SpliceAI Splicing Predictor (-D 10000) on $INPUT..."
if [ -n "$TIME_CMD" ]; then
    $TIME_CMD $PYTHON_BIN $SPLICEAI_SCRIPT \
        "$INPUT" \
        "$RAW_OUT" \
        "$FASTA" \
        -d 10000
    CMD_EXIT_CODE=$?
else
    $PYTHON_BIN $SPLICEAI_SCRIPT \
        "$INPUT" \
        "$RAW_OUT" \
        "$FASTA" \
        -d 10000
    CMD_EXIT_CODE=$?
fi

echo "Compressing and indexing output VCF..."
if [[ "$OUTPUT" == *.gz ]]; then
    RAW_OUT="${OUTPUT%.gz}"
    if [ -f "$RAW_OUT" ]; then
        $SPLICEAI_ENV/bin/bgzip -f "$RAW_OUT"
    fi
    FINAL_OUT="$OUTPUT"
else
    if [ -f "$OUTPUT" ]; then
        $SPLICEAI_ENV/bin/bgzip -f "$OUTPUT"
    fi
    FINAL_OUT="${OUTPUT}.gz"
fi
$SPLICEAI_ENV/bin/tabix -f "$FINAL_OUT"

echo "SpliceAI annotation completed. Final output: $FINAL_OUT"

# ----------------- Logger End -----------------
if [ -f "$LOGGER_SCRIPT" ] && [ -n "$LOG_DIR" ]; then
    log_step_end "spliceai" "${OUTPUT}.gz" "$CMD_EXIT_CODE" "$TIME_LOG"
fi
# ----------------------------------------------

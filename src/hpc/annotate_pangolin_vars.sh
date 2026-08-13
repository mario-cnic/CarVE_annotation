#!/bin/bash
#$ -P PGP
#$ -N PANGOLIN
#$ -A PGP
#$ -l thread=4
#$ -l h_vmem=20G
#$ -o _log/master/annotation_pangolin.stdout
#$ -e _log/master/annotation_pangolin.stderr

INPUT=$1
OUTPUT=$2

FASTA=/home/mruizp/data_references/genomes/Homo_sapiens/GATK_bundle/v0/Homo_sapiens_assembly38.fasta
if [ ! -f "$FASTA" ]; then
    FASTA=/references/genomes/Homo_sapiens/GATK_bundle/v0/Homo_sapiens_assembly38.fasta
fi

DB=/data_lab_PGP/shared/utils/pangolin_db/pangolin_grch38.db
if [ ! -f "$DB" ]; then
    DB=/home/mruizp/data_lab_PGP/shared/utils/pangolin_db/pangolin_grch38.db
fi

PANGOLIN_ENV=/data_lab_PGP/shared/utils/conda_envs/pangolin_env
if [ ! -d "$PANGOLIN_ENV" ]; then
    PANGOLIN_ENV=/home/mruizp/data_lab_PGP/shared/utils/conda_envs/pangolin_env
fi

export OMP_NUM_THREADS=${THREADS:-4}
export MKL_NUM_THREADS=${THREADS:-4}
export OPENBLAS_NUM_THREADS=${THREADS:-4}

# Decompress if input is gzipped
if [[ $INPUT == *.gz ]]; then
    gunzip -k -f $INPUT
    INPUT=${INPUT%.gz}
fi

# ----------------- Logger Setup -----------------
LOGGER_SCRIPT="src/hpc/logger.sh"
if [ -f "$LOGGER_SCRIPT" ] && [ -n "$LOG_DIR" ]; then
    source "$LOGGER_SCRIPT"
    log_step_start "pangolin" "$INPUT"
    TIME_LOG="$LOG_DIR/pangolin.time"
    TIME_CMD=$(get_time_cmd "$TIME_LOG")
else
    TIME_CMD=""
fi
# ------------------------------------------------

PYTHON_BIN=$PANGOLIN_ENV/bin/python3
PANGOLIN_REPO=src/external/Pangolin-main
BGZIP_BIN=$PANGOLIN_ENV/bin/bgzip
if [ ! -x "$BGZIP_BIN" ]; then
    BGZIP_BIN=bgzip
fi
TABIX_BIN=$PANGOLIN_ENV/bin/tabix
if [ ! -x "$TABIX_BIN" ]; then
    TABIX_BIN=tabix
fi

export PYTHONPATH="$PANGOLIN_REPO:$PANGOLIN_ENV/lib/python3.12/site-packages:$PANGOLIN_ENV/lib/python3.10/site-packages:${PYTHONPATH:-}"

echo "Running Pangolin Splicing Predictor on $INPUT..."
if [ -n "$TIME_CMD" ]; then
    $TIME_CMD $PYTHON_BIN -m pangolin.pangolin \
        "$INPUT" \
        "$FASTA" \
        "$DB" \
        "$OUTPUT" \
        -d 10000
    CMD_EXIT_CODE=$?
else
    $PYTHON_BIN -m pangolin.pangolin \
        "$INPUT" \
        "$FASTA" \
        "$DB" \
        "$OUTPUT" \
        -d 10000
    CMD_EXIT_CODE=$?
fi

echo "Compressing and indexing output VCF..."
$BGZIP_BIN -f "$OUTPUT"
$TABIX_BIN -f "${OUTPUT}.gz"

echo "Pangolin annotation completed. Final output: ${OUTPUT}.gz"

# ----------------- Logger End -----------------
if [ -f "$LOGGER_SCRIPT" ] && [ -n "$LOG_DIR" ]; then
    log_step_end "pangolin" "${OUTPUT}.gz" "$CMD_EXIT_CODE" "$TIME_LOG"
fi
# ----------------------------------------------

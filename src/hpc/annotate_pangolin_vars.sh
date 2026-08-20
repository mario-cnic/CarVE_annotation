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
if [ ! -d "$PANGOLIN_ENV" ] || [ ! -f "$PANGOLIN_ENV/lib/python3.12/site-packages/numpy/__init__.py" ]; then
    PANGOLIN_ENV=/home/mruizp/conda_envs/pangolin_env
fi
if [ ! -d "$PANGOLIN_ENV" ]; then
    PANGOLIN_ENV=/home/mruizp/data_lab_PGP/shared/utils/conda_envs/pangolin_env
fi

export OMP_NUM_THREADS=${THREADS:-4}
export MKL_NUM_THREADS=${THREADS:-4}
export OPENBLAS_NUM_THREADS=${THREADS:-4}

LOCAL_INPUT="$INPUT"
CLEAN_LOCAL_INPUT=false

# Decompress safely to a job-isolated temporary file if input is gzipped
if [[ $INPUT == *.gz ]]; then
    LOCAL_INPUT="${OUTPUT}.input_tmp.vcf"
    gunzip -c "$INPUT" > "$LOCAL_INPUT"
    CLEAN_LOCAL_INPUT=true
fi

# ----------------- Logger Setup -----------------
LOGGER_SCRIPT="src/hpc/logger.sh"
if [ -f "$LOGGER_SCRIPT" ] && [ -n "${LOG_DIR:-}" ]; then
    source "$LOGGER_SCRIPT"
    log_step_start "pangolin" "$INPUT"
    TIME_LOG="$LOG_DIR/pangolin.time"
    TIME_CMD=$(get_time_cmd "$TIME_LOG")
else
    TIME_CMD=""
fi
# ------------------------------------------------

PYTHON_BIN=$PANGOLIN_ENV/bin/python3
PANGOLIN_REPO="src/external/Pangolin-main"
BGZIP_BIN=$PANGOLIN_ENV/bin/bgzip
if [ ! -x "$BGZIP_BIN" ]; then
    BGZIP_BIN=bgzip
fi
TABIX_BIN=$PANGOLIN_ENV/bin/tabix
if [ ! -x "$TABIX_BIN" ]; then
    TABIX_BIN=tabix
fi

export PYTHONPATH="$PANGOLIN_ENV/lib/python3.12/site-packages:$PANGOLIN_REPO:${PYTHONPATH:-}"

# Ensure output directory exists
mkdir -p "$(dirname "$OUTPUT")"
RAW_OUT="${OUTPUT%.gz}"

echo "Running Pangolin Splicing Predictor on $LOCAL_INPUT..."
if [ -n "$TIME_CMD" ]; then
    $TIME_CMD $PYTHON_BIN -m pangolin.pangolin \
        "$LOCAL_INPUT" \
        "$FASTA" \
        "$DB" \
        "$RAW_OUT" \
        -d 10000
    CMD_EXIT_CODE=$?
else
    $PYTHON_BIN -m pangolin.pangolin \
        "$LOCAL_INPUT" \
        "$FASTA" \
        "$DB" \
        "$RAW_OUT" \
        -d 10000
    CMD_EXIT_CODE=$?
fi

# Clean isolated temporary uncompressed input
if [ "$CLEAN_LOCAL_INPUT" = true ] && [ -f "$LOCAL_INPUT" ]; then
    rm -f "$LOCAL_INPUT"
fi

echo "Compressing and indexing output VCF..."
if [ -f "$RAW_OUT" ]; then
    $BGZIP_BIN -f "$RAW_OUT"
    $TABIX_BIN -f "${RAW_OUT}.gz"
fi

echo "Pangolin annotation completed. Final output: ${RAW_OUT}.gz"

# ----------------- Logger End -----------------
if [ -f "$LOGGER_SCRIPT" ] && [ -n "${LOG_DIR:-}" ]; then
    log_step_end "pangolin" "${RAW_OUT}.gz" "$CMD_EXIT_CODE" "$TIME_LOG"
fi
# ----------------------------------------------

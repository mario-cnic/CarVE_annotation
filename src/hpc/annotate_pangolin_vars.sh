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

PANGOLIN_ENV=""
for candidate in "/home/mruizp/conda_envs/pangolin_env" "/data_lab_PGP/shared/utils/conda_envs/pangolin_env" "/home/mruizp/data_lab_PGP/shared/utils/conda_envs/pangolin_env"; do
    if [ -x "$candidate/bin/python3" ] && "$candidate/bin/python3" -c "import pysam" >/dev/null 2>&1; then
        PANGOLIN_ENV="$candidate"
        break
    fi
done

if [ -z "$PANGOLIN_ENV" ]; then
    PANGOLIN_ENV="/home/mruizp/conda_envs/pangolin_env"
fi

export TMPDIR=/tmp
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=${THREADS:-4}
export MKL_NUM_THREADS=${THREADS:-4}
export OPENBLAS_NUM_THREADS=${THREADS:-4}

LOCAL_INPUT="$INPUT"
CLEAN_LOCAL_INPUT=false

# Decompress safely to a job-isolated temporary file if input is gzipped
if [[ $INPUT == *.gz ]]; then
    LOCAL_INPUT="${OUTPUT}.input_tmp_$$.vcf"
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

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PANGOLIN_REPO="$SCRIPT_DIR/src/external/Pangolin-main"
PYTHON_BIN=$PANGOLIN_ENV/bin/python3
BGZIP_BIN="$PANGOLIN_ENV/bin/bgzip"
if [ ! -x "$BGZIP_BIN" ]; then
    BGZIP_BIN=$(which bgzip || echo "bgzip")
fi
TABIX_BIN="$PANGOLIN_ENV/bin/tabix"
if [ ! -x "$TABIX_BIN" ]; then
    TABIX_BIN=$(which tabix || echo "tabix")
fi

export PATH="$PANGOLIN_ENV/bin:$PATH"
export PYTHONPATH="$PANGOLIN_REPO:${PYTHONPATH:-}"

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
    $BGZIP_BIN -f -@ 4 "$RAW_OUT"
    $TABIX_BIN -f "${RAW_OUT}.gz"
fi

echo "Pangolin annotation completed. Final output: ${RAW_OUT}.gz"

# ----------------- Logger End -----------------
if [ -f "$LOGGER_SCRIPT" ] && [ -n "${LOG_DIR:-}" ]; then
    log_step_end "pangolin" "${RAW_OUT}.gz" "$CMD_EXIT_CODE" "$TIME_LOG"
fi
# ----------------------------------------------

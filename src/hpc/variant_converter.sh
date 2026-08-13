#!/bin/bash
#$ -P PGP
#$ -N variant_converter
#$ -A PGP
#$ -pe smp 1
#$ -l h_vmem=10G
#$ -o _log/variant_converter.stdout
#$ -e _log/variant_converter.stderr

set -e

# Default Parameters
INPUT=""
OUTPUT=""
BUILD="GRCh38"
BATCH_SIZE="2500"
COLUMN_PARAM=""

# Do not source ~/.bashrc in headless SGE cluster jobs

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export ARROW_IO_THREADS=1

# Conda Python Environment
PYTHON_BIN="/data_lab_PGP/shared/utils/conda_envs/liftover/bin/python"
if [ ! -x "$PYTHON_BIN" ]; then
    PYTHON_BIN="python3"
fi

POSITIONAL_COUNT=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --id_columns*|id_columns*|--id-columns*|id-columns*)
            if [[ "$1" == *" "* ]]; then
                COLUMN_PARAM="$1"
                if [[ "$COLUMN_PARAM" != --* ]]; then
                    COLUMN_PARAM="--$COLUMN_PARAM"
                fi
                shift 1
            else
                COLUMN_PARAM="--id_columns"
                shift 1
                while [[ $# -gt 0 && ! "$1" =~ ^- ]]; do
                    COLUMN_PARAM="$COLUMN_PARAM $1"
                    shift 1
                done
            fi
            ;;
        --cdna_column*|cdna_column*|--cdna-column*|cdna-column*)
            if [[ "$1" == *" "* ]]; then
                COLUMN_PARAM="$1"
                if [[ "$COLUMN_PARAM" != --* ]]; then
                    COLUMN_PARAM="--$COLUMN_PARAM"
                fi
                shift 1
            else
                COLUMN_PARAM="--cdna_column"
                shift 1
                while [[ $# -gt 0 && ! "$1" =~ ^- ]]; do
                    COLUMN_PARAM="$COLUMN_PARAM $1"
                    shift 1
                done
            fi
            ;;
        --build)
            BUILD="$2"
            shift 2
            ;;
        --batch_size)
            BATCH_SIZE="$2"
            shift 2
            ;;
        -*)
            echo "Unknown option $1"
            exit 1
            ;;
        *)
            if [[ $POSITIONAL_COUNT -eq 0 ]]; then
                INPUT="$1"
                POSITIONAL_COUNT=1
            elif [[ $POSITIONAL_COUNT -eq 1 ]]; then
                OUTPUT="$1"
                POSITIONAL_COUNT=2
            fi
            shift
            ;;
    esac
done

if [ -z "$INPUT" ] || [ -z "$OUTPUT" ]; then
    echo "Usage: $0 <input_file> <output_vcf> [options]"
    exit 1
fi

LOGGER_SCRIPT="src/hpc/logger.sh"
if [ -f "$LOGGER_SCRIPT" ] && [ -n "${LOG_DIR:-}" ]; then
    source "$LOGGER_SCRIPT"
    log_step_start "variant_converter" "$INPUT"
    TIME_LOG="$LOG_DIR/variant_converter.time"
    TIME_CMD=$(get_time_cmd "$TIME_LOG")
else
    TIME_CMD=""
fi

echo "Running universal_variant_converter.py on $INPUT -> $OUTPUT (Build: $BUILD)"

CONVERTER_SCRIPT="src/python/universal_variant_converter.py"

# Build execution flags
FLAGS="--input \"$INPUT\" --output \"$OUTPUT\" --build \"$BUILD\" --batch-size \"$BATCH_SIZE\" --overwrite"

if [[ "$COLUMN_PARAM" == *"--id_columns "* ]]; then
    ID_COLS="${COLUMN_PARAM/--id_columns /}"
    FLAGS="$FLAGS --id-columns $ID_COLS"
elif [[ "$COLUMN_PARAM" == *"--cdna_column "* ]]; then
    CDNA_COL="${COLUMN_PARAM/--cdna_column /}"
    FLAGS="$FLAGS --cdna-column $CDNA_COL"
fi

if [ -n "$TIME_CMD" ]; then
    eval $TIME_CMD $PYTHON_BIN $CONVERTER_SCRIPT $FLAGS
    CMD_EXIT_CODE=$?
else
    eval $PYTHON_BIN $CONVERTER_SCRIPT $FLAGS
    CMD_EXIT_CODE=$?
fi

# Sort and index output VCF if successful
if [[ "$OUTPUT" == *.vcf && $CMD_EXIT_CODE -eq 0 ]]; then
    BCFTOOLS="/data_lab_PGP/shared/utils/conda_envs/genomics/bin/bcftools"
    TABIX="/data_lab_PGP/shared/utils/conda_envs/genomics/bin/tabix"
    if [ -x "$BCFTOOLS" ] && [ -x "$TABIX" ]; then
        $BCFTOOLS sort -Oz -o "${OUTPUT}.gz" "$OUTPUT"
        $TABIX -f "${OUTPUT}.gz"
    fi
fi

if [ -f "$LOGGER_SCRIPT" ] && [ -n "${LOG_DIR:-}" ]; then
    log_step_end "variant_converter" "$OUTPUT" "$CMD_EXIT_CODE" "${TIME_LOG:-}"
fi

exit $CMD_EXIT_CODE
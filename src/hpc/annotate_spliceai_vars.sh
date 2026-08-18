#!/bin/bash
#$ -P PGP
#$ -N SPLICEAI
#$ -A PGP
#$ -l thread=4
#$ -l h_vmem=20G
#$ -o _log/master/annotation_spliceai.stdout
#$ -e _log/master/annotation_spliceai.stderr

set -euo pipefail

INPUT=$1
OUTPUT=$2

FASTA=/references/genomes/Homo_sapiens/GATK_bundle/v0/Homo_sapiens_assembly38.fasta
if [ ! -f "$FASTA" ]; then
    FASTA=/home/mruizp/data_references/genomes/Homo_sapiens/GATK_bundle/v0/Homo_sapiens_assembly38.fasta
fi

SPLICEAI_ENV=/data_lab_PGP/shared/utils/conda_envs/spliceai_env
if [ ! -d "$SPLICEAI_ENV" ]; then
    SPLICEAI_ENV=/home/mruizp/conda_envs/spliceai_env
fi

GENOMICS_ENV=/data_lab_PGP/shared/utils/conda_envs/genomics
if [ ! -d "$GENOMICS_ENV" ]; then
    GENOMICS_ENV=/home/mruizp/apps/miniforge3/envs/genomics
fi

PYTHON_BIN=$SPLICEAI_ENV/bin/python3
if [ ! -x "$PYTHON_BIN" ]; then
    PYTHON_BIN=python3
fi

BGZIP_BIN=$SPLICEAI_ENV/bin/bgzip
if [ ! -x "$BGZIP_BIN" ]; then
    BGZIP_BIN=$GENOMICS_ENV/bin/bgzip
fi
if [ ! -x "$BGZIP_BIN" ]; then
    BGZIP_BIN=bgzip
fi

TABIX_BIN=$SPLICEAI_ENV/bin/tabix
if [ ! -x "$TABIX_BIN" ]; then
    TABIX_BIN=$GENOMICS_ENV/bin/tabix
fi
if [ ! -x "$TABIX_BIN" ]; then
    TABIX_BIN=tabix
fi

BCFTOOLS_BIN=$GENOMICS_ENV/bin/bcftools
if [ ! -x "$BCFTOOLS_BIN" ]; then
    BCFTOOLS_BIN=bcftools
fi

SPLICEAI_SCRIPT=src/python/annotate_spliceai.py
CHUNK_SCRIPT=src/python/split_vcf_chunks.py
CHUNK_THRESHOLD=15000
CHUNK_SIZE=20000

export TF_NUM_INTEROP_THREADS=2
export TF_NUM_INTRAOP_THREADS=2
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2
export TF_ENABLE_ONEDNN_OPTS=0

# Determine final output paths
if [[ "$OUTPUT" == *.gz ]]; then
    FINAL_OUT="$OUTPUT"
    RAW_OUT="${OUTPUT%.gz}"
else
    FINAL_OUT="${OUTPUT}.gz"
    RAW_OUT="$OUTPUT"
fi

# ----------------- Logger Setup -----------------
LOGGER_SCRIPT="src/hpc/logger.sh"
if [ -f "$LOGGER_SCRIPT" ] && [ -n "${LOG_DIR:-}" ]; then
    source "$LOGGER_SCRIPT"
    log_step_start "spliceai" "$INPUT"
    TIME_LOG="$LOG_DIR/spliceai.time"
    TIME_CMD=$(get_time_cmd "$TIME_LOG")
else
    TIME_CMD=""
fi
# ------------------------------------------------

echo "============================================================================"
echo " Starting SpliceAI Neural Net Annotation (-D 10000) on: $INPUT"
echo "============================================================================"

# Count variants to decide between Single Worker vs Parallel Chunking
NUM_VARS=$($BCFTOOLS_BIN view -H "$INPUT" 2>/dev/null | wc -l || true)
if [ -z "$NUM_VARS" ]; then
    NUM_VARS=0
fi

echo "Total variant records detected: $NUM_VARS"

if [ "$NUM_VARS" -ge "$CHUNK_THRESHOLD" ]; then
    echo "High variant volume locus detected (>= $CHUNK_THRESHOLD). Activating Parallel VCF Chunker..."
    
    CHUNK_SCRATCH="${RAW_OUT}_spliceai_chunks"
    mkdir -p "$CHUNK_SCRATCH"
    
    # Check if raw chunks already exist or need splitting
    EXISTING_CHUNKS=$(find "$CHUNK_SCRATCH" -maxdepth 1 -name "chunk_*.vcf.gz" ! -name "*.annSpliceAI.*" | wc -l)
    if [ "$EXISTING_CHUNKS" -eq 0 ]; then
        echo "Splitting $INPUT into chunks of $CHUNK_SIZE records in $CHUNK_SCRATCH..."
        $PYTHON_BIN $CHUNK_SCRIPT --input "$INPUT" --output-dir "$CHUNK_SCRATCH" --chunk-size "$CHUNK_SIZE"
    fi
    
    MAX_WORKERS=${NSLOTS:-4}
    if [ "$MAX_WORKERS" -lt 1 ]; then
        MAX_WORKERS=4
    fi
    echo "Executing parallel inference with $MAX_WORKERS concurrent worker processes..."
    
    worker_pids=()
    for chunk_file in "$CHUNK_SCRATCH"/chunk_[0-9][0-9][0-9][0-9].vcf.gz; do
        [ -f "$chunk_file" ] || continue
        c_base=$(basename "$chunk_file" .vcf.gz)
        c_out_raw="$CHUNK_SCRATCH/${c_base}.annSpliceAI.vcf"
        c_out_gz="${c_out_raw}.gz"
        
        # Checkpoint: Skip pre-computed chunks
        if [ -s "$c_out_gz" ]; then
            echo "  [Chunk Checkpoint] Already completed: $c_base"
            continue
        fi
        
        (
            echo "  [Worker Pool] Starting SpliceAI on $c_base..."
            $PYTHON_BIN $SPLICEAI_SCRIPT "$chunk_file" "$c_out_raw" "$FASTA" -d 10000
            $BGZIP_BIN -f "$c_out_raw"
            $TABIX_BIN -f "$c_out_gz"
            echo "  [Worker Pool] Finished: $c_base"
        ) &
        
        worker_pids+=($!)
        if [ ${#worker_pids[@]} -ge "$MAX_WORKERS" ]; then
            wait "${worker_pids[0]}"
            worker_pids=("${worker_pids[@]:1}")
        fi
    done
    
    # Wait for all remaining background workers
    for pid in "${worker_pids[@]}"; do
        wait "$pid"
    done
    
    echo "All chunk workers finished successfully. Concatenating annotated chunks into $FINAL_OUT..."
    
    # Sort and collect all chunk outputs
    annotated_chunks=()
    while IFS= read -r f; do
        annotated_chunks+=("$f")
    done < <(find "$CHUNK_SCRATCH" -maxdepth 1 -name "chunk_*.annSpliceAI.vcf.gz" | sort)
    
    if [ ${#annotated_chunks[@]} -eq 0 ]; then
        echo "Error: No annotated chunk outputs found in $CHUNK_SCRATCH"
        exit 1
    fi
    
    $BCFTOOLS_BIN concat -a "${annotated_chunks[@]}" -Oz -o "$FINAL_OUT"
    $TABIX_BIN -f "$FINAL_OUT"
    
    # Verify concatenated output record count matches input
    FINAL_COUNT=$($BCFTOOLS_BIN view -H "$FINAL_OUT" 2>/dev/null | wc -l || true)
    echo "Concatenation complete. Verified $FINAL_COUNT / $NUM_VARS records in final output: $FINAL_OUT"
    
    # Cleanup chunk scratch directory
    rm -rf "$CHUNK_SCRATCH"
    CMD_EXIT_CODE=0

else
    echo "Running standard sequential SpliceAI inference..."
    if [ -n "$TIME_CMD" ]; then
        $TIME_CMD $PYTHON_BIN $SPLICEAI_SCRIPT "$INPUT" "$RAW_OUT" "$FASTA" -d 10000
        CMD_EXIT_CODE=$?
    else
        $PYTHON_BIN $SPLICEAI_SCRIPT "$INPUT" "$RAW_OUT" "$FASTA" -d 10000
        CMD_EXIT_CODE=$?
    fi
    
    echo "Compressing and indexing output VCF..."
    if [ -f "$RAW_OUT" ]; then
        $BGZIP_BIN -f "$RAW_OUT"
    fi
    $TABIX_BIN -f "$FINAL_OUT"
fi

echo "============================================================================"
echo " SpliceAI annotation completed successfully. Final output: $FINAL_OUT"
echo "============================================================================"

# ----------------- Logger End -----------------
if [ -f "$LOGGER_SCRIPT" ] && [ -n "${LOG_DIR:-}" ]; then
    log_step_end "spliceai" "$FINAL_OUT" "$CMD_EXIT_CODE" "$TIME_LOG"
fi
# ----------------------------------------------

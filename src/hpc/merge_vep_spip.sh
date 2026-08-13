#!/bin/bash
#$ -P PGP
#$ -N merge
#$ -A PGP
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -o _log/master/merge.stdout
#$ -e _log/master/merge.stderr

# This is the last script, run after VEP and SPIP files are already present

SPIP=$1
VEP=$2
OUTPUT=$3
LARGE_VCF=${4:-""}
PANGOLIN=${5:-""}
SPLICEAI=${6:-""}
BRANCHPOINT=${7:-""}

BCFTOOLS=/data_lab_PGP/shared/utils/conda_envs/genomics/bin/bcftools
if [ ! -x "$BCFTOOLS" ]; then
    BCFTOOLS=/home/mruizp/apps/miniforge3/envs/genomics/bin/bcftools
fi
if [ ! -x "$BCFTOOLS" ]; then
    BCFTOOLS=/home/mruizp/conda_envs/spliceai_env/bin/bcftools
fi
if [ ! -x "$BCFTOOLS" ]; then
    BCFTOOLS=bcftools
fi

TABIX=/data_lab_PGP/shared/utils/conda_envs/genomics/bin/tabix
if [ ! -x "$TABIX" ]; then
    TABIX=/home/mruizp/apps/miniforge3/envs/genomics/bin/tabix
fi
if [ ! -x "$TABIX" ]; then
    TABIX=/home/mruizp/conda_envs/spliceai_env/bin/tabix
fi
if [ ! -x "$TABIX" ]; then
    TABIX=tabix
fi

BGZIP=/data_lab_PGP/shared/utils/conda_envs/genomics/bin/bgzip
if [ ! -x "$BGZIP" ]; then
    BGZIP=/home/mruizp/apps/miniforge3/envs/genomics/bin/bgzip
fi
if [ ! -x "$BGZIP" ]; then
    BGZIP=/home/mruizp/conda_envs/spliceai_env/bin/bgzip
fi
if [ ! -x "$BGZIP" ]; then
    BGZIP=bgzip
fi

# ----------------- Logger Setup -----------------
LOGGER_SCRIPT="src/hpc/logger.sh"
if [ -f "$LOGGER_SCRIPT" ] && [ -n "$LOG_DIR" ]; then
    source "$LOGGER_SCRIPT"
    log_step_start "merge_ann" "$VEP"
    TIME_LOG="$LOG_DIR/merge_ann.time"
    TIME_CMD=$(get_time_cmd "$TIME_LOG")
else
    TIME_CMD=""
fi
if ! command -v run_command_timed &>/dev/null; then
    run_command_timed() { "$@"; }
fi
# ------------------------------------------------

CURRENT_VCF="$VEP"
if [ ! -f "${CURRENT_VCF}.tbi" ]; then
    echo "Indexing $CURRENT_VCF..."
    $TABIX "$CURRENT_VCF" 2>/dev/null || true
fi

TEMP_MERGE_FILES=()

merge_predictor() {
    local pred_file="$1"
    local pred_cols="$2"
    local pred_name="$3"
    
    if [ -z "$pred_file" ] || [ "$pred_file" = "NONE" ]; then
        return 0
    fi
    
    # Check if uncompressed version exists and compress if needed
    if [ ! -f "$pred_file" ] && [ -f "${pred_file%.gz}" ]; then
        $BGZIP -f "${pred_file%.gz}"
    fi
    
    if [ ! -s "$pred_file" ]; then
        echo "Warning: $pred_name file not found or empty ($pred_file). Skipping."
        return 0
    fi
    
    if [ ! -f "${pred_file}.tbi" ]; then
        $TABIX "$pred_file" 2>/dev/null || true
    fi
    
    echo "Merging $pred_name annotations from $pred_file..."
    local tmp_out="${OUTPUT}.tmp_${pred_name}.vcf.gz"
    if $BCFTOOLS annotate -a "$pred_file" -c "$pred_cols" "$CURRENT_VCF" -Oz -o "$tmp_out" && [ -s "$tmp_out" ]; then
        $TABIX -f "$tmp_out" 2>/dev/null || true
        CURRENT_VCF="$tmp_out"
        TEMP_MERGE_FILES+=("$tmp_out" "${tmp_out}.tbi")
        echo "Successfully merged $pred_name."
    else
        echo "Warning: bcftools annotate failed for $pred_name. Keeping current VCF intact."
        rm -f "$tmp_out" "${tmp_out}.tbi"
    fi
}

# 1. Merge SPiP if present
merge_predictor "$SPIP" "SPiP" "spip"

# 2. Merge Pangolin if present
merge_predictor "$PANGOLIN" "Pangolin" "pangolin"

# 3. Merge SpliceAI if present
merge_predictor "$SPLICEAI" "SpliceAI" "spliceai"

# 4. Merge Branchpointer / LaBranchoR if present
merge_predictor "$BRANCHPOINT" "Branchpointer_prob,Branchpointer_U2_energy,Branchpoint_disrupted,LaBranchoR_score,LaBranchoR_acc_dist" "branchpoint"

# Move final merged VCF to target OUTPUT
cp "$CURRENT_VCF" "$OUTPUT"
$TABIX -f "$OUTPUT"

# Concatenate large structural variants if present
if [ -n "$LARGE_VCF" ] && [ -f "$LARGE_VCF" ]; then
    NUM_RECORDS=$($BCFTOOLS view -H "$LARGE_VCF" | head -n 1 | wc -l)
    if [ "$NUM_RECORDS" -gt 0 ]; then
        echo "Concatenating $NUM_RECORDS large structural variants back into the final annotated VCF..."
        $BCFTOOLS concat -a "$OUTPUT" "$LARGE_VCF" -Oz -o "${OUTPUT}.tmp.gz"
        $BCFTOOLS sort -Oz -o "${OUTPUT}.sorted.gz" "${OUTPUT}.tmp.gz"
        mv "${OUTPUT}.sorted.gz" "$OUTPUT"
        rm -f "${OUTPUT}.tmp.gz"
        $TABIX -f "$OUTPUT"
    fi
fi

# Cleanup temporary intermediate merge files
for tmp_f in "${TEMP_MERGE_FILES[@]}"; do
    if [ -f "$tmp_f" ] && [ "$tmp_f" != "$OUTPUT" ]; then
        rm -f "$tmp_f"
    fi
done

CMD_EXIT_CODE=$?
echo "Multi-predictor merge completed successfully. Final output: $OUTPUT"

# ----------------- Logger End -----------------
if [ -f "$LOGGER_SCRIPT" ] && [ -n "$LOG_DIR" ]; then
    log_step_end "merge_ann" "$OUTPUT" "$CMD_EXIT_CODE" "$TIME_LOG"
fi
# ----------------------------------------------
#!/bin/bash
#$ -P PGP
#$ -N SPIP
#$ -A PGP
#$ -l thread=1
#$ -l h_vmem=60G
#$ -o _log/master/annotation_spip.stdout
#$ -e _log/master/annotation_spip.stderr


INPUT=$1
OUTPUT=$2

# Do not source ~/.bashrc in headless SGE cluster jobs
# decompress if input is compressed
if [[ $INPUT == *.gz ]]; then
    # Use standard system gunzip (always available in PATH) with force overwrite
    gunzip -k -f $INPUT
    INPUT=${INPUT%.gz}
fi
# [AUTO-UPDATE] Path updated to shared/utils
export PATH="/data_lab_PGP/shared/utils/conda_envs/spip_env/bin:$PATH"

# Prevent CPU thread contention across worker loops
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Determine CPU thread count dynamically (default to available cores, capped at 12 for optimal parallel efficiency)
NPROCS=$(nproc 2>/dev/null || echo 4)
NUM_THREADS=${THREADS:-$NPROCS}
if [ -z "$THREADS" ] && [ "$NUM_THREADS" -gt 12 ]; then
    NUM_THREADS=12
fi
echo "Running SPiP v2.1 using $NUM_THREADS threads..."

# ----------------- Logger Setup -----------------
LOGGER_SCRIPT="src/hpc/logger.sh"
if [ -f "$LOGGER_SCRIPT" ] && [ -n "$LOG_DIR" ]; then
    source "$LOGGER_SCRIPT"
    log_step_start "spip" "$INPUT"
    TIME_LOG="$LOG_DIR/spip.time"
    TIME_CMD=$(get_time_cmd "$TIME_LOG")
else
    TIME_CMD=""
fi
# ------------------------------------------------

if [ -n "$TIME_CMD" ]; then
    $TIME_CMD /data_lab_PGP/shared/utils/conda_envs/spip_env/bin/Rscript /data_lab_PGP/resources/annotation/SPiP/SPiPv2.1_main.r \
     --input $INPUT \
     --output $OUTPUT \
     -g hg38 \
     -t $NUM_THREADS \
     --maxLines 22000 \
     --VCF
    CMD_EXIT_CODE=$?
else
    /data_lab_PGP/shared/utils/conda_envs/spip_env/bin/Rscript /data_lab_PGP/resources/annotation/SPiP/SPiPv2.1_main.r \
     --input $INPUT \
     --output $OUTPUT \
     -g hg38 \
     -t $NUM_THREADS \
     --maxLines 22000 \
     --VCF
    CMD_EXIT_CODE=$?
fi
echo "SPiP annotation completed."

# remove bad format line
TMP_FILE=${OUTPUT%.vcf}.tmp.vcf
grep -v "^##SPiP output v2.1$" $OUTPUT > $TMP_FILE
mv $TMP_FILE $OUTPUT

echo "Compressing the output file..."
# compress the output file
# [AUTO-UPDATE] Path updated to shared/utils
/data_lab_PGP/shared/utils/conda_envs/genomics/bin/bgzip -f $OUTPUT
# [AUTO-UPDATE] Path updated to shared/utils
/data_lab_PGP/shared/utils/conda_envs/genomics/bin/tabix ${OUTPUT}.gz

echo "SPiP annotation and compression completed. Final output: ${OUTPUT}.gz"

# ----------------- Logger End -----------------
if [ -f "$LOGGER_SCRIPT" ] && [ -n "$LOG_DIR" ]; then
    log_step_end "spip" "${OUTPUT}.gz" "$CMD_EXIT_CODE" "$TIME_LOG"
fi
# ----------------------------------------------
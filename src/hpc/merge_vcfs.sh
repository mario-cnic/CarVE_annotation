#!/bin/bash
#$ -P PGP
#$ -N merge_vcfs
#$ -A PGP
#$ -pe smp 1
#$ -l h_vmem=40G
#$ -o _log/merge_vcfs.stdout
#$ -e _log/merge_vcfs.stderr

INPUT_VCF=$1
GNOMAD_VCF=$2
OUTPUT_VCF=$3

# [AUTO-UPDATE] Path updated to shared/utils
MERGE_PYHON_SCRIPT=/data_lab_PGP/shared/utils/src/merge_vcfs.py
source ~/.bashrc

# Limit threads to prevent resource/process limits exhaustion on concurrent cluster slots
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export ARROW_IO_THREADS=1
# ----------------- Logger Setup -----------------
LOGGER_SCRIPT="src/hpc/logger.sh"
if [ -f "$LOGGER_SCRIPT" ] && [ -n "$LOG_DIR" ]; then
    source "$LOGGER_SCRIPT"
    log_step_start "merge_gnomAD" "$INPUT_VCF"
    TIME_LOG="$LOG_DIR/merge_gnomAD.time"
    TIME_CMD=$(get_time_cmd "$TIME_LOG")
else
    TIME_CMD=""
fi
# ------------------------------------------------

echo "Merging $INPUT_VCF and $GNOMAD_VCF into $OUTPUT_VCF"
# [AUTO-UPDATE] Path updated to shared/utils
if [ -n "$TIME_CMD" ]; then
    $TIME_CMD /data_lab_PGP/shared/utils/conda_envs/vcf_parser/bin/python $MERGE_PYHON_SCRIPT \
        --vcf1 "$INPUT_VCF" \
        --vcf2 "$GNOMAD_VCF" \
        --output "$OUTPUT_VCF" \
        --loglevel INFO \
        --overwrite
    CMD_EXIT_CODE=$?
else
    /data_lab_PGP/shared/utils/conda_envs/vcf_parser/bin/python $MERGE_PYHON_SCRIPT \
        --vcf1 "$INPUT_VCF" \
        --vcf2 "$GNOMAD_VCF" \
        --output "$OUTPUT_VCF" \
        --loglevel INFO \
        --overwrite
    CMD_EXIT_CODE=$?
fi

# create a temporary sorted VCF file to ensure the final output is sorted
if [ $CMD_EXIT_CODE -eq 0 ]; then
    echo "Sorting and indexing the merged VCF file"
    # [AUTO-UPDATE] Path updated to shared/utils
    /data_lab_PGP/shared/utils/conda_envs/genomics/bin/bcftools sort -Oz -o "$OUTPUT_VCF.tmp" "$OUTPUT_VCF"
    mv "$OUTPUT_VCF.tmp" "$OUTPUT_VCF"
    /data_lab_PGP/shared/utils/conda_envs/genomics/bin/bcftools index -t "$OUTPUT_VCF"
    echo "Merged VCF file created at $OUTPUT_VCF and indexed successfully"
    
    # Extract large structural variants (span >= 1000 bp)
    LARGE_VCF="${OUTPUT_VCF%.merge.gnomAD.vcf.gz}.large_svs.vcf.gz"
    echo "Filtering out large structural variants (>=1000bp) to $LARGE_VCF..."
    /data_lab_PGP/shared/utils/conda_envs/genomics/bin/bcftools filter -i 'strlen(REF) >= 1000 || strlen(ALT) >= 1000' "$OUTPUT_VCF" -Oz -o "$LARGE_VCF"
    /data_lab_PGP/shared/utils/conda_envs/genomics/bin/bcftools index -t "$LARGE_VCF"
    
    # Log any large variants found and save them to a human-readable TSV file
    NUM_LARGE=$(/data_lab_PGP/shared/utils/conda_envs/genomics/bin/bcftools view -H "$LARGE_VCF" | wc -l)
    if [ "$NUM_LARGE" -gt 0 ]; then
        LARGE_TSV="${OUTPUT_VCF%.merge.gnomAD.vcf.gz}.large_svs.tsv"
        echo -e "CHROM\tPOS\tREF_LEN\tALT_LEN\tHGVS_input" > "$LARGE_TSV"
        /data_lab_PGP/shared/utils/conda_envs/genomics/bin/bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/HGVS_input\n' "$LARGE_VCF" | \
            awk -F'\t' '{print $1"\t"$2"\t"length($3)"\t"length($4)"\t"$5}' >> "$LARGE_TSV"
        
        echo "[LOG] Found $NUM_LARGE large structural variant(s) for $(basename $OUTPUT_VCF .merge.gnomAD.vcf.gz) and saved them to $(basename $LARGE_TSV):"
        cat "$LARGE_TSV"
    fi

    # Keep only small variants in OUTPUT_VCF
    echo "Filtering output VCF to keep only small variants (<1000bp)..."
    /data_lab_PGP/shared/utils/conda_envs/genomics/bin/bcftools filter -i 'strlen(REF) < 1000 && strlen(ALT) < 1000' "$OUTPUT_VCF" -Oz -o "$OUTPUT_VCF.small"
    mv "$OUTPUT_VCF.small" "$OUTPUT_VCF"
    /data_lab_PGP/shared/utils/conda_envs/genomics/bin/bcftools index -t -f "$OUTPUT_VCF"
fi

# ----------------- Logger End -----------------
if [ -f "$LOGGER_SCRIPT" ] && [ -n "$LOG_DIR" ]; then
    log_step_end "merge_gnomAD" "$OUTPUT_VCF" "$CMD_EXIT_CODE" "$TIME_LOG"
fi
# ----------------------------------------------
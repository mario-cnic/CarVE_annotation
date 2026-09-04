#!/bin/bash
#$ -P PGP
#$ -N subset_gnomAD
#$ -A PGP
#$ -pe smp 1
#$ -l h_vmem=10G
#$ -o _log/subset_gnomAD.stdout
#$ -e _log/subset_gnomAD.stderr

GENES_LIST=${1:-resources/cardio_genes_loc.bed}
OUTPUT=${2:-annotation/gnomAD_subset}
GENE=$3
# [AUTO-UPDATE] Path updated to shared/utils
CHR_MAPPING=/data_lab_PGP/shared/utils/data/map_chr.txt
GNOMAD_VCF=/references/genomes/Homo_sapiens/annotations/gnomAD/GRCh38/v4.1/gnomAD.v4.1.vcf.gz

# Run the loop only if a GENE was not specified
if [ -z "$GENE" ]; then
    mapfile -t cardio_genes < <(cut -f1,2,3,4 "$GENES_LIST")
else
    mapfile -t cardio_genes < <(grep "$GENE" "$GENES_LIST")
fi

source ~/.bashrc

for gene_element in "${cardio_genes[@]}"; do
    IFS=$'\t' read -r chr start end gene <<< "$gene_element"
    # extract first element sep by _ in gene
    gene=${gene%%_*}
    # Add +/- 20kb buffer to capture full intronic and regulatory boundaries
    buf_start=$(( start > 20000 ? start - 20000 : 0 ))
    buf_end=$(( end + 20000 ))
    echo "Buffered Coords: $chr:$buf_start-$buf_end"

    # [AUTO-UPDATE] Path updated to shared/utils
    if [ -n "$TIME_CMD" ]; then
        $TIME_CMD /data_lab_PGP/shared/utils/conda_envs/genomics/bin/bcftools annotate -r $chr:$buf_start-$buf_end --rename-chrs $CHR_MAPPING --write-index=tbi -x INFO -Oz -o $OUTPUT/gnomAD.v4.1.${gene}.vcf.gz $GNOMAD_VCF
        CMD_EXIT_CODE=$?
    else
        /data_lab_PGP/shared/utils/conda_envs/genomics/bin/bcftools annotate -r $chr:$buf_start-$buf_end --rename-chrs $CHR_MAPPING --write-index=tbi -x INFO -Oz -o $OUTPUT/gnomAD.v4.1.${gene}.vcf.gz $GNOMAD_VCF
        CMD_EXIT_CODE=$?
    fi

    # ----------------- Logger End -----------------
    if [ -f "$LOGGER_SCRIPT" ] && [ -n "$LOG_DIR" ]; then
        log_step_end "gnomAD_subset" "$OUTPUT/gnomAD.v4.1.${gene}.vcf.gz" "$CMD_EXIT_CODE" "$TIME_LOG"
    fi
    # ----------------------------------------------
done
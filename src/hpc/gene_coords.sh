#!/bin/bash
#$ -P PGP
#$ -N gene_coords
#$ -A PGP
#$ -pe smp 2
#$ -l h_vmem=4G
#$ -o _log/gene_coords.stdout
#$ -e _log/gene_coords.stderr

# [AUTO-UPDATE] Path updated to shared/utils
GENE_COORDS=/data_lab_PGP/shared/utils/src/gene_coords.py
MARGIN_5=0
MARGIN_3=0
GENE_TRANSCRIPT_MAPPING=${1:-resources/gene_transcript_mapping.txt}
BED_FILE=${2:-resources/cardio_genes_loc.bed}
GENE_NAME=${3:-}
# read gene_transcipt_mapping first and second column and save it to an array
mapfile -t gene_transcript_mapping < <(cut -d',' -f1,2 "$GENE_TRANSCRIPT_MAPPING")
source ~/.bashrc

# if GENE_NAME is provided, filter to only that gene and check if it's not already in BED_FILE
if [[ -n "$GENE_NAME" ]]; then
    if grep -q "^$GENE_NAME" "$BED_FILE" 2>/dev/null; then
        echo "Gene $GENE_NAME already in $BED_FILE"
        exit 0
    fi
    gene_transcript_mapping=("$(grep "^$GENE_NAME," <(printf '%s\n' "${gene_transcript_mapping[@]}"))")
fi

# loop through the values 
for gene_transcript in "${gene_transcript_mapping[@]}"; do
    [[ -z "$gene_transcript" ]] && continue
    IFS=',' read -r gene transcript <<< "$gene_transcript"
    # drop version number from transcript
    transcript=$(echo $transcript | cut -d'.' -f1)
    echo "Gene: $gene, Transcript: $transcript"
    # check if the gene is already in the BED_FILE, if it is skip to the next one
    if grep -q "^$gene" "$BED_FILE" 2>/dev/null; then
        echo "Gene $gene already in $BED_FILE, skipping..."
        continue
    fi
    # [AUTO-UPDATE] Path updated to shared/utils
    mamba run -p /data_lab_PGP/shared/utils/conda_envs/datasci python $GENE_COORDS \
        --gene $gene \
        --refseq $transcript \
        --margin5 $MARGIN_5 \
        --margin3 $MARGIN_3 \
        --mode a \
        --output $BED_FILE
        # --no-chr \
    sleep 3
done
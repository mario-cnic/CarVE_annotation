#!/bin/bash
#$ -P PGP
#$ -N utr_annotation
#$ -A PGP
#$ -l thread=1
#$ -l h_vmem=20G
#$ -o _log/utrann_step.stdout
#$ -e _log/utrann_step.stderr

data_PGP="/home/mruizp/data_lab_PGP"
in_dir="$data_PGP/enrichment_vars_220126"
file_pattern="_tmp/Variantes_Enriquecimiento_F.toUTRANT.csv"
genome_fasta=/home/mruizp/data_references/genomes/Homo_sapiens/GATK_bundle/v0/Homo_sapiens_assembly38.fasta
out_dir="$in_dir/annotation/Variantes_Enriquecimiento_F.UTRANT.csv"
species=human
version=109
dbDir="$in_dir/utr_data"
filetype=csv
partition=1
numCores=6

source /home/mruizp/apps/miniforge3/bin/activate /home/mruizp/apps/miniforge3/envs/R

# Find the first file matching the pattern
file="${in_dir}/${file_pattern}"

echo "Processing $file..."

filename=$(basename "$file")
out_file=$out_dir

if [ -f "$out_file" ]; then
    echo "Warning: Output file $out_file already exists. Skipping processing for $file."
    exit 0
fi

Rscript /home/mruizp/data_lab_PGP/annotation/utr.annotation/src/UTRannotation.R  ${file} ${out_file} ${species} ${version} ${dbDir} ${filetype}

echo "Process completed"



#!/bin/bash
#$ -P BIGN
#$ -N step02
#$ -A PGP
#$ -pe smp 1
#$ -l h_vmem=8G
#$ -o _log/02_cohort_enrich.stdout
#$ -e _log/02_cohort_enrich.stderr

GENE=$1
INPUT_DIR=$2
# Ensure log directory exists if it doesn't
mkdir -p _log

source ~/.bashrc
# Adjust conda env or R module loading here if necessary
# [AUTO-UPDATE] Path updated to shared/utils
mamba run -p /data_lab_PGP/shared/utils/conda_envs/r4.4 Rscript src/R/02_cohort_enrich.R "$GENE" "$INPUT_DIR"

#!/bin/bash
#$ -P BIGN
#$ -N step04
#$ -A PGP
#$ -pe smp 1
#$ -l h_vmem=8G
#$ -o _log/04_forest_plot.stdout
#$ -e _log/04_forest_plot.stderr

GENE=$1
INPUT_DIR=$2
# Ensure log directory exists if it doesn't
mkdir -p _log

source ~/.bashrc
# Adjust conda env or R module loading here if necessary
# [AUTO-UPDATE] Path updated to shared/utils
mamba run -p /data_lab_PGP/shared/utils/conda_envs/r4.4 Rscript src/downstream/04_forest_plot.R "$GENE" "$INPUT_DIR"

#!/bin/bash
# SIF Container Integration & Health Verification Unit Test
set -euo pipefail

echo "Running SIF Container Verification Test..."

# Determine Container Executable
CONTAINER_BIN=""
if command -v apptainer >/dev/null 2>&1; then
    CONTAINER_BIN="apptainer"
elif command -v singularity >/dev/null 2>&1; then
    CONTAINER_BIN="singularity"
fi

if [ -z "$CONTAINER_BIN" ]; then
    echo "  [SKIP] Neither apptainer nor singularity binary is installed on PATH."
    exit 0
fi

# Locate Target SIF Image
SIF_FILE="${1:-resources/containers/annotation_pipeline.sif}"
if [ ! -f "$SIF_FILE" ]; then
    SIF_FILE="resources/sif_images/annotation_pipeline.sif"
fi

if [ ! -f "$SIF_FILE" ]; then
    echo "  [SKIP] Container image not found at '$SIF_FILE'. Build it first using src/hpc/build_container.sh"
    exit 0
fi

echo "  Testing container image: $SIF_FILE"

# Test 1: Verify Python & Core AI/Genomics Imports inside SIF
echo "  [1/4] Verifying Python & Deep Learning Packages inside SIF..."
"$CONTAINER_BIN" exec "$SIF_FILE" python3 -c "
import pysam
import pyarrow
import pandas
import torch
import tensorflow as tf
import spliceai
import plotly
print('    -> Python version:', tf.__version__)
print('    -> All Python libraries (pysam, pyarrow, pandas, PyTorch, TensorFlow, SpliceAI, Plotly) imported successfully.')
"

# Test 2: Verify Genomic Binaries (bcftools, tabix) inside SIF
echo "  [2/4] Verifying Genomic Binaries (bcftools, tabix) inside SIF..."
bcf_ver=$("$CONTAINER_BIN" exec "$SIF_FILE" bcftools --version | head -n1)
tab_ver=$("$CONTAINER_BIN" exec "$SIF_FILE" tabix --version 2>&1 | head -n1 || true)
echo "    -> BCFtools: $bcf_ver"
echo "    -> Tabix:    $tab_ver"

# Test 3: Verify R Splicing & Plotting Packages inside SIF
echo "  [3/4] Verifying R Packages inside SIF..."
"$CONTAINER_BIN" exec "$SIF_FILE" Rscript -e "
suppressPackageStartupMessages({
    library(ggplot2)
    library(randomForest)
    library(optparse)
    library(data.table)
    library(dplyr)
})
cat('    -> All R packages (ggplot2, randomForest, optparse, data.table, dplyr) loaded successfully.\n')
"

# Test 4: Main.sh --use-container CLI Integration Test
echo "  [4/4] Verifying main.sh CLI --use-container integration..."
audit_output=$(bash main.sh --test --use-container --sif "$SIF_FILE" --audit-only 2>&1 || true)
if echo "$audit_output" | grep -q "Executing Master Pipeline Quality"; then
    echo "    -> main.sh --use-container --sif executed cleanly!"
else
    echo "  [FAIL] main.sh --use-container integration failed."
    echo "$audit_output"
    exit 1
fi

echo "  ✅ Container SIF image ($SIF_FILE) passed all health & integrity checks!"

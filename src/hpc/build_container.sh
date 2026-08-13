#!/bin/bash
#$ -P PGP
#$ -N build_sif
#$ -A PGP
#$ -pe smp 4
#$ -l h_vmem=16G
#$ -o _log/build_sif.stdout
#$ -e _log/build_sif.stderr

# Script to build the unified Annotation Pipeline Apptainer/Singularity container (.sif)
set -euo pipefail

DEF_FILE="resources/containers/annotation_pipeline.def"
OUTPUT_SIF="/data_lab_PGP/resources/sif_images/annotation_pipeline.sif"

if [ ! -f "$DEF_FILE" ]; then
    echo "Error: Definition file $DEF_FILE not found."
    exit 1
fi

mkdir -p "$(dirname "$OUTPUT_SIF")"

echo "Building Apptainer SIF container image..."
echo "Definition file: $DEF_FILE"
echo "Output image   : $OUTPUT_SIF"

# Build using apptainer (or singularity if apptainer is alias)
if command -v apptainer &>/dev/null; then
    apptainer build --fakeroot "$OUTPUT_SIF" "$DEF_FILE"
elif command -v singularity &>/dev/null; then
    singularity build --fakeroot "$OUTPUT_SIF" "$DEF_FILE"
else
    echo "Error: Neither apptainer nor singularity found in PATH."
    exit 1
fi

echo "Successfully built $OUTPUT_SIF!"

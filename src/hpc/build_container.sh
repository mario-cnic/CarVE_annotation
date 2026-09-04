#!/bin/bash
#$ -P PGP
#$ -N build_sif
#$ -A PGP
#$ -pe smp 4
#$ -l h_vmem=16G
#$ -o _log/build_sif.stdout
#$ -e _log/build_sif.stderr

# Master Script to build the unified Annotation Pipeline Apptainer/Singularity container (.sif)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

if [ -f "$PIPELINE_ROOT/config/env.sh" ]; then
    source "$PIPELINE_ROOT/config/env.sh"
fi

DEF_FILE="$PIPELINE_ROOT/resources/containers/annotation_pipeline.def"
OUTPUT_SIF="${ANNOTATION_SIF:-${DATA_LAB_PGP}/resources/sif_images/annotation_pipeline.sif}"

if [ ! -f "$DEF_FILE" ]; then
    echo "Error: Definition file $DEF_FILE not found."
    exit 1
fi

mkdir -p "$(dirname "$OUTPUT_SIF")"

echo "============================================================================"
echo " Building Apptainer SIF Container Image"
echo "============================================================================"
echo " Definition file: $DEF_FILE"
echo " Target image   : $OUTPUT_SIF"
echo "============================================================================"

if command -v apptainer &>/dev/null; then
    apptainer build --fakeroot "$OUTPUT_SIF" "$DEF_FILE"
elif command -v singularity &>/dev/null; then
    singularity build --fakeroot "$OUTPUT_SIF" "$DEF_FILE"
else
    echo "Error: Neither apptainer nor singularity binary found in PATH."
    exit 1
fi

echo "============================================================================"
echo " Successfully built container image: $OUTPUT_SIF"
echo "============================================================================"

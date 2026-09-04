#!/usr/bin/env bash
# Universal Genomic Variant Annotation Pipeline
# Centralized Configuration & Environment Resolver
# Auto-detects execution environment: Local Workstation Sandbox vs SGE HPC Cluster

# ----------------- Base System & Volume Paths -----------------
if [ -d "/data_lab_PGP" ]; then
    export DATA_LAB_PGP="${DATA_LAB_PGP:-/data_lab_PGP}"
else
    export DATA_LAB_PGP="${DATA_LAB_PGP:-/home/mruizp/data_lab_PGP}"
fi

if [ -d "/references" ]; then
    export REFERENCES_DIR="${REFERENCES_DIR:-/references}"
else
    export REFERENCES_DIR="${REFERENCES_DIR:-/home/mruizp/data_references}"
fi

export SHARED_UTILS_DIR="${SHARED_UTILS_DIR:-${DATA_LAB_PGP}/shared/utils}"
export SHARED_SRC_DIR="${SHARED_SRC_DIR:-${SHARED_UTILS_DIR}/src}"
export SHARED_DATA_DIR="${SHARED_DATA_DIR:-${SHARED_UTILS_DIR}/data}"

# ----------------- Conda Environments & Binary Registries -----------------
if [ -d "${SHARED_UTILS_DIR}/conda_envs" ]; then
    export CONDA_ENVS_DIR="${CONDA_ENVS_DIR:-${SHARED_UTILS_DIR}/conda_envs}"
elif [ -d "/home/mruizp/apps/miniforge3/envs" ]; then
    export CONDA_ENVS_DIR="${CONDA_ENVS_DIR:-/home/mruizp/apps/miniforge3/envs}"
else
    export CONDA_ENVS_DIR="${CONDA_ENVS_DIR:-/data_lab_PGP/shared/utils/conda_envs}"
fi

# Python Executables
if [ -x "/home/mruizp/apps/miniforge3/envs/datasci/bin/python3" ]; then
    export PYTHON_DATASCI="${PYTHON_DATASCI:-/home/mruizp/apps/miniforge3/envs/datasci/bin/python3}"
elif [ -x "${CONDA_ENVS_DIR}/datasci/bin/python3" ]; then
    export PYTHON_DATASCI="${PYTHON_DATASCI:-${CONDA_ENVS_DIR}/datasci/bin/python3}"
else
    export PYTHON_DATASCI="${PYTHON_DATASCI:-python3}"
fi

if [ -x "${CONDA_ENVS_DIR}/vcf_parser/bin/python" ]; then
    export PYTHON_VCF_PARSER="${PYTHON_VCF_PARSER:-${CONDA_ENVS_DIR}/vcf_parser/bin/python}"
elif [ -x "/home/mruizp/apps/miniforge3/envs/datasci/bin/python3" ]; then
    export PYTHON_VCF_PARSER="${PYTHON_VCF_PARSER:-/home/mruizp/apps/miniforge3/envs/datasci/bin/python3}"
else
    export PYTHON_VCF_PARSER="${PYTHON_VCF_PARSER:-python3}"
fi

if [ -x "${CONDA_ENVS_DIR}/pangolin_env/bin/python" ]; then
    export PYTHON_PANGOLIN="${PYTHON_PANGOLIN:-${CONDA_ENVS_DIR}/pangolin_env/bin/python}"
else
    export PYTHON_PANGOLIN="${PYTHON_PANGOLIN:-$PYTHON_DATASCI}"
fi

if [ -x "${CONDA_ENVS_DIR}/spliceai_env/bin/python" ]; then
    export PYTHON_SPLICEAI="${PYTHON_SPLICEAI:-${CONDA_ENVS_DIR}/spliceai_env/bin/python}"
else
    export PYTHON_SPLICEAI="${PYTHON_SPLICEAI:-$PYTHON_DATASCI}"
fi

# R Executables
if [ -x "${CONDA_ENVS_DIR}/datasci/bin/Rscript" ]; then
    export RSCRIPT_DATASCI="${RSCRIPT_DATASCI:-${CONDA_ENVS_DIR}/datasci/bin/Rscript}"
else
    export RSCRIPT_DATASCI="${RSCRIPT_DATASCI:-Rscript}"
fi

if [ -x "${CONDA_ENVS_DIR}/spip_env/bin/Rscript" ]; then
    export RSCRIPT_SPIP="${RSCRIPT_SPIP:-${CONDA_ENVS_DIR}/spip_env/bin/Rscript}"
else
    export RSCRIPT_SPIP="${RSCRIPT_SPIP:-$RSCRIPT_DATASCI}"
fi

# Genomic Binaries (BCFtools, Tabix, BGZIP)
if [ -x "/home/mruizp/apps/miniforge3/envs/genomics/bin/bcftools" ]; then
    export BCFTOOLS_BIN="${BCFTOOLS_BIN:-/home/mruizp/apps/miniforge3/envs/genomics/bin/bcftools}"
elif [ -x "${CONDA_ENVS_DIR}/genomics/bin/bcftools" ]; then
    export BCFTOOLS_BIN="${BCFTOOLS_BIN:-${CONDA_ENVS_DIR}/genomics/bin/bcftools}"
else
    export BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
fi

if [ -x "/home/mruizp/apps/miniforge3/envs/genomics/bin/tabix" ]; then
    export TABIX_BIN="${TABIX_BIN:-/home/mruizp/apps/miniforge3/envs/genomics/bin/tabix}"
elif [ -x "${CONDA_ENVS_DIR}/genomics/bin/tabix" ]; then
    export TABIX_BIN="${TABIX_BIN:-${CONDA_ENVS_DIR}/genomics/bin/tabix}"
else
    export TABIX_BIN="${TABIX_BIN:-tabix}"
fi

if [ -x "/home/mruizp/apps/miniforge3/envs/genomics/bin/bgzip" ]; then
    export BGZIP_BIN="${BGZIP_BIN:-/home/mruizp/apps/miniforge3/envs/genomics/bin/bgzip}"
elif [ -x "${CONDA_ENVS_DIR}/genomics/bin/bgzip" ]; then
    export BGZIP_BIN="${BGZIP_BIN:-${CONDA_ENVS_DIR}/genomics/bin/bgzip}"
else
    export BGZIP_BIN="${BGZIP_BIN:-bgzip}"
fi

# ----------------- Genomic References & Resources -----------------
export GRCH38_FASTA="${GRCH38_FASTA:-${REFERENCES_DIR}/genomes/Homo_sapiens/GATK_bundle/v0/Homo_sapiens_assembly38.fasta}"
export LABRANCHOR_BED="${LABRANCHOR_BED:-${DATA_LAB_PGP}/resources/annotation/labranchor/labranchor_grch38_top.bed.gz}"
export VCF_PARSER_SCRIPT="${VCF_PARSER_SCRIPT:-${SHARED_SRC_DIR}/vcf_parser_pysam.py}"
export FILTER_VARIANTS_SCRIPT="${FILTER_VARIANTS_SCRIPT:-${SHARED_SRC_DIR}/filter_variants.py}"
export VEP_COLUMNS_FILE="${VEP_COLUMNS_FILE:-${SHARED_DATA_DIR}/vcf_all_columns.txt}"

# ----------------- Container Execution Engine -----------------
export USE_CONTAINER="${USE_CONTAINER:-false}"
export ANNOTATION_SIF="${ANNOTATION_SIF:-${DATA_LAB_PGP}/resources/containers/annotation_pipeline.sif}"

# Apptainer / Singularity Executable
if command -v apptainer >/dev/null 2>&1; then
    export CONTAINER_BIN="apptainer"
elif command -v singularity >/dev/null 2>&1; then
    export CONTAINER_BIN="singularity"
else
    export CONTAINER_BIN=""
fi

# Function to run commands inside container if USE_CONTAINER is true
run_container_cmd() {
    if [ "${USE_CONTAINER:-false}" = "true" ] && [ -n "${CONTAINER_BIN:-}" ] && [ -f "${ANNOTATION_SIF:-}" ]; then
        "${CONTAINER_BIN}" exec -B /data_lab_PGP:/data_lab_PGP -B /references:/references -B /home/mruizp:/home/mruizp "${ANNOTATION_SIF}" "$@"
    else
        "$@"
    fi
}

# ----------------- Default Thresholds & Constants -----------------
export DEFAULT_MAX_AF=0.001
export DEFAULT_MIN_REVEL=0.75
export DEFAULT_MIN_ALPHAMISSENSE=0.80
export DEFAULT_MIN_SPIP=0.50
export DEFAULT_MIN_CADD=20.0


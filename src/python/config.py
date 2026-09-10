#!/usr/bin/env python3
"""
Centralized Python Configuration & Environment Resolver for Annotation Pipeline
Auto-detects workstation sandbox vs HPC cluster mounts and provides fallback defaults.
"""

import os
import sys

# 1. Base Directory Auto-Resolution
DATA_LAB_PGP = os.getenv("DATA_LAB_PGP", "/data_lab_PGP" if os.path.exists("/data_lab_PGP") else "/home/mruizp/data_lab_PGP")
REFERENCES_DIR = os.getenv("REFERENCES_DIR", "/references" if os.path.exists("/references") else "/home/mruizp/data_references")

# 2. Shared Utilities Directories
SHARED_UTILS_DIR = os.getenv("SHARED_UTILS_DIR", os.path.join(DATA_LAB_PGP, "shared", "utils"))
SHARED_SRC_DIR = os.getenv("SHARED_SRC_DIR", os.path.join(SHARED_UTILS_DIR, "src"))
SHARED_DATA_DIR = os.getenv("SHARED_DATA_DIR", os.path.join(SHARED_UTILS_DIR, "data"))

# Ensure SHARED_SRC_DIR is on sys.path for python modules
if os.path.exists(SHARED_SRC_DIR) and SHARED_SRC_DIR not in sys.path:
    sys.path.insert(0, SHARED_SRC_DIR)

# 3. Executable & Binary Paths
CONDA_ENVS_DIR = os.getenv("CONDA_ENVS_DIR", os.path.join(SHARED_UTILS_DIR, "conda_envs"))

PYTHON_DATASCI = os.getenv("PYTHON_DATASCI", sys.executable)
BCFTOOLS_BIN = os.getenv("BCFTOOLS_BIN", "/home/mruizp/apps/miniforge3/envs/genomics/bin/bcftools" if os.path.exists("/home/mruizp/apps/miniforge3/envs/genomics/bin/bcftools") else "bcftools")
TABIX_BIN = os.getenv("TABIX_BIN", "/home/mruizp/apps/miniforge3/envs/genomics/bin/tabix" if os.path.exists("/home/mruizp/apps/miniforge3/envs/genomics/bin/tabix") else "tabix")

# 4. Genomic References
GRCH38_FASTA = os.getenv("GRCH38_FASTA", os.path.join(REFERENCES_DIR, "genomes", "Homo_sapiens", "GATK_bundle", "v0", "Homo_sapiens_assembly38.fasta"))
LABRANCHOR_BED = os.getenv("LABRANCHOR_BED", os.path.join(DATA_LAB_PGP, "resources", "annotation", "labranchor", "labranchor_grch38_top.bed.gz"))

# 5. Default Threshold Constants
DEFAULT_MAX_AF = float(os.getenv("DEFAULT_MAX_AF", 0.001))
DEFAULT_MIN_REVEL = float(os.getenv("DEFAULT_MIN_REVEL", 0.75))
DEFAULT_MIN_ALPHAMISSENSE = float(os.getenv("DEFAULT_MIN_ALPHAMISSENSE", 0.80))
DEFAULT_MIN_SPIP = float(os.getenv("DEFAULT_MIN_SPIP", 0.50))
DEFAULT_MIN_CADD = float(os.getenv("DEFAULT_MIN_CADD", 20.0))

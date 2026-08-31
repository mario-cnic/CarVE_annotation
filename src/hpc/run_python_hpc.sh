#!/bin/bash
# HPC Python Execution Wrapper
# Ensures clean execution of Python scripts on SGE compute nodes
# preventing NFS mmap shared object memory mapping errors (failed to map segment).

set -euo pipefail

export TMPDIR=/tmp
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

PYTHON_EXE=$1
shift

if [ ! -x "$PYTHON_EXE" ]; then
    # Try local path fallback if HPC path differs
    ALT_PYTHON="/home/mruizp/data_lab_PGP/shared/utils/conda_envs/datasci/bin/python3"
    if [ -x "$ALT_PYTHON" ]; then
        PYTHON_EXE="$ALT_PYTHON"
    fi
fi

exec "$PYTHON_EXE" "$@"

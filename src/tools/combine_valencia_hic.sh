#!/usr/bin/env bash
#$ -P BIGN
#$ -N combine_valencia_hic
#$ -A PGP
#$ -pe smp 1
#$ -l h_vmem=30G
#$ -o _log/combine_valencia_hic.stdout
#$ -e _log/combine_valencia_hic.stderr

# Exit on error
set -e

# Source shell configuration
source ~/.bashrc

echo "=========================================================="
echo "Starting combine_valencia_hic execution on SGE HPC cluster"
echo "Date: $(date)"
echo "=========================================================="

# Inspect arguments without destroying "$@"
ORACLE_DIR="_RAW/data_Valencia_plus_hic/oracle"
NDG_DIR="_RAW/data_Valencia_plus_hic/ndg"
COMBINED_DIR=""
TEST_FLAG=""
ONLY_LIFTOVER_FLAG=""
OVERWRITE_ALL_FLAG=""

for ((i=1; i<=$#; i++)); do
    arg="${!i}"
    next_idx=$((i+1))
    next_arg="${!next_idx}"
    if [ "$arg" = "--test" ]; then
        TEST_FLAG="--test"
    elif [ "$arg" = "--oracle-dir" ]; then
        ORACLE_DIR="$next_arg"
    elif [ "$arg" = "--ndg-dir" ]; then
        NDG_DIR="$next_arg"
    elif [ "$arg" = "--out-dir" ]; then
        COMBINED_DIR="$next_arg"
    elif [ "$arg" = "--only-liftover" ]; then
        ONLY_LIFTOVER_FLAG="--only-liftover"
    elif [ "$arg" = "--overwrite-all" ]; then
        OVERWRITE_ALL_FLAG="--overwrite-all"
    fi
done

if [ -n "$ONLY_LIFTOVER_FLAG" ]; then
    echo "Mode: Liftover-Only (bypassing external GeneBe APIs)"
fi
if [ -n "$OVERWRITE_ALL_FLAG" ]; then
    echo "Mode: Overwrite all (re-combining all files)"
else
    echo "Mode: Skipping existing combined files (only merging missing ones by default)"
fi


# Step 1: Pre-merging Audit
echo "----------------------------------------------------------"
echo "Step 1: Running Pre-merging Audit..."
echo "----------------------------------------------------------"
/data_lab_PGP/shared/utils/conda_envs/liftover/bin/python3 src/exploratory/explore_valencia_hic.py \
    --oracle-dir "$ORACLE_DIR" \
    --ndg-dir "$NDG_DIR"

# Step 2: Combination/Integration and Liftover
echo "----------------------------------------------------------"
echo "Step 2: Running Valencia + HiC Combination/Integration..."
echo "----------------------------------------------------------"
/data_lab_PGP/shared/utils/conda_envs/liftover/bin/python3 -u src/python/combine_valencia_hic.py "$@"

# Step 3: Post-merging Audit
echo "----------------------------------------------------------"
echo "Step 3: Running Post-merging Audit..."
echo "----------------------------------------------------------"
AUDIT_ARGS=""
if [ -n "$COMBINED_DIR" ]; then
    AUDIT_ARGS="--combined-dir $COMBINED_DIR"
elif [ -n "$TEST_FLAG" ]; then
    AUDIT_ARGS="--test"
fi

/data_lab_PGP/shared/utils/conda_envs/liftover/bin/python3 src/exploratory/explore_combined_valencia_hic.py $AUDIT_ARGS

echo "=========================================================="
echo "Finished all orchestrated pipeline steps successfully"
echo "Date: $(date)"
echo "=========================================================="

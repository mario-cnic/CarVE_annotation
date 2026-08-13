#!/bin/bash
# Master Unit Test Suite Executor for Genomic Annotation Pipeline

set -euo pipefail

echo "============================================================================"
echo " 🧪 Master Pipeline Unit Testing Suite Runner "
echo "============================================================================"
echo "Timestamp: $(date '+%Y-%m-%d %H:%M:%S')"
echo "Host     : $(hostname)"
echo "Target   : GRCh38 Pipeline"
echo "============================================================================"
echo ""

FAILURES=0

# 1. Execute Python Unit Test Suite
echo "----------------------------------------------------------------------------"
echo " [1/3] Running Python Unit Tests (pytest)..."
echo "----------------------------------------------------------------------------"
if python3 -m pytest -v tests/python/; then
    echo "  ✅ Python unit test suite PASSED."
else
    echo "  ❌ Python unit test suite FAILED."
    FAILURES=$((FAILURES + 1))
fi
echo ""

# 2. Execute R Unit Test Suite
echo "----------------------------------------------------------------------------"
echo " [2/3] Running R Unit Tests (Rscript)..."
echo "----------------------------------------------------------------------------"
if Rscript tests/R/test_statsJPO.R && Rscript tests/R/test_plot_annotation_results.R; then
    echo "  ✅ R unit test suite PASSED."
else
    echo "  ❌ R unit test suite FAILED."
    FAILURES=$((FAILURES + 1))
fi
echo ""

# 3. Execute Bash Unit Test Suite
echo "----------------------------------------------------------------------------"
echo " [3/3] Running Bash Unit Tests (sh)..."
echo "----------------------------------------------------------------------------"
if bash tests/bash/test_logger.sh && bash tests/bash/test_cli_args_main.sh; then
    echo "  ✅ Bash unit test suite PASSED."
else
    echo "  ❌ Bash unit test suite FAILED."
    FAILURES=$((FAILURES + 1))
fi
echo ""

echo "============================================================================"
if [ "$FAILURES" -eq 0 ]; then
    echo " 🎉 ALL UNIT TEST SUITES PASSED CLEANLY! (100% SUCCESS)"
    echo "============================================================================"
    exit 0
else
    echo " ⚠️ $FAILURES TEST SUITE(S) FAILED. CHECK LOGS ABOVE."
    echo "============================================================================"
    exit 1
fi

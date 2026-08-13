#!/bin/bash
set -euo pipefail

echo "Running Bash unit test for logger.sh..."

TEST_DIR=$(mktemp -d)
export LOG_DIR="$TEST_DIR/_log"

source src/hpc/logger.sh

# Test log_init
RAW_DUMMY="$TEST_DIR/dummy_raw.csv"
echo "A,B,C" > "$RAW_DUMMY"
log_init "MYBPC3" "$RAW_DUMMY"

if [ ! -f "$REPORT_FILE" ]; then
    echo "  [FAIL] Report file was not created by log_init."
    exit 1
fi

# Test log_step_start
log_step_start "Variant_Converter" "$RAW_DUMMY"

# Test log_step_end
OUT_DUMMY="$TEST_DIR/out.vcf"
echo "##fileformat=VCFv4.2" > "$OUT_DUMMY"
log_step_end "Variant_Converter" "$OUT_DUMMY" 0 ""

# Test log_final
log_final

if grep -q "MYBPC3" "$REPORT_FILE" && grep -q "Variant_Converter" "$REPORT_FILE"; then
    echo "  [PASS] logger.sh unit test passed cleanly!"
else
    echo "  [FAIL] Report content verification failed."
    exit 1
fi

rm -rf "$TEST_DIR"

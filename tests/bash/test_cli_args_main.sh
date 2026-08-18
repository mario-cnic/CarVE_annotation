#!/bin/bash
set -euo pipefail

echo "Running Bash unit test for main.sh CLI options & modular checkpoint flags..."

# Test 1: main.sh without --raw-dir should print usage and exit with 1
output=$(bash main.sh 2>&1 || true)
if echo "$output" | grep -q "Usage: main.sh --raw-dir"; then
    echo "  [PASS] main.sh correctly rejected missing --raw-dir parameter."
else
    echo "  [FAIL] main.sh did not display expected usage notice."
    exit 1
fi

# Test 2: main.sh with unknown parameter should exit with error
unknown_output=$(bash main.sh --unknown-flag 2>&1 || true)
if echo "$unknown_output" | grep -q "Unknown option: --unknown-flag"; then
    echo "  [PASS] main.sh correctly caught unknown CLI parameter."
else
    echo "  [FAIL] main.sh failed to reject unknown parameter."
    exit 1
fi

# Test 3: main.sh help includes checkpointing and force flags
if echo "$output" | grep -q -- "--force-spliceai" && echo "$output" | grep -q -- "--overwrite-all"; then
    echo "  [PASS] main.sh help screen correctly includes modular checkpointing and force flags."
else
    echo "  [FAIL] main.sh help screen is missing expected checkpoint flags."
    exit 1
fi

echo "All main.sh CLI parameter unit tests passed cleanly!"

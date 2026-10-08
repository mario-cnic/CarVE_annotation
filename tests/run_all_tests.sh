#!/bin/bash
# Run each test module independently to isolate their Python import paths.
set -uo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT" || exit 1
PYTHON_BIN="${PYTHON_BIN:-python3}"
failures=0
for test_file in tests/python/test_*.py; do
    echo "Testing $test_file"
    if "$PYTHON_BIN" -m pytest -q "$test_file"; then
        continue
    fi
    failures=$((failures + 1))
done
if [[ "$failures" -ne 0 ]]; then
    echo "$failures test module(s) failed."
    exit 1
fi
echo "All active Python test modules passed."

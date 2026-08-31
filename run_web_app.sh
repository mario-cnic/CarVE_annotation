#!/bin/bash
# Standalone Clinical Reporting Web Application Launcher

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

PYTHON_BIN="/home/mruizp/apps/miniforge3/envs/datasci/bin/python3"
if [ ! -x "$PYTHON_BIN" ]; then
    PYTHON_BIN="/data_lab_PGP/shared/utils/conda_envs/datasci/bin/python3"
fi
if [ ! -x "$PYTHON_BIN" ]; then
    PYTHON_BIN="python3"
fi

PORT="${1:-8501}"

echo "=========================================================================="
echo " 🧬 Starting Standalone Clinical Reporting Web Application..."
echo "=========================================================================="
echo "  - Script: src/web_app/app.py"
echo "  - Python: $PYTHON_BIN"
echo "  - Port  : $PORT"
echo "=========================================================================="

"$PYTHON_BIN" -m streamlit run src/web_app/app.py --server.port "$PORT" --server.address 0.0.0.0

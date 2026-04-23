#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

ENV_NAME="DPPP"

if ! command -v conda >/dev/null 2>&1; then
  echo "[ERROR] conda was not found. Please install Miniconda or Mambaforge first."
  exit 1
fi

CONDA_BASE="$(conda info --base)"
# shellcheck disable=SC1091
source "${CONDA_BASE}/etc/profile.d/conda.sh"

if ! conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  echo "[INFO] Conda environment ${ENV_NAME} was not found. Creating it from environment.yml..."
  conda env create -f environment.yml
fi

conda activate "${ENV_NAME}"
export PATHOGEN_WEB_RUNTIME="${PATHOGEN_WEB_RUNTIME:-native}"

echo "[OK] Activated conda environment: ${ENV_NAME}"
echo "[OK] Open http://127.0.0.1:8080 after the server starts."
python tools/pipeline_web_app.py --host 127.0.0.1 --port 8080

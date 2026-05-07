#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

ENV_NAME="DPPP"
RCF_ENV_NAME="DPPP-rcf"
MAIN_ENV_FILE="environment.yml"
RCF_ENV_FILE="environment.recentrifuge.yml"

if ! command -v conda >/dev/null 2>&1; then
  echo "[ERROR] conda was not found. Please install Miniconda or Mambaforge first."
  exit 1
fi

CONDA_BASE="$(conda info --base)"
# shellcheck disable=SC1091
source "${CONDA_BASE}/etc/profile.d/conda.sh"

ensure_conda_env() {
  local env_name="$1"
  local env_file="$2"
  if ! conda env list | awk '{print $1}' | grep -qx "${env_name}"; then
    echo "[INFO] Conda environment ${env_name} was not found. Creating it from ${env_file}..."
    conda env create -f "${env_file}"
  fi
}

ensure_conda_env "${ENV_NAME}" "${MAIN_ENV_FILE}"
ensure_conda_env "${RCF_ENV_NAME}" "${RCF_ENV_FILE}"

conda activate "${ENV_NAME}"
export PATHOGEN_WEB_RUNTIME="${PATHOGEN_WEB_RUNTIME:-native}"
export PATHOGEN_RCF_ENV="${PATHOGEN_RCF_ENV:-${RCF_ENV_NAME}}"

echo "[OK] Activated conda environment: ${ENV_NAME}"
echo "[OK] Recentrifuge conda environment: ${PATHOGEN_RCF_ENV}"
echo "[OK] Open http://127.0.0.1:8080 after the server starts."
python tools/pipeline_web_app.py --host 127.0.0.1 --port 8080

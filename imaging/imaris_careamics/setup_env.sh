#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="${SCRIPT_DIR}/environment.yml"
ENV_NAME="${ENV_NAME:-imaris-careamics}"

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda is not on PATH. Load or initialize conda first." >&2
  exit 1
fi

CONDA_EXE="conda"
if command -v mamba >/dev/null 2>&1; then
  CONDA_EXE="mamba"
fi

if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  echo "Updating conda environment: ${ENV_NAME}"
  "${CONDA_EXE}" env update -n "${ENV_NAME}" -f "${ENV_FILE}" --prune
else
  echo "Creating conda environment: ${ENV_NAME}"
  "${CONDA_EXE}" env create -n "${ENV_NAME}" -f "${ENV_FILE}"
fi

echo "Environment ready."
echo "Activate with: conda activate ${ENV_NAME}"
echo "Check with: python ${SCRIPT_DIR}/check_env.py"

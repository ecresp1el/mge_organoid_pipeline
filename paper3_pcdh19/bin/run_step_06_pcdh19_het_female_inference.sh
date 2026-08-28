#!/usr/bin/env bash
# Supported inference-only entry point for formal Step 06. The exact full-fit
# Step 05 binary and count-informed models are applied unchanged to JZ-7--9.
set -Eeuo pipefail

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
source "${BUNDLE_DIR}/config/greatlakes.env"
: "${PAPER3_ROOT:?PAPER3_ROOT is required}"

PYTHON_BIN="${PAPER3_PYTHON_BIN:-python3}"
PYTHON_VERSION="$(${PYTHON_BIN} -c 'import sys; print("{}.{}".format(sys.version_info[0], sys.version_info[1]))')"
[[ "${PYTHON_VERSION}" == "3.6" ]] || {
  echo "Step 06 requires the Paper 3 Python 3.6 environment; observed ${PYTHON_VERSION}" >&2
  exit 2
}

LOCK="${BUNDLE_DIR}/config/step_06_pcdh19_het_female_inference.lock.json"
REQUIREMENTS="${BUNDLE_DIR}/config/step_06_pcdh19_het_female_inference.requirements.txt"
SCRIPT="${BUNDLE_DIR}/scripts/Step_06_PCDH19_HET_Female_Inference.py"
SOFTWARE_ROOT="${PAPER3_ROOT}/software/step_06_pcdh19_het_female_inference"
SITE_PACKAGES="${SOFTWARE_ROOT}/python3.6-site-packages"
LOG_ROOT="${PAPER3_ROOT}/logs/step_06_pcdh19_het_female_inference"
mkdir -p "${SOFTWARE_ROOT}" "${LOG_ROOT}"

if ! PYTHONPATH="${SITE_PACKAGES}" "${PYTHON_BIN}" -c 'import matplotlib, numpy; assert matplotlib.__version__ == "3.3.4"; assert numpy.__version__ == "1.19.5"' 2>/dev/null; then
  [[ ! -e "${SITE_PACKAGES}" ]] || {
    echo "Existing Step 06 dependency directory is invalid and will not be overwritten: ${SITE_PACKAGES}" >&2
    exit 3
  }
  BUILD_ROOT="$(mktemp -d "${SOFTWARE_ROOT}/.python-build.XXXXXX")"
  trap 'rm -rf "${BUILD_ROOT}"' EXIT
  mkdir -p "${BUILD_ROOT}/pip-tmp" "${BUILD_ROOT}/site-packages"
  TMPDIR="${BUILD_ROOT}/pip-tmp" "${PYTHON_BIN}" -m pip install --disable-pip-version-check --no-input --target "${BUILD_ROOT}/site-packages" -r "${REQUIREMENTS}"
  PYTHONPATH="${BUILD_ROOT}/site-packages" "${PYTHON_BIN}" -c 'import matplotlib, numpy; assert matplotlib.__version__ == "3.3.4"; assert numpy.__version__ == "1.19.5"'
  mv "${BUILD_ROOT}/site-packages" "${SITE_PACKAGES}"
  rm -rf "${BUILD_ROOT}"
  trap - EXIT
fi

STEP02_ROOT="${PAPER3_ROOT}/results/pcdh19_probe_audit"
STEP05_ROOT="${PAPER3_ROOT}/results/step_05_pcdh19_logistic_regression_baseline"
OUTPUT_ROOT="${PAPER3_ROOT}/results/step_06_pcdh19_het_female_inference"
SAMPLE_KEY="${BUNDLE_DIR}/config/sample_key.csv"
MPL_CACHE="${SOFTWARE_ROOT}/matplotlib-cache"
mkdir -p "${MPL_CACHE}"

RUN_STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
LOG_PATH="${LOG_ROOT}/run_${RUN_STAMP}_$$.log"
export PYTHONHASHSEED=0 LC_ALL=C TZ=UTC MPLCONFIGDIR="${MPL_CACHE}"

{
  echo "started_utc=$(date -u --iso-8601=seconds)"
  echo "host=$(hostname)"
  echo "step=06_pcdh19_het_female_inference"
  echo "pipeline=${SCRIPT}"
  echo "lock=${LOCK}"
  echo "step02_root=${STEP02_ROOT}"
  echo "step05_root=${STEP05_ROOT}"
  echo "output_root=${OUTPUT_ROOT}"
  PYTHONPATH="${SITE_PACKAGES}" "${PYTHON_BIN}" "${SCRIPT}" \
    --lock "${LOCK}" \
    --requirements "${REQUIREMENTS}" \
    --step02-root "${STEP02_ROOT}" \
    --step05-root "${STEP05_ROOT}" \
    --sample-key "${SAMPLE_KEY}" \
    --output-root "${OUTPUT_ROOT}" \
    "$@"
  echo "completed_utc=$(date -u --iso-8601=seconds)"
} 2>&1 | tee "${LOG_PATH}"

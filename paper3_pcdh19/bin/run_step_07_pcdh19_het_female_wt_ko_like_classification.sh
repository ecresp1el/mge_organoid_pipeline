#!/usr/bin/env bash
# Supported Step 07 entry point. Selects high-confidence thresholds from
# immutable held-out controls before loading Step 06 HET probabilities.
set -Eeuo pipefail

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
source "${BUNDLE_DIR}/config/greatlakes.env"
: "${PAPER3_ROOT:?PAPER3_ROOT is required}"

PYTHON_BIN="${PAPER3_PYTHON_BIN:-python3}"
PYTHON_VERSION="$(${PYTHON_BIN} -c 'import sys; print("{}.{}".format(sys.version_info[0], sys.version_info[1]))')"
[[ "${PYTHON_VERSION}" == "3.6" ]] || { echo "Step 07 requires Python 3.6; observed ${PYTHON_VERSION}" >&2; exit 2; }

LOCK="${BUNDLE_DIR}/config/step_07_pcdh19_het_female_wt_ko_like_classification.lock.json"
REQUIREMENTS="${BUNDLE_DIR}/config/step_07_pcdh19_het_female_wt_ko_like_classification.requirements.txt"
SCRIPT="${BUNDLE_DIR}/scripts/Step_07_PCDH19_HET_Female_WT_KO_Like_Classification.py"
SOFTWARE_ROOT="${PAPER3_ROOT}/software/step_07_pcdh19_het_female_wt_ko_like_classification"
SITE_PACKAGES="${SOFTWARE_ROOT}/python3.6-site-packages"
LOG_ROOT="${PAPER3_ROOT}/logs/step_07_pcdh19_het_female_wt_ko_like_classification"
mkdir -p "${SOFTWARE_ROOT}" "${LOG_ROOT}"

if ! PYTHONPATH="${SITE_PACKAGES}" "${PYTHON_BIN}" -c 'import matplotlib, numpy; assert matplotlib.__version__ == "3.3.4"; assert numpy.__version__ == "1.19.5"' 2>/dev/null; then
  [[ ! -e "${SITE_PACKAGES}" ]] || { echo "Existing Step 07 environment is invalid and will not be overwritten" >&2; exit 3; }
  BUILD_ROOT="$(mktemp -d "${SOFTWARE_ROOT}/.python-build.XXXXXX")"
  trap 'rm -rf "${BUILD_ROOT}"' EXIT
  mkdir -p "${BUILD_ROOT}/pip-tmp" "${BUILD_ROOT}/site-packages"
  TMPDIR="${BUILD_ROOT}/pip-tmp" "${PYTHON_BIN}" -m pip install --disable-pip-version-check --no-input --target "${BUILD_ROOT}/site-packages" -r "${REQUIREMENTS}"
  PYTHONPATH="${BUILD_ROOT}/site-packages" "${PYTHON_BIN}" -c 'import matplotlib, numpy; assert matplotlib.__version__ == "3.3.4"; assert numpy.__version__ == "1.19.5"'
  mv "${BUILD_ROOT}/site-packages" "${SITE_PACKAGES}"
  rm -rf "${BUILD_ROOT}"; trap - EXIT
fi

STEP05_ROOT="${PAPER3_ROOT}/results/step_05_pcdh19_logistic_regression_baseline"
STEP06_ROOT="${PAPER3_ROOT}/results/step_06_pcdh19_het_female_inference"
OUTPUT_ROOT="${PAPER3_ROOT}/results/step_07_pcdh19_het_female_wt_ko_like_classification"
MPL_CACHE="${SOFTWARE_ROOT}/matplotlib-cache"; mkdir -p "${MPL_CACHE}"
RUN_STAMP="$(date -u +%Y%m%dT%H%M%SZ)"; LOG_PATH="${LOG_ROOT}/run_${RUN_STAMP}_$$.log"
export PYTHONHASHSEED=0 LC_ALL=C TZ=UTC MPLCONFIGDIR="${MPL_CACHE}"

{
  echo "started_utc=$(date -u --iso-8601=seconds)"
  echo "host=$(hostname)"
  echo "step=07_pcdh19_het_female_wt_ko_like_classification"
  echo "pipeline=${SCRIPT}"
  echo "control_threshold_source=${STEP05_ROOT}"
  echo "het_probability_source=${STEP06_ROOT}"
  echo "output_root=${OUTPUT_ROOT}"
  PYTHONPATH="${SITE_PACKAGES}" "${PYTHON_BIN}" "${SCRIPT}" \
    --lock "${LOCK}" --requirements "${REQUIREMENTS}" \
    --step05-root "${STEP05_ROOT}" --step06-root "${STEP06_ROOT}" \
    --output-root "${OUTPUT_ROOT}" "$@"
  echo "completed_utc=$(date -u --iso-8601=seconds)"
} 2>&1 | tee "${LOG_PATH}"

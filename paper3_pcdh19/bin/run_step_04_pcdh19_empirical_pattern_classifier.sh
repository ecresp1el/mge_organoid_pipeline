#!/usr/bin/env bash
# PURPOSE
#   Supported Step 04 entry point for the empirical PCDH19 binary A/B/C
#   probe-pattern classifier.
#
# INPUTS
#   Reads only the checksum-locked Step 03 WT-male/KO-male classification-ready
#   table. It does not read HET cells or source matrices.
#
# OUTPUTS AND SCOPE
#   Publishes the eight-pattern empirical model, descriptive distribution,
#   validations, environment provenance, three diagnostic plots, and manifest
#   below results/step_04_pcdh19_empirical_pattern_classifier/. A timestamped
#   run log is written below the matching logs directory. No hard calls,
#   train/test split, confusion matrix, logistic regression, or HET prediction
#   is performed.
set -Eeuo pipefail

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
CONFIG="${BUNDLE_DIR}/config/greatlakes.env"
LOCK="${BUNDLE_DIR}/config/step_04_pcdh19_empirical_pattern_classifier.lock.json"
REQUIREMENTS="${BUNDLE_DIR}/config/step_04_pcdh19_empirical_pattern_classifier.requirements.txt"
SCRIPT="${BUNDLE_DIR}/scripts/Step_04_PCDH19_Empirical_Pattern_Classifier.py"

# shellcheck disable=SC1090
source "${CONFIG}"
: "${PAPER3_ROOT:?PAPER3_ROOT is required}"

PYTHON_BIN="${PAPER3_PYTHON_BIN:-python3}"
PYTHON_VERSION="$(${PYTHON_BIN} -c 'import sys; print("{}.{}".format(sys.version_info[0], sys.version_info[1]))')"
[[ "${PYTHON_VERSION}" == "3.6" ]] || {
  echo "Step 04 requires the Paper 3 Python 3.6 environment; observed ${PYTHON_VERSION}" >&2
  exit 2
}

SOFTWARE_ROOT="${PAPER3_ROOT}/software/step_04_pcdh19_empirical_pattern_classifier"
SITE_PACKAGES="${SOFTWARE_ROOT}/python3.6-site-packages"
LOG_ROOT="${PAPER3_ROOT}/logs/step_04_pcdh19_empirical_pattern_classifier"
mkdir -p "${SOFTWARE_ROOT}" "${LOG_ROOT}"

if ! PYTHONPATH="${SITE_PACKAGES}" "${PYTHON_BIN}" -c 'import matplotlib, numpy; assert matplotlib.__version__ == "3.3.4"; assert numpy.__version__ == "1.19.5"' 2>/dev/null; then
  [[ ! -e "${SITE_PACKAGES}" ]] || {
    echo "Existing Step 04 Python dependency directory is invalid and will not be overwritten: ${SITE_PACKAGES}" >&2
    exit 3
  }
  BUILD_ROOT="$(mktemp -d "${SOFTWARE_ROOT}/.python-build.XXXXXX")"
  trap 'rm -rf "${BUILD_ROOT}"' EXIT
  mkdir -p "${BUILD_ROOT}/pip-tmp" "${BUILD_ROOT}/site-packages"
  TMPDIR="${BUILD_ROOT}/pip-tmp" "${PYTHON_BIN}" -m pip install \
    --disable-pip-version-check --no-input --target "${BUILD_ROOT}/site-packages" \
    -r "${REQUIREMENTS}"
  PYTHONPATH="${BUILD_ROOT}/site-packages" "${PYTHON_BIN}" -c \
    'import matplotlib, numpy; assert matplotlib.__version__ == "3.3.4"; assert numpy.__version__ == "1.19.5"'
  mv "${BUILD_ROOT}/site-packages" "${SITE_PACKAGES}"
  rm -rf "${BUILD_ROOT}"
  trap - EXIT
fi

STEP03_ROOT="${PAPER3_ROOT}/results/step_03_pcdh19_genotype_classification_setup"
MPL_CACHE="${SOFTWARE_ROOT}/matplotlib-cache"
mkdir -p "${MPL_CACHE}"

RUN_STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
LOG_PATH="${LOG_ROOT}/run_${RUN_STAMP}_$$.log"
export PYTHONHASHSEED=0 LC_ALL=C TZ=UTC MPLCONFIGDIR="${MPL_CACHE}"

{
  echo "started_utc=$(date -u --iso-8601=seconds)"
  echo "host=$(hostname)"
  echo "step=04_pcdh19_empirical_pattern_classifier"
  echo "pipeline=${SCRIPT}"
  echo "lock=${LOCK}"
  echo "requirements=${REQUIREMENTS}"
  echo "step03_root=${STEP03_ROOT}"
  echo "paper3_root=${PAPER3_ROOT}"
  PYTHONPATH="${SITE_PACKAGES}" "${PYTHON_BIN}" "${SCRIPT}" \
    --lock "${LOCK}" \
    --requirements "${REQUIREMENTS}" \
    --step03-root "${STEP03_ROOT}" \
    --paper3-root "${PAPER3_ROOT}" \
    "$@"
  echo "completed_utc=$(date -u --iso-8601=seconds)"
} 2>&1 | tee "${LOG_PATH}"

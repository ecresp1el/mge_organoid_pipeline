#!/usr/bin/env bash
# PURPOSE
#   Supported Step 05 entry point for the unpenalized three-feature PCDH19
#   logistic-regression baseline.
#
# INPUTS AND SCOPE
#   Reads the validated Step 03 WT-male/KO-male table through the shared Step 04
#   encoder/reader and the manifested Step 04 empirical model for comparison.
#   Predictors are exactly A_detected, B_detected, and C_detected. WT=0 and
#   KO=1. The expanded Step 05 also performs leave-one-registered-sample-out
#   validation, always leaves 000 uncalled, and uses a fixed non-optimized 0.5
#   probability threshold for informative patterns. No HET cells,
#   interactions, nonlinear terms, UMI-count predictors, cell types, or
#   transcriptome features are used by the logistic model. Existing Step 05
#   fit outputs remain unchanged; held-out products are published in a
#   manifested validation subdirectory.
#   A second separately manifested Step 05 cohort uses WT males plus WT females
#   as WT ground truth and KO males as KO ground truth. HET females are excluded
#   from loading, fitting, validation, and comparison.
#   A descriptive A-negative module additionally summarizes unnormalized raw B,
#   C, and B+C UMI evidence and joint B/C counts. It fits no classifier and does
#   not alter any logistic-regression result or decision rule.
#   A paired count-informed implementation uses A_detected plus raw B/C
#   probe-level UMI/ligation evidence on the identical expanded-cohort LOSO
#   folds and compares it directly with the immutable binary model. Counts are
#   not transcript numbers. No weighting, interaction, threshold change, or
#   HET-cell loading is introduced.
#   A final read-only sample-level diagnostic verifies the immutable paired
#   cells against manifested raw probe rows and summarizes flow, distributions,
#   probe-state errors, sample pairs, and AUC orientation. It changes no model.
set -Eeuo pipefail

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
CONFIG="${BUNDLE_DIR}/config/greatlakes.env"
LOCK="${BUNDLE_DIR}/config/step_05_pcdh19_logistic_regression_baseline.lock.json"
REQUIREMENTS="${BUNDLE_DIR}/config/step_05_pcdh19_logistic_regression_baseline.requirements.txt"
STEP04_LOCK="${BUNDLE_DIR}/config/step_04_pcdh19_empirical_pattern_classifier.lock.json"
SCRIPT="${BUNDLE_DIR}/scripts/Step_05_PCDH19_Logistic_Regression_Baseline.py"

# shellcheck disable=SC1090
source "${CONFIG}"
: "${PAPER3_ROOT:?PAPER3_ROOT is required}"

PYTHON_BIN="${PAPER3_PYTHON_BIN:-python3}"
PYTHON_VERSION="$(${PYTHON_BIN} -c 'import sys; print("{}.{}".format(sys.version_info[0], sys.version_info[1]))')"
[[ "${PYTHON_VERSION}" == "3.6" ]] || {
  echo "Step 05 requires the Paper 3 Python 3.6 environment; observed ${PYTHON_VERSION}" >&2
  exit 2
}

SOFTWARE_ROOT="${PAPER3_ROOT}/software/step_05_pcdh19_logistic_regression_baseline"
SITE_PACKAGES="${SOFTWARE_ROOT}/python3.6-site-packages"
LOG_ROOT="${PAPER3_ROOT}/logs/step_05_pcdh19_logistic_regression_baseline"
mkdir -p "${SOFTWARE_ROOT}" "${LOG_ROOT}"

if ! PYTHONPATH="${SITE_PACKAGES}" "${PYTHON_BIN}" -c 'import matplotlib, numpy; assert matplotlib.__version__ == "3.3.4"; assert numpy.__version__ == "1.19.5"' 2>/dev/null; then
  [[ ! -e "${SITE_PACKAGES}" ]] || {
    echo "Existing Step 05 dependency directory is invalid and will not be overwritten: ${SITE_PACKAGES}" >&2
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
STEP04_ROOT="${PAPER3_ROOT}/results/step_04_pcdh19_empirical_pattern_classifier"
STEP02_ROOT="${PAPER3_ROOT}/results/pcdh19_probe_audit"
SAMPLE_KEY="${BUNDLE_DIR}/config/sample_key.csv"
MPL_CACHE="${SOFTWARE_ROOT}/matplotlib-cache"
mkdir -p "${MPL_CACHE}"

RUN_STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
LOG_PATH="${LOG_ROOT}/run_${RUN_STAMP}_$$.log"
export PYTHONHASHSEED=0 LC_ALL=C TZ=UTC MPLCONFIGDIR="${MPL_CACHE}"

{
  echo "started_utc=$(date -u --iso-8601=seconds)"
  echo "host=$(hostname)"
  echo "step=05_pcdh19_logistic_regression_baseline"
  echo "pipeline=${SCRIPT}"
  echo "lock=${LOCK}"
  echo "requirements=${REQUIREMENTS}"
  echo "step03_root=${STEP03_ROOT}"
  echo "step04_root=${STEP04_ROOT}"
  echo "step02_root=${STEP02_ROOT}"
  echo "sample_key=${SAMPLE_KEY}"
  echo "paper3_root=${PAPER3_ROOT}"
  PYTHONPATH="${SITE_PACKAGES}" "${PYTHON_BIN}" "${SCRIPT}" \
    --lock "${LOCK}" \
    --requirements "${REQUIREMENTS}" \
    --step04-lock "${STEP04_LOCK}" \
    --step03-root "${STEP03_ROOT}" \
    --step04-root "${STEP04_ROOT}" \
    --step02-root "${STEP02_ROOT}" \
    --sample-key "${SAMPLE_KEY}" \
    --paper3-root "${PAPER3_ROOT}" \
    --bundle-root "${BUNDLE_DIR}" \
    "$@"
  echo "completed_utc=$(date -u --iso-8601=seconds)"
} 2>&1 | tee "${LOG_PATH}"

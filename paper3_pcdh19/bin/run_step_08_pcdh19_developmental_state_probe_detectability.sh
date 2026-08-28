#!/usr/bin/env bash
# Supported Step 08 entry point. Joins frozen probe evidence to independent
# GSE94641 broad developmental-state labels; does not open or modify models.
set -Eeuo pipefail

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
source "${BUNDLE_DIR}/config/greatlakes.env"
: "${PAPER3_ROOT:?PAPER3_ROOT is required}"

PYTHON_BIN="${PAPER3_PYTHON_BIN:-python3}"
PYTHON_VERSION="$(${PYTHON_BIN} -c 'import sys; print("{}.{}".format(sys.version_info[0], sys.version_info[1]))')"
[[ "${PYTHON_VERSION}" == "3.6" ]] || { echo "Step 08 requires Python 3.6; observed ${PYTHON_VERSION}" >&2; exit 2; }

LOCK="${BUNDLE_DIR}/config/step_08_pcdh19_developmental_state_probe_detectability.lock.json"
REQUIREMENTS="${BUNDLE_DIR}/config/step_08_pcdh19_developmental_state_probe_detectability.requirements.txt"
SCRIPT="${BUNDLE_DIR}/scripts/Step_08_PCDH19_Developmental_State_Probe_Detectability.py"
SOFTWARE_ROOT="${PAPER3_ROOT}/software/step_08_pcdh19_developmental_state_probe_detectability"
SITE_PACKAGES="${SOFTWARE_ROOT}/python3.6-site-packages"
LOG_ROOT="${PAPER3_ROOT}/logs/step_08_pcdh19_developmental_state_probe_detectability"
mkdir -p "${SOFTWARE_ROOT}" "${LOG_ROOT}"

if ! PYTHONPATH="${SITE_PACKAGES}" "${PYTHON_BIN}" -c 'import matplotlib, numpy; assert matplotlib.__version__ == "3.3.4"; assert numpy.__version__ == "1.19.5"' 2>/dev/null; then
  [[ ! -e "${SITE_PACKAGES}" ]] || { echo "Existing Step 08 environment is invalid and will not be overwritten" >&2; exit 3; }
  BUILD_ROOT="$(mktemp -d "${SOFTWARE_ROOT}/.python-build.XXXXXX")"
  trap 'rm -rf "${BUILD_ROOT}"' EXIT
  mkdir -p "${BUILD_ROOT}/pip-tmp" "${BUILD_ROOT}/site-packages"
  TMPDIR="${BUILD_ROOT}/pip-tmp" "${PYTHON_BIN}" -m pip install --disable-pip-version-check --no-input --target "${BUILD_ROOT}/site-packages" -r "${REQUIREMENTS}"
  PYTHONPATH="${BUILD_ROOT}/site-packages" "${PYTHON_BIN}" -c 'import matplotlib, numpy; assert matplotlib.__version__ == "3.3.4"; assert numpy.__version__ == "1.19.5"'
  mv "${BUILD_ROOT}/site-packages" "${SITE_PACKAGES}"
  rm -rf "${BUILD_ROOT}"; trap - EXIT
fi

LABEL_ROOT="${PAPER3_ROOT}/results/mge_reference_mapping_gse94641/label_transfer_e15_5"
PROBE_ROOT="${PAPER3_ROOT}/results/pcdh19_probe_audit"
OUTPUT_ROOT="${PAPER3_ROOT}/results/step_08_pcdh19_developmental_state_probe_detectability"
MPL_CACHE="${SOFTWARE_ROOT}/matplotlib-cache"; mkdir -p "${MPL_CACHE}"
RUN_STAMP="$(date -u +%Y%m%dT%H%M%SZ)"; LOG_PATH="${LOG_ROOT}/run_${RUN_STAMP}_$$.log"
export PYTHONHASHSEED=0 LC_ALL=C TZ=UTC MPLCONFIGDIR="${MPL_CACHE}"

{
  echo "started_utc=$(date -u --iso-8601=seconds)"
  echo "host=$(hostname)"
  echo "step=08_pcdh19_developmental_state_probe_detectability"
  echo "pipeline=${SCRIPT}"
  echo "label_source=${LABEL_ROOT}"
  echo "probe_source=${PROBE_ROOT}"
  echo "output_root=${OUTPUT_ROOT}"
  PYTHONPATH="${SITE_PACKAGES}" "${PYTHON_BIN}" "${SCRIPT}" \
    --lock "${LOCK}" --requirements "${REQUIREMENTS}" \
    --sample-key "${BUNDLE_DIR}/config/sample_key.csv" \
    --label-root "${LABEL_ROOT}" --probe-root "${PROBE_ROOT}" \
    --output-root "${OUTPUT_ROOT}" "$@"
  echo "completed_utc=$(date -u --iso-8601=seconds)"
} 2>&1 | tee "${LOG_PATH}"


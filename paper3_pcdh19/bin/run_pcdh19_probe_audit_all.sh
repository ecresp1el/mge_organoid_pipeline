#!/usr/bin/env bash
# Supported Step 02a entry point for the locked 12-sample Pcdh19 probe audit.
# Exact scientific inputs, output tables, validation gates, restart behavior,
# and interpretation limits are defined in
# ../PIPELINE_IO_AND_BIOLOGICAL_SCOPE.md. Outputs are technical probe-pattern
# measurements and do not assign genotype, condition, or cell identity.
set -Eeuo pipefail

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
CONFIG="${BUNDLE_DIR}/config/greatlakes.env"
LOCK="${BUNDLE_DIR}/config/pcdh19_probe_audit.lock.json"
REQUIREMENTS="${BUNDLE_DIR}/config/pcdh19_probe_audit.requirements.txt"
SCRIPT="${BUNDLE_DIR}/scripts/pcdh19_probe_audit.py"

# shellcheck disable=SC1090
source "${CONFIG}"
: "${PAPER3_ROOT:?PAPER3_ROOT is required}"
: "${PAPER3_CELLRANGER_ROOT:?PAPER3_CELLRANGER_ROOT is required}"
: "${PAPER3_GRCM39_GTF:?PAPER3_GRCM39_GTF is required}"
: "${PAPER3_GRCM39_REFERENCE_JSON:?PAPER3_GRCM39_REFERENCE_JSON is required}"

PYTHON_BIN="${PAPER3_PYTHON_BIN:-python3}"
PYTHON_VERSION="$(${PYTHON_BIN} -c 'import sys; print("{}.{}".format(sys.version_info[0], sys.version_info[1]))')"
[[ "${PYTHON_VERSION}" == "3.6" ]] || {
  echo "Locked audit environment requires Python 3.6; observed ${PYTHON_VERSION}" >&2
  exit 2
}

SOFTWARE_ROOT="${PAPER3_ROOT}/software/pcdh19_probe_audit"
SITE_PACKAGES="${SOFTWARE_ROOT}/python3.6-site-packages"
mkdir -p "${SOFTWARE_ROOT}" "${PAPER3_ROOT}/logs/pcdh19_probe_audit"

if ! PYTHONPATH="${SITE_PACKAGES}" "${PYTHON_BIN}" -c 'import h5py, numpy; assert h5py.__version__ == "3.1.0"; assert numpy.__version__ == "1.19.5"' 2>/dev/null; then
  [[ ! -e "${SITE_PACKAGES}" ]] || {
    echo "Existing locked Python dependency directory is invalid and will not be overwritten: ${SITE_PACKAGES}" >&2
    exit 3
  }
  BUILD_ROOT="$(mktemp -d "${SOFTWARE_ROOT}/.python-build.XXXXXX")"
  trap 'rm -rf "${BUILD_ROOT}"' EXIT
  mkdir -p "${BUILD_ROOT}/pip-tmp" "${BUILD_ROOT}/site-packages"
  TMPDIR="${BUILD_ROOT}/pip-tmp" "${PYTHON_BIN}" -m pip install \
    --disable-pip-version-check --no-input --target "${BUILD_ROOT}/site-packages" \
    -r "${REQUIREMENTS}"
  PYTHONPATH="${BUILD_ROOT}/site-packages" "${PYTHON_BIN}" -c \
    'import h5py, numpy; assert h5py.__version__ == "3.1.0"; assert numpy.__version__ == "1.19.5"'
  mv "${BUILD_ROOT}/site-packages" "${SITE_PACKAGES}"
  rm -rf "${BUILD_ROOT}"
  trap - EXIT
fi

RUN_STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
LOG_PATH="${PAPER3_ROOT}/logs/pcdh19_probe_audit/run_${RUN_STAMP}_$$.log"
export PAPER3_GRCM39_GTF PAPER3_GRCM39_REFERENCE_JSON
export PYTHONHASHSEED=0 LC_ALL=C TZ=UTC

{
  echo "started_utc=$(date -u --iso-8601=seconds)"
  echo "host=$(hostname)"
  echo "pipeline=${SCRIPT}"
  echo "lock=${LOCK}"
  echo "cellranger_root=${PAPER3_CELLRANGER_ROOT}"
  echo "paper3_root=${PAPER3_ROOT}"
  PYTHONPATH="${SITE_PACKAGES}" "${PYTHON_BIN}" "${SCRIPT}" run-all \
    --lock "${LOCK}" \
    --cellranger-root "${PAPER3_CELLRANGER_ROOT}" \
    --paper3-root "${PAPER3_ROOT}" \
    "$@"
  echo "completed_utc=$(date -u --iso-8601=seconds)"
} 2>&1 | tee "${LOG_PATH}"

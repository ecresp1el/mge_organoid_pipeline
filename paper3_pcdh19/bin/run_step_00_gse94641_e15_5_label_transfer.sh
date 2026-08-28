#!/usr/bin/env bash
# E15.5-focused rapid MGE label transfer using the full validated GSE94641 PCA.
set -Eeuo pipefail

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
source "${BUNDLE_DIR}/config/greatlakes.env"
: "${PAPER3_ROOT:?PAPER3_ROOT is required}"
: "${PAPER3_CELLRANGER_ROOT:?PAPER3_CELLRANGER_ROOT is required}"

PYTHON_BIN="${PAPER3_PYTHON_BIN:-python3}"
PYTHON_VERSION="$(${PYTHON_BIN} -c 'import sys; print("{}.{}".format(sys.version_info[0], sys.version_info[1]))')"
[[ "${PYTHON_VERSION}" == "3.6" ]] || { echo "GSE94641 mapping requires Python 3.6; observed ${PYTHON_VERSION}" >&2; exit 2; }

SCRIPT="${BUNDLE_DIR}/scripts/reference_mapping/02_transfer_GSE94641_E15_5.py"
LOCK="${BUNDLE_DIR}/config/gse94641_e15_5_label_transfer.lock.json"
REQUIREMENTS="${BUNDLE_DIR}/config/gse94641_e15_5_label_transfer.requirements.txt"
REFERENCE_ROOT="${BUNDLE_DIR}/references/GSE94641"
VALIDATION_ROOT="${PAPER3_ROOT}/results/mge_reference_mapping_gse94641/reference_validation"
OUTPUT_ROOT="${PAPER3_ROOT}/results/mge_reference_mapping_gse94641/label_transfer_e15_5"
SOFTWARE_ROOT="${PAPER3_ROOT}/software/gse94641_e15_5_label_transfer"
SITE_PACKAGES="${SOFTWARE_ROOT}/python3.6-site-packages"
LOG_ROOT="${PAPER3_ROOT}/logs/mge_reference_mapping_gse94641"
mkdir -p "${SOFTWARE_ROOT}" "${LOG_ROOT}"

if ! PYTHONPATH="${SITE_PACKAGES}" "${PYTHON_BIN}" -c 'import h5py, matplotlib, numpy, scipy; assert h5py.__version__ == "3.1.0"; assert matplotlib.__version__ == "3.3.4"; assert numpy.__version__ == "1.19.5"; assert scipy.__version__ == "1.5.4"' 2>/dev/null; then
  [[ ! -e "${SITE_PACKAGES}" ]] || { echo "Existing mapping environment is invalid and will not be overwritten" >&2; exit 3; }
  BUILD_ROOT="$(mktemp -d "${SOFTWARE_ROOT}/.python-build.XXXXXX")"
  trap 'rm -rf "${BUILD_ROOT}"' EXIT
  mkdir -p "${BUILD_ROOT}/pip-tmp" "${BUILD_ROOT}/site-packages"
  TMPDIR="${BUILD_ROOT}/pip-tmp" "${PYTHON_BIN}" -m pip install --disable-pip-version-check --no-input --target "${BUILD_ROOT}/site-packages" -r "${REQUIREMENTS}"
  PYTHONPATH="${BUILD_ROOT}/site-packages" "${PYTHON_BIN}" -c 'import h5py, matplotlib, numpy, scipy'
  mv "${BUILD_ROOT}/site-packages" "${SITE_PACKAGES}"
  rm -rf "${BUILD_ROOT}"; trap - EXIT
fi

RUN_STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
LOG_PATH="${LOG_ROOT}/step_00_e15_5_label_transfer_${RUN_STAMP}_$$.log"
MPL_CACHE="${SOFTWARE_ROOT}/matplotlib-cache"; mkdir -p "${MPL_CACHE}"
export PYTHONHASHSEED=0 LC_ALL=C TZ=UTC MPLCONFIGDIR="${MPL_CACHE}"

{
  echo "started_utc=$(date -u --iso-8601=seconds)"
  echo "host=$(hostname)"
  echo "workflow=gse94641_e15_5_label_transfer"
  echo "reference_scope=full_age_pca_e15_5_primary_neighbors"
  echo "output_root=${OUTPUT_ROOT}"
  PYTHONPATH="${SITE_PACKAGES}" "${PYTHON_BIN}" "${SCRIPT}" \
    --lock "${LOCK}" \
    --reference-root "${REFERENCE_ROOT}" \
    --validation-root "${VALIDATION_ROOT}" \
    --query-root "${PAPER3_CELLRANGER_ROOT}" \
    --sample-key "${BUNDLE_DIR}/config/sample_key.csv" \
    --output-root "${OUTPUT_ROOT}" "$@"
  echo "completed_utc=$(date -u --iso-8601=seconds)"
} 2>&1 | tee "${LOG_PATH}"

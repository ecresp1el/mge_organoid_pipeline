#!/usr/bin/env bash
# Step 00 reference-only entry point. Downloads/validates GSE94641 and stops;
# it does not load query cells or run label transfer.
set -Eeuo pipefail

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
source "${BUNDLE_DIR}/config/greatlakes.env"
: "${PAPER3_ROOT:?PAPER3_ROOT is required}"

PYTHON_BIN="${PAPER3_PYTHON_BIN:-python3}"
SCRIPT="${BUNDLE_DIR}/scripts/reference_mapping/01_download_validate_GSE94641.py"
LOCK="${BUNDLE_DIR}/config/gse94641_reference_validation.lock.json"
REFERENCE_ROOT="${BUNDLE_DIR}/references/GSE94641"
OUTPUT_ROOT="${PAPER3_ROOT}/results/mge_reference_mapping_gse94641/reference_validation"
LOG_ROOT="${PAPER3_ROOT}/logs/mge_reference_mapping_gse94641"
mkdir -p "${LOG_ROOT}"

RUN_STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
LOG_PATH="${LOG_ROOT}/step_00_reference_validation_${RUN_STAMP}_$$.log"
export PYTHONHASHSEED=0 LC_ALL=C TZ=UTC

{
  echo "started_utc=$(date -u --iso-8601=seconds)"
  echo "host=$(hostname)"
  echo "step=00_gse94641_reference_download_validation"
  echo "reference_root=${REFERENCE_ROOT}"
  echo "output_root=${OUTPUT_ROOT}"
  "${PYTHON_BIN}" "${SCRIPT}" \
    --lock "${LOCK}" \
    --reference-root "${REFERENCE_ROOT}" \
    --output-root "${OUTPUT_ROOT}" "$@"
  echo "completed_utc=$(date -u --iso-8601=seconds)"
} 2>&1 | tee "${LOG_PATH}"

#!/usr/bin/env bash
# PURPOSE
#   Supported Step 03 entry point for creation of the PCDH19 WT-male/KO-male
#   genotype-classification-ready intermediate table.
#
# INPUTS
#   Reads the checksum-locked registered sample key and the existing validated
#   Step 02a one-row-per-filtered-barcode Pcdh19 probe tables. Source delivery
#   matrices and frozen Step 02a outputs remain read-only.
#
# OUTPUTS AND SCOPE
#   Publishes a restart-safe package beneath
#   results/step_03_pcdh19_genotype_classification_setup/ and writes a
#   timestamped log beneath logs/step_03_pcdh19_genotype_classification_setup/.
#   This step joins labels and represents raw probe counts/detection only. It
#   does not split data, fit a classifier, evaluate predictions, or score HET
#   cells.
set -Eeuo pipefail

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
CONFIG="${BUNDLE_DIR}/config/greatlakes.env"
LOCK="${BUNDLE_DIR}/config/step_03_pcdh19_genotype_classification_setup.lock.json"
SAMPLE_KEY="${BUNDLE_DIR}/config/sample_key.csv"
SCRIPT="${BUNDLE_DIR}/scripts/Step_03_PCDH19_Genotype_Classification_Setup.py"

# shellcheck disable=SC1090
source "${CONFIG}"
: "${PAPER3_ROOT:?PAPER3_ROOT is required}"

PYTHON_BIN="${PAPER3_PYTHON_BIN:-python3}"
PYTHON_VERSION="$(${PYTHON_BIN} -c 'import sys; print("{}.{}".format(sys.version_info[0], sys.version_info[1]))')"
[[ "${PYTHON_VERSION}" == "3.6" ]] || {
  echo "Step 03 requires the Paper 3 Python 3.6 environment; observed ${PYTHON_VERSION}" >&2
  exit 2
}

PROBE_AUDIT_ROOT="${PAPER3_ROOT}/results/pcdh19_probe_audit"
LOG_ROOT="${PAPER3_ROOT}/logs/step_03_pcdh19_genotype_classification_setup"
mkdir -p "${LOG_ROOT}"

RUN_STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
LOG_PATH="${LOG_ROOT}/run_${RUN_STAMP}_$$.log"
export PYTHONHASHSEED=0 LC_ALL=C TZ=UTC

{
  echo "started_utc=$(date -u --iso-8601=seconds)"
  echo "host=$(hostname)"
  echo "step=03_pcdh19_genotype_classification_setup"
  echo "pipeline=${SCRIPT}"
  echo "lock=${LOCK}"
  echo "sample_key=${SAMPLE_KEY}"
  echo "probe_audit_root=${PROBE_AUDIT_ROOT}"
  echo "paper3_root=${PAPER3_ROOT}"
  "${PYTHON_BIN}" "${SCRIPT}" \
    --lock "${LOCK}" \
    --sample-key "${SAMPLE_KEY}" \
    --probe-audit-root "${PROBE_AUDIT_ROOT}" \
    --paper3-root "${PAPER3_ROOT}" \
    "$@"
  echo "completed_utc=$(date -u --iso-8601=seconds)"
} 2>&1 | tee "${LOG_PATH}"

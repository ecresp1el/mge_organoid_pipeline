#!/usr/bin/env bash
# Step 02b: exact construct-level compatibility audit between the three custom
# Flex EGFP probes and the reporter sequence specified for D4/XEGFP.
set -Eeuo pipefail

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
CONFIG="${BUNDLE_DIR}/config/greatlakes.env"
LOCK="${BUNDLE_DIR}/config/xgfp_probe_audit.lock.json"
SCRIPT="${BUNDLE_DIR}/scripts/xgfp_probe_compatibility_audit.py"

# shellcheck disable=SC1090
source "${CONFIG}"
: "${PAPER3_ROOT:?PAPER3_ROOT is required}"
: "${PAPER3_CELLRANGER_ROOT:?PAPER3_CELLRANGER_ROOT is required}"

PYTHON_BIN="${PAPER3_PYTHON_BIN:-python3}"
PROBE_SET="${PAPER3_CELLRANGER_ROOT}/probe_set.csv"
mkdir -p "${PAPER3_ROOT}/logs/xgfp_probe_audit"
RUN_STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
LOG_PATH="${PAPER3_ROOT}/logs/xgfp_probe_audit/run_${RUN_STAMP}_$$.log"

{
  echo "started_utc=$(date -u --iso-8601=seconds)"
  echo "host=$(hostname)"
  echo "pipeline=${SCRIPT}"
  echo "lock=${LOCK}"
  echo "probe_set=${PROBE_SET}"
  "${PYTHON_BIN}" "${SCRIPT}" \
    --lock "${LOCK}" \
    --probe-set "${PROBE_SET}" \
    --paper3-root "${PAPER3_ROOT}"
  echo "completed_utc=$(date -u --iso-8601=seconds)"
} 2>&1 | tee "${LOG_PATH}"

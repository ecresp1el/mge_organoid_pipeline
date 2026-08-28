#!/usr/bin/env bash
# PURPOSE
#   Step 02b: test exact construct-level sequence compatibility between the
#   three delivered Flex EGFP probes and the EGFP CDS specified for D4/XEGFP.
#
# INPUTS
#   greatlakes.env supplies PAPER3_CELLRANGER_ROOT and PAPER3_ROOT. The script
#   reads the delivered probe_set.csv, xgfp_probe_audit.lock.json, and a cached
#   or downloaded checksum-locked NCBI U55762.1 FASTA.
#
# COMPUTATION AND OUTPUTS
#   The Python audit requires one unique 50/50 reverse-complement EGFP-CDS match
#   for every probe and atomically publishes alignments, validations, source
#   provenance, extracted CDS, conclusion, software versions, and checksums
#   below results/xgfp_probe_audit/. This runner adds a timestamped log below
#   logs/xgfp_probe_audit/. An existing validated package is verified/reused.
#
# SCIENTIFIC SCOPE
#   A PASS establishes theoretical sequence compatibility only. No count
#   matrix, barcode, sample label, expression state, or genotype is analyzed;
#   reporter transcription or functional GFP expression is not established.
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

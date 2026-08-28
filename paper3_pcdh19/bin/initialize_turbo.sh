#!/usr/bin/env bash
# PURPOSE
#   Initialize the Paper 3 output namespace configured in greatlakes.env.
#
# INPUTS
#   config/greatlakes.env supplies PAPER3_ROOT. TURBO_ROOT_README.md and
#   HANDOFF.md are the documentation copied into that root.
#
# OUTPUTS AND RESTART BEHAVIOR
#   Creates inputs/, results/, logs/, jobs/, and final_figures/ idempotently;
#   then refreshes PAPER3_ROOT/README.md and HANDOFF.md with repository copies.
#   Existing scientific assets below those directories are not enumerated,
#   read, moved, or removed. The two copied documentation files may be replaced.
#
# SCIENTIFIC SCOPE
#   Directory preparation only: no expression data, barcode, reference, sample
#   metadata, normalization, or biological inference is processed. See
#   ../PIPELINE_IO_AND_BIOLOGICAL_SCOPE.md for the complete contract.
set -Eeuo pipefail

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
CONFIG="${BUNDLE_DIR}/config/greatlakes.env"

# shellcheck disable=SC1090
source "${CONFIG}"
: "${PAPER3_ROOT:?PAPER3_ROOT is required}"

mkdir -p \
  "${PAPER3_ROOT}/inputs" \
  "${PAPER3_ROOT}/results" \
  "${PAPER3_ROOT}/logs" \
  "${PAPER3_ROOT}/jobs" \
  "${PAPER3_ROOT}/final_figures"

cp -p "${BUNDLE_DIR}/templates/TURBO_ROOT_README.md" "${PAPER3_ROOT}/README.md"
cp -p "${BUNDLE_DIR}/HANDOFF.md" "${PAPER3_ROOT}/HANDOFF.md"

echo "Paper 3 Turbo workstream initialized: ${PAPER3_ROOT}"
find "${PAPER3_ROOT}" -maxdepth 1 -mindepth 1 -printf '%f\n' | sort

#!/usr/bin/env bash
set -Eeuo pipefail

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
CONFIG="${BUNDLE_DIR}/config/greatlakes.env"

# shellcheck disable=SC1090
source "${CONFIG}"
: "${PAPER2_ROOT:?PAPER2_ROOT is required}"

mkdir -p \
  "${PAPER2_ROOT}/results" \
  "${PAPER2_ROOT}/logs" \
  "${PAPER2_ROOT}/jobs" \
  "${PAPER2_ROOT}/final_figures"

cp -p "${BUNDLE_DIR}/templates/TURBO_ROOT_README.md" "${PAPER2_ROOT}/README.md"
cp -p "${BUNDLE_DIR}/HANDOFF.md" "${PAPER2_ROOT}/HANDOFF.md"

echo "Paper 2 Turbo workstream initialized: ${PAPER2_ROOT}"
find "${PAPER2_ROOT}" -maxdepth 1 -mindepth 1 -printf '%f\n' | sort

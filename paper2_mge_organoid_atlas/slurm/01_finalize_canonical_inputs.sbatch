#!/usr/bin/env bash

set -Eeuo pipefail

: "${PAPER2_CANONICAL_STAGING:?canonical staging directory is required}"
: "${PAPER2_CANONICAL_DIR:?canonical destination is required}"
: "${PAPER2_CANONICAL_MANIFEST:?canonical manifest is required}"
: "${CONDA_ENV_BIN:?Python environment bin directory is required}"

FINALIZER="${PAPER2_CANONICAL_STAGING}/code/01_finalize_canonical_inputs.py"

echo "Canonical input finalization started: $(date --iso-8601=seconds)"
"${CONDA_ENV_BIN}/python" "${FINALIZER}" \
  --manifest "${PAPER2_CANONICAL_MANIFEST}" \
  --staging-dir "${PAPER2_CANONICAL_STAGING}" \
  --canonical-dir "${PAPER2_CANONICAL_DIR}"
echo "Canonical input finalization completed: $(date --iso-8601=seconds)"

#!/usr/bin/env bash
# Execute frozen input preparation and Python HiCAT pilot; never use live repo code.
set -Eeuo pipefail
: "${HICAT_RUN_DIR:?}" "${HICAT_PROJECT_ROOT:?}"
export PYTHONPATH="${HICAT_RUN_DIR}/code:${HICAT_RUN_DIR}/upstream"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export PYTHONHASHSEED=0 MPLCONFIGDIR="${HICAT_RUN_DIR}/provenance/matplotlib"
export TMPDIR="${HICAT_RUN_DIR}/scratch"
mkdir -p "${MPLCONFIGDIR}" "${TMPDIR}"
# Preserve failure evidence; no automatic retries or success-shaped fallbacks.
trap 'printf "FAILED line=%s status=%s\n" "${LINENO}" "$?" > "${HICAT_RUN_DIR}/COMPUTATION_FAILED.txt"' ERR
# Recheck every frozen code/config/reference asset before executing it.
/home/elcrespo/miniconda3/envs/pcdh19-hicat-py38/bin/python - "${HICAT_RUN_DIR}" <<'PY'
import hashlib
import json
from pathlib import Path
import sys
root = Path(sys.argv[1])
for row in json.loads((root/'provenance/frozen_source_manifest.json').read_text()):
    path = root/row['path']
    if path.stat().st_size != row['bytes'] or hashlib.sha256(path.read_bytes()).hexdigest() != row['sha256']:
        raise ValueError('Frozen asset identity changed: '+row['path'])
print('Frozen source/configuration/reference manifest: PASS', flush=True)
PY
/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python -m hicat.prepare \
  --run-dir "${HICAT_RUN_DIR}" --project-root "${HICAT_PROJECT_ROOT}"
/home/elcrespo/miniconda3/envs/pcdh19-hicat-py38/bin/python -m hicat.cli --run-dir "${HICAT_RUN_DIR}"
printf 'COMPLETED\nstatus=IN_REVIEW\nscope=technical_pilot_only\n' > "${HICAT_RUN_DIR}/COMPUTATION_SUCCESS.txt"

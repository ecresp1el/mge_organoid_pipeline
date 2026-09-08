#!/usr/bin/env bash
# Run frozen Step 07 code; full-data operations are scoring/display only.
set -Eeuo pipefail
: "${HICAT_RUN_DIR:?}"
export PYTHONPATH="${HICAT_RUN_DIR}/code:${HICAT_RUN_DIR}/upstream"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
export PYTHONHASHSEED=0 MPLCONFIGDIR="${HICAT_RUN_DIR}/provenance/matplotlib"
export TMPDIR="${HICAT_RUN_DIR}/scratch" PYTHONDONTWRITEBYTECODE=1
mkdir -p "${MPLCONFIGDIR}" "${TMPDIR}"
trap 'printf "FAILED line=%s status=%s\n" "${LINENO}" "$?" > "${HICAT_RUN_DIR}/COMPUTATION_FAILED.txt"' ERR
/home/elcrespo/miniconda3/envs/pcdh19-hicat-py38/bin/python - "${HICAT_RUN_DIR}" <<'PY'
import hashlib
import json
from pathlib import Path
import sys
root = Path(sys.argv[1])
for row in json.loads((root/'provenance/frozen_source_manifest.json').read_text()):
    path = root/row['path']
    if path.stat().st_size != row['bytes'] or hashlib.sha256(path.read_bytes()).hexdigest() != row['sha256']:
        raise ValueError('Frozen asset changed: '+row['path'])
print('Frozen source/config/reference manifest: PASS', flush=True)
PY
/home/elcrespo/miniconda3/envs/pcdh19-hicat-py38/bin/python -m unittest discover \
  -s "${HICAT_RUN_DIR}/code/tests" -p 'test_hicat_validation*.py' -v \
  > "${HICAT_RUN_DIR}/provenance/frozen_tests.log" 2>&1
/home/elcrespo/miniconda3/envs/pcdh19-hicat-py38/bin/python -m hicat.validation_workflow \
  --run-dir "${HICAT_RUN_DIR}"
printf 'COMPLETED\nstatus=IN_REVIEW\nscope=pilot_hierarchy_validation\n' > "${HICAT_RUN_DIR}/COMPUTATION_SUCCESS.txt"

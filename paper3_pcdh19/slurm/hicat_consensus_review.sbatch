#!/usr/bin/env bash
# Validate and report saved consensus checkpoints; no scientific refitting.
set -Eeuo pipefail
: "${HICAT_RUN_DIR:?}"
export PYTHONPATH="${HICAT_RUN_DIR}/code:${HICAT_RUN_DIR}/upstream"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 NUMBA_NUM_THREADS=1
export MPLCONFIGDIR="${HICAT_RUN_DIR}/provenance/matplotlib"
export TMPDIR="${HICAT_RUN_DIR}/scratch" PYTHONDONTWRITEBYTECODE=1
mkdir -p "${MPLCONFIGDIR}" "${TMPDIR}"
trap 'printf "FAILED line=%s status=%s\n" "${LINENO}" "$?" > "${HICAT_RUN_DIR}/COMPUTATION_FAILED.txt"' ERR
/home/elcrespo/miniconda3/envs/pcdh19-hicat-py38/bin/python - "${HICAT_RUN_DIR}" <<'PY'
from pathlib import Path
import json,sys,subprocess
from hicat.provenance import sha256
root=Path(sys.argv[1])
for row in json.loads((root/'provenance/frozen_source_manifest.json').read_text()):
    path=root/row['path']
    if path.stat().st_size!=row['bytes'] or sha256(path)!=row['sha256']:
        raise ValueError('Frozen source/input changed: '+row['path'])
if subprocess.check_output([sys.executable,'-m','pip','freeze'],text=True)!=(root/'config/python_pip_freeze.txt').read_text():
    raise ValueError('Environment differs from frozen dependency record')
print('Frozen source and environment verification passed',flush=True)
PY
/home/elcrespo/miniconda3/envs/pcdh19-hicat-py38/bin/python -m hicat.consensus_review_workflow --run-dir "${HICAT_RUN_DIR}"

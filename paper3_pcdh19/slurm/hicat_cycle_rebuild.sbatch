#!/bin/bash
#SBATCH --account=parent0
#SBATCH --partition=standard
#SBATCH --cpus-per-task=8
#SBATCH --mem=150G
#SBATCH --time=2-00:00:00
set -euo pipefail
: "${HICAT_CYCLE_RUN:?}"
: "${HICAT_CYCLE_REFERENCE:?}"
export PYTHONPATH="${HICAT_CYCLE_RUN}/code:${HICAT_CYCLE_REFERENCE}"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 NUMBA_NUM_THREADS=1
export MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1
export MPLCONFIGDIR="${SLURM_TMPDIR:-/tmp}/pcdh_hicat_mpl_${SLURM_JOB_ID}"
mkdir -p "${MPLCONFIGDIR}"
stage="${1:?stage}"; iteration="${2:--1}"
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then iteration="${SLURM_ARRAY_TASK_ID}"; fi
/home/elcrespo/miniconda3/envs/pcdh19-hicat-py38/bin/python -m hicat.cycle_rebuild --run-dir "${HICAT_CYCLE_RUN}" --stage "${stage}" --iteration "${iteration}"

#!/bin/bash
#SBATCH --job-name=pcdh_cycle_compare
#SBATCH --account=parent0
#SBATCH --partition=standard
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00
#SBATCH --output=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/logs/cycle_comparison_20261103-%j.log
set -euo pipefail
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 NUMBA_NUM_THREADS=16
export MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1
export MPLCONFIGDIR="${SLURM_TMPDIR:-/tmp}/pcdh_cycle_mpl_${SLURM_JOB_ID}"
/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline/paper3_pcdh19/scripts/compare_pcdh19_cell_cycle_umaps.py "$@"

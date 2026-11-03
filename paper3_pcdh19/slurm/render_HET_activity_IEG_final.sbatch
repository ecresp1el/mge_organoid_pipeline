#!/bin/bash -l
#SBATCH --job-name=pcdh_IEG_final
#SBATCH --account=parent0
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/logs/activity_IEG_final-%j.log
set -euo pipefail
module load R/4.4.3
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
export PYTHONPATH=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/software/gse94641_e15_5_label_transfer/python3.6-site-packages
export MPLCONFIGDIR=/tmp/pcdh_IEG_final_${SLURM_JOB_ID}
repo=/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
python3 "${repo}/paper3_pcdh19/scripts/render_HET_activity_IEG_final.py" "$@"

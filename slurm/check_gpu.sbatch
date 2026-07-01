#!/usr/bin/env bash
# Check the Imaris CAREamics conda environment and GPU visibility on Great Lakes.
#
# Submit:
#   mkdir -p /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs
#   sbatch slurm/check_gpu.sbatch

#SBATCH --job-name=imaris-gpu-check
#SBATCH --account=parent0
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/imaris-gpu-check-%j.out
#SBATCH --error=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/imaris-gpu-check-%j.err

set -euo pipefail

SAMPLE_CONFIG="${SAMPLE_CONFIG:-}"

if [[ -n "${SAMPLE_CONFIG}" ]]; then
  source "${SAMPLE_CONFIG}"
fi

REPO_ROOT="${REPO_ROOT:-/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline}"
PROJECT_ROOT="${PROJECT_ROOT:-/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder}"
OUTPUT_DIR="${OUTPUT_DIR:-${PROJECT_ROOT}/results/imaging/imaris_careamics}"
CONDA_BASE="${CONDA_BASE:-/home/elcrespo/miniconda3}"
CONDA_ENV="${CONDA_ENV:-imaris-careamics}"

mkdir -p "${OUTPUT_DIR}/logs" "${PROJECT_ROOT}/logs" "${PROJECT_ROOT}/jobs"
mkdir -p "${OUTPUT_DIR}/hdf5_plugins_empty"
export HDF5_PLUGIN_PATH="${HDF5_PLUGIN_PATH:-${OUTPUT_DIR}/hdf5_plugins_empty}"
exec > >(tee -a "${OUTPUT_DIR}/logs/check_gpu_${SLURM_JOB_ID:-manual}.log") 2>&1

echo "Imaris CAREamics GPU check started: $(date)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-}"
echo "REPO_ROOT=${REPO_ROOT}"
echo "PROJECT_ROOT=${PROJECT_ROOT}"
echo "OUTPUT_DIR=${OUTPUT_DIR}"
echo "CONDA_ENV=${CONDA_ENV}"

if command -v nvidia-smi >/dev/null 2>&1; then
  nvidia-smi
else
  echo "WARNING: nvidia-smi is not on PATH."
fi

source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"
python --version
which python
python "${REPO_ROOT}/imaging/imaris_careamics/check_env.py"

echo "Imaris CAREamics GPU check finished: $(date)"

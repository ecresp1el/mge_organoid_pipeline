#!/usr/bin/env bash
# Run the full Imaris .ims CAREamics denoising workflow on Great Lakes.
#
# Submit:
#   IMS_INPUT=/path/to/file.ims OUTPUT_DIR=/nfs/turbo/.../results/imaging/imaris_careamics/my_sample \
#     sbatch slurm/run_full_denoising.sbatch

#SBATCH --job-name=imaris-n2v-full
#SBATCH --account=parent0
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/imaris-n2v-full-%j.out
#SBATCH --error=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/imaris-n2v-full-%j.err

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
IMS_INPUT="${IMS_INPUT:-}"
EPOCHS="${EPOCHS:-50}"
NUM_STEPS="${NUM_STEPS:-200}"
BATCH_SIZE="${BATCH_SIZE:-1}"
PATCH_SIZE_Z="${PATCH_SIZE_Z:-16}"
PATCH_SIZE_YX="${PATCH_SIZE_YX:-64}"
REQUIRE_GPU="${REQUIRE_GPU:-true}"

if [[ -z "${IMS_INPUT}" ]]; then
  echo "ERROR: set IMS_INPUT=/path/to/file.ims before submitting." >&2
  exit 2
fi

mkdir -p "${OUTPUT_DIR}/logs" "${PROJECT_ROOT}/logs" "${PROJECT_ROOT}/jobs"
mkdir -p "${OUTPUT_DIR}/hdf5_plugins_empty"
export HDF5_PLUGIN_PATH="${HDF5_PLUGIN_PATH:-${OUTPUT_DIR}/hdf5_plugins_empty}"
export TMPDIR="${TMPDIR:-/tmp/${USER:-elcrespo}/imaris_careamics_${SLURM_JOB_ID:-manual}}"
mkdir -p "${TMPDIR}"
exec > >(tee -a "${OUTPUT_DIR}/logs/run_full_denoising_${SLURM_JOB_ID:-manual}.log") 2>&1

echo "Imaris CAREamics full denoising started: $(date)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-}"
echo "REPO_ROOT=${REPO_ROOT}"
echo "PROJECT_ROOT=${PROJECT_ROOT}"
echo "OUTPUT_DIR=${OUTPUT_DIR}"
echo "IMS_INPUT=${IMS_INPUT}"
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

gpu_arg=()
if [[ "${REQUIRE_GPU}" == "true" ]]; then
  gpu_arg+=(--require-gpu)
fi

python "${REPO_ROOT}/imaging/imaris_careamics/denoise_ims_careamics.py" \
  --input "${IMS_INPUT}" \
  --output-dir "${OUTPUT_DIR}" \
  --epochs "${EPOCHS}" \
  --num-steps "${NUM_STEPS}" \
  --batch-size "${BATCH_SIZE}" \
  --patch-size-z "${PATCH_SIZE_Z}" \
  --patch-size-yx "${PATCH_SIZE_YX}" \
  "${gpu_arg[@]}"

echo "Imaris CAREamics full denoising finished: $(date)"
find "${OUTPUT_DIR}/full_run" -maxdepth 3 -type f | sort

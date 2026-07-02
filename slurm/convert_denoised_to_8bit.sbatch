#!/usr/bin/env bash
# Convert CAREamics denoised OME-TIFFs to fixed-scale 8-bit display outputs.
#
# Submit:
#   FULL_RUN_DIR=/nfs/turbo/.../results/imaging/imaris_careamics/sample/full_run \
#     sbatch slurm/convert_denoised_to_8bit.sbatch

#SBATCH --job-name=imaris-8bit
#SBATCH --account=parent0
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/imaris-8bit-%j.out
#SBATCH --error=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/imaris-8bit-%j.err

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline}"
PROJECT_ROOT="${PROJECT_ROOT:-/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder}"
CONDA_BASE="${CONDA_BASE:-/home/elcrespo/miniconda3}"
CONDA_ENV="${CONDA_ENV:-imaris-careamics}"
FULL_RUN_DIR="${FULL_RUN_DIR:-}"
OUTPUT_DIR_8BIT="${OUTPUT_DIR_8BIT:-}"
LOWER_PERCENTILE="${LOWER_PERCENTILE:-0.5}"
UPPER_PERCENTILE="${UPPER_PERCENTILE:-99.8}"
PERCENTILE_Z_STEP="${PERCENTILE_Z_STEP:-2}"
PERCENTILE_YX_STEP="${PERCENTILE_YX_STEP:-2}"
OVERWRITE="${OVERWRITE:-true}"

if [[ -z "${FULL_RUN_DIR}" ]]; then
  echo "ERROR: set FULL_RUN_DIR=/path/to/sample/full_run before submitting." >&2
  exit 2
fi

mkdir -p "${PROJECT_ROOT}/logs" "${FULL_RUN_DIR}/logs"
exec > >(tee -a "${FULL_RUN_DIR}/logs/convert_denoised_to_8bit_${SLURM_JOB_ID:-manual}.log") 2>&1

echo "Imaris 8-bit display conversion started: $(date)"
echo "Host: $(hostname)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-}"
echo "REPO_ROOT=${REPO_ROOT}"
echo "FULL_RUN_DIR=${FULL_RUN_DIR}"
echo "OUTPUT_DIR_8BIT=${OUTPUT_DIR_8BIT}"
echo "CONDA_ENV=${CONDA_ENV}"
echo "LOWER_PERCENTILE=${LOWER_PERCENTILE}"
echo "UPPER_PERCENTILE=${UPPER_PERCENTILE}"
echo "PERCENTILE_Z_STEP=${PERCENTILE_Z_STEP}"
echo "PERCENTILE_YX_STEP=${PERCENTILE_YX_STEP}"

source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"
python --version
which python

overwrite_arg=()
if [[ "${OVERWRITE}" == "true" ]]; then
  overwrite_arg+=(--overwrite)
fi

output_arg=()
if [[ -n "${OUTPUT_DIR_8BIT}" ]]; then
  output_arg+=(--output-dir "${OUTPUT_DIR_8BIT}")
fi

python "${REPO_ROOT}/imaging/imaris_careamics/convert_denoised_to_8bit.py" \
  --input-dir "${FULL_RUN_DIR}" \
  "${output_arg[@]}" \
  --lower-percentile "${LOWER_PERCENTILE}" \
  --upper-percentile "${UPPER_PERCENTILE}" \
  --percentile-z-step "${PERCENTILE_Z_STEP}" \
  --percentile-yx-step "${PERCENTILE_YX_STEP}" \
  "${overwrite_arg[@]}"

echo "Imaris 8-bit display conversion finished: $(date)"
find "${FULL_RUN_DIR}/display_8bit" -maxdepth 1 -type f | sort

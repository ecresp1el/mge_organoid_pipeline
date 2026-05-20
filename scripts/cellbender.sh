#!/usr/bin/env bash
#SBATCH --job-name=cellbender
#SBATCH --account=parent0
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --output=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/cellbender-%A_%a.out
#SBATCH --error=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/cellbender-%A_%a.err

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  sbatch --array=0-5 scripts/cellbender.sh
  sbatch scripts/cellbender.sh --input /path/to/sample.h5ad

Environment alternatives:
  INPUT_DIR   Directory of input .h5/.h5ad files for Slurm array mode.
              Default: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/raw_adata
  INPUT_H5    Single input file if --input is not supplied.

Outputs:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/clean_adata/<input_basename>_cellbender_denoised.<input_extension>

Notes:
  The output and log directories are created only if absent. Existing files are
  not deleted. CellBender is run from the Great Lakes Bioinformatics module.
EOF
}

project_root="/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"
input_dir="${INPUT_DIR:-${project_root}/raw_adata}"
output_dir="${project_root}/clean_adata"
input_h5="${INPUT_H5:-}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input|-i)
      input_h5="${2:-}"
      shift 2
      ;;
    --input-dir)
      input_dir="${2:-}"
      shift 2
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ -z "${input_h5}" ]]; then
  if [[ ! -d "${input_dir}" ]]; then
    echo "Error: input directory does not exist: ${input_dir}" >&2
    exit 1
  fi

  mapfile -t input_files < <(find "${input_dir}" -maxdepth 1 -type f \( -name '*.h5' -o -name '*.h5ad' -o -name '*.hdf5' \) | sort)

  if [[ "${#input_files[@]}" -eq 0 ]]; then
    echo "Error: no .h5/.h5ad/.hdf5 files found in ${input_dir}" >&2
    exit 1
  fi

  if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "Error: no --input was supplied and SLURM_ARRAY_TASK_ID is unset." >&2
    echo "Submit directory mode with: sbatch --array=0-$(( ${#input_files[@]} - 1 )) scripts/cellbender.sh" >&2
    exit 2
  fi

  if (( SLURM_ARRAY_TASK_ID < 0 || SLURM_ARRAY_TASK_ID >= ${#input_files[@]} )); then
    echo "Error: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} is outside 0-$(( ${#input_files[@]} - 1 ))" >&2
    exit 1
  fi

  input_h5="${input_files[${SLURM_ARRAY_TASK_ID}]}"
fi

if [[ ! -f "${input_h5}" ]]; then
  echo "Error: input file does not exist: ${input_h5}" >&2
  exit 1
fi

input_filename="$(basename "${input_h5}")"
case "${input_filename}" in
  *.h5ad)
    sample_name="${input_filename%.h5ad}"
    output_extension="h5ad"
    ;;
  *.hdf5)
    sample_name="${input_filename%.hdf5}"
    output_extension="hdf5"
    ;;
  *.h5)
    sample_name="${input_filename%.h5}"
    output_extension="h5"
    ;;
  *)
    echo "Error: unsupported input extension: ${input_filename}" >&2
    exit 1
    ;;
esac

mkdir -p "${output_dir}" "${project_root}/logs"

module purge
module load Bioinformatics
module load cellbender/0.3.0

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Input: ${input_h5}"
echo "Output directory: ${output_dir}"
echo "CellBender: $(command -v cellbender)"
echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES:-unset}"
echo "SLURM_JOB_GPUS: ${SLURM_JOB_GPUS:-unset}"

if ! command -v nvidia-smi >/dev/null 2>&1; then
  echo "Error: nvidia-smi is not available. This job is not on a GPU-ready node." >&2
  exit 1
fi

nvidia-smi

output_h5="${output_dir}/${sample_name}_cellbender_denoised.${output_extension}"

if [[ -e "${output_h5}" ]]; then
  echo "Error: output already exists and will not be overwritten: ${output_h5}" >&2
  exit 1
fi

cellbender remove-background \
  --input "${input_h5}" \
  --output "${output_h5}" \
  --total-droplets-included 500000 \
  --cuda

echo "Job finished: $(date)"
echo "Wrote: ${output_h5}"

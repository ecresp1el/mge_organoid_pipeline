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
  sbatch --array=1,5 scripts/cellbender.sh
  sbatch scripts/cellbender.sh --input /path/to/sample.h5ad

Environment alternatives:
  INPUT_DIR   Directory of input .h5/.h5ad files for Slurm array mode.
              Default: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/raw_adata
  INPUT_H5    Single input file if --input is not supplied.
  CELLBENDER_WORK_DIR
              Base directory for per-sample CellBender working directories.
              Default: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/cellbender_work
  INCLUDE_DUPLICATE_BASENAMES
              Set to 1 to include files with "__" in the basename.
              Default: 0 (skip duplicate-suffixed files)
  TOTAL_DROPLETS_INCLUDED
              Optional override for CellBender --total-droplets-included.
              Default: unset (let CellBender use its internal default)

Inputs:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/raw_adata/*.h5ad

Outputs and runtime files:
  Slurm logs:
    /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/cellbender-%A_%a.out
    /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/cellbender-%A_%a.err
  CellBender outputs:
    /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/clean_adata/<input_basename>_cellbender_denoised.h5
    plus CellBender sidecar files such as *_filtered.h5, *_metrics.csv,
    *_cell_barcodes.csv, *_report.html, *.pdf, and *.log.
  CellBender working files and checkpoints:
    /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/cellbender_work/<input_basename>/

Nothing from the Slurm job should be written to the GitHub checkout. The job
changes into the per-sample NFS work directory before running CellBender, so
files such as ckpt.tar.gz and tmp.report.* are kept out of the repository.

Notes:
  The output, log, and work directories are created only if absent. Existing
  primary .h5 outputs are not overwritten; samples with an existing primary
  output are treated as already complete and skipped. CellBender is run from
  the Great Lakes Bioinformatics module.
EOF
}

project_root="/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"
input_dir="${INPUT_DIR:-${project_root}/raw_adata}"
output_dir="${project_root}/clean_adata"
logs_dir="${project_root}/logs"
work_base_dir="${CELLBENDER_WORK_DIR:-${project_root}/cellbender_work}"
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

  if [[ "${INCLUDE_DUPLICATE_BASENAMES:-0}" != "1" ]]; then
    mapfile -t input_files < <(printf '%s\n' "${input_files[@]}" | grep -v '/[^/]*__[^/]*\.[^.]*$' || true)
  fi

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

input_h5="$(realpath "${input_h5}")"
input_filename="$(basename "${input_h5}")"
case "${input_filename}" in
  *.h5ad)
    sample_name="${input_filename%.h5ad}"
    ;;
  *.hdf5)
    sample_name="${input_filename%.hdf5}"
    ;;
  *.h5)
    sample_name="${input_filename%.h5}"
    ;;
  *)
    echo "Error: unsupported input extension: ${input_filename}" >&2
    exit 1
    ;;
esac

output_h5="${output_dir}/${sample_name}_cellbender_denoised.h5"
work_dir="${work_base_dir}/${sample_name}"
tmp_dir="${work_dir}/tmp"
total_droplets_included="${TOTAL_DROPLETS_INCLUDED:-}"

mkdir -p "${output_dir}" "${logs_dir}" "${work_dir}" "${tmp_dir}"

if [[ -e "${output_h5}" ]]; then
  echo "Output already exists; skipping without overwriting: ${output_h5}"
  exit 0
fi

cd "${work_dir}"
export TMPDIR="${tmp_dir}"

module purge
module load Bioinformatics
module load cellbender/0.3.0

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo "Submit directory: ${SLURM_SUBMIT_DIR:-unset}"
echo "Input: ${input_h5}"
echo "Output directory: ${output_dir}"
echo "Work directory: ${work_dir}"
echo "TMPDIR: ${TMPDIR}"
echo "CellBender: $(command -v cellbender)"
echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES:-unset}"
echo "SLURM_JOB_GPUS: ${SLURM_JOB_GPUS:-unset}"

if ! command -v nvidia-smi >/dev/null 2>&1; then
  echo "Error: nvidia-smi is not available. This job is not on a GPU-ready node." >&2
  exit 1
fi

nvidia-smi

cellbender_cmd=(
  cellbender remove-background
  --input "${input_h5}"
  --output "${output_h5}"
  --cuda
)

if [[ -n "${total_droplets_included}" ]]; then
  cellbender_cmd+=(--total-droplets-included "${total_droplets_included}")
fi

echo "Command: ${cellbender_cmd[*]}"
"${cellbender_cmd[@]}"

echo "Job finished: $(date)"
echo "Wrote: ${output_h5}"

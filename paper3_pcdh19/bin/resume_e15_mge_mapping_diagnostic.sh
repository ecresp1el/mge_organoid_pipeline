#!/usr/bin/env bash
set -Eeuo pipefail

[[ $# -eq 1 ]] || { echo "Usage: resume_e15_mge_mapping_diagnostic.sh RUN_ID" >&2; exit 2; }
BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
# shellcheck disable=SC1090
source "${BUNDLE_DIR}/config/greatlakes.env"
# shellcheck disable=SC1090
source "${BUNDLE_DIR}/config/e15_mge_mapping_diagnostic.env"
RUN_ID="$1"
STEP=01_e15_mge_mapping_diagnostic
if [[ "${RUN_ID}" == *[!A-Za-z0-9._-]* || "${RUN_ID}" != "${STEP}"_* ]]; then
  echo "Unsafe run ID: ${RUN_ID}" >&2
  exit 2
fi
RUN_DIR="${PAPER3_ROOT}/results/${STEP}/${RUN_ID}"
[[ -d "${RUN_DIR}" && -f "${RUN_DIR}/provenance/job_ids.tsv" ]] || { echo "Run not found" >&2; exit 2; }
PREPARE_JOB_ID="$(awk -F '\t' '$1=="prepare_query" {print $2}' "${RUN_DIR}/provenance/job_ids.tsv")"
[[ -n "${PREPARE_JOB_ID}" ]] || { echo "Missing prepare job" >&2; exit 2; }
for stage in bandler_mind mapmycells report; do
  [[ -z "$(awk -F '\t' -v stage="${stage}" '$1==stage {print $2}' "${RUN_DIR}/provenance/job_ids.tsv")" ]] || {
    echo "Stage ${stage} is already recorded; refusing duplicate submission" >&2
    exit 2
  }
done

ASSET_ROOT="${PAPER3_ROOT}/inputs/mapmycells/WMB-10X/${MAPMYCELLS_ASSET_VERSION}"
MARKERS="${ASSET_ROOT}/mouse_markers_230821.json"
STATS="${ASSET_ROOT}/precomputed_stats_ABC_revision_230821.h5"
BANDLER_MATRIX="${PAPER3_ROOT}/inputs/developing_mouse_mge/Bandler2022/source/GSM5684876_CA301_filtered_RNA_counts.RDS.gz"
for required in "${MARKERS}" "${STATS}" "${BANDLER_MATRIX}" \
  "${RUN_DIR}/code/01b_transfer_bandler_mind_labels.sbatch" \
  "${RUN_DIR}/code/01c_run_mapmycells_comparator.sbatch" \
  "${RUN_DIR}/code/01d_build_e15_mapping_diagnostic_report.sbatch"; do
  [[ -f "${required}" ]] || { echo "Missing continuation input: ${required}" >&2; exit 2; }
done

cp -p "${BASH_SOURCE[0]}" "${RUN_DIR}/code/"
cp -p "${BUNDLE_DIR}/config/e15_mge_mapping_diagnostic.env" "${RUN_DIR}/config/submitted_mapping.env"
{
  echo "continued_at=$(date --iso-8601=seconds)"
  echo "reason=initial branch submission rejected because 192G exceeded the 180G standard-node capacity"
  printf 'command='; printf '%q ' "$0" "$@"; printf '\n'
} > "${RUN_DIR}/provenance/submission_continuation.txt"

EXPORTS="ALL,PAPER3_MAPPING_RUN_DIR=${RUN_DIR},PAPER3_CELLRANGER_ROOT=${PAPER3_CELLRANGER_ROOT},PAPER3_MAPMYCELLS_PYTHON=${MAPMYCELLS_PYTHON_BIN},PAPER3_MAPPING_R_MODULE=${MAPPING_R_MODULE},PAPER3_MAPPING_COMPILER_MODULE=${MAPPING_COMPILER_MODULE},PAPER3_MAPPING_HDF5_MODULE=${MAPPING_HDF5_MODULE},PAPER3_MAPPING_R_ADDON_LIB=${MAPPING_R_ADDON_LIB},PAPER3_MIND_THRESHOLD=${MIND_UNASSIGNED_THRESHOLD},PAPER3_MAPPING_SEED=${MAPPING_RANDOM_SEED},PAPER3_BANDLER_CA301_MATRIX=${BANDLER_MATRIX},PAPER3_MAPMYCELLS_MARKERS=${MARKERS},PAPER3_MAPMYCELLS_STATS=${STATS},PAPER3_MAPMYCELLS_DROP_LEVEL=${MAPMYCELLS_DROP_LEVEL}"
MIND_JOB_ID="$(sbatch --parsable --job-name=p3-map-mind --dependency="afterok:${PREPARE_JOB_ID}" --account="${ACCOUNT}" --partition="${PARTITION}" \
  --nodes=1 --ntasks=1 --cpus-per-task="${MAPPING_MIND_CPUS}" --mem="${MAPPING_MIND_MEMORY}" --time="${MAPPING_MIND_WALLTIME}" \
  --output="${RUN_DIR}/logs/mind-%j.out" --error="${RUN_DIR}/logs/mind-%j.err" --export="${EXPORTS}" "${RUN_DIR}/code/01b_transfer_bandler_mind_labels.sbatch")"
printf 'bandler_mind\t%s\n' "${MIND_JOB_ID}" >> "${RUN_DIR}/provenance/job_ids.tsv"
MMC_JOB_ID="$(sbatch --parsable --job-name=p3-map-mmc --dependency="afterok:${PREPARE_JOB_ID}" --account="${ACCOUNT}" --partition="${PARTITION}" \
  --nodes=1 --ntasks=1 --cpus-per-task="${MAPPING_MMC_CPUS}" --mem="${MAPPING_MMC_MEMORY}" --time="${MAPPING_MMC_WALLTIME}" \
  --output="${RUN_DIR}/logs/mapmycells-%j.out" --error="${RUN_DIR}/logs/mapmycells-%j.err" --export="${EXPORTS}" "${RUN_DIR}/code/01c_run_mapmycells_comparator.sbatch")"
printf 'mapmycells\t%s\n' "${MMC_JOB_ID}" >> "${RUN_DIR}/provenance/job_ids.tsv"
REPORT_JOB_ID="$(sbatch --parsable --job-name=p3-map-report --dependency="afterok:${MIND_JOB_ID}:${MMC_JOB_ID}" --account="${ACCOUNT}" --partition="${PARTITION}" \
  --nodes=1 --ntasks=1 --cpus-per-task="${MAPPING_REPORT_CPUS}" --mem="${MAPPING_REPORT_MEMORY}" --time="${MAPPING_REPORT_WALLTIME}" \
  --output="${RUN_DIR}/logs/report-%j.out" --error="${RUN_DIR}/logs/report-%j.err" --export="${EXPORTS}" "${RUN_DIR}/code/01d_build_e15_mapping_diagnostic_report.sbatch")"
printf 'report\t%s\n' "${REPORT_JOB_ID}" >> "${RUN_DIR}/provenance/job_ids.tsv"

echo "Continued run ${RUN_ID}"
echo "Prepare dependency: ${PREPARE_JOB_ID}"
echo "Bandler/MIND: ${MIND_JOB_ID}"
echo "MapMyCells: ${MMC_JOB_ID}"
echo "Report: ${REPORT_JOB_ID}"

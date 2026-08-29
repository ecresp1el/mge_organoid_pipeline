#!/usr/bin/env bash
set -Eeuo pipefail

[[ $# -eq 1 ]] || { echo "Usage: retry_e15_mge_mind_branch.sh RUN_ID" >&2; exit 2; }
BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
# shellcheck disable=SC1090
source "${BUNDLE_DIR}/config/greatlakes.env"
# shellcheck disable=SC1090
source "${BUNDLE_DIR}/config/e15_mge_mapping_diagnostic.env"
RUN_ID="$1"
STEP=01_e15_mge_mapping_diagnostic
[[ "${RUN_ID}" != *[!A-Za-z0-9._-]* && "${RUN_ID}" == "${STEP}"_* ]] || { echo "Unsafe run ID" >&2; exit 2; }
RUN_DIR="${PAPER3_ROOT}/results/${STEP}/${RUN_ID}"
[[ -f "${RUN_DIR}/PREPARE_SUCCESS.txt" ]] || { echo "Prepare stage is not complete" >&2; exit 2; }
FAILED_JOB_ID="$(awk -F '\t' '$1 ~ /^bandler_mind/ {value=$2} END {print value}' "${RUN_DIR}/provenance/job_ids.tsv")"
MMC_JOB_ID="$(awk -F '\t' '$1=="mapmycells" {print $2}' "${RUN_DIR}/provenance/job_ids.tsv")"
[[ "$(sacct -n -X -j "${FAILED_JOB_ID}" --format=State -P | head -n 1)" == FAILED* ]] || { echo "Original MIND job is not failed" >&2; exit 2; }
RETRY_NUMBER="$(awk -F '\t' '$1 ~ /^bandler_mind/ {n++} END {print n+0}' "${RUN_DIR}/provenance/job_ids.tsv")"
MIND_STAGE=bandler_mind_retry
REPORT_STAGE=report_retry
[[ "${RETRY_NUMBER}" -eq 1 ]] || { MIND_STAGE="bandler_mind_retry${RETRY_NUMBER}"; REPORT_STAGE="report_retry${RETRY_NUMBER}"; }

cp -p "${RUN_DIR}/code/transfer_bandler_mind_labels.R" "${RUN_DIR}/code/failed_${FAILED_JOB_ID}_transfer_bandler_mind_labels.R"
cp -p "${RUN_DIR}/code/01b_transfer_bandler_mind_labels.sbatch" "${RUN_DIR}/code/failed_${FAILED_JOB_ID}_01b_transfer_bandler_mind_labels.sbatch"
if [[ -f "${RUN_DIR}/code/retry_e15_mge_mind_branch.sh" ]]; then
  cp -p "${RUN_DIR}/code/retry_e15_mge_mind_branch.sh" "${RUN_DIR}/code/failed_${FAILED_JOB_ID}_retry_e15_mge_mind_branch.sh"
fi
cp -p "${BUNDLE_DIR}/scripts/reference_mapping/transfer_bandler_mind_labels.R" "${RUN_DIR}/code/transfer_bandler_mind_labels.R"
cp -p "${BUNDLE_DIR}/slurm/01b_transfer_bandler_mind_labels.sbatch" "${RUN_DIR}/code/01b_transfer_bandler_mind_labels.sbatch"
cp -p "${BASH_SOURCE[0]}" "${RUN_DIR}/code/"
cp -p "${BUNDLE_DIR}/config/e15_mge_mapping_diagnostic.env" "${RUN_DIR}/config/submitted_mapping.env"
module purge
module load "${MAPPING_COMPILER_MODULE}" "${MAPPING_HDF5_MODULE}" Bioinformatics "${MAPPING_R_MODULE}"
R_LIBS_USER="${MAPPING_R_ADDON_LIB}" Rscript -e '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())); packages <- installed.packages(lib.loc=Sys.getenv("R_LIBS_USER")); write.table(data.frame(package=rownames(packages),version=packages[,"Version"]), file=commandArgs(TRUE)[1], sep="\t", quote=FALSE, row.names=FALSE)' "${RUN_DIR}/config/R_addon_packages.tsv"

BANDLER_MATRIX="${PAPER3_ROOT}/inputs/developing_mouse_mge/Bandler2022/source/GSM5684876_CA301_filtered_RNA_counts.RDS.gz"
EXPORTS="ALL,PAPER3_MAPPING_RUN_DIR=${RUN_DIR},PAPER3_CELLRANGER_ROOT=${PAPER3_CELLRANGER_ROOT},PAPER3_MAPMYCELLS_PYTHON=${MAPMYCELLS_PYTHON_BIN},PAPER3_MAPPING_R_MODULE=${MAPPING_R_MODULE},PAPER3_MAPPING_COMPILER_MODULE=${MAPPING_COMPILER_MODULE},PAPER3_MAPPING_HDF5_MODULE=${MAPPING_HDF5_MODULE},PAPER3_MAPPING_R_ADDON_LIB=${MAPPING_R_ADDON_LIB},PAPER3_MIND_THRESHOLD=${MIND_UNASSIGNED_THRESHOLD},PAPER3_MAPPING_SEED=${MAPPING_RANDOM_SEED},PAPER3_BANDLER_CA301_MATRIX=${BANDLER_MATRIX}"
MIND_JOB_ID="$(sbatch --parsable --job-name="p3-map-mind-r${RETRY_NUMBER}" --account="${ACCOUNT}" --partition="${PARTITION}" \
  --nodes=1 --ntasks=1 --cpus-per-task="${MAPPING_MIND_CPUS}" --mem="${MAPPING_MIND_MEMORY}" --time="${MAPPING_MIND_WALLTIME}" \
  --output="${RUN_DIR}/logs/mind-retry${RETRY_NUMBER}-%j.out" --error="${RUN_DIR}/logs/mind-retry${RETRY_NUMBER}-%j.err" --export="${EXPORTS}" "${RUN_DIR}/code/01b_transfer_bandler_mind_labels.sbatch")"
printf '%s\t%s\n' "${MIND_STAGE}" "${MIND_JOB_ID}" >> "${RUN_DIR}/provenance/job_ids.tsv"
REPORT_JOB_ID="$(sbatch --parsable --job-name="p3-map-report-r${RETRY_NUMBER}" --dependency="afterok:${MIND_JOB_ID}:${MMC_JOB_ID}" --account="${ACCOUNT}" --partition="${PARTITION}" \
  --nodes=1 --ntasks=1 --cpus-per-task="${MAPPING_REPORT_CPUS}" --mem="${MAPPING_REPORT_MEMORY}" --time="${MAPPING_REPORT_WALLTIME}" \
  --output="${RUN_DIR}/logs/report-retry${RETRY_NUMBER}-%j.out" --error="${RUN_DIR}/logs/report-retry${RETRY_NUMBER}-%j.err" --export="${EXPORTS},PAPER3_MIND_THRESHOLD=${MIND_UNASSIGNED_THRESHOLD}" "${RUN_DIR}/code/01d_build_e15_mapping_diagnostic_report.sbatch")"
printf '%s\t%s\n' "${REPORT_STAGE}" "${REPORT_JOB_ID}" >> "${RUN_DIR}/provenance/job_ids.tsv"
{
  echo "retried_at=$(date --iso-8601=seconds)"
  echo "failed_job=${FAILED_JOB_ID}"
  echo "reason=Seurat module lacked optional hdf5r; installed pinned hdf5r 1.3.12 in isolated addon library"
  echo "mind_retry_job=${MIND_JOB_ID}"
  echo "report_retry_job=${REPORT_JOB_ID}"
} > "${RUN_DIR}/provenance/mind_retry${RETRY_NUMBER}.txt"
echo "Submitted MIND retry ${MIND_JOB_ID}; replacement report ${REPORT_JOB_ID} also awaits MapMyCells ${MMC_JOB_ID}"

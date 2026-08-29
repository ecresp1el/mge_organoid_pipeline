#!/usr/bin/env bash
set -Eeuo pipefail

[[ $# -eq 1 ]] || { echo "Usage: replace_e15_mapping_report_job.sh RUN_ID" >&2; exit 2; }
BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
# shellcheck disable=SC1090
source "${BUNDLE_DIR}/config/greatlakes.env"
# shellcheck disable=SC1090
source "${BUNDLE_DIR}/config/e15_mge_mapping_diagnostic.env"
RUN_ID="$1"; STEP=01_e15_mge_mapping_diagnostic
[[ "${RUN_ID}" != *[!A-Za-z0-9._-]* && "${RUN_ID}" == "${STEP}"_* ]] || { echo "Unsafe run ID" >&2; exit 2; }
RUN_DIR="${PAPER3_ROOT}/results/${STEP}/${RUN_ID}"
[[ -f "${RUN_DIR}/provenance/job_ids.tsv" ]] || { echo "Run not found" >&2; exit 2; }
MIND_JOB_ID="$(awk -F '\t' '$1 ~ /^bandler_mind/ {value=$2} END {print value}' "${RUN_DIR}/provenance/job_ids.tsv")"
MMC_JOB_ID="$(awk -F '\t' '$1=="mapmycells" {print $2}' "${RUN_DIR}/provenance/job_ids.tsv")"
OLD_REPORT_ID="$(awk -F '\t' '$1 ~ /^report/ {value=$2} END {print value}' "${RUN_DIR}/provenance/job_ids.tsv")"
OLD_STATE="$(sacct -n -X -j "${OLD_REPORT_ID}" --format=State -P | head -n 1)"
if [[ "${OLD_STATE}" == PENDING* ]]; then scancel "${OLD_REPORT_ID}"; fi
[[ "${OLD_STATE}" == PENDING* || "${OLD_STATE}" == CANCELLED* ]] || { echo "Latest report is not safely replaceable: ${OLD_STATE}" >&2; exit 2; }

cp -p "${RUN_DIR}/code/01d_build_e15_mapping_diagnostic_report.sbatch" "${RUN_DIR}/code/replaced_${OLD_REPORT_ID}_01d_build_e15_mapping_diagnostic_report.sbatch"
cp -p "${BUNDLE_DIR}/slurm/01d_build_e15_mapping_diagnostic_report.sbatch" "${RUN_DIR}/code/01d_build_e15_mapping_diagnostic_report.sbatch"
cp -p "${BUNDLE_DIR}/scripts/reference_mapping/attach_mapmycells_to_seurat.R" "${RUN_DIR}/code/attach_mapmycells_to_seurat.R"
cp -p "${BUNDLE_DIR}/scripts/reference_mapping/build_e15_mapping_diagnostic_report.py" "${RUN_DIR}/code/build_e15_mapping_diagnostic_report.py"
cp -p "${BASH_SOURCE[0]}" "${RUN_DIR}/code/"
cp -p "${BUNDLE_DIR}/config/e15_mge_mapping_diagnostic.env" "${RUN_DIR}/config/submitted_mapping.env"
EXPORTS="ALL,PAPER3_MAPPING_RUN_DIR=${RUN_DIR},PAPER3_MAPMYCELLS_PYTHON=${MAPMYCELLS_PYTHON_BIN},PAPER3_MAPPING_R_MODULE=${MAPPING_R_MODULE},PAPER3_MIND_THRESHOLD=${MIND_UNASSIGNED_THRESHOLD}"
REPORT_JOB_ID="$(sbatch --parsable --job-name=p3-map-report-final --dependency="afterok:${MIND_JOB_ID}:${MMC_JOB_ID}" --account="${ACCOUNT}" --partition="${PARTITION}" \
  --nodes=1 --ntasks=1 --cpus-per-task="${MAPPING_REPORT_CPUS}" --mem="${MAPPING_REPORT_MEMORY}" --time="${MAPPING_REPORT_WALLTIME}" \
  --output="${RUN_DIR}/logs/report-final-%j.out" --error="${RUN_DIR}/logs/report-final-%j.err" --export="${EXPORTS}" "${RUN_DIR}/code/01d_build_e15_mapping_diagnostic_report.sbatch")"
printf 'report_final\t%s\n' "${REPORT_JOB_ID}" >> "${RUN_DIR}/provenance/job_ids.tsv"
{
  echo "replaced_at=$(date --iso-8601=seconds)"
  echo "old_report_job=${OLD_REPORT_ID}"
  echo "old_state=${OLD_STATE}"
  echo "new_report_job=${REPORT_JOB_ID}"
  echo "reason=attach all MapMyCells fields to one validated combined Seurat RDS before reporting"
} > "${RUN_DIR}/provenance/report_replacement.txt"
echo "Replaced pending/cancelled report ${OLD_REPORT_ID} with combined-object/report job ${REPORT_JOB_ID}"

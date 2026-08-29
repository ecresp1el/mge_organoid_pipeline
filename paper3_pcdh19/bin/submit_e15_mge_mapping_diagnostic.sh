#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
  cat >&2 <<'EOF'
Usage: submit_e15_mge_mapping_diagnostic.sh [--dry-run] [--replace-run RUN_ID]

Create or intentionally replace one inactive run of the first diagnostic
reference-mapping step. The DAG is prepare -> {MIND, MapMyCells} -> report.
EOF
}

ORIGINAL_ARGS=("$@")
DRY_RUN=false
REPLACE_RUN_ID=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run) DRY_RUN=true; shift ;;
    --replace-run) [[ $# -ge 2 ]] || { usage; exit 2; }; REPLACE_RUN_ID="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
GREATLAKES_CONFIG="${BUNDLE_DIR}/config/greatlakes.env"
MAPPING_CONFIG="${BUNDLE_DIR}/config/e15_mge_mapping_diagnostic.env"
REQUIREMENTS="${BUNDLE_DIR}/config/e15_mge_mapping_diagnostic.requirements.txt"
SAMPLE_KEY="${BUNDLE_DIR}/config/sample_key.csv"
PREPARE_SCRIPT="${BUNDLE_DIR}/scripts/reference_mapping/prepare_e15_mge_mapping_query.py"
MIND_SCRIPT="${BUNDLE_DIR}/scripts/reference_mapping/transfer_bandler_mind_labels.R"
MMC_SCRIPT="${BUNDLE_DIR}/scripts/reference_mapping/run_mapmycells_comparator.py"
REPORT_SCRIPT="${BUNDLE_DIR}/scripts/reference_mapping/build_e15_mapping_diagnostic_report.py"
ATTACH_SCRIPT="${BUNDLE_DIR}/scripts/reference_mapping/attach_mapmycells_to_seurat.R"
PREPARE_SBATCH="${BUNDLE_DIR}/slurm/01a_prepare_e15_mge_mapping_query.sbatch"
MIND_SBATCH="${BUNDLE_DIR}/slurm/01b_transfer_bandler_mind_labels.sbatch"
MMC_SBATCH="${BUNDLE_DIR}/slurm/01c_run_mapmycells_comparator.sbatch"
REPORT_SBATCH="${BUNDLE_DIR}/slurm/01d_build_e15_mapping_diagnostic_report.sbatch"

# shellcheck disable=SC1090
source "${GREATLAKES_CONFIG}"
# shellcheck disable=SC1090
source "${MAPPING_CONFIG}"
: "${REPO_ROOT:?}" "${PAPER3_ROOT:?}" "${PAPER3_CELLRANGER_ROOT:?}"
: "${ACCOUNT:?}" "${PARTITION:?}" "${MAPMYCELLS_PYTHON_BIN:?}"
: "${MAPPING_CURATION_RUN_ID:?}" "${MIND_UNASSIGNED_THRESHOLD:?}"

for required in "${GREATLAKES_CONFIG}" "${MAPPING_CONFIG}" "${REQUIREMENTS}" "${SAMPLE_KEY}" \
  "${PREPARE_SCRIPT}" "${MIND_SCRIPT}" "${MMC_SCRIPT}" "${REPORT_SCRIPT}" "${ATTACH_SCRIPT}" \
  "${PREPARE_SBATCH}" "${MIND_SBATCH}" "${MMC_SBATCH}" "${REPORT_SBATCH}" \
  "${MAPMYCELLS_PYTHON_BIN}"; do
  [[ -f "${required}" ]] || { echo "Missing required file: ${required}" >&2; exit 2; }
done

CURATION_RUN="${PAPER3_ROOT}/results/00_developing_mouse_mge_reference_curation/${MAPPING_CURATION_RUN_ID}"
BANDLER_JOIN="${CURATION_RUN}/Bandler2022/interactive_atlas/barcode_recovery/metadata/CA301_later_atlas_barcode_join.tsv"
BANDLER_MATRIX="${PAPER3_ROOT}/inputs/developing_mouse_mge/Bandler2022/source/GSM5684876_CA301_filtered_RNA_counts.RDS.gz"
[[ -f "${CURATION_RUN}/BANDLER_E15_BARCODE_RECOVERY_SUCCESS.txt" ]] || { echo "Missing successful barcode-recovery source run" >&2; exit 2; }
[[ -f "${BANDLER_JOIN}" && -f "${BANDLER_MATRIX}" ]] || { echo "Missing Bandler mapping input" >&2; exit 2; }

STEP=01_e15_mge_mapping_diagnostic
STEP_ROOT="${PAPER3_ROOT}/results/${STEP}"
STAMP="$(date +%Y%m%d_%H%M%S)"
GIT_SHORT="$(git -C "${REPO_ROOT}" rev-parse --short HEAD)"
RUN_ID="${STEP}_${STAMP}_${GIT_SHORT}"
OUTPUT_MODE=versioned
if [[ -n "${REPLACE_RUN_ID}" ]]; then
  OUTPUT_MODE=replace
  RUN_ID="${REPLACE_RUN_ID}"
  if [[ "${RUN_ID}" == *[!A-Za-z0-9._-]* || "${RUN_ID}" != "${STEP}"_* ]]; then
    echo "Unsafe replacement run ID: ${RUN_ID}" >&2
    exit 2
  fi
fi
RUN_DIR="${STEP_ROOT}/${RUN_ID}"
if [[ "$(dirname "$(realpath -m "${RUN_DIR}")")" != "$(realpath -m "${STEP_ROOT}")" ]]; then
  echo "Run escaped its exact step root" >&2
  exit 2
fi

if [[ "${OUTPUT_MODE}" == replace && -f "${RUN_DIR}/provenance/job_ids.tsv" ]]; then
  while IFS=$'\t' read -r stage job_id; do
    [[ "${stage}" == stage || -z "${job_id}" ]] && continue
    state="$(squeue -h -j "${job_id}" -o '%T' | head -n 1)"
    [[ -z "${state}" ]] || { echo "Refusing replacement: ${job_id} is ${state}" >&2; exit 2; }
  done < "${RUN_DIR}/provenance/job_ids.tsv"
fi

ASSET_ROOT="${PAPER3_ROOT}/inputs/mapmycells/WMB-10X/${MAPMYCELLS_ASSET_VERSION}"
MARKERS="${ASSET_ROOT}/mouse_markers_230821.json"
STATS="${ASSET_ROOT}/precomputed_stats_ABC_revision_230821.h5"
echo "Step: ${STEP}"
echo "Output mode: ${OUTPUT_MODE}"
echo "Run directory: ${RUN_DIR}"
echo "Query: 12 original Cell Ranger samples; original counts, clusters, and per-sample UMAPs"
echo "Bandler/MIND: 4,481 CA301 cells; fixed Unassigned threshold ${MIND_UNASSIGNED_THRESHOLD}"
echo "MapMyCells: pinned local WMB hierarchical mapping with bootstrap outputs"
echo "DAG: prepare -> Bandler/MIND + MapMyCells in parallel -> concordance/report"

"${MAPMYCELLS_PYTHON_BIN}" -m py_compile "${PREPARE_SCRIPT}" "${MMC_SCRIPT}" "${REPORT_SCRIPT}"
"${MAPMYCELLS_PYTHON_BIN}" -m cell_type_mapper.cli.from_specified_markers --help >/dev/null
bash -n "${PREPARE_SBATCH}" "${MIND_SBATCH}" "${MMC_SBATCH}" "${REPORT_SBATCH}" "${BASH_SOURCE[0]}"
if [[ "${DRY_RUN}" == true ]]; then
  echo "Dry run passed; nothing was downloaded, copied, replaced, or submitted."
  exit 0
fi

mkdir -p "${STEP_ROOT}" "${ASSET_ROOT}"
download_asset() {
  local url="$1" target="$2" sidecar="${2}.sha256" actual temporary
  if [[ -f "${target}" && -f "${sidecar}" ]]; then
    actual="$(sha256sum "${target}" | awk '{print $1}')"
    [[ "${actual}" == "$(awk '{print $1}' "${sidecar}")" ]] || { echo "Cached asset checksum failed: ${target}" >&2; exit 2; }
    return
  fi
  [[ ! -e "${target}" ]] || { echo "Unverified cached asset exists: ${target}" >&2; exit 2; }
  temporary="$(mktemp "${ASSET_ROOT}/.download.XXXXXX")"
  trap 'rm -f "${temporary}"' RETURN
  curl --fail --location --retry 4 --output "${temporary}" "${url}"
  [[ -s "${temporary}" ]] || { echo "Empty asset download: ${url}" >&2; exit 2; }
  mv "${temporary}" "${target}"
  sha256sum "${target}" > "${sidecar}"
  trap - RETURN
}
download_asset "${MAPMYCELLS_MARKERS_URL}" "${MARKERS}"
download_asset "${MAPMYCELLS_STATS_URL}" "${STATS}"

if [[ "${OUTPUT_MODE}" == versioned ]]; then
  [[ ! -e "${RUN_DIR}" ]] || { echo "Versioned run already exists" >&2; exit 2; }
elif [[ -d "${RUN_DIR}" ]]; then
  echo "Intentionally replacing inactive run within ${STEP}: ${RUN_DIR}"
  find "${RUN_DIR}" -mindepth 1 -depth -delete
elif [[ -e "${RUN_DIR}" ]]; then
  echo "Replacement target is not a directory" >&2
  exit 2
fi
mkdir -p "${RUN_DIR}/code" "${RUN_DIR}/config" "${RUN_DIR}/inputs" "${RUN_DIR}/query" \
  "${RUN_DIR}/Bandler_MIND" "${RUN_DIR}/MapMyCells" "${RUN_DIR}/tables" \
  "${RUN_DIR}/figures" "${RUN_DIR}/logs" "${RUN_DIR}/provenance"
cp -p "${PREPARE_SCRIPT}" "${MIND_SCRIPT}" "${MMC_SCRIPT}" "${REPORT_SCRIPT}" "${ATTACH_SCRIPT}" \
  "${PREPARE_SBATCH}" "${MIND_SBATCH}" "${MMC_SBATCH}" "${REPORT_SBATCH}" \
  "${BASH_SOURCE[0]}" "${RUN_DIR}/code/"
cp -p "${GREATLAKES_CONFIG}" "${RUN_DIR}/config/submitted_greatlakes.env"
cp -p "${MAPPING_CONFIG}" "${RUN_DIR}/config/submitted_mapping.env"
cp -p "${REQUIREMENTS}" "${RUN_DIR}/config/declared_requirements.txt"
cp -p "${SAMPLE_KEY}" "${RUN_DIR}/config/sample_key.csv"
cp -p "${BANDLER_JOIN}" "${RUN_DIR}/inputs/CA301_later_atlas_barcode_join.tsv"
"${MAPMYCELLS_PYTHON_BIN}" -m pip freeze > "${RUN_DIR}/config/mapmycells_python_pip_freeze.txt"

cat > "${RUN_DIR}/README.md" <<EOF
# First E15 MGE mapping diagnostic

This derived run preserves the original 12 Cell Ranger count matrices,
sample-qualified graph clusters, and per-sample UMAPs. It adds Bandler/MIND
Seurat label transfer and an independent Allen WMB MapMyCells comparator.
No clustering, integration, manual annotation, or embedding is recomputed.

Exact executable copies are in \`code/\`. The small CA301 barcode/label join is
copied into \`inputs/\`; large immutable source matrices and Allen assets remain
in their recorded shared input caches. Use the submission wrapper's
\`--replace-run ${RUN_ID}\` only to intentionally rerun and overwrite this exact
inactive step package.
EOF
{
  echo "RUN_ID=${RUN_ID}"
  echo "RUN_DIR=${RUN_DIR}"
  echo "OUTPUT_MODE=${OUTPUT_MODE}"
  echo "PAPER3_CELLRANGER_ROOT=${PAPER3_CELLRANGER_ROOT}"
  echo "BANDLER_CA301_MATRIX=${BANDLER_MATRIX}"
  echo "BANDLER_JOIN_SOURCE=${BANDLER_JOIN}"
  echo "MAPMYCELLS_MARKERS=${MARKERS}"
  echo "MAPMYCELLS_STATS=${STATS}"
  echo "MAPMYCELLS_PYTHON=${MAPMYCELLS_PYTHON_BIN}"
  echo "MIND_UNASSIGNED_THRESHOLD=${MIND_UNASSIGNED_THRESHOLD}"
  echo "MAPPING_RANDOM_SEED=${MAPPING_RANDOM_SEED}"
} > "${RUN_DIR}/config/resolved.env"
{
  echo "submitted_at=$(date --iso-8601=seconds)"
  echo "submit_host=$(hostname)"
  printf 'command='; printf '%q ' "$0" "${ORIGINAL_ARGS[@]}"; printf '\n'
} > "${RUN_DIR}/provenance/submission.txt"
{
  printf 'asset\tpath\tbytes\tsha256\n'
  printf 'markers\t%s\t%s\t%s\n' "${MARKERS}" "$(stat -c '%s' "${MARKERS}")" "$(sha256sum "${MARKERS}" | awk '{print $1}')"
  printf 'precomputed_stats\t%s\t%s\t%s\n' "${STATS}" "$(stat -c '%s' "${STATS}")" "$(sha256sum "${STATS}" | awk '{print $1}')"
  printf 'bandler_CA301_matrix\t%s\t%s\t%s\n' "${BANDLER_MATRIX}" "$(stat -c '%s' "${BANDLER_MATRIX}")" "$(sha256sum "${BANDLER_MATRIX}" | awk '{print $1}')"
  printf 'bandler_CA301_join_copy\t%s\t%s\t%s\n' "${RUN_DIR}/inputs/CA301_later_atlas_barcode_join.tsv" "$(stat -c '%s' "${RUN_DIR}/inputs/CA301_later_atlas_barcode_join.tsv")" "$(sha256sum "${RUN_DIR}/inputs/CA301_later_atlas_barcode_join.tsv" | awk '{print $1}')"
} > "${RUN_DIR}/provenance/input_manifest.tsv"

EXPORTS="ALL,PAPER3_MAPPING_RUN_DIR=${RUN_DIR},PAPER3_CELLRANGER_ROOT=${PAPER3_CELLRANGER_ROOT},PAPER3_MAPMYCELLS_PYTHON=${MAPMYCELLS_PYTHON_BIN},PAPER3_MAPPING_R_MODULE=${MAPPING_R_MODULE},PAPER3_MAPPING_COMPILER_MODULE=${MAPPING_COMPILER_MODULE},PAPER3_MAPPING_HDF5_MODULE=${MAPPING_HDF5_MODULE},PAPER3_MAPPING_R_ADDON_LIB=${MAPPING_R_ADDON_LIB},PAPER3_MIND_THRESHOLD=${MIND_UNASSIGNED_THRESHOLD},PAPER3_MAPPING_SEED=${MAPPING_RANDOM_SEED},PAPER3_BANDLER_CA301_MATRIX=${BANDLER_MATRIX},PAPER3_MAPMYCELLS_MARKERS=${MARKERS},PAPER3_MAPMYCELLS_STATS=${STATS},PAPER3_MAPMYCELLS_DROP_LEVEL=${MAPMYCELLS_DROP_LEVEL}"
printf 'stage\tjob_id\n' > "${RUN_DIR}/provenance/job_ids.tsv"
PREPARE_JOB_ID="$(sbatch --parsable --job-name=p3-map-prep --account="${ACCOUNT}" --partition="${PARTITION}" \
  --nodes=1 --ntasks=1 --cpus-per-task="${MAPPING_PREPARE_CPUS}" --mem="${MAPPING_PREPARE_MEMORY}" --time="${MAPPING_PREPARE_WALLTIME}" \
  --output="${RUN_DIR}/logs/prepare-%j.out" --error="${RUN_DIR}/logs/prepare-%j.err" --export="${EXPORTS}" "${RUN_DIR}/code/01a_prepare_e15_mge_mapping_query.sbatch")"
printf 'prepare_query\t%s\n' "${PREPARE_JOB_ID}" >> "${RUN_DIR}/provenance/job_ids.tsv"
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

echo "Submitted prepare job ${PREPARE_JOB_ID}"
echo "Submitted Bandler/MIND job ${MIND_JOB_ID} and MapMyCells job ${MMC_JOB_ID} in parallel"
echo "Submitted concordance/report job ${REPORT_JOB_ID}"
echo "Run directory: ${RUN_DIR}"

#!/usr/bin/env bash

set -Eeuo pipefail

usage() {
  echo "Usage: $0 [--dry-run]" >&2
}

DRY_RUN=false
ORIGINAL_ARGS=("$@")
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run)
      DRY_RUN=true
      shift
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 2
      ;;
  esac
done

BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
GREATLAKES_CONFIG="${BUNDLE_DIR}/config/greatlakes.env"
RECOVERY_CONFIG="${BUNDLE_DIR}/config/legacy_gene_id_recovery.env"
PY_SCRIPT="${BUNDLE_DIR}/scripts/02b_recover_legacy_gene_ids.py"
SBATCH_TEMPLATE="${BUNDLE_DIR}/slurm/02b_recover_legacy_gene_ids.sbatch"

# shellcheck disable=SC1090
source "${GREATLAKES_CONFIG}"
# shellcheck disable=SC1090
source "${RECOVERY_CONFIG}"
: "${REPO_ROOT:?REPO_ROOT is required}"
: "${PAPER2_ROOT:?PAPER2_ROOT is required}"
: "${ACCOUNT:?ACCOUNT is required}"
: "${PARTITION:?PARTITION is required}"
: "${CONDA_ENV_BIN:?CONDA_ENV_BIN is required}"

STEP02_DIR="${PAPER2_ROOT}/results/02_harmonize_genes"
STEP02_RUN_DIR="$(find "${STEP02_DIR}" -mindepth 1 -maxdepth 1 -type d -name '02_harmonize_genes_*' | sort | tail -n 1)"
[[ -n "${STEP02_RUN_DIR}" && -s "${STEP02_RUN_DIR}/SUCCESS.txt" ]] || {
  echo "A completed Step 02 package is required." >&2
  exit 2
}
grep -q '^review_stop=YES$' "${STEP02_RUN_DIR}/SUCCESS.txt" || {
  echo "Step 02 package does not retain its review stop." >&2
  exit 2
}

for required in \
  "${RECOVERY_CONFIG}" "${PY_SCRIPT}" "${SBATCH_TEMPLATE}" \
  "${STEP02_RUN_DIR}/tables/feature_mapping_long.tsv.gz" \
  "${STEP02_RUN_DIR}/reference/gencode.v50.annotation.gtf.gz" \
  "${VARELA_DIV30_FEATURES}" "${VARELA_DIV90_FEATURES}" \
  "${LEGACY_GENCODE32_GTF}" "${LEGACY_GENCODE35_GTF}" "${LEGACY_GENCODE44_GTF}" \
  "${WALSH_GEO_SOFT}" "${BERSHTEYN_2025_GEO_SOFT}" \
  "${BERSHTEYN_2023_GEO_SOFT}" "${SIEBERT_2026_METADATA}"; do
  [[ -s "${required}" ]] || { echo "Missing required evidence file: ${required}" >&2; exit 2; }
done

echo "Step: 02b_legacy_gene_id_recovery (report-only review extension)"
echo "Parent Step 02: ${STEP02_RUN_DIR}"
echo "Confirmed source feature tables: Varela DIV30 and DIV90"
echo "Historical reference comparison: GENCODE 27, 32, 35, 44, and 50"
echo "No canonical object or Step 02 table will be modified."
if [[ "${DRY_RUN}" == "true" ]]; then
  echo "Dry run passed; nothing was downloaded, created, or submitted."
  exit 0
fi

"${BUNDLE_DIR}/bin/initialize_turbo.sh" >/dev/null
STAMP="$(date +%Y%m%d_%H%M%S)"
GIT_SHORT="$(git -C "${REPO_ROOT}" rev-parse --short HEAD)"
RUN_ID="02b_legacy_gene_id_recovery_${STAMP}_${GIT_SHORT}"
STEP_DIR="${PAPER2_ROOT}/results/02b_legacy_gene_id_recovery"
RUN_DIR="${STEP_DIR}/${RUN_ID}"
JOB_FILE="${PAPER2_ROOT}/jobs/${RUN_ID}.sbatch"
[[ ! -e "${RUN_DIR}" ]] || { echo "Run already exists: ${RUN_DIR}" >&2; exit 2; }

mkdir -p \
  "${RUN_DIR}/code" \
  "${RUN_DIR}/config" \
  "${RUN_DIR}/reference" \
  "${RUN_DIR}/tables" \
  "${RUN_DIR}/logs" \
  "${RUN_DIR}/provenance"

cp -p "${PY_SCRIPT}" "${SBATCH_TEMPLATE}" "${RUN_DIR}/code/"
cp -p "${GREATLAKES_CONFIG}" "${RUN_DIR}/config/submitted_greatlakes.env"
cp -p "${RECOVERY_CONFIG}" "${RUN_DIR}/config/legacy_gene_id_recovery.env"
cp -p "${SBATCH_TEMPLATE}" "${JOB_FILE}"

echo "Downloading frozen GENCODE 27 primary-assembly annotation..."
curl --fail --location --retry 3 --output \
  "${RUN_DIR}/reference/gencode.v27.primary_assembly.annotation.gtf.gz" \
  "${LEGACY_GENCODE27_URL}"
printf '%s  %s\n' "${LEGACY_GENCODE27_MD5}" \
  "${RUN_DIR}/reference/gencode.v27.primary_assembly.annotation.gtf.gz" | md5sum --check --status

{
  echo -e "evidence_id\tevidence_type\tstudy_id\treference_release\tpath\tsize_bytes\tsha256\tclaim"
  add_evidence() {
    local evidence_id="$1" evidence_type="$2" study_id="$3" release="$4" path="$5" claim="$6"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${evidence_id}" "${evidence_type}" "${study_id}" "${release}" "${path}" \
      "$(stat -c %s "${path}")" "$(sha256sum "${path}" | cut -d ' ' -f 1)" "${claim}"
  }
  add_evidence step02_mapping parent_mapping all 50 \
    "${STEP02_RUN_DIR}/tables/feature_mapping_long.tsv.gz" \
    "completed Step 02 feature mapping; read only"
  add_evidence varela_div30_features source_feature_table varela_div30 original \
    "${VARELA_DIV30_FEATURES}" "confirmed project Cell Ranger feature table"
  add_evidence varela_div90_features source_feature_table varela_div90 original \
    "${VARELA_DIV90_FEATURES}" "confirmed project Cell Ranger feature table"
  add_evidence gencode27 historical_annotation all 27 \
    "${RUN_DIR}/reference/gencode.v27.primary_assembly.annotation.gtf.gz" \
    "historical identity reference; not proof of study reference bundle"
  add_evidence gencode32 historical_annotation all 32 \
    "${LEGACY_GENCODE32_GTF}" "historical identity reference; not proof of study reference bundle"
  add_evidence gencode35 historical_annotation all 35 \
    "${LEGACY_GENCODE35_GTF}" "historical identity reference; not proof of study reference bundle"
  add_evidence gencode44 historical_annotation all 44 \
    "${LEGACY_GENCODE44_GTF}" "historical identity reference; not proof of study reference bundle"
  add_evidence gencode50 current_annotation all 50 \
    "${STEP02_RUN_DIR}/reference/gencode.v50.annotation.gtf.gz" "current mapping identity reference"
  add_evidence walsh_geo_soft provenance_metadata walsh NA \
    "${WALSH_GEO_SOFT}" "Cell Ranger/assembly provenance; exact reference bundle not specified"
  add_evidence bershteyn2025_geo_soft provenance_metadata bershteyn_2025 NA \
    "${BERSHTEYN_2025_GEO_SOFT}" "Cell Ranger/assembly provenance; exact reference bundle not specified"
  add_evidence bershteyn2023_geo_soft provenance_metadata bershteyn_2023 NA \
    "${BERSHTEYN_2023_GEO_SOFT}" "Cell Ranger/assembly provenance; exact reference bundle not specified"
  add_evidence siebert2026_metadata provenance_metadata siebert_2026 NA \
    "${SIEBERT_2026_METADATA}" "local NeMO metadata; exact feature table/reference bundle absent"
} > "${RUN_DIR}/provenance/source_evidence_registry.tsv"

{
  echo "REPO_ROOT=${REPO_ROOT}"
  echo "PAPER2_ROOT=${PAPER2_ROOT}"
  echo "RUN_ID=${RUN_ID}"
  echo "RUN_DIR=${RUN_DIR}"
  echo "STEP02_RUN_DIR=${STEP02_RUN_DIR}"
  echo "CONDA_ENV_BIN=${CONDA_ENV_BIN}"
} > "${RUN_DIR}/config/resolved.env"

{
  echo "submitted_at=$(date --iso-8601=seconds)"
  echo "submit_host=$(hostname)"
  echo "submit_user=${USER:-unknown}"
  printf 'command='
  printf '%q ' "$0" "${ORIGINAL_ARGS[@]}"
  printf '\n'
} > "${RUN_DIR}/provenance/submission.txt"
git -C "${REPO_ROOT}" rev-parse HEAD > "${RUN_DIR}/provenance/git_commit.txt"
git -C "${REPO_ROOT}" status --short > "${RUN_DIR}/provenance/git_status.txt"
git -C "${REPO_ROOT}" diff --binary > "${RUN_DIR}/provenance/git_working_tree.patch"

JOB_ID="$(sbatch --parsable \
  --job-name="p2-gene-recover" \
  --account="${ACCOUNT}" \
  --partition="${PARTITION}" \
  --nodes=1 \
  --ntasks=1 \
  --cpus-per-task="${LEGACY_RECOVERY_CPUS}" \
  --mem="${LEGACY_RECOVERY_MEMORY}" \
  --time="${LEGACY_RECOVERY_WALLTIME}" \
  --output="${RUN_DIR}/logs/slurm-%j.out" \
  --error="${RUN_DIR}/logs/slurm-%j.err" \
  --export="ALL,PAPER2_LEGACY_RUN_DIR=${RUN_DIR},CONDA_ENV_BIN=${CONDA_ENV_BIN}" \
  "${JOB_FILE}")"
printf '%s\n' "${JOB_ID}" > "${RUN_DIR}/provenance/slurm_job_id.txt"
cp -p "${RUN_DIR}/provenance/slurm_job_id.txt" "${RUN_DIR}/provenance/job_id.txt"

echo "Submitted job ${JOB_ID}"
echo "Run directory: ${RUN_DIR}"
echo "Monitor: ${BUNDLE_DIR}/bin/status.sh ${RUN_DIR}"

#!/usr/bin/env bash
# PURPOSE
#   Freeze and submit approved Step 06 unintegrated diagnostics.
set -Eeuo pipefail

# usage: write CLI syntax to stderr; no submission or filesystem mutation.
usage() {
  cat >&2 <<'EOF'
Usage: submit_primary_processing_step_06.sh [--dry-run] [--replace-run RUN_ID]

The default creates a new versioned run. The job consumes approved Step 02
directly, uses no Step 03-05 output, and stops IN_REVIEW without integration.
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

# Resolve editable repository assets. See STEP06_CODE_AND_TUNING_GUIDE.md for their roles.
BUNDLE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
GREATLAKES_CONFIG="${BUNDLE_DIR}/config/greatlakes.env"
STEP_CONFIG="${BUNDLE_DIR}/config/primary_processing_step06.env"
REQUIREMENTS="${BUNDLE_DIR}/config/primary_processing_step06.requirements.txt"
PYTHON_PACKAGE="${BUNDLE_DIR}/scripts/primary_processing"
SBATCH_SOURCE="${BUNDLE_DIR}/slurm/primary_processing_06_diagnostics.sbatch"
PACKAGE_README="${BUNDLE_DIR}/templates/PRIMARY_PROCESSING_STEP06_OUTPUT_PACKAGE_README.md"
HANDOFF="${BUNDLE_DIR}/PCDH19_PRIMARY_PROCESSING_HANDOFF.md"
BYPASS_DECISION="${BUNDLE_DIR}/PCDH19_PRIMARY_PROCESSING_STEP04_STEP05_BYPASS_DECISION.md"
SPECIFICATION="${BUNDLE_DIR}/STEP06_DIAGNOSTIC_SPECIFICATION.txt"
SAMPLE_KEY="${BUNDLE_DIR}/config/sample_key.csv"
TECHNICAL_MANIFEST="${BUNDLE_DIR}/config/sample_manifest_draft.tsv"

# shellcheck disable=SC1090
source "${GREATLAKES_CONFIG}"
# shellcheck disable=SC1090
source "${STEP_CONFIG}"
: "${REPO_ROOT:?}" "${PAPER3_ROOT:?}" "${ACCOUNT:?}"
: "${PRIMARY_PROCESSING_PYTHON_BIN:?}" "${PRIMARY_PROCESSING_STEP06_STEP02_RUN_ID:?}"

STEP02_RUN_DIR="${PAPER3_ROOT}/results/primary_processing/02_qc_filtering/${PRIMARY_PROCESSING_STEP06_STEP02_RUN_ID}"
INPUT_H5AD="${STEP02_RUN_DIR}/${PRIMARY_PROCESSING_STEP06_STEP02_CHECKPOINT}"
STEP02_STATUS="${STEP02_RUN_DIR}/STEP_STATUS.tsv"
STEP02_MANIFEST="${STEP02_RUN_DIR}/tables/output_manifest.tsv"
APPROVAL_LEDGER="${PAPER3_ROOT}/results/primary_processing/APPROVAL_LEDGER.tsv"
for required in "${GREATLAKES_CONFIG}" "${STEP_CONFIG}" "${REQUIREMENTS}" "${SBATCH_SOURCE}" "${PACKAGE_README}" "${HANDOFF}" "${BYPASS_DECISION}" "${SPECIFICATION}" "${SAMPLE_KEY}" "${TECHNICAL_MANIFEST}" "${INPUT_H5AD}" "${STEP02_STATUS}" "${STEP02_MANIFEST}" "${APPROVAL_LEDGER}" "${PRIMARY_PROCESSING_PYTHON_BIN}"; do
  [[ -f "${required}" ]] || { echo "Missing required Step 06 asset: ${required}" >&2; exit 2; }
done
for module in step06_cli.py step06_models.py step06_analysis.py step06_metrics.py step06_plots.py step06_progress.py step06_progress_cli.py step06_publishing.py step06_validation.py step06_workflow.py; do
  [[ -f "${PYTHON_PACKAGE}/${module}" ]] || { echo "Missing Step 06 module: ${module}" >&2; exit 2; }
done

# Verify approval/bypass lineage before creating a run. Identity checks are not tuning knobs.
awk -F '\t' -v run="${PRIMARY_PROCESSING_STEP06_STEP02_RUN_ID}" 'NR==1 {for(i=1;i<=NF;i++) h[$i]=i; next} $(h["run_id"])==run && $(h["status"])=="APPROVED" && $(h["approved_run_id"])==run {ok++} END {exit(ok==1?0:1)}' "${APPROVAL_LEDGER}" || { echo "Exact Step 02 run is not uniquely APPROVED" >&2; exit 2; }
for pair in "03_scdblfinder:REJECTED" "04_ambient_rna_contamination_assessment:SKIPPED" "05_broad_biological_contaminant_assessment:SKIPPED"; do
  ledger_step="${pair%%:*}"
  ledger_status="${pair##*:}"
  awk -F '\t' -v step="${ledger_step}" -v status="${ledger_status}" 'NR==1 {for(i=1;i<=NF;i++) h[$i]=i; next} $(h["step"])==step && $(h["status"])==status {ok++} END {exit(ok==1?0:1)}' "${APPROVAL_LEDGER}" || { echo "Missing unique bypass ledger record: ${pair}" >&2; exit 2; }
done
[[ "$(stat -c '%s' "${INPUT_H5AD}")" == "${PRIMARY_PROCESSING_STEP06_INPUT_BYTES}" ]] || { echo "Step 02 H5AD size mismatch" >&2; exit 2; }

# Verify installed versions against exact pins; this does not install or update packages.
"${PRIMARY_PROCESSING_PYTHON_BIN}" - "${REQUIREMENTS}" <<'PY'
from importlib import metadata
from pathlib import Path
import sys
expected = {}
for line in Path(sys.argv[1]).read_text().splitlines():
    if line and not line.startswith("#"):
        name, version = line.split("==", 1)
        expected[name] = version
observed = {name: metadata.version(name) for name in expected}
if observed != expected:
    raise SystemExit(f"Step 06 environment mismatch: {observed!r}")
PY
PYTHONPATH="${BUNDLE_DIR}/scripts" "${PRIMARY_PROCESSING_PYTHON_BIN}" -m compileall -q "${PYTHON_PACKAGE}"
PYTHONPATH="${BUNDLE_DIR}/scripts" "${PRIMARY_PROCESSING_PYTHON_BIN}" -m primary_processing.step06_cli --help >/dev/null
PYTHONPATH="${BUNDLE_DIR}/scripts" "${PRIMARY_PROCESSING_PYTHON_BIN}" -m primary_processing.step06_progress_cli --help >/dev/null

STEP="06_technical_sample_batch_diagnostics"
WORKFLOW_ROOT="${PAPER3_ROOT}/results/primary_processing"
STEP_ROOT="${WORKFLOW_ROOT}/${STEP}"
LOG_ROOT="${PAPER3_ROOT}/logs/primary_processing"
STAMP="$(date +%Y%m%d_%H%M%S)"
GIT_SHORT="$(git -C "${REPO_ROOT}" rev-parse --short HEAD)"
OUTPUT_MODE=versioned
RUN_ID="${STEP}_${STAMP}_${GIT_SHORT}"
if [[ -n "${REPLACE_RUN_ID}" ]]; then
  OUTPUT_MODE=replace
  RUN_ID="${REPLACE_RUN_ID}"
  [[ "${RUN_ID}" != *[!A-Za-z0-9._-]* && "${RUN_ID}" == "${STEP}"_* ]] || { echo "Unsafe replacement run ID" >&2; exit 2; }
fi
RUN_DIR="${STEP_ROOT}/${RUN_ID}"
[[ "$(dirname "$(realpath -m "${RUN_DIR}")")" == "$(realpath -m "${STEP_ROOT}")" ]] || { echo "Run escaped Step 06 root" >&2; exit 2; }

echo "Step: ${STEP}"
echo "Approved input: ${PRIMARY_PROCESSING_STEP06_STEP02_RUN_ID}"
echo "Lineage: Step 03 rejected; Steps 04-05 skipped"
echo "Output: ${RUN_DIR}"
echo "Resources: ${PRIMARY_PROCESSING_STEP06_CPUS} CPU, ${PRIMARY_PROCESSING_STEP06_MEMORY}, ${PRIMARY_PROCESSING_STEP06_WALLTIME}"
echo "Boundary: unintegrated A-L diagnostics; zero removals"
# Dry-run exits before run creation/submission, but earlier compileall may write Python bytecode caches.
if [[ "${DRY_RUN}" == true ]]; then
  echo "Dry run passed; no directory was created and no job was submitted."
  exit 0
fi

# Actual submission requires a clean commit so the frozen run has an auditable source identity.
[[ -z "$(git -C "${REPO_ROOT}" status --porcelain)" ]] || { echo "Refusing submission from a dirty repository; commit exact Step 06 assets first" >&2; exit 2; }
"${BUNDLE_DIR}/bin/initialize_turbo.sh" >/dev/null
mkdir -p "${STEP_ROOT}" "${LOG_ROOT}"
JOB_FILE="${PAPER3_ROOT}/jobs/${RUN_ID}.sbatch"
if [[ "${OUTPUT_MODE}" == versioned ]]; then
  [[ ! -e "${RUN_DIR}" && ! -e "${JOB_FILE}" ]] || { echo "Refusing existing versioned target" >&2; exit 2; }
elif [[ -d "${RUN_DIR}" ]]; then
  while IFS= read -r prior_job; do
    [[ -n "${prior_job}" ]] || continue
    [[ -z "$(squeue -h -j "${prior_job}" -o '%T' | head -n 1)" ]] || { echo "Refusing active replacement" >&2; exit 2; }
  done < <(awk -F '\t' 'NR>1 {print $2}' "${RUN_DIR}/provenance/job_ids.tsv" 2>/dev/null || true)
  # Explicit replacement is destructive within this named inactive run; default versioning avoids it.
  find "${RUN_DIR}" -mindepth 1 -depth -delete
elif [[ -e "${RUN_DIR}" ]]; then
  echo "Replacement target is not a directory" >&2
  exit 2
fi

mkdir -p "${RUN_DIR}/code/primary_processing" "${RUN_DIR}/config" "${RUN_DIR}/logs" "${RUN_DIR}/provenance"
# Freeze source/config/lineage now. Jobs use these snapshots, not later repository edits.
cp -p "${PYTHON_PACKAGE}"/*.py "${RUN_DIR}/code/primary_processing/"
cp -p "${BASH_SOURCE[0]}" "${SBATCH_SOURCE}" "${RUN_DIR}/code/"
cp -p "${GREATLAKES_CONFIG}" "${RUN_DIR}/config/submitted_greatlakes.env"
cp -p "${STEP_CONFIG}" "${RUN_DIR}/config/submitted_step06.env"
cp -p "${REQUIREMENTS}" "${RUN_DIR}/config/requirements.txt"
cp -p "${HANDOFF}" "${RUN_DIR}/config/authoritative_handoff.md"
cp -p "${BYPASS_DECISION}" "${RUN_DIR}/config/step04_step05_bypass_decision.md"
cp -p "${SPECIFICATION}" "${RUN_DIR}/config/STEP06_DIAGNOSTIC_SPECIFICATION.txt"
cp -p "${SAMPLE_KEY}" "${RUN_DIR}/config/sample_key.csv"
cp -p "${TECHNICAL_MANIFEST}" "${RUN_DIR}/config/technical_manifest.tsv"
cp -p "${STEP02_STATUS}" "${RUN_DIR}/config/input_step02_STEP_STATUS.tsv"
cp -p "${STEP02_MANIFEST}" "${RUN_DIR}/config/input_step02_output_manifest.tsv"
cp -p "${APPROVAL_LEDGER}" "${RUN_DIR}/config/approval_ledger_at_submission.tsv"
cp -p "${PACKAGE_README}" "${RUN_DIR}/README.md"
cp -p "${RUN_DIR}/code/$(basename "${SBATCH_SOURCE}")" "${JOB_FILE}"

{
  printf 'REPO_ROOT=%q\n' "${REPO_ROOT}"
  printf 'PAPER3_ROOT=%q\n' "${PAPER3_ROOT}"
  printf 'INPUT_H5AD=%q\n' "${INPUT_H5AD}"
  echo "STEP=${STEP}"
  echo "RUN_ID=${RUN_ID}"
  printf 'RUN_DIR=%q\n' "${RUN_DIR}"
  printf 'WORKFLOW_ROOT=%q\n' "${WORKFLOW_ROOT}"
  echo "OUTPUT_MODE=${OUTPUT_MODE}"
  echo "PRIMARY_PROCESSING_PYTHON_BIN=${PRIMARY_PROCESSING_PYTHON_BIN}"
  grep '^PRIMARY_PROCESSING_STEP06_' "${STEP_CONFIG}"
} > "${RUN_DIR}/config/resolved.env"
{
  echo "submitted_at=$(date --iso-8601=seconds)"
  echo "submit_host=$(hostname)"
  echo "submit_user=${USER:-unknown}"
  printf 'command='; printf '%q ' "$0" "${ORIGINAL_ARGS[@]}"; printf '\n'
} > "${RUN_DIR}/provenance/submission.txt"
{
  echo "git_commit=$(git -C "${REPO_ROOT}" rev-parse HEAD)"
  echo "git_short=${GIT_SHORT}"
  echo "git_status_begin"; git -C "${REPO_ROOT}" status --short; echo "git_status_end"
} > "${RUN_DIR}/provenance/repository_state.txt"
"${PRIMARY_PROCESSING_PYTHON_BIN}" -m pip freeze > "${RUN_DIR}/config/python_pip_freeze.txt"

# Only this final section calls sbatch. Scheduler resources come from the Step 06 env file.
EXPORTS="ALL,PAPER3_PRIMARY_RUN_DIR=${RUN_DIR},PAPER3_PRIMARY_WORKFLOW_ROOT=${WORKFLOW_ROOT},PAPER3_PRIMARY_PYTHON_BIN=${PRIMARY_PROCESSING_PYTHON_BIN}"
JOB_ID="$(sbatch --parsable --job-name=pcdh19-primary-06 --account="${ACCOUNT}" --partition="${PRIMARY_PROCESSING_STEP06_PARTITION}" --nodes=1 --ntasks=1 --cpus-per-task="${PRIMARY_PROCESSING_STEP06_CPUS}" --mem="${PRIMARY_PROCESSING_STEP06_MEMORY}" --time="${PRIMARY_PROCESSING_STEP06_WALLTIME}" --output="${LOG_ROOT}/${RUN_ID}-%j.out" --error="${LOG_ROOT}/${RUN_ID}-%j.err" --export="${EXPORTS}" "${JOB_FILE}")"
printf 'stage\tjob_id\nstep06\t%s\n' "${JOB_ID}" > "${RUN_DIR}/provenance/job_ids.tsv"
ln -s "${LOG_ROOT}/${RUN_ID}-${JOB_ID}.out" "${RUN_DIR}/logs/scheduler.out"
ln -s "${LOG_ROOT}/${RUN_ID}-${JOB_ID}.err" "${RUN_DIR}/logs/scheduler.err"
echo "Submitted Step 06 job ${JOB_ID}"
echo "Run ID: ${RUN_ID}"
echo "Run directory: ${RUN_DIR}"
echo "Successful computation remains IN_REVIEW."

#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  python_notebooks/scripts/start_vscode_compute_job.sh [options]

Purpose:
  Request an interactive Slurm allocation for VS Code + notebook work.
  Keep the terminal running while VS Code is connected to the allocated node.

Options:
  --account <name>      Slurm account. Default: parent0
  --partition <name>    Slurm partition. Default: standard
  --mem <amount>        Memory per node. Default: 128GB
  --cpus <int>          CPUs for the notebook/kernel. Default: 4
  --time <HH:MM:SS>     Walltime. Default: 04:00:00
  --job-name <name>     Job name. Default: mge-py-notebook
  --help                Show this help.

Examples:
  python_notebooks/scripts/start_vscode_compute_job.sh
  python_notebooks/scripts/start_vscode_compute_job.sh --account louisdan0 --mem 160GB --time 06:00:00
  python_notebooks/scripts/start_vscode_compute_job.sh --partition largemem --mem 512GB --time 04:00:00
EOF
}

account="${SLURM_ACCOUNT:-parent0}"
partition="${SLURM_PARTITION:-standard}"
mem="${SLURM_MEM:-128GB}"
cpus="${SLURM_CPUS_PER_TASK:-4}"
time_limit="${SLURM_TIME:-04:00:00}"
job_name="${SLURM_JOB_NAME:-mge-py-notebook}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --account)
      account="${2:-}"
      shift 2
      ;;
    --partition)
      partition="${2:-}"
      shift 2
      ;;
    --mem)
      mem="${2:-}"
      shift 2
      ;;
    --cpus)
      cpus="${2:-}"
      shift 2
      ;;
    --time)
      time_limit="${2:-}"
      shift 2
      ;;
    --job-name)
      job_name="${2:-}"
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

cat <<EOF
Requesting interactive VS Code notebook allocation:
  account:   ${account}
  partition: ${partition}
  mem:       ${mem}
  cpus:      ${cpus}
  time:      ${time_limit}
  job name:  ${job_name}

Keep this terminal open while you use the compute node.
After Slurm grants the allocation, note the node name in the salloc output.
EOF

exec salloc \
  --account="${account}" \
  --partition="${partition}" \
  --nodes=1 \
  --ntasks-per-node=1 \
  --cpus-per-task="${cpus}" \
  --mem="${mem}" \
  --time="${time_limit}" \
  --job-name="${job_name}"

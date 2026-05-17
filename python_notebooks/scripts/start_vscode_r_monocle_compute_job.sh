#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  python_notebooks/scripts/start_vscode_r_monocle_compute_job.sh [options]

Purpose:
  Request an interactive Slurm allocation for VS Code + Monocle3 R work.
  Keep this terminal running while VS Code is connected to the allocated node.

Options:
  --account <name>      Slurm account. Default: parent0
  --partition <name>    Slurm partition. Default: standard
  --mem <amount>        Memory per node. Default: 160GB
  --cpus <int>          CPUs for the R session. Default: 8
  --time <HH:MM:SS>     Walltime. Default: 08:00:00
  --job-name <name>     Job name. Default: mge-r-monocle
  --help                Show this help.

Examples:
  python_notebooks/scripts/start_vscode_r_monocle_compute_job.sh
  python_notebooks/scripts/start_vscode_r_monocle_compute_job.sh --mem 240GB --time 12:00:00
EOF
}

account="${SLURM_ACCOUNT:-parent0}"
partition="${SLURM_PARTITION:-standard}"
mem="${SLURM_MEM:-160GB}"
cpus="${SLURM_CPUS_PER_TASK:-8}"
time_limit="${SLURM_TIME:-08:00:00}"
job_name="${SLURM_JOB_NAME:-mge-r-monocle}"

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
Requesting interactive VS Code Monocle3/R allocation:
  account:   ${account}
  partition: ${partition}
  mem:       ${mem}
  cpus:      ${cpus}
  time:      ${time_limit}
  job name:  ${job_name}

Keep this terminal open while you use the compute node.
After Slurm grants the allocation:
  1. Note the compute node name in the salloc output.
  2. Connect VS Code to that compute node using your Great Lakes SSH setup.
  3. In the compute-node VS Code terminal, run:
       cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
       export PROJECT_ROOT=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
       module load Bioinformatics
       module load Rmonocle3/1.3.7
       R
  4. Open and run section-by-section:
       exploration/09_mgeo_rgc_ipc_monocle3_interactive.R

This is an interactive R script workflow in VS Code. It uses # %% sections like
a notebook, but avoids needing a Jupyter R kernel/IRkernel.
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

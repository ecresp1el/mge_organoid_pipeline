#!/usr/bin/env bash
#SBATCH --job-name=gse94641_validate
#SBATCH --account=parent0
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/logs/%x-%j.out

set -Eeuo pipefail
REPO_ROOT=/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
"${REPO_ROOT}/paper3_pcdh19/bin/run_step_00_gse94641_reference_validation.sh"

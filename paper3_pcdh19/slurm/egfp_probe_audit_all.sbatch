#!/usr/bin/env bash
# PURPOSE
#   Great Lakes batch wrapper for Step 02c. It only requests resources and
#   invokes bin/run_egfp_probe_audit_all.sh; all scientific logic remains in
#   the version-controlled runner and Python implementation.
#
# SCHEDULER CONTRACT
#   Account parent0, standard partition, one task using two CPUs, 8 GiB RAM,
#   and a four-hour wall-time limit. SLURM expands %x/%j in the scheduler log.
#
# INPUTS AND OUTPUTS
#   The fixed checkout supplies code/configuration. The runner reads the locked
#   Step 02b, panel, sample-key, and Cell Ranger inputs and publishes Step 02c
#   assets below PAPER3_ROOT/results/egfp_probe_audit/. It also creates its own
#   timestamped run log. Scheduler stdout/stderr goes to the declared
#   PAPER3_ROOT/logs/%x-%j.out path, whose parent must exist at submission.
#
# FAILURE/RESTART AND SCIENTIFIC SCOPE
#   Strict shell mode propagates failures to SLURM. Reruns reuse only validated
#   identical outputs. No normalization, cell calling, cell typing, genotype
#   inference, or reporter-positive classification occurs. See
#   ../PIPELINE_IO_AND_BIOLOGICAL_SCOPE.md for the full scientific contract.
#SBATCH --job-name=egfp_probe_audit
#SBATCH --account=parent0
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/logs/%x-%j.out

set -Eeuo pipefail

REPO_ROOT=/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
"${REPO_ROOT}/paper3_pcdh19/bin/run_egfp_probe_audit_all.sh"

# Overnight URD Run Check

Date submitted: 2026-06-12

Two separate checkpointed/full URD jobs were submitted so they can run overnight and be checked tomorrow morning.

## Jobs

| Dataset | Job ID | Run label | Status at last check | Node | Purpose |
|---|---:|---|---|---|---|
| DIV30 | 51691451 | `div30_first_urd_paper_radial_glia_30k_checkpoint_v2` | Running | `gl3307` | Checkpointed 30k rerun with new stage/resource logs before any all-cell DIV30 run |
| DIV90 | 51691452 | `div90_urd_jia_lineage_full_v4_glia_tips` | Running | `gl3112` | Full retained-cell DIV90 run using frozen v4 Jia-lineage/glia-tip logic |

## Morning Check

Checked: 2026-06-12 09:18 EDT

| Dataset | Job ID | Final state | Exit code | Elapsed | MaxRSS | Notes |
|---|---:|---|---|---:|---:|---|
| DIV30 | 51691451 | COMPLETED | 0:0 | 06:37:08 | 71304172K | Wrote `div30_first_urd_object.rds`, stage/resource logs, pseudotime table, plots, and lineage decision report. Pseudotime is finite for 30000/30000 cells. |
| DIV90 | 51691452 | FAILED | 1:0 | 02:13:13 | 27271968K | PCA/diffusion and flood pseudotime completed, then `write_tables` failed because all exported pseudotime values were `NA` for 20049/20049 cells. |

DIV90 failed in:

```text
tables/div30_first_urd_stage_timings.tsv
write_tables	failed	2026-06-12 02:25:50	2026-06-12 02:25:52	error=no rows to aggregate
```

The relevant log tail is:

```text
[2026-06-12 02:25:46] Processing flood pseudotime max.frac.NA=0.4 stability.div=10
[2026-06-12 02:25:47] STAGE_END flood_pseudotime elapsed_seconds=882.82
[2026-06-12 02:25:50] STAGE_START write_tables
[2026-06-12 02:25:52] STAGE_FAIL write_tables elapsed_seconds=1.36 error=no rows to aggregate
Error in aggregate.data.frame(lhs, mf[-1L], FUN = FUN, ...) :
  no rows to aggregate
```

DIV90 root/export checks:

```text
n_selected_cells = 20049
n_root_cells = 8
root_cluster = 12
root_top_percent = 2.0
root_min_cells = 8
non_na_urd_pseudotime = 0 / 20049
```

Recommended DIV90 next step: do not continue to tree building from this run. First rerun a smaller/cheaper DIV90 diagnostic with a more permissive root definition and/or flood processing setting, then confirm `non_na_urd_pseudotime > 0` before launching the full tree workflow. The current failure is upstream of the summary table in pseudotime/root processing; the `no rows to aggregate` error is only the first downstream symptom.

## Root Convention Going Forward

Use a consistent top-percent-within-designated-root-pool rule, not a top-percent over all cells.

The current intended diagnostic root rule is:

```text
1. Define the dataset-specific root pool.
2. Rank cells inside that pool by the dataset's Jia/progenitor RootScore.
3. Select the top 10% within that pool as URD roots.
4. Record the root pool column/value, pool size, selected percent, selected root cells, RootScore table, and program/marker comparison table in the run output.
```

Current pool definitions and expected top-10% root counts:

| Dataset | Root pool | Pool size | Selected percent | Expected root cells |
|---|---|---:|---:|---:|
| DIV30 | `paper_cluster_annotation == "Radial glia"` | 12246 | 10% | 1225 |
| DIV90 | `cluster_id_numeric == 12` | 358 | 10% | 36 |

This intentionally replaces the inconsistent previous pairing:

```text
DIV30 previous root: all Radial glia cells = 12246 / 30000
DIV90 failed root: top 2% of cluster 12 = 8 / 20049
```

Implementation note: `scripts/19_div30_jia_rootscore_candidates.R` now accepts `--pool-col` and `--pool-value` so DIV30 score-defined root candidates can be selected as top X% within Radial glia while preserving the same output columns (`jia_rootscore_selected_root`, `jia_rootscore_selected_top_percent`, and top-percent candidate flags).

The required root comparison table is:

```text
root_score_program_marker_summary.tsv
```

It reports selected roots versus the full designated root pool, all scored/retained cells, and cluster/pool-level comparison rows. Required mean columns include:

```text
mean_jia_score_RGC1
mean_jia_score_RGC2
mean_jia_score_IPC
mean_logupx_HES1
mean_logupx_VIM
mean_logupx_NES
mean_logupx_DLX1
mean_logupx_DLX2
mean_logupx_ASCL1
```

DIV30 implementation:

```text
scripts/19_div30_jia_rootscore_candidates.R
slurm_templates/35_div30_jia_rootscore_root10_reflood.sbatch.template
```

DIV90 implementation:

```text
python_notebooks/scripts/export_div90_jia_lineage_urd_inputs.py
```

Planned top-10% rerun labels:

```text
DIV30: div30_first_urd_paper_radial_glia_30k_checkpoint_v2/jia_rootscore_root10_radial_glia_pool_v1
DIV90: div90_urd_jia_lineage_full_v4_glia_tips_root10_v1
```

## Top-10% Root Reruns Submitted

Submitted: 2026-06-12 09:54 EDT

| Dataset | Job ID | Run/root label | Status at submit check | Purpose |
|---|---:|---|---|---|
| DIV30 | 51699299 | `jia_rootscore_root10_radial_glia_pool_v1` under `div30_first_urd_paper_radial_glia_30k_checkpoint_v2` | PENDING `(Priority)` | Score Radial glia pool, select top 10% RootScore roots, reflood existing 30k URD object, write decision report. |
| DIV90 | 51699300 | `div90_urd_jia_lineage_full_v4_glia_tips_root10_v1` | FAILED `1:0` after 00:00:43 | Exporter completed with 36 root cells and wrote `root_score_program_marker_summary.tsv`; R step failed immediately with `Error: unexpected end of input` while repo files were being merged locally. Superseded by job 51699477. |
| DIV90 | 51699477 | `div90_urd_jia_lineage_full_v4_glia_tips_root10_v1` | PENDING `(None)` | Active resubmission of full retained DIV90 top-10% cluster-12 root run from stable merged checkout. |

Logs:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/35_div30_jia_rootscore_root10_reflood_51699299.log
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/34_div90_jia_lineage_urd_smoke_51699300.log
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/34_div90_jia_lineage_urd_smoke_51699477.log
```

Output roots:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_30k_checkpoint_v2/jia_rootscore_root10_radial_glia_pool_v1/
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_urd_jia_lineage_full_v4_glia_tips_root10_v1/
```

Check:

```bash
squeue -j 51699299,51699477
sacct -j 51699299,51699300,51699477 --format=JobID,JobName%32,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS,NodeList
```

## Logs

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/31_div30_first_urd_51691451.log
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/34_div90_jia_lineage_urd_smoke_51691452.log
```

## Output Roots

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_30k_checkpoint_v2/
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_urd_jia_lineage_full_v4_glia_tips/
```

## First Things To Check

```bash
squeue -j 51691451,51691452
sacct -j 51691451,51691452 --format=JobID,JobName%32,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS,NodeList
```

Check these files if the jobs are still running or if a run stops unexpectedly:

```text
tables/slurm_resource_heartbeat.tsv
tables/div30_first_urd_stage_timings.tsv
lineage_tree_jia_endpoint_tips_v1/tables/lineage_tree_stage_timings.tsv
lineage_tree_jia_endpoint_tips_v1/tables/slurm_resource_heartbeat.tsv
```

## Interpretation

DIV90 full is the intended production-scale v4 run.

DIV30 is intentionally not all-cells yet. It is a checkpointed 30k rerun because the previous 30k job consumed substantial resources but did not leave the expected final `div30_first_urd_object.rds`. If this checkpointed 30k run completes cleanly and writes the expected stage/resource logs, then DIV30 can be promoted to all-cells next.

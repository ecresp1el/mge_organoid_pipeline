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

Live check: 2026-06-12 09:57 EDT

| Dataset | Job ID | State | Elapsed at check | Node | Note |
|---|---:|---|---:|---|---|
| DIV30 | 51699299 | RUNNING | 00:02:41 | gl3252 | `root_score_program_marker_summary.tsv` already written under the top-10% root run `tables/` directory; candidate scoring completed/underway before reflood/report stages. |
| DIV90 | 51699477 | RUNNING | 00:01:09 | gl3262 | Resubmission started successfully and reached the R URD build step; exporter already wrote `inputs/root_score_program_marker_summary.tsv`. |

Live check: 2026-06-12 11:53 EDT

| Dataset | Job ID | State | Current stage signal | Output audit status |
|---|---:|---|---|---|
| DIV30 | 51699299 | RUNNING | Reflood in progress; log reached flood 16/20 from 1225 root cells. | `root_score` stage complete in manifest: 16/16 expected artifacts present. `reflood` and reflood `lineage_decision_report` artifacts not present yet because job is still running. |
| DIV90 | 51699477 | RUNNING | Initial URD build in `calc_pca_and_diffusion_map`; exporter/input stage complete. | `input_export` stage complete in manifest: 9/9 expected artifacts present. Initial URD/post-URD plots and tree/marker reports are pending until the R stages complete. |

## Output Reproducibility Audit

Added:

```text
python_notebooks/scripts/audit_urd_run_outputs.py
```

Purpose: write a stable manifest for each URD run so expected files, generated plots, missing artifacts, sizes, and modification times are explicit instead of inferred from directory browsing.

Manifest outputs:

```text
tables/urd_output_manifest.tsv
tables/urd_output_manifest_summary.tsv
tables/urd_output_discovered_files.tsv
```

The manifest script is now called at the end of these templates:

```text
slurm_templates/34_div90_jia_lineage_urd_smoke.sbatch.template
slurm_templates/35_div30_jia_rootscore_root10_reflood.sbatch.template
```

Important: jobs `51699299` and `51699477` were submitted before the template audit hooks were added. The manual manifests below capture the current in-progress state; rerun the audit commands after the jobs finish to refresh them as final completion manifests.

Refresh commands after completion:

```bash
source /home/elcrespo/miniconda3/etc/profile.d/conda.sh
conda activate mge-organoid-python

python python_notebooks/scripts/audit_urd_run_outputs.py \
  --run-root /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_30k_checkpoint_v2/jia_rootscore_root10_radial_glia_pool_v1 \
  --mode div30_root10 \
  --selected-pct 10 \
  --outdir /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_30k_checkpoint_v2/jia_rootscore_root10_radial_glia_pool_v1/tables

python python_notebooks/scripts/audit_urd_run_outputs.py \
  --run-root /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_urd_jia_lineage_full_v4_glia_tips_root10_v1 \
  --mode div90_jia \
  --outdir /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_urd_jia_lineage_full_v4_glia_tips_root10_v1/tables
```

Manual manifests were also generated for the active runs:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_30k_checkpoint_v2/jia_rootscore_root10_radial_glia_pool_v1/tables/urd_output_manifest.tsv
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_urd_jia_lineage_full_v4_glia_tips_root10_v1/tables/urd_output_manifest.tsv
```

Current manifest summaries while jobs are still running:

| Dataset | Manifest mode | Expected | Present | Missing required | Interpretation |
|---|---|---:|---:|---:|---|
| DIV30 | `div30_root10` | 35 | 16 | 19 | RootScore layer exists; reflood object and reflood decision-report plots/tables are pending. |
| DIV90 | `div90_jia` | 77 | 12 | 64 | Input/export layer and Slurm heartbeat exist; initial URD, decision report, lineage-tree, Jia marker-validation, and candidate-marker outputs are pending. |

## Plot Refresh: Annotation Grids And Marker Overlays

Submitted: 2026-06-12 14:34 EDT

| Job ID | Job name | Purpose |
|---:|---|---|
| 51718989 | `urd_plot_refresh` | Regenerate post-hoc plots/reports from completed saved RDS objects after plotting-code updates. Does not rerun PCA, diffusion maps, floods, random walks, or `buildTree()`. |
| 51742318 | `urd_plot_refresh` | Corrected replot after visual QC: larger aspect-aware annotation-grid panels, and multi-marker grids with visible per-gene autoscaled logUPX colorbars. |

Template:

```text
slurm_templates/36_refresh_urd_posthoc_plots.sbatch.template
```

Log:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/36_refresh_urd_posthoc_plots_51718989.log
```

New annotation-grid contract:

```text
lineage_decision_report/plots/diffusion_map_annotation.png
lineage_decision_report/plots/diffusion_map_annotation_grid.png
```

The first file remains the original all-groups-colored diffusion map. The new `*_grid.png` file facets by annotation group; each panel shows all cells in grey and highlights only one group in red.

```text
lineage_tree_jia_endpoint_tips_v1/plots/urd_tree_annotation.png
lineage_tree_jia_endpoint_tips_v1/plots/urd_tree_annotation_grid.png
lineage_tree_cluster_number_name_v1/plots/urd_tree_annotation.png
lineage_tree_cluster_number_name_v1/plots/urd_tree_annotation_grid.png
```

The first tree file remains the original all-groups-colored tree. The new `*_grid.png` file facets by annotation group; each panel shows all tree branches/cells in grey and highlights only one group in red.

New individual marker-intensity overlay contract:

```text
jia_fig_s11_style_marker_validation_v1/plots/jia_fig_s11_marker_tree_overlay_<GENE>.png
jia_fig_s11_style_marker_validation_v1/plots/jia_fig_s11_marker_tree_overlay_<GENE>.pdf
candidate_pv_marker_projection_v1/plots/div90_candidate_marker_tree_overlay_<GENE>.png
candidate_pv_marker_projection_v1/plots/div90_candidate_marker_tree_overlay_<GENE>.pdf
```

Each individual marker overlay shows the full URD tree/branches with marker intensity as logUPX and a visible colorbar. This is separate from the compact multi-gene panel files, which are useful for overview but harder to audit quantitatively.

Correction after visual QC: the first annotation-grid implementation produced panels that were too compressed, especially for `urd_tree_annotation_grid.png`. The plotting code now uses aspect-aware dimensions and fewer columns for tree grids. Multi-gene marker grids now keep a visible logUPX colorbar inside each gene panel, and each gene panel still autoscales to that gene's own max expression, matching the individual per-gene files.

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

## 2026-06-14 Final All-Cell URD Rerun Plan

This section supersedes the older note above about waiting before promoting DIV30 to all-cells. The final reruns below are intentionally configured with `MAX_CELLS=0`, which means export all available/retained cells rather than a stratified cap.

Code and logging changes made before submission:

```text
local commit 7723daf Log planned URD outputs before final runs
local commit 28da3f6 Record final all-cell URD job submissions

python_notebooks/scripts/audit_urd_run_outputs.py
  Adds --planned-only and --artifact-types so Slurm logs can print a plot/report plan before heavy work starts.
  DIV30 first-pass expected artifacts now match the actual template: initial URD + lineage decision report, not a missing tree stage.

slurm_templates/31_div30_first_urd.sbatch.template
  Prints planned DIV30 first-pass plots/reports at startup.
  Runs the output audit at completion.

slurm_templates/35_div30_jia_rootscore_root10_reflood.sbatch.template
  Prints planned RootScore/reflood plots/reports at startup.

slurm_templates/34_div90_jia_lineage_urd_smoke.sbatch.template
  Prints planned DIV90 plots/reports at startup.
```

GitHub push status:

```text
Push to origin/main is currently blocked by the local HTTPS credential prompt:
fatal: could not read Username for 'https://github.com': No such device or address
```

The startup log block is bracketed by:

```text
=== URD EXPECTED ARTIFACT PLAN BEGIN ===
=== URD EXPECTED ARTIFACT PLAN END ===
```

Dry-run planned manifests were generated before submission:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_allcells_final_v1/tables/urd_expected_artifacts.tsv
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_allcells_final_v1/jia_rootscore_root10_radial_glia_pool_allcells_final_v1/tables/urd_expected_artifacts.tsv
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_urd_jia_lineage_allcells_root10_final_v1/tables/urd_expected_artifacts.tsv
```

Planned plot/report counts from those manifests:

```text
DIV30 first-pass all-cell URD: 11 plot/report artifacts
DIV30 Jia RootScore top10 reflood: 14 plot/report artifacts
DIV90 all-cell Jia lineage/root10 URD: 62 plot/report artifacts
```

### Final Root Definitions

DIV30 is a two-step chain:

```text
Step 1, all-cell geometry:
  MAX_CELLS=0
  initial URD root label = paper_cluster_annotation == "Radial glia"

Step 2, final reflood root:
  root pool = paper_cluster_annotation == "Radial glia"
  root score = Jia/progenitor RootScore
  selected root candidates = top 10% by RootScore within the Radial glia pool
  count rule = ceiling(n_radial_glia_pool * 10 / 100)
```

DIV90 final all-cell run:

```text
MAX_CELLS=0
retained clusters = 0,1,2,3,4,5,8,9,10,11,12
root pool = cluster_id_numeric == 12
root score = div90_jia_rootscore
selected root candidates = top 10% by RootScore within cluster 12
count rule = ceiling(n_cluster12_pool * 10 / 100)
DIV90_ROOT_MIN_CELLS=1, so the old 8-cell smoke-test floor cannot override the 10% definition
```

Root/program summary tables to inspect after completion:

```text
DIV30:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_allcells_final_v1/jia_rootscore_root10_radial_glia_pool_allcells_final_v1/tables/root_score_program_marker_summary.tsv

DIV90:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_urd_jia_lineage_allcells_root10_final_v1/inputs/root_score_program_marker_summary.tsv
```

These tables compare selected roots against the root pool and all cells, including mean RGC1, mean RGC2, mean IPC, mean HES1, mean VIM, mean NES, mean DLX1, mean DLX2, mean ASCL1, RootScore, and group percentages/metadata.

### Plot Contract For Final Runs

The final reruns inherit the corrected plot formatting from `a74f54f`:

```text
diffusion_map_annotation.png
  all annotation groups colored together

diffusion_map_annotation_grid.png
  every facet shows all cells in grey and overlays one group in red
  aspect-aware dimensions keep panels from being squished

urd_tree_annotation.png
  all tree groups colored together

urd_tree_annotation_grid.png
  every facet shows all tree branches/cells in grey and overlays one group in red
  tree grids use wider panels/fewer columns to keep the tree readable

marker tree overlays
  individual per-gene PNG/PDF outputs include visible logUPX colorbars
  multi-gene grids keep per-panel colorbars and per-gene autoscaling
```

Validation before submission:

```text
python -m py_compile python_notebooks/scripts/audit_urd_run_outputs.py: OK
bash -n Slurm templates 31/34/35: OK
R parse scripts 15/17/25/27: OK
git diff --check: OK
```

### Submission Commands

Use `160G` on `standard`; `200G` exceeds the practical memory available on these nodes.

```bash
RUN_LABEL="div30_first_urd_paper_radial_glia_allcells_final_v1" \
ROOT_LABEL="Radial glia" \
MAX_CELLS=0 \
HEARTBEAT_INTERVAL_SECONDS=300 \
sbatch --parsable --export=ALL --job-name=div30_all_urd1 --time=72:00:00 --cpus-per-task=8 --mem=160G \
  slurm_templates/31_div30_first_urd.sbatch.template
```

```bash
RUN_LABEL="div30_first_urd_paper_radial_glia_allcells_final_v1" \
ROOT_RUN_LABEL="jia_rootscore_root10_radial_glia_pool_allcells_final_v1" \
SELECTED_PCT=10 \
POOL_COL="paper_cluster_annotation" \
POOL_VALUE="Radial glia" \
PSEUDOTIME_NAME="jia_rootscore_top10pct_radial_glia_root" \
sbatch --parsable --export=ALL --dependency=afterok:<DIV30_FIRST_JOB_ID> --job-name=div30_all_root10 --time=48:00:00 --cpus-per-task=4 --mem=160G \
  slurm_templates/35_div30_jia_rootscore_root10_reflood.sbatch.template
```

```bash
RUN_LABEL="div90_urd_jia_lineage_allcells_root10_final_v1" \
MAX_CELLS=0 \
DIV90_ROOT_TOP_PERCENT=10 \
DIV90_ROOT_MIN_CELLS=1 \
HEARTBEAT_INTERVAL_SECONDS=300 \
sbatch --parsable --export=ALL --job-name=div90_all_root10 --time=72:00:00 --cpus-per-task=8 --mem=160G \
  slurm_templates/34_div90_jia_lineage_urd_smoke.sbatch.template
```

### Submitted Jobs

Submitted:

```text
2026-06-14 12:24 EDT
```

```text
51772269 div30_all_urd1
  State at submission check: PENDING
  Reason: ReqNodeNotAvail,_Reserved_for_maintenance
  Dependency: none
  Time/mem/cpus: 72:00:00, 160G, 8 CPUs
  Log: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/31_div30_first_urd_51772269.log
  Output root: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_allcells_final_v1/

51772270 div30_all_root10
  State at submission check: PENDING
  Reason: Dependency
  Dependency: afterok:51772269
  Time/mem/cpus: 48:00:00, 160G, 4 CPUs
  Log: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/35_div30_jia_rootscore_root10_reflood_51772270.log
  Output root: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_allcells_final_v1/jia_rootscore_root10_radial_glia_pool_allcells_final_v1/

51772271 div90_all_root10
  State at submission check: PENDING
  Reason: ReqNodeNotAvail,_Reserved_for_maintenance
  Dependency: none
  Time/mem/cpus: 72:00:00, 160G, 8 CPUs
  Log: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/34_div90_jia_lineage_urd_smoke_51772271.log
  Output root: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_urd_jia_lineage_allcells_root10_final_v1/
```

Immediate check commands:

```bash
squeue -j 51772269,51772270,51772271 -o "%.18i %.26j %.10T %.12M %.9l %.6D %R"
sacct -j 51772269,51772270,51772271 --format=JobID,JobName%32,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS,NodeList
tail -f /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/31_div30_first_urd_51772269.log
tail -f /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/34_div90_jia_lineage_urd_smoke_51772271.log
```

The log files will appear once Slurm allocates a node and starts the scripts. The first visible content should include the Slurm context echo block and the planned plot/report artifact plan.

### Backfill Time-Limit Update

At 2026-06-14 17:03 EDT, `standard` had idle nodes but Slurm would not start the 72h/48h jobs because the `SM2026_Maintenance` reservation begins at 2026-06-15 04:00 EDT. To let the final all-cell jobs run in the available backfill window, the existing pending jobs were updated in place rather than resubmitted:

```bash
scontrol update JobId=51772269 TimeLimit=10:00:00
scontrol update JobId=51772271 TimeLimit=10:00:00
scontrol update JobId=51772270 TimeLimit=04:00:00
```

Current live status after the update:

```text
51772269 div30_all_urd1   RUNNING 10:00:00 gl3206
51772270 div30_all_root10 PENDING 04:00:00 Dependency
51772271 div90_all_root10 RUNNING 10:00:00 gl3299
```

The startup logs confirm the planned plot/report blocks printed correctly. Early all-cell export counts:

```text
DIV30 selected_cells=90631
DIV30 Radial glia root pool=36991
DIV90 retained selected_cells=20049
DIV90 cluster 12 root pool=358
DIV90 selected top10 root_cells=36
```

### Max-Backfill Replacement

At 2026-06-14 17:05 EDT, the 10h backfill attempt above was intentionally stopped so the final rerun could use the largest schedulable walltime before the 2026-06-15 04:00 EDT maintenance reservation.

Canceled superseded jobs:

```text
51772269 div30_all_urd1   CANCELLED after ~2 minutes
51772270 div30_all_root10 CANCELLED before start
51772271 div90_all_root10 CANCELLED after ~2 minutes
```

Fresh max-backfill jobs use new run labels to avoid mixing partial output from the canceled attempt:

```text
DIV30 run label:
  div30_first_urd_paper_radial_glia_allcells_backfillmax_20260614_v1

DIV30 RootScore/reflood label:
  jia_rootscore_root10_radial_glia_pool_backfillmax_20260614_v1

DIV90 run label:
  div90_urd_jia_lineage_allcells_root10_backfillmax_20260614_v1
```

Submitted max-backfill jobs:

```text
51778136 div30_max_urd1
  State after submission check: RUNNING on gl3206
  Time/mem/cpus: 10:50:00, 160G, 8 CPUs
  Log: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/31_div30_first_urd_51778136.log
  Output root: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_allcells_backfillmax_20260614_v1/

51778137 div30_max_root10
  State after submission check: PENDING
  Dependency: afterok:51778136
  Time/mem/cpus: 04:00:00, 160G, 4 CPUs
  Log: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/35_div30_jia_rootscore_root10_reflood_51778137.log
  Output root: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_allcells_backfillmax_20260614_v1/jia_rootscore_root10_radial_glia_pool_backfillmax_20260614_v1/

51778138 div90_max_root10
  State after submission check: RUNNING on gl3299
  Time/mem/cpus: 10:50:00, 160G, 8 CPUs
  Log: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/34_div90_jia_lineage_urd_smoke_51778138.log
  Output root: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_urd_jia_lineage_allcells_root10_backfillmax_20260614_v1/
```

Early log checks confirm the planned plot/report blocks printed and the all-cell exports started correctly:

```text
DIV30 selected_cells=90631
DIV30 Radial glia root pool=36991
DIV90 retained selected_cells=20049
DIV90 cluster 12 root pool=358
DIV90 selected top10 root_cells=36
```

Current check commands:

```bash
squeue -j 51778136,51778137,51778138 -o "%.18i %.26j %.10T %.12M %.9l %.6D %R"
sacct -j 51778136,51778137,51778138 --format=JobID,JobName%32,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS,NodeList
tail -f /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/31_div30_first_urd_51778136.log
tail -f /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/34_div90_jia_lineage_urd_smoke_51778138.log
```

## Resumable Production Refactor

Updated after cleanup/refactor. The previous canceled, superseded, and smoke
URD result directories were removed to reduce confusion and storage use. Kept
directories:

```text
DIV30 completed baseline:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_30k_checkpoint_v2

DIV90 completed baselines:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_urd_jia_lineage_allcells_root10_backfillmax_20260614_v1
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_urd_jia_lineage_allcells_root10_neuron_s9_7tips_v1
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_urd_jia_lineage_full_v4_glia_tips_root10_v1
```

Deleted directories:

```text
div30_first_urd_paper_radial_glia_allcells_backfillmax_20260614_v1
div30_first_urd_paper_radial_glia_allcells_final_v1
div30_first_urd_paper_radial_glia_smoke5k_knn100_v1
div30_first_urd_paper_radial_glia_v1
div90_urd_jia_lineage_allcells_root10_final_v1
div90_urd_jia_lineage_full_v4_glia_tips
div90_urd_jia_lineage_smoke5k_knn100_v1
div90_urd_jia_lineage_smoke5k_knn100_v2_crabp1_tip
div90_urd_jia_lineage_smoke5k_knn100_v3_glia_cells
div90_urd_jia_lineage_smoke5k_knn100_v4_glia_tips
```

The shared URD geometry runner `scripts/14_div30_first_urd.R` is now resumable.
It preserves the same biological/statistical parameters and only changes
execution mechanics. It writes:

```text
checkpoints/urd_after_filter.rds
checkpoints/urd_after_variable_genes.rds
checkpoints/urd_after_pca.rds
checkpoints/urd_after_diffusion_map.rds
checkpoints/urd_after_flood_pseudotime.rds
tables/div30_first_urd_checkpoint_manifest.tsv
```

If `URD_RESUME=true` and `URD_FORCE_RECOMPUTE=false`, resubmitting the same run
label resumes from the latest available checkpoint. The final plots/reports are
unchanged: the same posthoc decision report, tree reports, marker validation
figures, candidate marker projections, and output manifests are still generated
after the URD object exists.

New clear run labels for the next production reruns:

```text
DIV30 first-pass all-cell URD:
div30_allcells_radial_glia_firstpass_urd_resumable_v1

DIV30 Jia RootScore top10 reflood:
div30_allcells_radial_glia_jia_rootscore_top10_reflood_v1

DIV90 all-cell Jia root10 neuron-S9 seven-tip URD:
div90_allcells_jia_root10_neuron_s9_7tips_urd_resumable_v1
```

Parameter consistency is intentional:

```text
MAX_CELLS=0
URD_KNN=100
URD_N_FLOODS=20
URD_NUM_VARIABLE_GENES=3000
URD_PCA_MP_FACTOR=2
URD_SIGMA=local
URD_MIN_GENES=500
URD_MIN_CELLS=3
URD_MIN_COUNTS=10
URD_FLOOD_MINIMUM_CELLS=2
URD_FLOOD_MAX_FRAC_NA=0.4
```

Template defaults updated:

```text
slurm_templates/31_div30_first_urd.sbatch.template
slurm_templates/35_div30_jia_rootscore_root10_reflood.sbatch.template
slurm_templates/34_div90_jia_lineage_urd_smoke.sbatch.template
python_notebooks/scripts/audit_urd_run_outputs.py
```

The DIV90 template filename still has historical `smoke` in its name, but the
default output directory and new Slurm log prefix are production/all-cell names:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_allcells_jia_root10_neuron_s9_7tips_urd_resumable_v1
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/34_div90_jia_lineage_urd_allcells_%j.log
```

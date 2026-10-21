# DIV90 Jia-Lineage-Driven URD Plan

Date: 2026-06-11

Current status (2026-10-21): the canonical seven-tip DIV90 production tree is saved. The downstream branchpoint-gene analysis at the end of this handoff has been implemented; the user authorized execution after planning review. Technical validation precedes timestamped carbon-copy submission. Earlier smoke/planning sections below are historical, not the current run definition.

This handoff records the DIV90 URD plan after the DIV90 UMAP/metadata audit. It is intentionally aligned with the DIV30 URD workflow, but the DIV90 biological question is different because the DIV90 object has real terminal cluster labels.

## Primary Question

The DIV90 URD tree should ask:

```text
Which Jia fetal MGE lineage programs are represented by the terminal DIV90 states?
```

It should not primarily ask:

```text
How do progenitors generate SST, PV, and MGE neurons?
```

Therefore, DIV90 tips should approximate Jia lineage endpoints, not broad SST/PV/MGE paper labels.

## Inputs

Source Seurat object:

```text
/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/umap_props_output/clustered_day90_with_cluster_names_2.rds
```

Cached AnnData:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/varela_div90.h5ad
```

Exported AnnData audit inputs:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/varela_div90_srznf5gg/obs.tsv
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/varela_div90_srznf5gg/umap.tsv
```

Annotation column:

```text
cluster_number_name
```

Cluster ID column:

```text
cluster_id
```

UMAP audit outputs:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_umap_cluster_label_audit/div90_umap_cluster_label_audit_v1/plots/div90_umap_cluster_number_name_labeled.png
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_umap_cluster_label_audit/div90_umap_cluster_label_audit_v1/plots/div90_umap_cluster_overlay_grid.png
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_umap_cluster_label_audit/div90_umap_cluster_label_audit_v1/tables/div90_cluster_number_name_to_biology_mapping.tsv
```

Machine-readable local planning files:

```text
metadata/div90_jia_lineage_urd_plan.tsv
metadata/div90_jia_rootscore_markers.tsv
metadata/jia_lineage_modules.tsv
```

Run scripts added for the first smoke execution:

```text
python_notebooks/scripts/export_div90_jia_lineage_urd_inputs.py
slurm_templates/34_div90_jia_lineage_urd_smoke.sbatch.template
```

The Slurm run reuses the existing DIV30 URD runner/reporting scripts where the computational logic should match DIV30:

```text
scripts/14_urd_first_pass_pseudotime.R
scripts/15_div30_urd_lineage_decision_report.R
scripts/16_div30_urd_build_lineage_tree.R
scripts/17_div30_urd_finalize_lineage_tree_report.R
scripts/25_div30_urd_jia_fig_s11_marker_validation.R
scripts/27_div90_urd_project_candidate_lineage_markers.R
```

Compatibility note: the DIV90 exporter writes input-bundle files with the legacy `div30_first_urd_*` names because `scripts/14_urd_first_pass_pseudotime.R` consumes that bundle format. The run directory, metadata, root column, pseudotime column, reports, and manifest identify this as DIV90.

## DIV90 Cluster Mapping

| Cluster | Exact metadata name | n cells | First URD role |
|---:|---|---:|---|
| 0 | `0 - MGE Striatal/GP Fated` | 3,601 | LHX8/ISL1-like tip |
| 1 | `1 - SST+, NPY +, Cortical Fated` | 3,548 | LHX6/NFIA-like tip |
| 2 | `2 - CRABP1+/PV Precursors` | 3,503 | v2 CRABP1/ANGPT2-like tip |
| 3 | `3 - PV precursors/Migrating cells/Cortical-fated` | 2,283 | retain in graph, not a v2 tip |
| 4 | `4 - Pre-Astrocytes/Astrocytes 1` | 1,987 | v4 combined astrocyte tip with cluster 10 |
| 5 | `5 - LHX8+ vMGE GABergic Striatal/GP fated 1` | 1,924 | LHX8/ISL1-like tip |
| 6 | `6 - Stressed Cells` | 1,226 | exclude/track separately |
| 7 | `7 - Stressed Cells` | 1,063 | exclude/track separately |
| 8 | `8 - LHX8+ vMGE GABergic Striatal/GP fated 2` | 922 | LHX8/ISL1-like tip |
| 9 | `9 - Pre-OPCs/OPCs` | 915 | v4 OPC tip |
| 10 | `10 - Pre-Astrocytes/Astrocytes 2` | 583 | v4 combined astrocyte tip with cluster 4 |
| 11 | `11 - PV Precursors` | 425 | retain in graph, not a v2 tip |
| 12 | `12 - Dividing cells` | 358 | root candidate pool |

The metadata contains 13 cluster IDs. These reconcile to approximately 10 broader paper-level biology groups because three biology categories are split across two Seurat clusters:

- `4` + `10`: Pre-Astrocytes/Astrocytes
- `5` + `8`: LHX8+ vMGE GABAergic Striatal/GP-fated
- `6` + `7`: Stressed Cells

## Root Definition

Use cluster `12 - Dividing cells` as the candidate root pool.

Within cluster 12, compute:

```text
RootScore =
  z(RGC1_score)
+ z(RGC2_score)
+ z(HES1)
+ z(VIM)
+ z(NES)
- z(IPC_score)
- z(DLX1)
- z(DLX2)
- z(ASCL1)
```

Use the top 1-2% RootScore cells as the first DIV90 Jia-style progenitor root.

Interpretation: the desired root is a proliferative VZ/SVZ-like progenitor compartment: RGC1-high, RGC2-high, HES1/VIM/NES-high, IPC-low, and DLX1/DLX2/ASCL1-low.

Dependency: if DIV90 does not already contain `RGC1_score`, `RGC2_score`, and `IPC_score`, score them using the same Jia program scoring logic used for DIV30 before root selection.

## Tip Definitions

| Tip ID | Jia lineage target | Clusters | Validation genes | Expected biology |
|---|---|---|---|---|
| `tip_lhx8_isl1` | LHX8/ISL1-like lineage | `0`, `5`, `8` | `LHX8`, `ISL1`, `GBX2`, `TAC1` | subpallial cholinergic/striatal-GP-fated inhibitory lineage |
| `tip_lhx6_nfia` | LHX6/NFIA-like lineage | `1` | `LHX6`, `SST`, `NPY`, `ERBB4`, `CXCR4`, `ARX` | cortical interneuron lineage |
| `tip_crabp1_angpt2` | CRABP1/ANGPT2-like lineage | `2` | `CRABP1`, `ANGPT2`, `SPOCK1`, `DSCAM`, `RIPOR2`, `RBP1` | CRABP1/ANGPT2-like PV-precursor/subpallial GABAergic endpoint |
| `tip_astrocytes` | astrocyte endpoint | `4`, `10` | `AQP4`, `GFAP`, `SLC1A3`, `ALDH1L1` | combined pre-astrocyte/astrocyte endpoint |
| `tip_opc` | OPC endpoint | `9` | `PDGFRA`, `OLIG1`, `OLIG2`, `SOX10` | pre-OPC/OPC endpoint |

Do not define the neuronal tips as broad SST/PV/MGE categories. Those labels can be retained as metadata, but the neuronal URD tip construction should be Jia-lineage driven. In v4, glial endpoints are intentionally added as non-Jia terminal states so we can test whether the DIV90 manifold separates neuronal and glial terminal fates cleanly.

Clusters `3` and `11` are retained in the v2 manifold but are not used as tips. Their assignment is tested only after URD tree construction by projecting:

```text
MEF2C, EPHA5, LHX6, CRABP1, LHX8, NR2F1, NR2F2
```

onto the tree and comparing cluster 3/11 marker profiles to the v2 tip profiles.

## Updated Neuron-S9 Seven-Tip Plan

Later DIV90 neuron-only reclustering and Science Data S9 module scoring refined
the neuronal endpoint model. The old v4 tip definition above is still preserved
for reproducibility, but the next all-cell DIV90 run should use this newer
seven-tip setup.

### Science Data S9 Module Scoring Pointer

The detailed, canonical record for the Science Data S9 full-module and TF-only
module scoring workflow is:

```text
python_notebooks/HANDOFF_div90_parent_cluster_subclustering_audit.md
section: Science Data S9 / Five-Class Reference Module Scoring
```

That workflow scores the user-provided workbook:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/reference/science.adw1803_data_s9.xlsx
```

using:

```text
python_notebooks/scripts/div90_neuron_only_reference_module_scoring.py
slurm_templates/44_div90_neuron_only_reference_module_scoring.sbatch.template
```

and writes:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_neuron_only_reference_modules/div90_neuron_only_reference_modules_v1
```

Key outputs to inspect before changing DIV90 tips are:

```text
div90_neuron_only_s9_reference_module_report.md
div90_neuron_only_s9_gene_set_coverage.tsv
div90_neuron_only_s9_best_module_calls.tsv
div90_neuron_only_s9_module_scores_by_neuron_cluster.tsv
div90_neuron_only_s9_full_module_heatmap.png
div90_neuron_only_s9_tf_module_heatmap.png
```

Important interpretation rule: the workbook modules and the five biological
anchor classes are scored separately. Use the anchor calls as the clean
five-class labels, and use the full-module/TF-module heatmaps as supporting
evidence.

Summary of intermediate work:

```text
1. Per-parent subclustering audit showed clusters 1, 3, and 11 have interesting
   marker-supported splits, but those splits were sample/cell-line or depth
   associated under conservative gates.

2. Targeted validation preserved these as biological hypotheses:
   cluster 1 = EDIL3/CRABP1 vs ACKR3/ZEB2
   cluster 3 = BEX3/CRABP1/DCX neuroblast-like state
   cluster 11 = ACKR3/ZEB2/MAF vs CPE/TUBB2A/CNTN1

3. Neuron-only reclustering of original neuronal clusters 0,1,2,3,5,8,11
   selected resolution 0.6, giving 9 neuron-only states.

4. Science Data S9 scoring mapped those 9 states onto five anchor classes:
   EPHA5/MEF2C
   LHX6/NFIA
   CRABP1/ANGPT2
   NR2F1/NR2F2
   LHX8/ISL1
```

The current working neuronal branch model is:

```text
1. LHX8/ISL1-like ventral MGE branch
   neuron-only clusters 4 + 7
   neuron-only cluster 0 retained as possible upstream/intermediate state

2. CRABP1/ANGPT2-like fetal precursor branch
   neuron-only cluster 3
   caveat: ANGPT2 is missing from the H5AD, so this is CRABP1-driven

3. LHX6/NFIA / EPHA5/MEF2C cortical interneuron branch
   neuron-only clusters 1 + 2

4. NR2F1/NR2F2-like branch
   neuron-only clusters 5 + 6 + 8
```

Use all non-stressed DIV90 cells, as in the prior all-cell/glia-tip DIV90 run:

```text
retain original clusters: 0,1,2,3,4,5,8,9,10,11,12
exclude original clusters: 6,7
root pool: original cluster 12 - Dividing cells
root selection: top 10% RootScore within cluster 12
```

Use seven tips:

| Tip ID | Source cells | Expected biology |
|---|---|---|
| `tip_lhx8_isl1_state1` | neuron-only cluster 4 | LHX8/ISL1-like ventral MGE branch state 1 |
| `tip_lhx8_isl1_state2` | neuron-only cluster 7 | LHX8/ISL1-like ventral MGE branch state 2 |
| `tip_crabp1_angpt2_fetal_precursor` | neuron-only cluster 3 | CRABP1/ANGPT2-like fetal precursor branch |
| `tip_lhx6_nfia_epha5_mef2c_cortical` | neuron-only clusters 1 + 2 | LHX6/NFIA plus EPHA5/MEF2C cortical interneuron branch |
| `tip_nr2f1_nr2f2` | neuron-only clusters 5 + 6 + 8 | NR2F1/NR2F2-like branch |
| `tip_astrocytes` | original clusters 4 + 10 | combined astrocyte endpoint, unchanged from prior v4 |
| `tip_opc` | original cluster 9 | OPC endpoint, unchanged from prior v4 |

Neuron-only cluster `0` is deliberately retained but not used as a tip because
it looks immature/upstream by the marker and S9 scoring.

Machine-readable plan:

```text
metadata/div90_jia_lineage_urd_plan_neuron_s9_7tips.tsv
```

Scripts/templates updated for this mode:

```text
python_notebooks/scripts/div90_assign_neuron_s9_7tip_groups.py
python_notebooks/scripts/export_div90_jia_lineage_urd_inputs.py
slurm_templates/34_div90_jia_lineage_urd_smoke.sbatch.template
```

Run with:

```bash
DIV90_TIP_MODE=neuron_s9_7tip
RUN_LABEL=div90_urd_jia_lineage_allcells_root10_neuron_s9_7tips_v1
MAX_CELLS=0
DIV90_ROOT_TOP_PERCENT=10
sbatch slurm_templates/34_div90_jia_lineage_urd_smoke.sbatch.template
```

Submitted all-cell seven-tip run:

```text
initial job 51793603:
  cancelled because 12h walltime crossed a maintenance reservation

active job 51796286:
  RUN_LABEL = div90_urd_jia_lineage_allcells_root10_neuron_s9_7tips_v1
  DIV90_TIP_MODE = neuron_s9_7tip
  MAX_CELLS = 0, meaning all retained cells
  DIV90_ROOT_TOP_PERCENT = 10
  requested walltime = 04:00:00
  requested memory = 128G
  status at 2026-06-14 23:32 EDT = RUNNING on gl3156
```

Expected run root:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_urd_jia_lineage_allcells_root10_neuron_s9_7tips_v1
```

Check command:

```bash
squeue -j 51796286
sacct -j 51796286 --format=JobID,JobName%32,State,ExitCode,Elapsed,AllocCPUS,ReqMem,MaxRSS,NodeList
tail -n 120 /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/34_div90_jia_lineage_urd_smoke_51796286.log
```

Early checkpoint:

```text
The neuron-S9 assignment and input export completed successfully.
The run retained 20,049 cells:
  all non-stressed cells from original clusters 0,1,2,3,4,5,8,9,10,11,12

Root pool:
  original cluster 12, n = 358
  selected top 10% RootScore roots, n = 36

R URD stage reached:
  scripts/14_urd_first_pass_pseudotime.R
  create_filtered_urd / normalization stage
```

All-cell seven-tip assignment counts:

| Role | Tip group | n cells |
|---|---|---:|
| retained upstream | `retained_upstream_lhx8_isl1_neuron0` | 3,479 |
| root pool | `root_pool_cluster12_dividing_cells` | 358 |
| tip | `tip_lhx8_isl1_state1` | 1,866 |
| tip | `tip_lhx8_isl1_state2` | 926 |
| tip | `tip_crabp1_angpt2_fetal_precursor` | 1,948 |
| tip | `tip_lhx6_nfia_epha5_mef2c_cortical` | 4,747 |
| tip | `tip_nr2f1_nr2f2` | 3,240 |
| tip | `tip_astrocytes` | 2,570 |
| tip | `tip_opc` | 915 |
| excluded | `excluded_stressed` | 2,289 |

## First Smoke Run

Retain:

```text
0, 1, 2, 3, 5, 8, 11, 12
```

Exclude or track separately for the first neuronal lineage smoke test:

```text
4, 6, 7, 9, 10
```

Corrected v2 retained roles:

```text
root pool: 12
tip_lhx8_isl1: 0, 5, 8
tip_lhx6_nfia: 1
tip_crabp1_angpt2: 2
retained unassigned candidates: 3, 11
```

Run the same URD major stages used for DIV30:

```text
create URD object
calculate PCA/diffusion space if needed
flood pseudotime from RootScore root
pseudotimeDetermineLogistic()
pseudotimeWeightTransitionMatrix()
simulateRandomWalksFromTips()
processRandomWalksFromTips()
buildTree()
```

### Historical v1 Smoke Run

Submitted: 2026-06-11

Slurm job:

```text
51685516
```

Template:

```text
slurm_templates/34_div90_jia_lineage_urd_smoke.sbatch.template
```

Run label:

```text
div90_urd_jia_lineage_smoke5k_knn100_v1
```

Output root:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_urd_jia_lineage_smoke5k_knn100_v1/
```

Log:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/34_div90_jia_lineage_urd_smoke_51685516.log
```

Submitted parameters:

| Parameter | Value |
|---|---|
| `MAX_CELLS` | `5000` |
| retained clusters | `0,1,2,3,5,8,11,12` |
| excluded first-smoke clusters | `4,6,7,9,10` |
| root cluster | `12 - Dividing cells` |
| root score top percent | `2` |
| root minimum cells | `8` |
| pseudotime name | `div90_jia_rootscore_root` |
| URD KNN | `100` |
| URD variable genes | `3000` |
| URD floods | `20` |
| tree tips | `tip_lhx8_isl1`, `tip_epha5_mef2c`, `tip_lhx6_nfia` |
| random walks per tip | `5000` |

This v1 run is now historical/superseded for DIV90 tip logic because it grouped clusters `2`, `3`, and `11` together as `tip_epha5_mef2c`. The corrected v2 run keeps clusters `3` and `11` unassigned and uses cluster `2` as the explicit `tip_crabp1_angpt2` endpoint.

Pre-submit exporter smoke check:

```text
/tmp/div90_jia_urd_export_test/
```

The exporter smoke check selected 800 cells with all 358 cluster-12 dividing cells retained, and selected 8 RootScore root cells. R read the exported `urd_root_candidate` column as logical, so the root column format is valid for URD.

### Completed Smoke Run Summary

Completed: 2026-06-11 21:12:49 EDT

Slurm accounting:

| Job | State | Exit code | Elapsed | Max RSS | CPUs |
|---|---|---|---|---|---:|
| `51685516.batch` | `COMPLETED` | `0:0` | `00:20:25` | `4007404K` | 8 |

Root selection:

| Metric | Value |
|---|---:|
| root cluster | 12 |
| root top percent | 2.0 |
| root pool cells retained | 358 |
| root cells selected | 8 |
| selected root RootScore min | 24.4202 |
| selected root RootScore median | 24.7901 |
| selected root RootScore max | 27.1786 |

Pseudotime smoke result:

| Annotation | n cells | Median pseudotime |
|---|---:|---:|
| `12 - Dividing cells` | 358 | 0.244 |
| `2 - CRABP1+/PV Precursors` | 1003 | 0.383 |
| `8 - LHX8+ vMGE GABergic Striatal/GP fated 2` | 264 | 0.401 |
| `11 - PV Precursors` | 122 | 0.408 |
| `1 - SST+, NPY +, Cortical Fated` | 1016 | 0.412 |
| `5 - LHX8+ vMGE GABergic Striatal/GP fated 1` | 551 | 0.413 |
| `0 - MGE Striatal/GP Fated` | 1032 | 0.421 |
| `3 - PV precursors/Migrating cells/Cortical-fated` | 654 | 0.440 |

Flood stability reached 1.0 Spearman correlation to final pseudotime at 20 walks/cell. Early stability improved monotonically from 0.606 at 2 walks/cell to 0.984 at 18 walks/cell.

Tree result:

| Metric | Value |
|---|---:|
| requested tips | 3 |
| final segment joins | 2 |
| final segments | 3 |
| distinct branching detected | TRUE |
| tree-layout cells | 4512 |

Final tree topology:

```text
tip_epha5_mef2c + tip_lhx6_nfia fuse first.
The combined branch then splits from tip_lhx8_isl1 at pseudotime ~0.267.
```

URD reported that the difference between `tip_epha5_mef2c` and `tip_lhx6_nfia` was always false in the divergence scan, so this smoke tree does not support a strong EPHA5/MEF2C-vs-LHX6/NFIA split at 5k cells. It does support a stronger LHX8/ISL1-like branch versus the combined EPHA5/MEF2C/LHX6/NFIA-like branch.

Branch-specific gene examples:

| Tip | Top branch-enriched genes |
|---|---|
| `tip_lhx8_isl1` | `LHX8`, `ISLR2`, `TSHZ2`, `ECEL1`, `GAP43`, `SPOCK2`, `NRP1`, `VAT1L` |
| `tip_epha5_mef2c` | `ERBB4`, `NXPH1`, `GRIA3`, `TCF4`, `ZEB2`, `SST`, `QKI`, `NFIB`, `NPY`, `ARX`, `EPHA5`, `LHX6` |
| `tip_lhx6_nfia` | `VWC2`, `ERBB4`, `NFIB`, `RAB3IP`, `EDIL3`, `NXPH1`, `GRIA3`, `CXCR4`, `LHX6`, `SST`, `NR2F1`, `CRABP1` |

S11 marker validation:

All requested Jia Fig. S11-style markers were present:

```text
HES1, CACNA1E, DLX2, DCX, LHX8, NR2F1, EPHA5, MEF2C, CRABP1
```

### Corrected v2 Smoke Run

Purpose: rerun DIV90 URD with cluster `2 - CRABP1+/PV Precursors` as the CRABP1/ANGPT2-like tip and clusters `3`/`11` retained but not supplied as tips.

Original v2 submission: 2026-06-11 21:43 EDT

Slurm job:

```text
51687415
```

Status:

```text
CANCELLED before start because the initial resource request was overlarge for a 5k smoke run.
Initial request: 8 CPUs, 160G RAM, 48h.
```

Lean v2 resubmission: 2026-06-11 21:57 EDT

Active Slurm job:

```text
51687695
```

Active request:

```text
2 CPUs, 32G RAM, 4h wall time
```

Active-job status at 2026-06-11 21:57 EDT:

```text
RUNNING on gl3018
```

Completed: 2026-06-11 22:18:32 EDT

Slurm accounting:

| Job | State | Exit code | Elapsed | Max RSS | CPUs | Request |
|---|---|---|---|---|---:|---|
| `51687695.batch` | `COMPLETED` | `0:0` | `00:20:56` | `4401788K` | 2 | `32G` |

Corrected v2 tree result:

| Metric | Value |
|---|---:|
| requested tips | 3 |
| final segment joins | 2 |
| final segments | 3 |
| distinct branching detected | TRUE |

Tip mapping:

| Tip ID | Tip label | n tip cells | median pseudotime |
|---:|---|---:|---:|
| 1 | `tip_lhx8_isl1` | 1847 | 0.417 |
| 2 | `tip_lhx6_nfia` | 1016 | 0.412 |
| 3 | `tip_crabp1_angpt2` | 1003 | 0.383 |

Topology note:

```text
tip_lhx6_nfia and tip_crabp1_angpt2 join together first.
The combined segment then splits from tip_lhx8_isl1 at pseudotime ~0.267.
```

Candidate clusters retained but not used as tips:

| Candidate cluster | Tree cells | Best marker-profile match | Pearson | Highest lineage proxy |
|---:|---:|---|---:|---|
| 3 | 52 | `tip_crabp1_angpt2` | 0.962 | `LHX6_NFIA_proxy` |
| 11 | 35 | `tip_crabp1_angpt2` | 0.607 | `LHX6_NFIA_proxy` |

Interpretation guardrail: these are candidate assignments from marker projection only. Inspect `candidate_pv_marker_projection_v1/plots/div90_candidate_marker_tree_overlays.png` and the profile heatmap before promoting clusters `3`/`11` into explicit tips.

Run label:

```text
div90_urd_jia_lineage_smoke5k_knn100_v2_crabp1_tip
```

Output root:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_urd_jia_lineage_smoke5k_knn100_v2_crabp1_tip/
```

Log:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/34_div90_jia_lineage_urd_smoke_51687695.log
```

Corrected v2 submitted parameters:

| Parameter | Value |
|---|---|
| retained clusters | `0,1,2,3,5,8,11,12` |
| excluded clusters | `4,6,7,9,10` |
| root pool | `12 - Dividing cells` |
| root score | `z(RGC1)+z(RGC2)+z(HES1)+z(VIM)+z(NES)-z(IPC)-z(DLX1)-z(DLX2)-z(ASCL1)` |
| tree tips | `tip_lhx8_isl1`, `tip_lhx6_nfia`, `tip_crabp1_angpt2` |
| retained unassigned candidates | `3`, `11` |
| post-tree candidate marker genes | `MEF2C`, `EPHA5`, `LHX6`, `CRABP1`, `LHX8`, `NR2F1`, `NR2F2` |

Corrected v2 candidate-assignment outputs:

```text
candidate_pv_marker_projection_v1/plots/div90_candidate_marker_tree_overlays.png
candidate_pv_marker_projection_v1/plots/div90_candidate_and_tip_marker_profile_heatmap.png
candidate_pv_marker_projection_v1/tables/div90_candidate_cluster_lineage_assignment.tsv
candidate_pv_marker_projection_v1/tables/div90_candidate_cluster_to_tip_profile_correlations.tsv
candidate_pv_marker_projection_v1/tables/div90_candidate_and_tip_marker_profiles.tsv
candidate_pv_marker_projection_v1/tables/div90_marker_expression_summary_by_cluster.tsv
candidate_pv_marker_projection_v1/tables/div90_marker_expression_summary_by_tree_segment.tsv
```

### Corrected v3 Smoke Run: Glia Retained As Cells

Purpose: rerun the corrected v2 logic but add glial/OPC clusters back into the manifold as ordinary cells. They are not tips.

Change from v2:

| Role | v2 | v3 |
|---|---|---|
| retained clusters | `0,1,2,3,5,8,11,12` | `0,1,2,3,4,5,8,9,10,11,12` |
| glia/OPC clusters | excluded | retained as non-tip cells |
| stressed clusters | excluded | excluded |
| tips | `0+5+8`, `1`, `2` | unchanged |

v3 retained non-tip glia:

```text
4  = Pre-Astrocytes/Astrocytes 1
9  = Pre-OPCs/OPCs
10 = Pre-Astrocytes/Astrocytes 2
```

v3 excluded:

```text
6, 7 = Stressed Cells
```

Run label:

```text
div90_urd_jia_lineage_smoke5k_knn100_v3_glia_cells
```

Submitted: 2026-06-11 22:42 EDT

Slurm job:

```text
51688064
```

Final status:

```text
COMPLETED, exit code 0
```

Request:

```text
2 CPUs, 32G RAM, 4h wall time
```

Runtime/resource use:

```text
Elapsed: 00:20:44
MaxRSS: 4,124,912K
Node: gl3019
```

Output root:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_urd_jia_lineage_smoke5k_knn100_v3_glia_cells/
```

Log:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/34_div90_jia_lineage_urd_smoke_51688064.log
```

Expected extra tree plot for direct comparison to the DIV90 UMAP labels:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_urd_jia_lineage_smoke5k_knn100_v3_glia_cells/lineage_tree_cluster_number_name_v1/plots/urd_tree_annotation.png
```

This plot colors the URD tree by `cluster_number_name`, so it shows the original DIV90 biological cluster names, including retained glial cells.

v3 verified input manifest:

```text
retained_clusters: 0,1,2,3,4,5,8,9,10,11,12
excluded_clusters_first_smoke: 6,7
retained_glia_non_tip_clusters: 4,9,10
tip_lhx8_isl1_clusters: 0,5,8
tip_lhx6_nfia_clusters: 1
tip_crabp1_angpt2_clusters: 2
retained_unassigned_candidate_clusters: 3,11
n_selected_cells: 5000
n_root_cells: 8
```

v3 sampled cluster composition:

| Cluster | Exact metadata name | sampled cells | URD role |
|---:|---|---:|---|
| 0 | `0 - MGE Striatal/GP Fated` | 850 | `tip_lhx8_isl1` |
| 1 | `1 - SST+, NPY +, Cortical Fated` | 836 | `tip_lhx6_nfia` |
| 2 | `2 - CRABP1+/PV Precursors` | 826 | `tip_crabp1_angpt2` |
| 3 | `3 - PV precursors/Migrating cells/Cortical-fated` | 538 | retained candidate, not a tip |
| 4 | `4 - Pre-Astrocytes/Astrocytes 1` | 468 | retained glia, not a tip |
| 5 | `5 - LHX8+ vMGE GABergic Striatal/GP fated 1` | 454 | `tip_lhx8_isl1` |
| 8 | `8 - LHX8+ vMGE GABergic Striatal/GP fated 2` | 217 | `tip_lhx8_isl1` |
| 9 | `9 - Pre-OPCs/OPCs` | 216 | retained glia, not a tip |
| 10 | `10 - Pre-Astrocytes/Astrocytes 2` | 137 | retained glia, not a tip |
| 11 | `11 - PV Precursors` | 100 | retained candidate, not a tip |
| 12 | `12 - Dividing cells` | 358 | root candidate pool |

v3 tree result:

```text
tree_slot_length: 18
n_requested_tips: 3
n_segment_joins: 2
n_segments: 3
has_distinct_branching: TRUE
tips: 1=tip_lhx8_isl1; 2=tip_lhx6_nfia; 3=tip_crabp1_angpt2
```

URD joined the LHX6/NFIA-like and CRABP1/ANGPT2-like tips first, then split the LHX8/ISL1-like tip from that combined side:

```text
2 + 3 -> segment 4 at pseudotime 0.586
1 + 4 -> segment 5 at pseudotime 0.336
```

The retained glial/OPC clusters were not used as tips. Their median URD pseudotimes were:

| Cluster | Median URD pseudotime |
|---|---:|
| `10 - Pre-Astrocytes/Astrocytes 2` | 0.269 |
| `4 - Pre-Astrocytes/Astrocytes 1` | 0.275 |
| `9 - Pre-OPCs/OPCs` | 0.337 |

Tip composition remained unchanged from corrected v2:

| Tip | n tip cells | median pseudotime | composition |
|---|---:|---:|---|
| `tip_lhx8_isl1` | 1521 | 0.474 | clusters `0`, `5`, `8` |
| `tip_lhx6_nfia` | 836 | 0.472 | cluster `1` |
| `tip_crabp1_angpt2` | 826 | 0.442 | cluster `2` |

Clusters `3` and `11` remained non-tip candidate clusters. Post-tree marker-profile projection still matched both most closely to `tip_crabp1_angpt2`, but the guardrail remains: inspect tree position, marker overlays, and profile correlations before promoting either as a tip.

| Candidate cluster | n tree cells | best marker-profile match | Pearson | highest lineage proxy |
|---:|---:|---|---:|---|
| 3 | 46 | `tip_crabp1_angpt2` | 0.974 | `LHX6_NFIA_proxy` |
| 11 | 36 | `tip_crabp1_angpt2` | 0.588 | `LHX6_NFIA_proxy` |

Key v3 output paths:

```text
lineage_decision_report/plots/umap_pseudotime.png
lineage_decision_report/plots/diffusion_map_pseudotime.png
lineage_decision_report/plots/diffusion_map_annotation.png
lineage_decision_report/plots/flood_stability.png
lineage_decision_report/plots/gene_cascade_heatmap.png
lineage_tree_jia_endpoint_tips_v1/plots/urd_tree_annotation.png
lineage_tree_jia_endpoint_tips_v1/plots/urd_tree_pseudotime.png
lineage_tree_cluster_number_name_v1/plots/urd_tree_annotation.png
jia_fig_s11_style_marker_validation_v1/plots/jia_fig_s11_style_urd_marker_validation.png
candidate_pv_marker_projection_v1/plots/div90_candidate_marker_tree_overlays.png
candidate_pv_marker_projection_v1/plots/div90_candidate_and_tip_marker_profile_heatmap.png
```

### Corrected v4 Smoke Run: Glial Endpoints As Tips

Purpose: rerun the v3 inclusion logic but now promote the retained glial populations to explicit terminal tips. This is the final inclusion-logic smoke before moving toward the production/full-scale DIV90 run.

Change from v3:

| Role | v3 | v4 |
|---|---|---|
| retained clusters | `0,1,2,3,4,5,8,9,10,11,12` | unchanged |
| stressed clusters | excluded `6,7` | unchanged |
| neuronal tips | `0+5+8`, `1`, `2` | unchanged |
| clusters `3`/`11` | retained candidate cells, not tips | unchanged |
| astrocyte clusters `4`/`10` | retained non-tip cells | combined into `tip_astrocytes` |
| OPC cluster `9` | retained non-tip cells | promoted to `tip_opc` |

v4 tips:

```text
tip_lhx8_isl1      = clusters 0,5,8
tip_lhx6_nfia      = cluster 1
tip_crabp1_angpt2  = cluster 2
tip_astrocytes     = clusters 4,10
tip_opc            = cluster 9
```

v4 retained non-tip candidates:

```text
3  = PV precursors/Migrating cells/Cortical-fated
11 = PV Precursors
```

v4 excluded:

```text
6, 7 = Stressed Cells
```

Run label:

```text
div90_urd_jia_lineage_smoke5k_knn100_v4_glia_tips
```

Submitted: 2026-06-11

Slurm job:

```text
51690159
```

Final status:

```text
COMPLETED, exit code 0
```

Runtime/resource use:

```text
Elapsed: 00:23:55
MaxRSS: 4,331,164K
Node: gl3219
```

Log:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/34_div90_jia_lineage_urd_smoke_51690159.log
```

Output root:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_urd_jia_lineage_smoke5k_knn100_v4_glia_tips/
```

Primary comparison plots expected from v4:

```text
lineage_tree_jia_endpoint_tips_v1/plots/urd_tree_annotation.png
lineage_tree_cluster_number_name_v1/plots/urd_tree_annotation.png
lineage_tree_jia_endpoint_tips_v1/plots/urd_tree_pseudotime.png
lineage_decision_report/plots/umap_pseudotime.png
lineage_decision_report/plots/diffusion_map_annotation.png
lineage_decision_report/plots/diffusion_map_pseudotime.png
jia_fig_s11_style_marker_validation_v1/plots/jia_fig_s11_style_urd_marker_validation.png
```

Decision focus for v4:

```text
Does the DIV90 URD tree separate neuronal Jia-lineage endpoints from glial endpoints?
Do astrocytes and OPCs form distinct terminal branches?
Do clusters 3 and 11 remain aligned with neuronal lineages or drift toward a glial branch?
```

v4 verified input manifest:

```text
retained_clusters: 0,1,2,3,4,5,8,9,10,11,12
excluded_clusters_first_smoke: 6,7
tip_lhx8_isl1_clusters: 0,5,8
tip_lhx6_nfia_clusters: 1
tip_crabp1_angpt2_clusters: 2
tip_astrocytes_clusters: 4,10
tip_opc_clusters: 9
retained_unassigned_candidate_clusters: 3,11
retained_glia_non_tip_clusters: empty
n_selected_cells: 5000
n_root_cells: 8
```

v4 tree result:

```text
tree_slot_length: 18
n_requested_tips: 5
n_segment_joins: 5
n_segments: 6
has_distinct_branching: TRUE
tips: 1=tip_lhx8_isl1; 2=tip_lhx6_nfia; 3=tip_crabp1_angpt2; 4=tip_astrocytes; 5=tip_opc
```

Final segment-join table:

```text
parent 8 -> child 1 at pseudotime 0.314
parent 8 -> child 6 at pseudotime 0.314
parent 8 -> child 5 at pseudotime 0.314
parent 9 -> child 4 at pseudotime 0
parent 9 -> child 8 at pseudotime 0
```

The detailed tree-build log shows the neuronal LHX6/NFIA-like and CRABP1/ANGPT2-like tips joined first, matching v2/v3:

```text
2 + 3 -> segment 6 at pseudotime 0.586
1 + 6 -> segment 7/8 side at pseudotime 0.314
5 + neuronal side joins at pseudotime 0.314
4 + all other branches joins at pseudotime 0
```

v4 tip composition:

| Tip | n tip cells | median pseudotime | composition |
|---|---:|---:|---|
| `tip_lhx8_isl1` | 1521 | 0.474 | clusters `0`, `5`, `8` |
| `tip_lhx6_nfia` | 836 | 0.472 | cluster `1` |
| `tip_crabp1_angpt2` | 826 | 0.442 | cluster `2` |
| `tip_astrocytes` | 605 | 0.274 | clusters `4`, `10` |
| `tip_opc` | 216 | 0.337 | cluster `9` |

v4 candidate-cluster projection:

| Candidate cluster | n tree cells | best marker-profile match | Pearson | highest lineage proxy |
|---:|---:|---|---:|---|
| 3 | 47 | `tip_crabp1_angpt2` | 0.971 | `LHX6_NFIA_proxy` |
| 11 | 36 | `tip_opc` | 0.719 | `LHX6_NFIA_proxy` |

Interpretation guardrail: cluster `3` remains more CRABP1/ANGPT2-like by the current marker-profile test. Cluster `11` now correlates most strongly with the OPC tip after OPCs are included as a terminal state, so it should remain a candidate until marker overlays and full-scale behavior are reviewed.

Key v4 output paths:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_urd_jia_lineage_smoke5k_knn100_v4_glia_tips/lineage_tree_cluster_number_name_v1/plots/urd_tree_annotation.png
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_urd_jia_lineage_smoke5k_knn100_v4_glia_tips/lineage_tree_jia_endpoint_tips_v1/plots/urd_tree_annotation.png
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_urd_jia_lineage_smoke5k_knn100_v4_glia_tips/lineage_tree_jia_endpoint_tips_v1/plots/urd_tree_pseudotime.png
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_urd_jia_lineage_smoke5k_knn100_v4_glia_tips/lineage_decision_report/plots/umap_pseudotime.png
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_urd_jia_lineage_smoke5k_knn100_v4_glia_tips/lineage_decision_report/plots/diffusion_map_annotation.png
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_urd_jia_lineage_smoke5k_knn100_v4_glia_tips/jia_fig_s11_style_marker_validation_v1/plots/jia_fig_s11_style_urd_marker_validation.png
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_urd_jia_lineage_smoke5k_knn100_v4_glia_tips/candidate_pv_marker_projection_v1/plots/div90_candidate_marker_tree_overlays.png
```

### Production Scaling Decision

Date: 2026-06-12

Freeze this v4 inclusion logic for the next DIV90 full-cell run:

```text
root pool: cluster 12 Dividing cells
root cells: top 2% Jia RootScore within cluster 12, minimum 8 cells
retained clusters: 0,1,2,3,4,5,8,9,10,11,12
excluded clusters: 6,7 stressed cells
tips:
  tip_lhx8_isl1      = 0,5,8
  tip_lhx6_nfia      = 1
  tip_crabp1_angpt2  = 2
  tip_astrocytes     = 4,10
  tip_opc            = 9
retained non-tip candidates: 3,11
```

Why DIV90 full is reasonable to submit next:

```text
v4 smoke selected cells: 5,000
full retained DIV90 cells: approximately 20,248 after excluding stressed clusters 6/7
v4 smoke runtime: 00:23:55
v4 smoke MaxRSS: 4,331,164K
```

Prediction: full DIV90 is a moderate-risk scale-up, not a blind leap. It is roughly 4x the v4 smoke by cell count, and the 5k run had a wide memory margin. Diffusion-map and random-walk steps may scale nonlinearly, so use production resources and the new stage/resource logs.

Recommended full submission:

```bash
sbatch \
  --time=48:00:00 \
  --cpus-per-task=8 \
  --mem=160G \
  --export=ALL,RUN_LABEL=div90_urd_jia_lineage_full_v4_glia_tips,MAX_CELLS=0,SEED=7,URD_SEED=7,DIV90_ROOT_TOP_PERCENT=2,HEARTBEAT_INTERVAL_SECONDS=300 \
  slurm_templates/34_div90_jia_lineage_urd_smoke.sbatch.template
```

`MAX_CELLS=0` means all retained cells in the DIV90 exporter. The template name still says `smoke` for historical continuity, but this command promotes the same frozen v4 logic to the full retained DIV90 dataset.

Submitted full DIV90 run:

```text
Date submitted: 2026-06-12
Slurm job: 51691452
State at submission check: PENDING (Priority)
Run label: div90_urd_jia_lineage_full_v4_glia_tips
MAX_CELLS: 0
Resources requested: 8 CPUs, 160G, 48:00:00
Log:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/34_div90_jia_lineage_urd_smoke_51691452.log
Output root:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_urd_jia_lineage_full_v4_glia_tips/
```

New monitoring/checkpoint outputs for this and later runs:

```text
tables/slurm_resource_heartbeat.tsv
tables/div30_first_urd_stage_timings.tsv
lineage_tree_jia_endpoint_tips_v1/tables/lineage_tree_stage_timings.tsv
lineage_tree_jia_endpoint_tips_v1/tables/slurm_resource_heartbeat.tsv
```

Interpretation rule for full DIV90: if the full run changes cluster `11` away from the OPC-like behavior seen in the 5k tree, report that as a scale-dependent candidate reassignment rather than forcing the 5k interpretation.

## Required Outputs

Every DIV90 smoke/full run should write:

```text
URD object RDS
root candidate table
tip composition table
tree tip mapping table
tree segment table
run parameter table
```

Expected first-URD outputs:

```text
div30_first_urd_object.rds
tables/div30_first_urd_parameters.tsv
tables/div30_first_urd_pseudotime.tsv
tables/div30_first_urd_summary.tsv
tables/div30_first_urd_marker_correlations.tsv
plots/div30_first_urd_umap_pseudotime.png
lineage_decision_report/
```

Expected smoke-panel outputs matching the DIV30 smoke panel:

```text
lineage_decision_report/plots/umap_pseudotime.png
lineage_decision_report/plots/diffusion_map_pseudotime.png
lineage_decision_report/plots/diffusion_map_annotation.png
lineage_decision_report/plots/flood_stability.png
lineage_decision_report/plots/gene_cascade_heatmap.png
lineage_decision_report/plots/tree_visualization.png
lineage_decision_report/plots/lineage_decision_tree.png
lineage_decision_report/tables/top_positive_pseudotime_genes.tsv
lineage_decision_report/tables/top_negative_pseudotime_genes.tsv
lineage_decision_report/tables/pseudotime_ordering_by_annotation.tsv
lineage_decision_report/tables/root_annotation_composition.tsv
```

Expected tree outputs:

```text
lineage_tree_jia_endpoint_tips_v1/div30_urd_lineage_tree_object.rds
lineage_tree_jia_endpoint_tips_v1/tables/tree_tip_mapping.tsv
lineage_tree_jia_endpoint_tips_v1/tables/tree_tip_composition.tsv
lineage_tree_jia_endpoint_tips_v1/tables/tree_segment_joins.tsv
lineage_tree_jia_endpoint_tips_v1/tables/branch_specific_genes.tsv
lineage_tree_jia_endpoint_tips_v1/plots/pseudotime_logistic.png
lineage_tree_jia_endpoint_tips_v1/plots/urd_tree_annotation.png
lineage_tree_jia_endpoint_tips_v1/plots/urd_tree_pseudotime.png
```

Expected Jia Fig. S11-style validation outputs:

```text
jia_fig_s11_style_marker_validation_v1/plots/jia_fig_s11_style_urd_marker_validation.png
jia_fig_s11_style_marker_validation_v1/tables/jia_fig_s11_marker_order.tsv
jia_fig_s11_style_marker_validation_v1/tables/jia_fig_s11_marker_expression_long.tsv.gz
jia_fig_s11_style_marker_validation_v1/tables/jia_fig_s11_marker_expression_summary.tsv
```

Required smoke validation figures:

```text
UMAP pseudotime
diffusion map pseudotime
diffusion map annotation
flood stability
URD tree
top pseudotime gene cascade
top positive pseudotime genes table
top negative pseudotime genes table
```

Required Jia-lineage figures after tree construction:

```text
segment x Jia lineage score heatmap
segment x Jia lineage z-score heatmap
Jia lineage score tree overlays
branchpoint decision-gene heatmaps
branchpoint decision-gene volcano plots
ranked branchpoint marker tables
```

Required Jia Fig. S11-style marker-expression figure:

```text
Panel A: HES1 | CACNA1E | DLX2 | DCX
Panel B: LHX8 | NR2F1 | EPHA5 | MEF2C | CRABP1
```

This marker figure must use expression only, no z-score scaling across genes, and the same URD tree coordinates for every panel.

## Post-Tree Interpretation

Score the actual Jia lineage modules:

```text
LHX8/ISL1
NR2F1/NR2F2
EPHA5/MEF2C
LHX6/NFIA
CRABP1/ANGPT2
```

Then determine whether the chosen tips actually correspond to the expected Jia lineages.

Important separation of analyses:

- Jia lineage scores answer: which Jia programs localize to each branch or segment?
- Branchpoint DE answers: which genes define lineage divergence?
- Jia Fig. S11-style markers answer: whether the organoid URD tree recapitulates canonical developmental and lineage marker localization.

## Open Checks

- Confirm DIV90 Jia program score availability before root scoring.
- Confirm cluster 12 connects to the neuronal manifold.
- Confirm cluster 2 behavior after scoring because `CRABP1+/PV Precursors` may bridge EPHA5/MEF2C and CRABP1/ANGPT2 biology.
- Decide after the first tree whether NR2F1/NR2F2 or CRABP1/ANGPT2 require their own explicit terminal tips.
- Keep the Stiletti UMAP comparison as a separate upstream process using the same DIV30/DIV90 cluster mapping files.

## Current Production Rerun Name

The next DIV90 all-cell run should be treated as the new production/resumable
version, not as another smoke/backfill/final scratch run.

Use this run label:

```text
div90_allcells_jia_root10_neuron_s9_7tips_urd_resumable_v1
```

Full output path:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_allcells_jia_root10_neuron_s9_7tips_urd_resumable_v1
```

The biological setup is unchanged from the all-cell neuron-S9 seven-tip plan:

```text
MAX_CELLS=0
DIV90_TIP_MODE=neuron_s9_7tip
DIV90_ROOT_TOP_PERCENT=10
DIV90_ROOT_MIN_CELLS=1
TREE_TIP_LABELS=tip_lhx8_isl1_state1,tip_lhx8_isl1_state2,tip_crabp1_angpt2_fetal_precursor,tip_lhx6_nfia_epha5_mef2c_cortical,tip_nr2f1_nr2f2,tip_astrocytes,tip_opc
```

The computational URD parameters also remain unchanged for consistency with the
prior DIV30/DIV90 runs:

```text
URD_KNN=100
URD_N_FLOODS=20
URD_NUM_VARIABLE_GENES=3000
URD_PCA_MP_FACTOR=2
URD_SIGMA=local
```

Only execution mechanics changed. The shared R runner now writes resumable
checkpoints under:

```text
checkpoints/urd_after_filter.rds
checkpoints/urd_after_variable_genes.rds
checkpoints/urd_after_pca.rds
checkpoints/urd_after_diffusion_map.rds
checkpoints/urd_after_flood_pseudotime.rds
tables/div30_first_urd_checkpoint_manifest.tsv
```

Resubmitting the same run label with `URD_RESUME=true` will continue from the
latest checkpoint and then regenerate the same downstream decision reports,
lineage-tree reports, Jia Fig. S11 marker validation, candidate marker
projection, and output manifests.

## Proposed downstream step: genes associated with each DIV90 divergence (2026-10-21)

Planning revision timestamp: `2026-10-21T14:10:54-04:00`
(`America/Detroit`, EDT; `2026-10-21T18:10:54Z`). This records the planning
revision, not an analysis execution or submission time.

### Scope and review boundary

Analyze every divergence in the saved final seven-tip tree, including neuronal
and glial branches. Preserve topology, named pseudotime, terminal identities,
cell assignments, UMAP coordinates, existing results, and final figure packages.
Identify expression associated with divergence and candidate transcriptional
regulators; expression differences do not establish causal fate determinants.
An inferred URD split is not necessarily the biological time of commitment.

This section began as the implementation specification reviewed with the user.
The user subsequently explicitly authorized implementation and execution
("proceed ... run the analysis ... now"). That instruction supersedes the
earlier hold on execution. Read-only object inspection informed the design;
synthetic validation precedes biological execution. Submission/run records
below distinguish authorization, validation, and actual scheduler state.

### Storage and job conventions

Follow `WORKFLOW.md`: the checkout contains code/configuration/documentation;
runtime data, generated figures, tables, caches, job files, and logs live on Turbo.
Only this existing handoff was edited at the planning stage. Implementation
code and configuration are now present in the repository; runtime artifacts
remain on NFS, including clearly named synthetic validation packages.

Protected input (read-only):

```text
PROJECT_ROOT=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
SOURCE_RUN=${PROJECT_ROOT}/results/div90_jia_lineage_urd/div90_allcells_jia_root10_neuron_s9_7tips_urd_resumable_v1
TREE_RDS=${SOURCE_RUN}/lineage_tree_jia_endpoint_tips_v1/div30_urd_lineage_tree_object.rds
```

Proposed independent runtime destination, outside the protected source run:

```text
${PROJECT_ROOT}/results/div90_urd_branchpoint_genes/div90_final_7tips_branchpoint_genes_v1_<YYYYMMDD_HHMMSS_EDT-or-EST>_<git-short>/
  code/          # carbon copy of the exact reviewed code executed
  config/        # carbon copy of resolved settings and reference paths
  tables/
  plots/
  provenance/
  checkpoints/   # small per-comparison results only, not copies of the URD object
```

Implementation files:

- `scripts/40_div90_urd_branchpoint_genes.R`
- `config/div90_urd_branchpoint_genes.tsv`
- `slurm_templates/64_div90_urd_branchpoint_genes.sbatch.template`
- `scripts/40_submit_div90_urd_branchpoint_genes.py` (dry run, carbon copies,
  reference provenance, immutable job rendering and submission).

After implementation review and explicit execution approval, copy/render the
template to `${PROJECT_ROOT}/jobs/<run_id>.sbatch`
and submit that NFS job file with `sbatch`. Use absolute Slurm log paths under
`${PROJECT_ROOT}/logs/64_div90_urd_branchpoint_genes_%j.log`; Slurm does not
expand shell variables in `#SBATCH` paths. Use the established `standard`
partition, `parent0` account, and module-managed R/4.5.1 with the existing 4.4
URD library and 4.5 `svglite` library, matching the final URD figure's library
combination. R/4.4.3 alone lacks `svglite` here. No package installations in jobs.
Initial resources are one
node, 4 CPUs, 64G, 12 hours; these are planning estimates to review after the
implementation's memory strategy is known. Keep sparse matrices sparse and
process genes/comparisons in blocks. Do not call the original full-pipeline or
posthoc-refresh templates: they regenerate source outputs and have historical
defaults. The submission record below gives the actual job ID and paths.

#### Timestamping and carbon copies at submission

This follows the existing snapshot/submission convention documented in
`paper2_mge_organoid_atlas/HANDOFF.md` and implemented in its `bin/submit.sh`;
frozen-code execution is also used in `paper3_pcdh19/slurm/hicat_method_pilot.sbatch`.
It does not route this DIV90 analysis into either paper's data directories.

- Generate one run ID before submission, using `America/Detroit` wall time,
  an explicit EDT/EST suffix, and the Git short revision. Store the exact
  ISO-8601 timestamp with numeric UTC offset and a UTC counterpart in the
  manifest. Refuse a run-directory collision instead of overwriting it.
- Before `sbatch`, carbon-copy the reviewed R script, all sourced helpers,
  resolved configuration, and the relevant handoff into that NFS run package.
  Preserve relative helper paths. Save Git commit, working-tree status and
  relevant uncommitted diff: the commit alone does not capture edited code.
- Keep the exact rendered `.sbatch` at `${PROJECT_ROOT}/jobs/<run_id>.sbatch`
  and a byte-identical copy at `<run_dir>/provenance/submitted.sbatch`.
  Record file sizes and SHA256 hashes for the code/config/job snapshots.
- The job executes `<run_dir>/code/...` with `<run_dir>/config/...`, verifies
  those hashes before computation, and never reads mutable repo scripts or
  configuration while queued/running. The protected source RDS is read by
  its fixed NFS path and checksum, not duplicated into the package.
- Capture `sbatch --parsable` output, exact submission command, job ID,
  absolute scheduler log paths, submission/start/end timestamps, exit status,
  session/package versions, and stage timings. A failed submission/run is
  recorded as failed; do not create a success-shaped results package.
- A later plotting-only job also gets its own timestamped carbon-copy bundle
  and job record linked to the source analysis run. It consumes saved gene
  tables/trajectories, does not rerun DE, and does not overwrite reviewed plots.

Carbon-copy packages are created only by the explicit submission wrapper.

### Read-only inspection findings

The saved RDS was loaded using R/4.4.3 and installed URD 1.1.1. It is an `URD`
object with named `count.data` and `logupx.data` sparse matrices, each
15,239 genes by 20,049 cells; `meta` has 20,049 rows and 55 columns. Named
pseudotime `div90_jia_rootscore_root` is present, alongside seven tip-walk
columns. Tree segment membership, joins, pseudotime limits, segment layout,
and cell coordinates are present. No dimname reconstruction is needed for
these two expression matrices. Runtime implementation must still validate
cell/gene identity and order before testing.

The actual joins identify four branchpoints, not six assumed binary splits:

| Proposed branchpoint ID | Parent segment | Split pseudotime | Immediate daughters | Descendant terminal identities |
| --- | --- | --- | --- | --- |
| BP13 | 13 | 0 | 6 versus 12 | Astrocytes versus the subtree containing tips 1–5 and OPC |
| BP12 | 12 | 0.173664130168718 | 7 versus 11 | OPC versus neuronal tips 1–5 |
| BP11 | 11 | 0.235455231461149 | 5, 9, 3 | NR2F1/NR2F2; subtree containing tips 1, 2, 4; CRABP1/ANGPT2 precursor |
| BP9 | 9 | 0.272363130443725 | 1, 2, 4 | LHX8/ISL1 state 1; LHX8/ISL1 state 2; LHX6/NFIA/EPHA5/MEF2C cortical |

These are native run identities, not revised biological annotations. Preserve
them even where final-figure display labels differ. Reconstruct descendants
from the object rather than hard-code this table into the analysis.

Native segment membership counts are 1:1002, 2:1027, 3:2048, 4:343, 5:463,
6:1727, 7:985, 9:310, 11:301, 12:165, 13:0. Do not assume all 20,049 cells
are assigned to the plotted tree. Enumerate unique membership, duplicates,
unassigned cells, finite pseudotime, and sample representation in the preflight.
Segment 13 has zero width and no cells: an upstream expression baseline at
BP13 cannot be estimated. Do not invent one or shift the split to obtain it.

Direct inspection found 8,371 unique segment-assigned cells and 11,678 cells
without segment membership. The named and tree pseudotimes agree exactly
(maximum absolute difference zero), both ranging from 0 to approximately
0.5531137. All six `orig.ident` labels (`10496-MW-1` through `10496-MW-6`)
occur in every nonempty segment, but counts are markedly unequal; a sample
label's presence does not establish independent biological replication.

### Native URD methods and their roles

The installed URD 1.1.1 signatures and marker-function bodies were inspected.
Use installed behavior as authoritative and record package/session versions.

- `segChildren()` and `segChildrenAll()` identify daughters and descendants;
  verify these against `tree$segment.joins` and terminal mapping.
- `markersAUCPR(cells.1=..., cells.2=..., genes.use=...)` is the native
  marker-specificity comparison for explicitly selected cells. Call both
  directions. Report AUCPR, its random-classifier baseline/ratio, native
  expression fold change, and detection fractions. AUCPR is not a p-value.
- `markersBinom()` supplies native detection-frequency tests and p-values/FDR
  for the same cell sets. Its `log.effect` is the natural log of the ratio of
  detection fractions, not mean-expression log2 fold change. Label both
  quantities separately. This does not test mean expression among positive
  cells; mark mean-only candidates without detection support accordingly.
- `deTTest()` was inspected as a possible complementary native mean-expression
  test. It calls Welch's `t.test()` but has no built-in FDR correction, uses a
  ratio of input means for its effect field, and can fail on constant data.
  It is not selected as the primary significance method in this proposal.
  Adding a guarded mean-expression test with a complete p-value family is a
  specific implementation-review decision; detection significance must not be
  relabeled as significance of mean-expression change.
- `aucprTestAlongTree(..., report.debug=TRUE)` can provide a supplementary
  whole-lineage comparison and raw `marker.chain` entries. It traverses
  segment-defined ranges; it is not the primary engine for arbitrary nested
  windows around a split. `binomTestAlongTree()` also summarizes along paths
  rather than exposing the desired window-specific comparisons. Do not
  substitute their pooled summaries for branchpoint tables.
- `pseudotimeMovingWindow()` supports descriptive trajectories;
  `geneCascadeProcess()` and `geneCascadeHeatmap()` are optional native
  fitting/visualization tools for selected genes. Cascade onset is expression
  onset within a lineage, not the onset of a between-lineage difference.
  Failed/poor fits must not replace the observed binned data.

Native helper edge cases must be addressed in the implementation review:
empty or single-gene results, absent positive cells, nonfinite fold changes,
and `markersBinom()` output naming with explicit cell lists. Provide stable
population labels (or an in-memory grouping) and map output columns by
validated semantics. Never patch or overwrite the installed URD package.

Primary source references:
[URD differential-expression implementation](https://github.com/farrellja/URD/blob/master/R/diff-exp.R)
and [URD gene-cascade implementation](https://github.com/farrellja/URD/blob/master/R/gene-cascade.R).

### Define comparisons from actual daughters

At each parent, compare every unordered pair of its immediate daughters and
retain native AUCPR results in both directions. This produces eight unordered
pairwise contrasts across the four branchpoints (1 + 1 + 3 + 3), before window
variants. At a three-way split, a daughter-specific gene must beat both sibling
branches to receive an "all-siblings" label; genes beating only one sibling
remain explicitly labeled pair-specific. Do not impose an artificial binary
ordering on either three-way split.

For local tests, start with exclusive immediate-daughter segment cells. Audit
and exclude overlapping memberships between contrasted groups. Shared parent
cells provide upstream context only; never assign those cells a future fate
from a terminal label. For downstream trajectories, include descendants but
separate the paths when a daughter branches again. Never average later
descendants into one curve and interpret a changing mixture as specification.
Preserve unassigned cells as unassigned; do not recruit them using expression
of candidate genes. All cell selections must be exported with their rationale.

### Multiple windows and pseudotime support

Use the fixed named pseudotime on the source object's original scale. First
check that tree join/segment times and the named pseudotime are on the same
scale; disagreement is a blocking preflight issue, not permission to reflood.
Let `R` be the finite named-pseudotime maximum minus minimum across the source
object. Proposed nested half-widths are `0.02*R`, `0.05*R`, and `0.10*R`.
These are prespecified sensitivity settings, not estimates of split confidence.
For the inspected range these widths are approximately 0.01106, 0.02766, and
0.05531 pseudotime units; the implementation must calculate from the input
instead of hard-coding rounded values.

For a split at `b` and each candidate half-width `w`:

1. Show upstream context `[b-w,b)` from the native parent segment, clipped to
   its available interval. Show daughter expression in `[b,b+w]`.
2. Restrict a local pairwise test to the shared supported downstream interval
   of both daughters, before either daughter's next split or terminal end.
   Use identical time limits for both groups. Save requested versus realized
   limits, clipping reasons, counts, occupied bins, and pseudotime quantiles.
3. Bin this common interval into five equal-width pseudotime bins. Proposed
   minimum support: 30 cells per daughter overall, at least three jointly
   occupied bins, and at least five cells per daughter in each retained bin.
   Test on cells from jointly supported bins; expose the native-all-cells
   result and a deterministic seed-7, equal-count-per-bin sensitivity result.
   Report sample composition before and after balancing. If support fails,
   mark the comparison/window insufficient; do not silently widen or move it.
4. If clipping makes requested windows identical, test once and record that
   identity. They are not independent confirmations. Keep all four branchpoints
   in the report even if a particular contrast has no supported local window.

Also examine split-location sensitivity at `b-0.02*R`, `b`, and `b+0.02*R`
using original memberships and the same clipping/support rules. These are
alternative anchors for expression summaries, not modified tree branchpoints
or a biological confidence interval. Upstream cells cannot acquire a known
daughter identity merely because an anchor moves earlier. Never extrapolate
past available support. BP13 has no upstream parent data and must be flagged.

Use disjoint downstream bins after the widest supported local window for the
later-expression analysis. Extend along each actual daughter-to-tip path,
retain original pseudotime and later branchpoint markers, and stop pairwise
interpretation where common pseudotime support ends. If no later matched
interval exists, classify timing as unresolved. Quantile-based moving bins
may be a separately labeled descriptive sensitivity view, never a substitute
for the matched original-pseudotime comparisons.

### Gene universe, significance, and ranking

Test expressed genes from the full saved expression matrix, not just the
3,000 variable genes or the markers used to name tips. Proposed expression
filter: detection in at least 10% of one daughter for a contrast/window;
record the entire tested universe before outcome-based fold-change filtering.
Retain complete gene-level results, not only the plotted top genes.

Proposed display/candidate thresholds for review: absolute mean-expression
log2 fold change at least 0.25, detection-fraction difference at least 0.10,
AUCPR ratio at least 1.2 in the enriched direction, and detection-test FDR
at most 0.05. An AUCPR filter must never be described as significance.
Call native binomial tests with permissive return filters, retain raw p-values
for all eligible genes/directions, and recompute BH over the full intended
gene family before applying candidate thresholds. Explicitly distinguish
native within-comparison FDR from a study-wide BH family spanning all tested
branchpoint/pair/window/anchor/direction hypotheses. Do not report a gene as
globally significant merely because one of many windows passed local FDR.
Any native helper that prevents complete p-value recovery requires an explicit
reviewed adjustment; do not silently correct only an AUCPR-selected subset.

Implementation resolution of native helper constraints: `markersBinom()` is
called on every eligible gene with unequal detection, with stable in-memory
population labels, minimal effect filtering and permissive p/FDR return limits.
It returns the enriched direction's raw p-value. Multiply that p-value by two
(capped at one) to conservatively account for choosing the enriched direction,
then calculate local/global BH before selecting genes for AUCPR displays.
Equal-detection genes have no enriched detection direction: retain raw p as
missing, explicitly flag that status, and use one as a conservative BH input.
Mirrored A/B reporting rows share this one corrected gene test; duplicate cell
selections also share one test. Global BH spans all unique selection/gene tests,
including shifted windows, balanced variants and later-path comparisons.
No AUCPR preselection defines the detection p-value family. Single-gene native
dimension-dropping is handled with a temporary duplicated padding row, removed
after the call; native p-values are retained and BH is recomputed externally.
Native FDR values from padded calls are flagged, not used for candidate calls.

All native p-values are exploratory cell-level detection tests on the same
data that produced the tree. They are not donor-level validation and do not
account for tree-selection uncertainty. Inspect `orig.ident` representation;
show per-sample effects and support where possible. A replicate-level count
model would be an additional, separately reviewed method if independent
biological replication is adequate. No confirmatory claims from cell counts.

Classify genes by observed support:

- **Early-associated:** consistent direction/effect in at least two distinct
  supported local widths, including the smallest supported width, with
  candidate thresholds and significance flags reported independently.
- **Early and sustained:** early-associated plus consistent effects in at
  least two supported later bins.
- **Late-emerging:** supported effects in at least two later bins with no
  threshold-size early effect; report early estimates/uncertainty. Failure of
  an early significance test alone is not evidence of absence.
- **Transient, window-sensitive, or unresolved:** inconsistent direction,
  limited support, anchor sensitivity, or inability to separate timings.

Report first supported between-daughter difference as a bounded pseudotime
interval, not a precise commitment time. Rank each enriched daughter using
distinct-window robustness first, then AUCPR ratio and mean-expression effect;
display detection differences, FDR, sample consistency, and timing class.
Preserve a separate late list so local prioritization does not hide maturation.

### Candidate transcription factors

Annotate the tested gene universe with a versioned human DNA-binding TF list
from the [Lambert et al. author resource](https://humantfs.ccbr.utoronto.ca/download.php).
Record release/download URL, retrieval date, checksum, species, TF-inclusion
field, symbol/Ensembl mappings, and unmapped/ambiguous identifiers. Do not infer
TF status from gene-name prefixes. No suitable local TF list was identified
in the targeted metadata/configuration search. The submission wrapper acquires
the author v1.01 database, checks its schema, selects `Is TF? == Yes`, and
stores a content-addressed reference plus the exact run copy and manifest.

The reference data and mapping outputs belong under Turbo `reference/` and
the new run's `provenance/`, respectively; only paths/settings belong in repo
configuration. Subset the same DE results for separate TF rankings and TF
heatmaps, keeping the original all-gene FDR family. Do not retest TFs with a
smaller correction family to improve significance. Summarize the favored
daughter, specificity versus all siblings, timing, window stability, and
supporting expression. TF expression alone does not establish TF activity,
binding, target genes, or causality. Non-TF regulators may remain annotated
in the all-gene results but must not be mislabeled DNA-binding TFs.

### Planned outputs and plotting conventions

All generated files go to the independent NFS runtime destination above:

- **Single entry-point CSV: `<run_dir>/DIV90_URD_gene_results.csv`.** Use plain,
  uncompressed UTF-8 CSV with proper quoting, one header row and stable column
  names, for opening directly in a spreadsheet. One row represents one gene,
  branchpoint, daughter-pair/direction, window, anchor and analysis variant.
  Preserve all tested genes and nonsignificant results, not only top hits.
  Include native daughter IDs and readable descendant labels; split time and
  realized window; expression means/log2FC; detection fractions/difference;
  AUCPR/baseline/ratio; p-values and local/global FDR; TF flag/family/reference;
  early/late class and robustness counts; all-gene and TF ranks; sample/support
  flags; run ID, generated-at timestamp and method versions; and relative paths
  to the corresponding tree, trajectory and heatmap files. Mark untestable
  comparisons explicitly with status rows and missing statistics, never zeros.
  Supplementary per-cell membership and per-bin data remain in detailed tables;
  the master CSV contains their paths and scientific summaries without
  multiplying every gene row by every cell/bin. Record its row count/file size
  and spreadsheet row-limit warning in `analysis_summary.txt` if applicable.
- **Spreadsheet summary: `<run_dir>/DIV90_URD_gene_summary.csv`.** One ranked
  representative per gene/branchpoint/daughter contrast, including TF status,
  timing, distinct-window robustness, FDR and paths to the plots/full results.
  This smaller CSV is the convenient browsing entry point if the complete
  per-window master exceeds a worksheet's row limit; all tests remain retained.
- **Visual entry point: `<run_dir>/DIV90_URD_gene_overview.pdf`.** Start with
  the annotated unchanged tree and a daughter-label key. Give each of the four
  branchpoints a readable section showing leading early-associated genes,
  separately ranked TFs, observed daughter-expression trajectories, early/late
  heatmaps and window robustness. Include insufficient-support splits rather
  than dropping them. Mark inferred split positions and tested intervals on
  pseudotime plots; explain shared upstream cells. Put run ID, generation
  timestamp and source-tree identity in a caption/footer without crowding gene
  panels. Export the individual panels as PNG/PDF/editable SVG in `plots/`.
- `tables/branchpoints.tsv`: parent/daughter IDs, descendant tip labels,
  inferred split times, and root-boundary flags.
- `tables/window_support.tsv` and `tables/comparison_cells.tsv.gz`: requested
  and realized windows, matched bins, sample composition, exclusions, and
  exact cell IDs; include insufficient comparisons.
- `tables/branchpoint_gene_statistics.tsv.gz`: every tested gene with
  branchpoint/pair/direction/window/anchor, cell counts, expression means,
  mean-expression log2FC, detection fractions/difference, AUCPR/baseline/ratio,
  binomial p-values, native/local/global FDR, and diagnostic flags.
- `tables/ranked_divergence_genes.tsv`, `tables/ranked_late_genes.tsv`,
  `tables/ranked_candidate_tfs.tsv`, and `tables/timing_window_stability.tsv`.
- `tables/expression_trajectories.tsv.gz`: observed binned expression,
  detection, pseudotime bounds, counts and per-sample summaries; optional fit
  diagnostics and bootstrap intervals labeled by their sampling unit.
- `plots/urd_tree_branchpoints.{png,pdf,svg}`: unchanged native layout,
  explicit BP13/BP12/BP11/BP9 labels and daughter-to-terminal key.
- `plots/<branchpoint>/`: top-gene and TF trajectories, matched-pseudotime
  heatmaps, and window-sensitivity panels; all four branchpoints represented.
- `provenance/`: source hash, configuration, exact command, session/package
  versions, reference manifest, stage status/timing, and output manifest.
  `analysis_summary.txt` describes supported candidate programs and limitations.

Render PNG/PDF/SVG using existing R/ggplot conventions. Use 300 dpi PNGs for
the URD panels, Cairo PDF and `svglite` SVG with editable text, consistent
daughter colors across tree/curves/heatmaps, gene labels, cell counts, and
explicit pseudotime axes. Reuse final Figure 4 display mappings only as visual
annotations; native tip labels remain in tables. Heatmaps use a gene-wise scale
shared across compared daughters (never independently scaled per daughter),
with the transformation and clipping documented. Show unsmoothed binned data
with fitted curves; share the upstream trunk without implying daughter-specific
fates before the split. Plot failures must not trigger any URD recomputation.

### Implementation review checklist before any execution

Review the new script, configuration, and standalone template together, with
exact line references for topology extraction, membership selection, window
clipping/matching, native tests, p-value families, TF annotation, timing rules,
and plotting. Verify helper edge cases with small synthetic inputs rather than
running DE on this biological object. Confirm source and output paths cannot
overlap and that existing output files will not be silently overwritten.
Record the source RDS hash and compare after the eventual run. The existing
tree, topology, pseudotime, identities, source checkpoints, and figure packages
must remain unchanged. The user's explicit execution instruction now authorizes
submission after technical validation; it does not authorize source reruns.

### Execution implementation record

Implementation prepared at `2026-10-21T14:20:23-04:00` (`America/Detroit`).
No biological DE was run on the login node. Native-helper tests use synthetic
matrices; the small end-to-end validation uses a synthetic seven-tip tree and
saves its clearly identified validation artifacts only on NFS.

Code responsibilities:

- `scripts/40_div90_urd_branchpoint_genes.R`: window support/matching;
  unmodified native URD tests; actual topology extraction; exclusive local
  memberships and separate later daughter-to-tip paths; multiplicity handling;
  timing/TF rankings; consolidated CSV; observed expression trajectories;
  annotated tree, TF/gene panels, heatmaps and multi-page PDF overview.
- `config/div90_urd_branchpoint_genes.tsv`: all prespecified thresholds,
  matching windows, support limits, plot gene counts and random seed.
- `scripts/40_submit_div90_urd_branchpoint_genes.py`: timestamped versioned
  run ID; source/mapping checksums; code/config/handoff/reference carbon copies;
  exact submitted job; Git state/diff; `sbatch --parsable` and durable job ID.
- `slurm_templates/64_div90_urd_branchpoint_genes.sbatch.template`: frozen
  hash verification, module/library environment, compute-node execution,
  scheduler start/end/exit records, source hash recheck and output manifest.

The implementation uses observed binned trajectories rather than optional
impulse fitting. Genes can be visualized even if significance/support fails;
the CSV retains `candidate`, `global_fdr`, timing and support flags. Where no
local contrast is supported, later-path genes may still be plotted, explicitly
remaining unresolved for early-divergence timing. No absent local evidence is
presented as biological absence or fate commitment.

Validation before biological submission:

- R/Python parsing, shell syntax and `git diff --check`: passed.
- Synthetic native-method checks: matched-window support, detection direction,
  equal-detection genes, empty directional AUCPR results, and padded
  single-gene p-value recovery: passed.
- Synthetic seven-tip end-to-end check: passed at
  `2026-10-21T14:22:59-04:00`, including a 5,600-row CSV and the tree,
  gene/TF trajectories, heatmaps, sample curves and overview PDF. Validation
  artifacts: `PROJECT_ROOT/results/div90_urd_branchpoint_genes/synthetic_validation_20261021_142031/`.
  The synthetic plot was visually inspected. These are validation data, not
  DIV90 biological findings. No biological inference was run on the login node.

### Biological job submission attempt 1: 2026-10-21 14:23 EDT

- Slurm job: **64478969**, submitted at `2026-10-21T14:23:17-04:00`.
- Run ID: `div90_final_7tips_branchpoint_genes_v1_20261021_142315_EDT_2e9ba1e`.
- Run directory: `PROJECT_ROOT/results/div90_urd_branchpoint_genes/div90_final_7tips_branchpoint_genes_v1_20261021_142315_EDT_2e9ba1e/`.
- Submitted carbon copy: `PROJECT_ROOT/jobs/div90_final_7tips_branchpoint_genes_v1_20261021_142315_EDT_2e9ba1e.sbatch`;
  byte-identical copy in the run's `provenance/submitted.sbatch`.
- Scheduler log: `PROJECT_ROOT/logs/64_div90_urd_branchpoint_genes_64478969.log`.
- Durable scheduler record: `<run_dir>/provenance/submission.json` and
  `provenance/job_id.txt`; prepared/submitted timestamps include EDT offset
  and UTC. The frozen handoff records pre-submission state; this live section
  adds the returned job ID without modifying any frozen assets.
- TF reference: Lambert HumanTFs v1.01, 2,765 annotated entries including
  1,639 TF entries; the run stores the exact reference and its SHA256.
- Source RDS SHA256: `30322446a2f1013496f6bcbfdd14924b8f5ad931406a333a6c46ae935fa9f984`.
  Source/mapping and frozen-file checks are performed before computation;
  source checks repeat before the job writes `COMPLETED.txt`.
- State at submission verification: **PENDING (Priority)**, with no predicted
  start time supplied by Slurm. This is not a completed analysis. Use `squeue`
  and `sacct` for live state; `ANALYSIS_COMPLETE.txt` alone is not final success.
- Expected entry points after successful execution:
  `DIV90_URD_gene_summary.csv`, `DIV90_URD_gene_results.csv`,
  `DIV90_URD_gene_overview.pdf`, `analysis_summary.txt` and
  `provenance/output_manifest.tsv` inside the run directory.

Attempt 1 outcome: **FAILED 1:0 after 26 seconds** at
`2026-10-21T14:24:05-04:00`, before differential-expression comparisons. A
negative-anchor/narrow-width combination produced a floating-point interval
only about one numerical rounding unit wide; `cut()` rejected repeated bin
edges. The source and frozen hashes had passed. The implementation now treats
intervals narrower than 64 machine epsilons on this scale as unsupported,
with a regression check. No window is widened or moved. Preserve this failed
package and its log; the corrected attempt gets a fresh timestamp/carbon copy.

Monitor without resubmitting or regenerating upstream URD:

```bash
squeue -j 64478969
sacct -j 64478969 --format=JobID,State,ExitCode,Elapsed,MaxRSS
tail -n 40 /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/64_div90_urd_branchpoint_genes_64478969.log
```

### Corrected biological submission: 2026-10-21 14:25 EDT

- Current job: **64479256**, submitted at `2026-10-21T14:25:16-04:00`.
- Run ID: `div90_final_7tips_branchpoint_genes_v1_20261021_142513_EDT_2e9ba1e`.
- Run directory: `PROJECT_ROOT/results/div90_urd_branchpoint_genes/div90_final_7tips_branchpoint_genes_v1_20261021_142513_EDT_2e9ba1e/`.
- Job carbon copy: `PROJECT_ROOT/jobs/div90_final_7tips_branchpoint_genes_v1_20261021_142513_EDT_2e9ba1e.sbatch`.
- Scheduler log: `PROJECT_ROOT/logs/64_div90_urd_branchpoint_genes_64479256.log`.
- Includes the zero-width floating-point regression fix; native-method
  synthetic tests passed again. Source/reference identities and configured
  scientific windows are unchanged. The earlier failed package is preserved.
- Initial scheduler state: PENDING (Priority); inspect live Slurm state and
  this package's completion/exit records rather than assume completion.

At `2026-10-21T14:26:15-04:00`, job 64479256 was RUNNING on `gl3183`,
had passed all frozen/source checksum checks, loaded the unchanged saved
tree, and completed topology/window setup: 214 supported requests collapsed
to 182 unique native comparisons. Gene statistics then started in four
workers. This is a progress observation, not completion.

Support audit: BP13 has no supported local comparison under the prespecified
rules (zero-time parent boundary and insufficient nearby daughter cells);
24 later-path requests are supported. BP12 has 10 local and 40 later supported
requests, BP11 has 48 local and 36 later, and BP9 has 48 local and 8 later.
These counts include balanced variants and shifted anchors, not independent
biological replications. BP13 early-divergence timing must remain unresolved;
later markers cannot be relabeled as early fate regulators.

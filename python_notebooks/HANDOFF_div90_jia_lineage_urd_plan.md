# DIV90 Jia-Lineage-Driven URD Plan

Date: 2026-06-11

Status: DIV90 v4 glia-tip smoke completed. This is the final inclusion-logic swing before scaling: stressed cells remain excluded, clusters `3`/`11` remain retained non-tip candidates, astrocytes are combined as one terminal tip, and OPCs are their own terminal tip.

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
scripts/14_div30_first_urd.R
scripts/15_div30_urd_lineage_decision_report.R
scripts/16_div30_urd_build_lineage_tree.R
scripts/17_div30_urd_finalize_lineage_tree_report.R
scripts/25_div30_urd_jia_fig_s11_marker_validation.R
scripts/27_div90_urd_project_candidate_lineage_markers.R
```

Compatibility note: the DIV90 exporter writes input-bundle files with the legacy `div30_first_urd_*` names because `scripts/14_div30_first_urd.R` consumes that bundle format. The run directory, metadata, root column, pseudotime column, reports, and manifest identify this as DIV90.

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

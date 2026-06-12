# DIV90 Jia-Lineage-Driven URD Plan

Date: 2026-06-11

Status: first DIV90 Jia-lineage smoke run submitted.

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
```

Compatibility note: the DIV90 exporter writes input-bundle files with the legacy `div30_first_urd_*` names because `scripts/14_div30_first_urd.R` consumes that bundle format. The run directory, metadata, root column, pseudotime column, reports, and manifest identify this as DIV90.

## DIV90 Cluster Mapping

| Cluster | Exact metadata name | n cells | First URD role |
|---:|---|---:|---|
| 0 | `0 - MGE Striatal/GP Fated` | 3,601 | LHX8/ISL1-like tip |
| 1 | `1 - SST+, NPY +, Cortical Fated` | 3,548 | LHX6/NFIA-like tip |
| 2 | `2 - CRABP1+/PV Precursors` | 3,503 | EPHA5/MEF2C-like tip candidate; check CRABP1/ANGPT2 after scoring |
| 3 | `3 - PV precursors/Migrating cells/Cortical-fated` | 2,283 | EPHA5/MEF2C-like tip |
| 4 | `4 - Pre-Astrocytes/Astrocytes 1` | 1,987 | exclude/track separately |
| 5 | `5 - LHX8+ vMGE GABergic Striatal/GP fated 1` | 1,924 | LHX8/ISL1-like tip |
| 6 | `6 - Stressed Cells` | 1,226 | exclude/track separately |
| 7 | `7 - Stressed Cells` | 1,063 | exclude/track separately |
| 8 | `8 - LHX8+ vMGE GABergic Striatal/GP fated 2` | 922 | LHX8/ISL1-like tip |
| 9 | `9 - Pre-OPCs/OPCs` | 915 | exclude/track separately |
| 10 | `10 - Pre-Astrocytes/Astrocytes 2` | 583 | exclude/track separately |
| 11 | `11 - PV Precursors` | 425 | EPHA5/MEF2C-like tip |
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
| `tip_epha5_mef2c` | EPHA5/MEF2C-like lineage | `2`, `3`, `11` | `MEF2C`, `MAFB`, `ETV1`, `ERBB4` | cortical-fated GABAergic/PV precursor-associated lineage |
| `tip_lhx6_nfia` | LHX6/NFIA-like lineage | `1` | `LHX6`, `SST`, `NPY`, `ERBB4`, `CXCR4`, `ARX` | cortical interneuron lineage |

Do not define these tips as broad SST/PV/MGE categories. Those labels can be retained as metadata, but the URD tip construction should be Jia-lineage driven.

## First Smoke Run

Retain:

```text
0, 1, 2, 3, 5, 8, 11, 12
```

Exclude or track separately for the first neuronal lineage smoke test:

```text
4, 6, 7, 9, 10
```

Proposed retained roles:

```text
root pool: 12
tip_lhx8_isl1: 0, 5, 8
tip_epha5_mef2c: 2, 3, 11
tip_lhx6_nfia: 1
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

### Submitted Smoke Run

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

Pre-submit exporter smoke check:

```text
/tmp/div90_jia_urd_export_test/
```

The exporter smoke check selected 800 cells with all 358 cluster-12 dividing cells retained, and selected 8 RootScore root cells. R read the exported `urd_root_candidate` column as logical, so the root column format is valid for URD.

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

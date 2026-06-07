# Handoff: Shi Reference DIV30 Label Transfer

Standalone workflow for mapping DIV30 AnnData cells to Shi et al. reference
cell-level labels.

This is not Notebook 02. It is a separate reference-mapping workflow.

## Core Rule

DIV30 `seurat_clusters` are not used for prediction.

Predictions are made per DIV30 cell using expression similarity to Shi reference
cells. Existing DIV30 `seurat_clusters` are used only afterward for summaries,
stacked bars, heatmaps, and interpretation.

## Current Run Target

Default run label:

```text
shi_reference_div30_label_transfer_v2
```

Default run directory:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_reference_div30_label_transfer/shi_reference_div30_label_transfer_v2
```

Executed notebook path:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_reference_div30_label_transfer/executed/shi_reference_div30_label_transfer.shi_reference_div30_label_transfer_v2.executed.ipynb
```

## Code Added

Reusable module:

```text
python_notebooks/src/mge_organoid_python/shi_label_transfer.py
```

Notebook:

```text
python_notebooks/notebooks/shi_reference_div30_label_transfer.ipynb
```

Slurm template:

```text
slurm_templates/20_execute_shi_reference_div30_label_transfer.sbatch.template
```

## Inputs

DIV30 query AnnData:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/varela_div30.h5ad
```

Shi reference AnnData:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/shi_2019_paper_qc.h5ad
```

Shi Table S2 cell labels:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/reference/shi_2021_tables_s2_to_s9/science.abj6641_table_s2.xlsx
```

Table S2 provides cell-level `Major types` labels. The workflow joins those
labels to the cached Shi AnnData by normalized barcode.

Observed Table S2 labels during setup:

```text
MGE
LGE
CGE
progenitor
Excitatory IPC
Excitatory neuron
Thalamic neurons
Microglia
OPC
Endothelial
```

## Comparisons

The workflow runs two independent per-cell label-transfer comparisons.

### A. Full Relevant Developmental Reference

As of `shi_reference_div30_label_transfer_v2`, the default full reference
includes all 10 observed Shi Table S2 major-type labels:

```text
MGE
LGE
CGE
progenitor
Excitatory IPC
Excitatory neuron
Thalamic neurons
Microglia
OPC
Endothelial
```

### B. Restricted MGE/LGE/CGE Reference

Default included labels:

```text
MGE, LGE, CGE
```

This is intentionally compared against the full-reference result to identify
cells whose MGE/LGE/CGE label is potentially being forced by removing other
developmental states.

## Legacy v2 Method

This describes the completed `shi_reference_div30_label_transfer_v2` run only.
Do not use this Scanpy-side kNN vote-fraction method for the next scoring
target.

1. Load DIV30 query AnnData and Shi reference AnnData.
2. Read Shi Table S2 `Major types` labels.
3. Join Table S2 labels to Shi reference cells by barcode.
4. Harmonize query/reference genes with exact and unambiguous case-insensitive
   gene matching.
5. For each comparison:
   - subset Shi reference cells to the relevant label set
   - select top variable genes from the Shi reference among shared genes
   - fit sparse StandardScaler and TruncatedSVD on the Shi reference
   - project DIV30 cells into the same SVD space
   - run kNN from each DIV30 cell to Shi reference cells
   - predict labels by neighbor-vote at the single-cell level
   - predict Shi gestational week from the same neighbor set
6. Add predictions to DIV30 `adata.obs`.
7. Summarize predictions by DIV30 `seurat_clusters`.
8. Summarize Shi gestational-week predictions by DIV30 sample.

Prediction score is the fraction of k nearest Shi reference neighbors voting for
the winning label. Uncertainty is `1 - prediction_score`. Entropy is normalized
neighbor-label entropy. The gestational-week score uses the same vote-fraction
logic, but over Shi neighbor `week_label`; mean/median neighbor week use the
numeric Shi gestational week.

## Next Target: Seurat Prediction Scores

Default next run label:

```text
shi_reference_div30_seurat_label_transfer_v1
```

This next target should replace Scanpy-side prediction scoring with Seurat
anchor-based label transfer. Scanpy/AnnData should remain the plotting and
summary environment after Seurat predictions are exported and joined back onto
DIV30 cells.

Core rule for the next target:

```text
Do not compute prediction_score from Scanpy kNN neighbor-vote fractions.
Use Seurat's transferred prediction scores as the prediction-score source of
truth.
```

Recommended first runnable target:

1. Run one Seurat label-transfer comparison using the full Shi reference label
   set.
2. Export per-cell Seurat predictions and the full label-score matrix.
3. Import those tables into the existing DIV30 AnnData.
4. Generate a small smoke-test plot set:
   - UMAP colored by predicted Shi label
   - UMAP colored by maximum Seurat prediction score
   - UMAP panels for MGE, LGE, and CGE score columns
   - ridge/density-style score distributions by DIV30 sample
   - stacked bar of predicted labels by DIV30 `seurat_clusters`
5. Only after that runs, add the MGE/LGE/CGE-restricted comparison and
   gestational-week transfer.

Expected new output roots:

```text
RUN_DIR/seurat/
RUN_DIR/tables/
RUN_DIR/plots/
RUN_DIR/h5ad/
```

Expected Seurat-side files:

```text
seurat/div30_shi_seurat_full_predictions.tsv.gz
seurat/div30_shi_seurat_full_prediction_scores.tsv.gz
seurat/div30_shi_seurat_full_transfer_diagnostics.tsv
```

Expected AnnData-side files:

```text
tables/div30_shi_seurat_label_transfer_obs.tsv.gz
tables/div30_shi_seurat_label_transfer_label_scores_long.tsv.gz
h5ad/div30_shi_seurat_label_transfer_predictions.h5ad
```

Required new `adata.obs` columns for the first target:

```text
shi_seurat_full_predicted_shi_label
shi_seurat_full_prediction_score
shi_seurat_full_uncertainty_score
shi_seurat_full_broad_region_class
shi_seurat_full_developmental_class
```

Required per-label score columns in `adata.obs` for plotting:

```text
shi_seurat_full_prediction_score_MGE
shi_seurat_full_prediction_score_LGE
shi_seurat_full_prediction_score_CGE
shi_seurat_full_prediction_score_progenitor
shi_seurat_full_prediction_score_Excitatory_IPC
shi_seurat_full_prediction_score_Excitatory_neuron
shi_seurat_full_prediction_score_Thalamic_neurons
shi_seurat_full_prediction_score_Microglia
shi_seurat_full_prediction_score_OPC
shi_seurat_full_prediction_score_Endothelial
```

Column-name sanitization rule: preserve the original Shi label in a companion
long table, but replace spaces and punctuation with underscores for wide
`adata.obs` score columns.

Implementation added for the first smoke target:

```text
python_notebooks/scripts/seurat_shi_label_transfer_export.R
python_notebooks/scripts/run_shi_seurat_label_transfer_smoke.py
slurm_templates/21_run_shi_seurat_label_transfer_smoke.sbatch.template
```

Submitted smoke run:

```text
Slurm job: 51405243
Job script:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/21_run_shi_seurat_label_transfer_smoke.sbatch
Expected run directory:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_reference_div30_seurat_label_transfer/shi_reference_div30_seurat_label_transfer_v1
Logs:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/shi-seurat-label-transfer-shi-seurat-xfer-51405243.out
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/shi-seurat-label-transfer-shi-seurat-xfer-51405243.err
```

Completed Seurat transfer/plot runs:

```text
DIV30 refreshed plot run:
  Slurm job: 51410249
  State: COMPLETED
  ExitCode: 0:0
  Output root:
    /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_reference_div30_seurat_label_transfer/shi_reference_div30_seurat_label_transfer_v1
  Completion:
    n_query_cells = 90631
    n_score_columns = 10
    n_plots = 8
  Notes:
    - reused existing Seurat prediction tables
    - regenerated plots with non-squished cluster UMAP panel layout
    - subtype overlay panels now keep all Shi labels, including zero-count CGE/Endothelial

DIV30 unexpected-label marker summary:
  Slurm job: 51410250
  State: COMPLETED
  ExitCode: 0:0
  Output table:
    tables/markers/shi_seurat_unexpected_label_marker_summary.tsv
  Result summary:
    Thalamic neurons: 556 cells, mean max score 0.604
    Excitatory neuron: 54 cells, mean max score 0.534
    Excitatory IPC: 1 cell, mean max score 0.473
  Interpretation note:
    These are low-confidence Seurat edge calls. The "Excitatory neuron" cells
    are high for SST/GAD1, so do not interpret that label literally without
    checking markers and score thresholds.

DIV90 Seurat transfer:
  Slurm job: 51410263
  State: COMPLETED
  ExitCode: 0:0
  Elapsed: 00:06:08
  MaxRSS: 22621688K
  Output root:
    /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_reference_div90_seurat_label_transfer/shi_reference_div90_seurat_label_transfer_v1
  Query inputs:
    AnnData: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/varela_div90.h5ad
    Seurat: /nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/umap_props_output/clustered_day90_with_cluster_names_2.rds
  Completion:
    n_query_cells = 22338
    n_score_columns = 10
    n_plots = 8
  Diagnostics:
    reference_cells_labelled = 55704
    shared_features_used = 2542
    anchors = 3249
```

Completed Seurat plot inventory for each timepoint:

```text
<slug>_umap_shi_seurat_full_predicted_shi_label.png
<slug>_umap_shi_seurat_full_predicted_shi_label_with_subtype_overlays.png
<slug>_umap_shi_seurat_full_prediction_score.png
<slug>_umap_shi_seurat_full_mge_lge_cge_score_panel.png
<slug>_shi_seurat_full_prediction_score_density_by_sample.png
<slug>_shi_seurat_full_expected_shi_gw_ridge_by_sample.png
<slug>_shi_seurat_full_shi_label_stacked_bar_by_seurat_clusters.png
<slug>_shi_seurat_full_shi_label_stacked_bar_by_seurat_clusters_with_cluster_umaps.png
```

Rsync completed Seurat plot folders from a Mac:

```bash
mkdir -p /Users/ecrespo/Downloads/shi_reference_div30_seurat_label_transfer_v1/plots
mkdir -p /Users/ecrespo/Downloads/shi_reference_div90_seurat_label_transfer_v1/plots

rsync -avh --progress \
  elcrespo@greatlakes.arc-ts.umich.edu:/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_reference_div30_seurat_label_transfer/shi_reference_div30_seurat_label_transfer_v1/plots/ \
  /Users/ecrespo/Downloads/shi_reference_div30_seurat_label_transfer_v1/plots/

rsync -avh --progress \
  elcrespo@greatlakes.arc-ts.umich.edu:/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_reference_div90_seurat_label_transfer/shi_reference_div90_seurat_label_transfer_v1/plots/ \
  /Users/ecrespo/Downloads/shi_reference_div90_seurat_label_transfer_v1/plots/
```

Completed DIV30/DIV90 shared-axis comparison plots:

```text
Script:
  python_notebooks/scripts/plot_shi_seurat_timepoint_comparisons.py
Output root:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_reference_div30_div90_seurat_label_transfer_comparison/shi_reference_div30_div90_seurat_comparison_v1
Inputs:
  DIV30 tables/div30_shi_seurat_label_transfer_obs.tsv.gz
  DIV90 tables/div90_shi_seurat_label_transfer_obs.tsv.gz
Completion:
  n_cells = 112969
  n_plots = 8
  numeric_gw_range = [9, 18]
  even_gw_range = [0, 4]
  kde_bw_method = 0.08
  kde_min_density = 0.01
Plots:
  plots/div30_div90_shi_seurat_full_expected_shi_gw_ridge_by_sample_shared_axis.png
  plots/div30_div90_shi_seurat_full_expected_shi_gw_density_by_timepoint_shared_axis.png
  plots/div30_div90_shi_seurat_full_expected_shi_gw_ridge_by_sample_shared_axis_even_gw_spacing.png
  plots/div30_div90_shi_seurat_full_expected_shi_gw_density_by_timepoint_shared_axis_even_gw_spacing.png
  plots/div30_div90_shi_seurat_full_expected_shi_gw_ridge_by_sample_shared_axis_predicted_mge.png
  plots/div30_div90_shi_seurat_full_expected_shi_gw_density_by_timepoint_shared_axis_predicted_mge.png
  plots/div30_div90_shi_seurat_full_expected_shi_gw_ridge_by_sample_shared_axis_even_gw_spacing_predicted_mge.png
  plots/div30_div90_shi_seurat_full_expected_shi_gw_density_by_timepoint_shared_axis_even_gw_spacing_predicted_mge.png
Tables:
  tables/div30_div90_expected_shi_gw_summary_by_sample.tsv
  tables/div30_div90_expected_shi_gw_summary_by_sample_predicted_mge.tsv
  tables/div30_div90_shi_label_counts_by_timepoint.tsv
  tables/div30_div90_shi_label_counts_by_timepoint_predicted_mge.tsv
```

Rsync the shared-axis comparison plots from a Mac:

```bash
mkdir -p /Users/ecrespo/Downloads/shi_reference_div30_div90_seurat_comparison_v1/plots

rsync -avh --progress \
  elcrespo@greatlakes.arc-ts.umich.edu:/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_reference_div30_div90_seurat_label_transfer_comparison/shi_reference_div30_div90_seurat_comparison_v1/plots/ \
  /Users/ecrespo/Downloads/shi_reference_div30_div90_seurat_comparison_v1/plots/
```

## Cross-Study Prediction Inventory

This inventory was run before bringing the cross-study UMAP objects into the
notebook and before extending Shi-style Seurat prediction scores beyond the
currently plotted Varela DIV30/DIV90 objects. The corrected inventory includes
all studies from the cross-study Panel B config: Varela DIV30, Varela DIV90,
Walsh, Bershteyn 2025, Bershteyn 2023, Xiang, Samarasinghe, and Siebert 2026.

Script:

```text
scripts/07_cross_study_prediction_inventory.R
```

Slurm template:

```text
slurm_templates/22_cross_study_prediction_inventory.sbatch.template
```

Submitted job:

```text
Slurm job: 51414106
State: COMPLETED
ExitCode: 0:0
Elapsed: 00:03:39
MaxRSS: 37875628K
Node: gl3420
Job script:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/22_cross_study_prediction_inventory.sbatch
Logs:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/cross-study-inventory-cross-study-inv-51414106.out
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/cross-study-inventory-cross-study-inv-51414106.err
```

Output table directory:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_prediction_inventory/cross_study_prediction_inventory_v1/tables
```

Inventory tables:

```text
cross_study_prediction_inventory_complete.tsv
cross_study_prediction_object_summary.tsv
cross_study_prediction_readiness.tsv
cross_study_prediction_metadata_columns.tsv
cross_study_prediction_metadata_value_counts.tsv
cross_study_prediction_primary_sample_counts.tsv
cross_study_prediction_primary_cluster_counts.tsv
div30_div90_existing_shi_prediction_table_summary.tsv
div30_div90_existing_shi_prediction_sample_counts.tsv
div30_div90_existing_shi_prediction_cluster_counts.tsv
div30_div90_existing_shi_prediction_label_counts.tsv
```

Objects inspected:

| study_id | object | cells | RNA features | UMAP for plotting | shared Shi genes |
| --- | --- | ---: | ---: | --- | ---: |
| `varela_div30` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/varela_this_paper/varela_this_paper_seurat.rds` | 90,631 | 18,082 | `umap` | 14,354 |
| `varela_div90` | `/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/umap_props_output/clustered_day90_with_cluster_names_2.rds` | 22,338 | 18,082 | `umap` | 14,354 |
| `walsh` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/walsh_day75/walsh_day75_final_annotated.rds` | 4,519 | 20,194 | `umap_sel` | 14,945 |
| `bershteyn_2025` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/bershteyn_2025/bershteyn_2025_seurat.rds` | 124,583 | 45,068 | `umap` | 20,021 |
| `bershteyn_2023` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/bershteyn_2023/bershteyn_2023_seurat.rds` | 98,042 | 45,068 | `umap` | 20,021 |
| `xiang_2018` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/xiang_2018/xiang_2018_seurat.rds` | 58,950 | 23,287 | `umap` | 20,484 after Ensembl-to-symbol mapping |
| `samarasinghe_2021` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/samarasinghe_2021/samarasinghe_2021_seurat.rds` | missing object | missing object | `umap` requested | missing object |
| `siebert_2026` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siebert_2026/siebert_2026_seurat.rds` | 64,676 | 32,131 | `umap` | 15,305 |

Readiness conclusion:

```text
Seven of eight configured cross-study objects are ready for Seurat label
transfer, sample-level score plots, and cluster summaries. Varela DIV30,
Varela DIV90, Walsh, Bershteyn 2025, Bershteyn 2023, Xiang, and Siebert 2026
have an RNA assay, usable UMAP coordinates, sample metadata, cluster metadata,
and enough shared features with the Shi reference for transfer-anchor testing.

Samarasinghe is not ready because the canonical Seurat object is missing:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/samarasinghe_2021/samarasinghe_2021_seurat.rds

Xiang is ready only if the GSE98201 feature map is used. The raw Xiang rownames
are Ensembl IDs and have 0 direct overlaps with the Shi reference symbols. With
the configured feature map, Xiang has 20,484 mapped shared Shi genes:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/data/raw/xiang_2018_geo_files/suppl/GSE98201_genes.tsv.gz

Varela DIV30 and DIV90 already have exported Shi Seurat prediction obs tables.
Walsh, Bershteyn, Xiang, and Siebert do not yet have Shi prediction metadata
columns, so they need new Seurat TransferData runs before Scanpy-side plotting.
```

Primary metadata columns to use:

| study_id | sample column | cluster column | sample columns present | cluster columns present |
| --- | --- | --- | --- | --- |
| `varela_div30` | `orig.ident` | `seurat_clusters` | `orig.ident` | `seurat_clusters`, `RNA_snn_res.0.2` |
| `varela_div90` | `orig.ident` | `seurat_clusters` | `orig.ident` | `seurat_clusters`, `cluster_id`, `cluster_number_name`, `RNA_snn_res.0.5` |
| `walsh` | `sample_id` | `seurat_clusters` | `sample_id`, `orig.ident` | `seurat_clusters`, `RNA_snn_res.2` |
| `bershteyn_2025` | `sample` | `seurat_clusters` | `sample` | `seurat_clusters`, `predicted.GEcluster`, `predicted.GEtype`, `predicted.GEgws` |
| `bershteyn_2023` | `orig.ident` | `seurat_clusters` | `orig.ident`, `samples` | `seurat_clusters`, `celltype`, `process` |
| `xiang_2018` | `orig.ident` | `seurat_clusters` | `orig.ident` | `seurat_clusters`, `RNA_snn_res.0.5` |
| `samarasinghe_2021` | missing object | missing object | missing object | missing object |
| `siebert_2026` | `orig.ident` | `seurat_clusters` | `orig.ident`, `sample` | `seurat_clusters`, `SCT_snn_res.1`, `SCT_snn_res.0.8` |

Primary sample inventory:

| study_id | samples from primary metadata column |
| --- | --- |
| `varela_div30` | `9583-MW-6` 22,528; `9583-MW-5` 22,230; `9583-MW-3` 15,773; `9583-MW-1` 14,564; `9583-MW-4` 11,004; `9583-MW-2` 4,532 |
| `varela_div90` | `10496-MW-4` 6,314; `10496-MW-6` 4,850; `10496-MW-2` 3,533; `10496-MW-1` 3,095; `10496-MW-3` 2,714; `10496-MW-5` 1,832 |
| `walsh` | `GSM7979671` / `MEL1_dFB_d75` 2,273; `GSM7979672` / `MEL1_vFB_d75` 2,246 |
| `bershteyn_2025` | `010720S` 11,834; `200520S2` 11,123; `010519S1` 10,778; `280120S` 10,776; `010519S2` 10,661; `070120S` 10,537; `220720S1` 10,345; `100620S` 9,885; `220720S2` 8,802; `150120S` 8,294; `251219S` 6,929; `200520S1` 6,756; `111219S` 5,722; `200319S` 2,141 |
| `bershteyn_2023` | `MB279` 10,009; `MS35r41` 9,208; `r41v2ym` 8,722; `MS35mock` 8,403; `mockv2ym` 8,208; `mockv2dw` 8,206; `D0` 8,118; `MB528` 7,127; `r41v2dw` 6,656; `MB460` 6,447; `MB527` 4,933; `D14` 4,851; `MB280` 4,205; `MB461` 2,949 |
| `xiang_2018` | `Xiang2018` 58,950 |
| `samarasinghe_2021` | missing object |
| `siebert_2026` | `Old_1` 16,606; `Young_2` 16,377; `Old_2` 16,073; `Young_1` 15,620 |

Sample alias and metadata notes:

- Walsh `sample_id` values are GEO accessions. Use `MEL1_dFB_d75` as the
  readable alias for `GSM7979671` and `MEL1_vFB_d75` as the readable alias for
  `GSM7979672`. Both are day-75/DIV75 units; the object also carries
  `orig.ident` values `dFB` and `vFB`, plus `domain` values `dFB_domain` and
  `vFB_domain`.
- Siebert 2026 has four primary sample values in both `orig.ident` and
  `sample`: `Old_1`, `Young_2`, `Old_2`, and `Young_1`. Additional informative
  metadata columns include `cellLine`, `donor_batch`, and author annotation
  columns `predictions` and `predictions_class`. Use `orig.ident`/`sample` for
  primary sample-level plots; use `cellLine` or `donor_batch` only when the
  analysis needs line- or donor/batch-level stratification.

Siebert 2026 metadata value inventory was extracted by Slurm job `51479753`
using:

```text
scripts/08_extract_siebert_metadata_values.R
slurm_templates/23_extract_siebert_metadata_values.sbatch.template
```

Output tables:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siebert_2026/metadata_value_inventory
```

Key Siebert 2026 metadata values:

- `cellLine`: `CellLine1` 28,080; `CellLine3` 20,719; `CellLine2` 14,993;
  `CellLine4` 884.
- `donor_batch`: `Old_1_CellLine1` 7,554; `Old_2_CellLine1` 7,429;
  `Young_2_CellLine1` 6,792; `Young_1_CellLine1` 6,305;
  `Young_2_CellLine3` 5,707; `Young_1_CellLine3` 5,484;
  `Old_1_CellLine3` 4,880; `Old_2_CellLine3` 4,648;
  `Old_1_CellLine2` 4,172; `Old_2_CellLine2` 3,996;
  `Young_1_CellLine2` 3,413; `Young_2_CellLine2` 3,412;
  `Young_2_CellLine4` 466; `Young_1_CellLine4` 418.
- `predictions`: 16 author labels. Top labels are `progenitors` 44,680,
  `MGE_IN` 10,200, `CGE_NR2F2/PROX1` 3,775, and `MGE_SST` 1,673.
- `predictions_class`: 21 author labels. Top labels are `IPC_late` 21,578,
  `MGE_IN` 10,899, `Radial_glia` 10,386, `Neuroblast` 7,185, and
  `CGE_NR2F2/PROX1` 3,850.

Primary cluster inventory:

| study_id | cluster summary |
| --- | --- |
| `varela_div30` | `seurat_clusters`, 7 clusters: 0, 1, 2, 3, 4, 6, 7 |
| `varela_div90` | `seurat_clusters`, 13 clusters: 0 through 12 |
| `walsh` | `seurat_clusters`, 24 clusters: 0 through 23 |
| `bershteyn_2025` | `seurat_clusters`, 9 clusters: 0 through 8 |
| `bershteyn_2023` | `seurat_clusters`, 6 clusters: 0 through 5 |
| `xiang_2018` | `seurat_clusters`, 28 clusters: 0 through 27 |
| `samarasinghe_2021` | missing object |
| `siebert_2026` | `seurat_clusters`, 28 clusters: 0 through 27 |

Existing Varela DIV30/DIV90 Shi Seurat prediction tables:

| timepoint | cells | samples | clusters | label score columns | week score columns | predicted label column | predicted week column |
| --- | ---: | ---: | ---: | ---: | ---: | --- | --- |
| DIV30 | 90,631 | 6 | 7 | 11 | 6 | `shi_seurat_full_predicted_shi_label` | `shi_seurat_full_predicted_shi_week_label` |
| DIV90 | 22,338 | 6 | 13 | 11 | 6 | `shi_seurat_full_predicted_shi_label` | `shi_seurat_full_predicted_shi_week_label` |

Current Varela predicted-label inventory:

```text
DIV30:
  MGE 49,251
  progenitor 40,476
  Thalamic neurons 556
  OPC 266
  Excitatory neuron 54
  LGE 25
  Microglia 2
  Excitatory IPC 1

DIV90:
  MGE 17,525
  progenitor 3,247
  LGE 590
  CGE 490
  OPC 269
  Thalamic neurons 180
  Endothelial 20
  Microglia 13
  Excitatory neuron 4
```

What is needed to extend prediction scores across studies:

```text
1. Keep using the Shi reference Seurat object as the reference:
   /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_2019_paper_qc/shi_2019_seurat.rds

2. For each target object, run Seurat FindTransferAnchors/TransferData using
   RNA features shared with Shi. Do not use sample IDs or clusters as prediction
   inputs.

   For Xiang, first harmonize Ensembl rownames to gene symbols using:
   /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/data/raw/xiang_2018_geo_files/suppl/GSE98201_genes.tsv.gz

   For Samarasinghe, first generate or register the missing Seurat object at:
   /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/samarasinghe_2021/samarasinghe_2021_seurat.rds

3. Export per-cell predicted Shi label, max label score, full per-label score
   matrix, predicted Shi gestational-week label, max week score, and full
   per-week score matrix.

4. Bring the exported tables back into AnnData/Scanpy only for plotting and
   summaries, using the primary sample and cluster columns listed above.
```

Rsync inventory tables from a Mac:

```bash
mkdir -p /Users/ecrespo/Downloads/cross_study_prediction_inventory_v1/tables

rsync -avh --progress \
  elcrespo@greatlakes.arc-ts.umich.edu:/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_prediction_inventory/cross_study_prediction_inventory_v1/tables/ \
  /Users/ecrespo/Downloads/cross_study_prediction_inventory_v1/tables/
```

## HTML Reference For This Target

Local inspiration file from the Mac:

```text
file:///Users/ecrespo/Downloads/ms2024.html
```

Public Dataverse source:

```text
https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/ZQSPKW
```

Dataverse API metadata endpoint:

```text
https://dataverse.harvard.edu/api/datasets/:persistentId/?persistentId=doi:10.7910/DVN/ZQSPKW
```

The API lists one unrestricted HTML file:

```text
filename: ms2024.html
dataFile id: 11713377
contentType: text/html
filesize: 20,925,576 bytes
license: CC0 1.0
```

Download command for the analysis environment:

```bash
curl -L \
  "https://dataverse.harvard.edu/api/access/datafile/11713377" \
  -o /tmp/ms2024.html
```

The rendered Dataverse landing page may require JavaScript/browser
verification, but the API and direct file endpoint are readable from the Linux
workspace.

Use this HTML specifically for the Seurat prediction-score plotting pattern.
Do not copy the biological interpretation directly, because the context differs.
The useful transferable pattern is:

1. Treat Seurat `prediction.score.*` assays as score matrices.
2. Plot selected score rows on UMAP.
3. Melt full score matrices to long tables for density/ridge-style plots.
4. Use metadata `predicted.*` and `predicted.*.score` columns for winner-label
   summaries and confidence distributions.

```text
Observed relevant pattern in ms2024.html:
  prediction.score.adSupercluster
  prediction.score.adROIGroup
  prediction.score.mtg_type
  FeaturePlot(..., features = rownames(prediction.score.* assay), ...)
  reshape2::melt(t(GetAssayData(prediction.score.* assay)))
  FetchData(..., vars = c("predicted.*.score", "predicted.*", "mpt"))
```

The Seurat route should preserve the same core rule: DIV30 `seurat_clusters`
must not be used as prediction input. They are only used after prediction for
cluster summaries, stacked bars, heatmaps, and interpretation.

Candidate Seurat-side workflow for this target:

```r
anchors <- FindTransferAnchors(
  reference = shi_reference_seurat,
  query = div30_query_seurat,
  normalization.method = "LogNormalize",
  reference.assay = "RNA",
  query.assay = "RNA",
  dims = 1:50
)

label_predictions <- TransferData(
  anchorset = anchors,
  refdata = shi_reference_seurat$shi_label,
  dims = 1:50
)

div30_query_seurat <- AddMetaData(div30_query_seurat, label_predictions)
```

Expected Seurat prediction fields to export:

```text
cell_id
predicted.id or predicted.<transfer_name>
prediction.score.max or predicted.<transfer_name>.score
prediction.score.<label> columns or prediction.score.<transfer_name> assay rows
```

Suggested AnnData column names after import:

```text
shi_seurat_full_predicted_shi_label
shi_seurat_full_prediction_score
shi_seurat_full_uncertainty_score
shi_seurat_full_prediction_score_MGE
shi_seurat_full_prediction_score_LGE
shi_seurat_full_prediction_score_CGE
...
```

For gestational-week transfer, run a second `TransferData` call with
`refdata = shi_reference_seurat$week_label` and import analogous columns:

```text
shi_seurat_full_predicted_shi_week_label
shi_seurat_full_week_prediction_score
shi_seurat_full_week_uncertainty_score
```

After import, Scanpy plotting should remain unchanged in spirit: color the DIV30
Seurat UMAP by `shi_seurat_full_predicted_shi_label`,
`shi_seurat_full_prediction_score`, individual
`shi_seurat_full_prediction_score_*` label-score columns, broad/developmental
classes derived from the predicted label, and later week prediction columns.
Cluster-level summaries should use DIV30 `seurat_clusters` only after these
per-cell predictions have been joined.

## Legacy v2 Output Contract

This output contract belongs to the completed Scanpy-side
`shi_reference_div30_label_transfer_v2` run. For the next Seurat score target,
use the expected output roots and files listed in
`Next Target: Seurat Prediction Scores`.

Run outputs:

```text
RUN_DIR/tables/
RUN_DIR/plots/
RUN_DIR/h5ad/
```

For `shi_reference_div30_label_transfer_v2`, the completed run contains:

```text
plots/*.png: 14 files
tables/*: 34 files
```

Key tables:

```text
tables/shi_reference_table_s2_label_join_summary.tsv
tables/shi_reference_table_s2_label_counts.tsv
tables/shi_reference_comparison_label_sets.tsv
tables/shi_label_transfer_gene_harmonization_summary.tsv
tables/shi_label_transfer_gene_harmonization_detail.tsv
tables/div30_shi_label_transfer_obs.tsv.gz
tables/div30_shi_label_transfer_predictions_long.tsv.gz
tables/div30_shi_label_transfer_cluster_summaries.tsv
tables/div30_shi_label_transfer_label_fractions_by_cluster.tsv
tables/div30_shi_full_vs_mge_lge_cge_predictions.tsv.gz
tables/div30_shi_full_vs_mge_lge_cge_summary_by_seurat_clusters.tsv
tables/div30_shi_week_prediction_counts_by_sample.tsv
tables/div30_shi_week_prediction_summary_by_sample.tsv
tables/shi_reference_week_counts.tsv
tables/shi_reference_week_metadata.tsv
tables/shi_label_transfer_output_manifest.tsv
tables/shi_label_transfer_complete.tsv
```

Annotated AnnData:

```text
h5ad/div30_shi_label_transfer_predictions.h5ad
```

Prediction columns added to `adata.obs`:

```text
shi_full_predicted_shi_label
shi_full_prediction_score
shi_full_uncertainty_score
shi_full_prediction_entropy
shi_full_broad_region_class
shi_full_developmental_class
shi_full_predicted_shi_week_label
shi_full_week_prediction_score
shi_full_week_uncertainty_score
shi_full_week_prediction_entropy
shi_full_mean_neighbor_week_numeric
shi_full_median_neighbor_week_numeric
shi_full_std_neighbor_week_numeric

shi_mge_lge_cge_predicted_shi_label
shi_mge_lge_cge_prediction_score
shi_mge_lge_cge_uncertainty_score
shi_mge_lge_cge_prediction_entropy
shi_mge_lge_cge_broad_region_class
shi_mge_lge_cge_developmental_class
shi_mge_lge_cge_predicted_shi_week_label
shi_mge_lge_cge_week_prediction_score
shi_mge_lge_cge_week_uncertainty_score
shi_mge_lge_cge_week_prediction_entropy
shi_mge_lge_cge_mean_neighbor_week_numeric
shi_mge_lge_cge_median_neighbor_week_numeric
shi_mge_lge_cge_std_neighbor_week_numeric
```

## Legacy v2 Plot Inventory

`shi_reference_div30_label_transfer_v2` generated exactly 14 PNG plots:
7 plots for `full_relevant` and 7 plots for `mge_lge_cge_only`.

Run plot directory:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_reference_div30_label_transfer/shi_reference_div30_label_transfer_v2/plots
```

Shared plot inputs:

```text
DIV30 query AnnData:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/varela_div30.h5ad
  UMAP coordinates from .obsm["X_umap_seurat"]
  sample groups from .obs["orig.ident"]
  post hoc cluster groups from .obs["seurat_clusters"]

Shi reference AnnData:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/shi_2019_paper_qc.h5ad

Shi Table S2:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/reference/shi_2021_tables_s2_to_s9/science.abj6641_table_s2.xlsx
  cell-level Major types and gestational-week-derived labels
```

Plot helper functions are defined in:

```text
python_notebooks/src/mge_organoid_python/shi_label_transfer.py
```

Full-reference plots:

```text
1. div30_umap_shi_full_predicted_shi_label.png
   Input: query .obsm["X_umap_seurat"] + query.obs["shi_full_predicted_shi_label"]
   Function: plot_umap_categorical

2. div30_umap_shi_full_prediction_score.png
   Input: query .obsm["X_umap_seurat"] + query.obs["shi_full_prediction_score"]
   Function: plot_umap_continuous

3. div30_umap_shi_full_broad_region_class.png
   Input: query .obsm["X_umap_seurat"] + query.obs["shi_full_broad_region_class"]
   Function: plot_umap_categorical

4. div30_shi_full_shi_label_stacked_bar_by_seurat_clusters.png
   Input table: tables/div30_shi_full_label_counts.tsv
   Source columns: seurat_clusters, shi_full_predicted_shi_label, fraction
   Function: plot_stacked_bar

5. div30_shi_full_shi_label_heatmap_by_seurat_clusters.png
   Input table: tables/div30_shi_full_label_counts.tsv
   Source columns: seurat_clusters, shi_full_predicted_shi_label, fraction
   Function: plot_heatmap

6. div30_shi_full_mean_neighbor_gw_age_density_by_sample.png
   Input: query.obs["orig.ident"] + query.obs["shi_full_mean_neighbor_week_numeric"]
   Function: plot_overlaid_density_by_group

7. div30_shi_full_gw_prediction_score_density_by_sample.png
   Input: query.obs["orig.ident"] + query.obs["shi_full_week_prediction_score"]
   Function: plot_overlaid_density_by_group
```

MGE/LGE/CGE-restricted plots:

```text
8. div30_umap_shi_mge_lge_cge_predicted_shi_label.png
   Input: query .obsm["X_umap_seurat"] + query.obs["shi_mge_lge_cge_predicted_shi_label"]
   Function: plot_umap_categorical

9. div30_umap_shi_mge_lge_cge_prediction_score.png
   Input: query .obsm["X_umap_seurat"] + query.obs["shi_mge_lge_cge_prediction_score"]
   Function: plot_umap_continuous

10. div30_umap_shi_mge_lge_cge_broad_region_class.png
    Input: query .obsm["X_umap_seurat"] + query.obs["shi_mge_lge_cge_broad_region_class"]
    Function: plot_umap_categorical

11. div30_shi_mge_lge_cge_shi_label_stacked_bar_by_seurat_clusters.png
    Input table: tables/div30_shi_mge_lge_cge_label_counts.tsv
    Source columns: seurat_clusters, shi_mge_lge_cge_predicted_shi_label, fraction
    Function: plot_stacked_bar

12. div30_shi_mge_lge_cge_shi_label_heatmap_by_seurat_clusters.png
    Input table: tables/div30_shi_mge_lge_cge_label_counts.tsv
    Source columns: seurat_clusters, shi_mge_lge_cge_predicted_shi_label, fraction
    Function: plot_heatmap

13. div30_shi_mge_lge_cge_mean_neighbor_gw_age_density_by_sample.png
    Input: query.obs["orig.ident"] + query.obs["shi_mge_lge_cge_mean_neighbor_week_numeric"]
    Function: plot_overlaid_density_by_group

14. div30_shi_mge_lge_cge_gw_prediction_score_density_by_sample.png
    Input: query.obs["orig.ident"] + query.obs["shi_mge_lge_cge_week_prediction_score"]
    Function: plot_overlaid_density_by_group
```

Notes for tomorrow:

```text
tables/shi_reference_comparison_label_sets.tsv confirms full_relevant has all
10 Shi labels, and mge_lge_cge_only has MGE, LGE, and CGE.

tables/div30_shi_full_label_counts.tsv keeps zero-count label rows, so CGE,
Microglia, OPC, and Endothelial do not silently disappear when their predicted
fraction is 0.

The overlaid GW-age plots are not based on DIV30 collection age directly. They
are inferred from the Shi reference neighbors for each DIV30 cell and then
overlaid by DIV30 sample (`orig.ident`).
```

## Slurm Command

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
cp slurm_templates/20_execute_shi_reference_div30_label_transfer.sbatch.template \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/20_execute_shi_reference_div30_label_transfer.sbatch
sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/20_execute_shi_reference_div30_label_transfer.sbatch
```

## Rsync From Mac

From `/Users/ecrespo/Downloads`:

```bash
mkdir -p ./shi_reference_div30_label_transfer_v2

rsync -avh --progress \
  elcrespo@greatlakes.arc-ts.umich.edu:/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_reference_div30_label_transfer/shi_reference_div30_label_transfer_v2/plots/ \
  ./shi_reference_div30_label_transfer_v2/plots/

rsync -avh --progress \
  elcrespo@greatlakes.arc-ts.umich.edu:/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_reference_div30_label_transfer/shi_reference_div30_label_transfer_v2/tables/ \
  ./shi_reference_div30_label_transfer_v2/tables/
```

## Slurm Status

Latest corrected run:

```text
51359356:
  Run label: shi_reference_div30_label_transfer_v2
  State: COMPLETED
  ExitCode: 0:0
  Elapsed: 00:06:22
  MaxRSS: 24756096K
  Node: gl3076
  Change from v1:
    - full_relevant now includes all 10 Shi Table S2 major-type labels
    - zero-count labels are retained in label-count tables and plot legends
    - added Shi gestational-week prediction columns and sample-overlaid density plots
  Output log: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/shi-label-transfer-shi-div30-xfer-51359356.out
  Error log: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/shi-label-transfer-shi-div30-xfer-51359356.err
```

Setup/debug attempt:

```text
51358066:
  State: FAILED
  Reason: Table S2 join initially matched only 9.8% of Shi reference cells.
  Fix: updated label join to infer Table S2 numeric sample suffixes and use
       per-cell unique barcode keys before suffix-level fallback.
```

Corrected submitted run:

```text
51358138:
  State: COMPLETED
  ExitCode: 0:0
  Elapsed: 00:06:14
  MaxRSS: 24623288K
  Node: gl3118
  CPUs: 8
  Memory: 160G
  Output log: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/shi-label-transfer-shi-div30-xfer-51358138.out
  Error log: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/shi-label-transfer-shi-div30-xfer-51358138.err
  Pre-submit reference-label join check: 55,704 / 56,136 matched (99.23%)
  Completion marker: present
```

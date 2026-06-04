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

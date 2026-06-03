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
shi_reference_div30_label_transfer_v1
```

Default run directory:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_reference_div30_label_transfer/shi_reference_div30_label_transfer_v1
```

Executed notebook path:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_reference_div30_label_transfer/executed/shi_reference_div30_label_transfer.shi_reference_div30_label_transfer_v1.executed.ipynb
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

Default included labels are all observed Table S2 labels except:

```text
Microglia, OPC, Endothelial
```

This keeps:

```text
MGE
LGE
CGE
progenitor
Excitatory IPC
Excitatory neuron
Thalamic neurons
```

### B. Restricted MGE/LGE/CGE Reference

Default included labels:

```text
MGE, LGE, CGE
```

This is intentionally compared against the full-reference result to identify
cells whose MGE/LGE/CGE label is potentially being forced by removing other
developmental states.

## Method

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
6. Add predictions to DIV30 `adata.obs`.
7. Summarize predictions by DIV30 `seurat_clusters`.

Prediction score is the fraction of k nearest Shi reference neighbors voting for
the winning label. Uncertainty is `1 - prediction_score`. Entropy is normalized
neighbor-label entropy.

## Output Contract

Run outputs:

```text
RUN_DIR/tables/
RUN_DIR/plots/
RUN_DIR/h5ad/
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

shi_mge_lge_cge_predicted_shi_label
shi_mge_lge_cge_prediction_score
shi_mge_lge_cge_uncertainty_score
shi_mge_lge_cge_prediction_entropy
shi_mge_lge_cge_broad_region_class
shi_mge_lge_cge_developmental_class
```

Plots are generated for each comparison:

```text
UMAP colored by Shi predicted label
UMAP colored by prediction score
UMAP colored by broad MGE/LGE/CGE class
stacked barplot of Shi labels by seurat_clusters
heatmap of seurat_clusters x Shi predicted labels
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
mkdir -p ./shi_reference_div30_label_transfer_v1

rsync -avh --progress \
  elcrespo@greatlakes.arc-ts.umich.edu:/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_reference_div30_label_transfer/shi_reference_div30_label_transfer_v1/ \
  ./shi_reference_div30_label_transfer_v1/
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

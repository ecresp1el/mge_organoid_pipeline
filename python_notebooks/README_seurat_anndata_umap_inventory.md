# Seurat And AnnData UMAP Inventory

This is a standalone inspection workflow for the source Seurat `.rds` objects
and their cached AnnData `.h5ad` conversions. It is not Notebook 02 and does
not continue the Notebook 00/01 pipeline.

## Goal

The notebook answers two narrow questions:

1. What information is present in each original Seurat object?
2. What information was preserved in the cached AnnData files, and can we
   remake the exported Seurat UMAPs from those AnnData files?

It does not rerun Seurat-to-AnnData conversion, Notebook 00, Notebook 01,
Scanpy clustering, or new UMAP computation.

## Primary Files

Notebook:

```text
python_notebooks/notebooks/seurat_anndata_umap_inventory.ipynb
```

R Seurat inventory helper:

```text
python_notebooks/scripts/inspect_seurat_object.R
```

Slurm template:

```text
slurm_templates/17_execute_seurat_anndata_umap_inventory.sbatch.template
```

Handoff:

```text
python_notebooks/HANDOFF_seurat_anndata_umap_inventory.md
```

## Inputs

The notebook uses the canonical study definitions from
`mge_organoid_python.default_studies()`:

```text
shi_2019_paper_qc
varela_div30
varela_div90
```

Source Seurat files:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_2019_paper_qc/shi_2019_seurat.rds
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/varela_this_paper/varela_this_paper_seurat.rds
/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/umap_props_output/clustered_day90_with_cluster_names_2.rds
```

Cached AnnData files:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/shi_2019_paper_qc.h5ad
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/varela_div30.h5ad
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/varela_div90.h5ad
```

## Output Location

Default run label:

```text
seurat_anndata_umap_inventory_v1
```

Default run directory:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/seurat_anndata_umap_inventory/seurat_anndata_umap_inventory_v1
```

Main subdirectories:

```text
RUN_DIR/tables/
RUN_DIR/tables/seurat_object_inventory/<study_id>/
RUN_DIR/tables/anndata_inventory/
RUN_DIR/plots/umaps/<study_id>/
RUN_DIR/logs/
```

Executed notebooks are saved separately:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/seurat_anndata_umap_inventory/executed/
```

## Run With Slurm

Use Slurm for the full run because the R-side inventory loads full Seurat
objects, including the 19 GB DIV30 `.rds`.

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
cp slurm_templates/17_execute_seurat_anndata_umap_inventory.sbatch.template \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/17_execute_seurat_anndata_umap_inventory.sbatch
sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/17_execute_seurat_anndata_umap_inventory.sbatch
```

Optional overrides:

```bash
SEURAT_INVENTORY_RUN_LABEL=seurat_anndata_umap_inventory_test_v1 \
SEURAT_INVENTORY_STUDIES=varela_div90 \
sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/17_execute_seurat_anndata_umap_inventory.sbatch
```

```bash
SEURAT_INVENTORY_OVERWRITE_R=1 \
sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/17_execute_seurat_anndata_umap_inventory.sbatch
```


## Parallel Slurm Option

For future reruns, the expensive Seurat `.rds` inventory can run as a Slurm
array, one study per task. Then a dependent merge job can build the executed
notebook, combined tables, AnnData inventory, and UMAP plots.

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
cp slurm_templates/17a_run_seurat_inventory_array.sbatch.template \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/17a_run_seurat_inventory_array.sbatch
cp slurm_templates/17b_execute_seurat_anndata_umap_inventory_merge.sbatch.template \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/17b_execute_seurat_anndata_umap_inventory_merge.sbatch

array_job=$(sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/17a_run_seurat_inventory_array.sbatch | awk '{print $4}')
sbatch --dependency=afterok:${array_job} \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/17b_execute_seurat_anndata_umap_inventory_merge.sbatch
```

The array template calls:

```text
python_notebooks/scripts/run_seurat_inventory_for_study.py
```

The merge template runs the notebook with:

```text
SEURAT_INVENTORY_RUN_R=0
SEURAT_INVENTORY_REQUIRE_R_MARKERS=1
```

That means the merge job fails if any per-study `seurat_inventory_complete.tsv`
marker is missing.

## Key Tables

Run and input metadata:

```text
RUN_DIR/tables/seurat_anndata_umap_inventory_run_parameters.tsv
RUN_DIR/tables/seurat_anndata_umap_inventory_runtime_diagnostics.tsv
RUN_DIR/tables/seurat_anndata_umap_inventory_input_paths.tsv
```

Seurat object inventory:

```text
RUN_DIR/tables/combined_seurat_object_summary.tsv
RUN_DIR/tables/combined_seurat_assays.tsv
RUN_DIR/tables/combined_seurat_assay_layers.tsv
RUN_DIR/tables/combined_seurat_reductions.tsv
RUN_DIR/tables/combined_seurat_metadata_columns.tsv
RUN_DIR/tables/combined_seurat_metadata_top_values.tsv
RUN_DIR/tables/combined_seurat_identity_counts.tsv
RUN_DIR/tables/combined_seurat_commands.tsv
RUN_DIR/tables/seurat_inventory_file_manifest.tsv
```

AnnData inventory:

```text
RUN_DIR/tables/anndata_inventory/anndata_load_reports.tsv
RUN_DIR/tables/anndata_inventory/anndata_object_inventory.tsv
RUN_DIR/tables/anndata_inventory/anndata_obs_column_inventory.tsv
RUN_DIR/tables/anndata_inventory/anndata_var_column_inventory.tsv
RUN_DIR/tables/anndata_inventory/anndata_key_inventory.tsv
RUN_DIR/tables/anndata_inventory/anndata_categorical_obs_summary.tsv
RUN_DIR/tables/anndata_inventory/anndata_numeric_obs_summary.tsv
RUN_DIR/tables/anndata_inventory/anndata_conversion_manifest.tsv
```

UMAP outputs:

```text
RUN_DIR/tables/umap_color_column_selection.tsv
RUN_DIR/tables/umap_plot_manifest.tsv
RUN_DIR/plots/umaps/<study_id>/*.png
```

## Matrix And Object Meaning

The cached AnnData files come from the existing converter:

```text
.X = Seurat data-like layer selected by StudySpec.expression_layer
.layers["counts"] = Seurat counts layer
.obs = Seurat obj@meta.data plus cell_id
.var["feature_id"] = converted assay feature names
.obsm["X_umap_seurat"] = exported Seurat UMAP coordinates
.obsm["X_umap"] = compatibility alias for X_umap_seurat
.uns["conversion_manifest"] = converter assay/layer/reduction metadata
```

The inventory tables are the source of truth for existing cache state. In the
`seurat_anndata_umap_inventory_v1` run, DIV30 and DIV90 had
`.layers["counts"]` and `.obsm["X_umap_seurat"]`, while the current Shi cache
had `.obsm["X_umap"]` only and no counts layer. The Shi source Seurat object
does contain RNA counts/data/scale.data, so rebuild the cache if parity with
the current converter is required.

The current converter does not transfer all Seurat object internals. The
inventory records what exists in the source `.rds`, including other reductions,
graphs, command history, `misc`, `tools`, active identities, and variable
features.

# Handoff: Seurat And AnnData UMAP Inventory

This standalone notebook inspects the source Seurat `.rds` objects and cached
AnnData `.h5ad` files, then remakes UMAP plots from the exported Seurat
coordinates stored in AnnData. It is intentionally not named Notebook 02.

## Current Run Target

Default run label:

```text
seurat_anndata_umap_inventory_v1
```

Default run directory:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/seurat_anndata_umap_inventory/seurat_anndata_umap_inventory_v1
```

Executed notebook path:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/seurat_anndata_umap_inventory/executed/seurat_anndata_umap_inventory.seurat_anndata_umap_inventory_v1.executed.ipynb
```

Completed successful Slurm jobs:

```text
51347052:
  State: COMPLETED
  Purpose: full wrapper rerun after R inventories existed
  Elapsed: 00:00:47
  MaxRSS: 11123016K

51347110:
  State: COMPLETED
  Purpose: final merge-only validation of the array-compatible source notebook
  Elapsed: 00:00:46
  MaxRSS: 11311368K
  Executed notebook written: yes
  Completion marker lines: 5
  Tables reported by notebook: 61
  Plots reported by notebook: 33
  R inventory logs reported by notebook: 3
```

Earlier failed setup/debug jobs:

```text
51346967: failed while serializing an R function/closure in Seurat command params
51346990: failed on second serializer fallback edge case
51347008: completed all three R inventories, then failed merging TSVs with an existing study_id column
```

## Code Added

Notebook:

```text
python_notebooks/notebooks/seurat_anndata_umap_inventory.ipynb
```

R inventory helper:

```text
python_notebooks/scripts/inspect_seurat_object.R
```

Slurm template:

```text
slurm_templates/17_execute_seurat_anndata_umap_inventory.sbatch.template
```

README:

```text
python_notebooks/README_seurat_anndata_umap_inventory.md
```

## Inputs

Studies come from `mge_organoid_python.default_studies()`:

```text
shi_2019_paper_qc
varela_div30
varela_div90
```

Source Seurat `.rds` files:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_2019_paper_qc/shi_2019_seurat.rds
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/varela_this_paper/varela_this_paper_seurat.rds
/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/umap_props_output/clustered_day90_with_cluster_names_2.rds
```

Cached AnnData `.h5ad` files:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/shi_2019_paper_qc.h5ad
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/varela_div30.h5ad
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/varela_div90.h5ad
```

These files were confirmed present before creating the workflow:

```text
shi_2019_paper_qc.h5ad: 431M
varela_div30.h5ad: 4.4G
varela_div90.h5ad: 520M
shi_2019_seurat.rds: 724M
varela_this_paper_seurat.rds: 19G
clustered_day90_with_cluster_names_2.rds: 3.4G
```

## Output Contract

Run outputs:

```text
RUN_DIR/tables/
RUN_DIR/tables/seurat_object_inventory/<study_id>/
RUN_DIR/tables/anndata_inventory/
RUN_DIR/plots/umaps/<study_id>/
RUN_DIR/logs/
```

Executed notebooks:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/seurat_anndata_umap_inventory/executed/
```

Slurm logs:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/execute-seurat-anndata-inventory-seurat-ann-inv-<jobid>.out
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/execute-seurat-anndata-inventory-seurat-ann-inv-<jobid>.err
```

## Current Workflow

Notebook steps:

```text
1. Resolve repo root and PROJECT_ROOT.
2. Create RUN_DIR, TABLE_DIR, PLOT_DIR, LOG_DIR from run label.
3. Validate Seurat source paths and cached H5AD paths.
4. Run R Seurat inventory for each study unless existing tables are reused.
5. Combine Seurat inventory tables across studies.
6. Load cached AnnData objects with backed="r".
7. Write AnnData structure, obs, var, layer, obsm, obsp, varm, and uns summaries.
8. Write conversion-scope table explaining what the converter preserved.
9. Replot UMAPs from .obsm["X_umap_seurat"] or .obsm["X_umap"].
10. Write output manifest and completion marker.
```

The UMAP plots replot exported Seurat UMAP coordinates. They do not recompute
Scanpy PCA, neighbors, or UMAP.


## Parallel Slurm Option

The original full notebook Slurm wrapper runs study inventories sequentially in
one notebook kernel. The expensive part can now be parallelized with an optional
array/dependent-merge workflow.

Array inventory template:

```text
slurm_templates/17a_run_seurat_inventory_array.sbatch.template
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/17a_run_seurat_inventory_array.sbatch
```

Merge notebook template:

```text
slurm_templates/17b_execute_seurat_anndata_umap_inventory_merge.sbatch.template
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/17b_execute_seurat_anndata_umap_inventory_merge.sbatch
```

One-study inventory wrapper:

```text
python_notebooks/scripts/run_seurat_inventory_for_study.py
```

Submit pattern:

```bash
array_job=$(sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/17a_run_seurat_inventory_array.sbatch | awk '{print $4}')
sbatch --dependency=afterok:${array_job}   /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/17b_execute_seurat_anndata_umap_inventory_merge.sbatch
```

The merge job uses `SEURAT_INVENTORY_RUN_R=0` and
`SEURAT_INVENTORY_REQUIRE_R_MARKERS=1`, so it does not rerun `.rds` inspection
and fails if any array task did not write `seurat_inventory_complete.tsv`.

## Slurm Command

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
cp slurm_templates/17_execute_seurat_anndata_umap_inventory.sbatch.template \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/17_execute_seurat_anndata_umap_inventory.sbatch
sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/17_execute_seurat_anndata_umap_inventory.sbatch
```

## Important Interpretation

The existing converter preserves these AnnData pieces:

```text
.X = Seurat data-like layer requested by StudySpec.expression_layer
.layers["counts"] = Seurat counts layer
.obs = Seurat obj@meta.data plus cell_id
.var["feature_id"] = converted assay feature names
.obsm["X_umap_seurat"] = exported Seurat UMAP coordinates
.obsm["X_umap"] = compatibility alias for X_umap_seurat
.uns["source_seurat_path"], .uns["seurat_assay"], .uns["seurat_reduction"], .uns["conversion_manifest"]
```

The completed `seurat_anndata_umap_inventory_v1` run found this current cache
state:

```text
shi_2019_paper_qc:
  AnnData shape: 56136 x 21191
  layers["counts"]: absent
  obsm["X_umap"]: present
  obsm["X_umap_seurat"]: absent
  source Seurat object: RNA assay with counts/data/scale.data, pca and umap

varela_div30:
  AnnData shape: 90631 x 18082
  layers["counts"]: present
  obsm["X_umap"]: present
  obsm["X_umap_seurat"]: present
  source Seurat object: RNA assay, pca/integrated.cca/umap reductions

varela_div90:
  AnnData shape: 22338 x 18082
  layers["counts"]: present
  obsm["X_umap"]: present
  obsm["X_umap_seurat"]: present
  source Seurat object: RNA assay, pca/integrated.cca/umap reductions
```

This means the current Shi cache appears older or less complete than the
DIV30/DIV90 caches with respect to converter parity. The source Seurat object
does contain the missing count/data layers, so this can be fixed later by
rebuilding the Shi cache if needed.

The source Seurat object may contain more than this, including additional
assays, reductions, graphs, command history, active identities, `misc`,
`tools`, and variable feature lists. The standalone notebook records those
pieces in TSV tables so future analysis can decide whether any should be
exported into AnnData.

## Pickup Notes

After Slurm completes, check:

```bash
sacct -j <jobid> --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS,NodeList
```

Then inspect:

```text
RUN_DIR/tables/seurat_anndata_umap_inventory_output_manifest.tsv
RUN_DIR/tables/seurat_to_anndata_conversion_scope.tsv
RUN_DIR/tables/umap_plot_manifest.tsv
RUN_DIR/tables/combined_seurat_object_summary.tsv
RUN_DIR/tables/combined_seurat_metadata_columns.tsv
RUN_DIR/tables/anndata_inventory/anndata_object_inventory.tsv
```

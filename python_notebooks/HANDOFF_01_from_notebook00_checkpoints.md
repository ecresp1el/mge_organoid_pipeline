# Handoff: Notebook 01 From Notebook 00 Checkpoints

Branch: `codex/notebook01-from-nb00-checkpoints`

Start point: `main` after the verified Notebook 00 freeze merge.

Notebook 00 is now frozen as a loading, QC, filtering, and checkpoint-creation
notebook. Notebook 01 should start from the `.h5ad` checkpoints written by
Notebook 00 and should not rerun matrix loading, source comparison setup, or
manual_ec filtering.

## Current Handoff For Pickup On 2026-06-02

Use this section first. Older sections below include planning notes from before
the CCDifference run and before the Seurat-ordering correction; where they
conflict with this section, treat them as historical context.

### Seurat-Aligned CCDifference Rerun Completed

Completed run:

```text
NOTEBOOK00_RUN_LABEL=cellranger_filtered_manual_ec_div30_core_samples_freeze
RUN_LABEL=cellranger_filtered_manual_ec_div30_core_samples_freeze_ccdifference_seurat_order_v1
RUN_DIR=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook01/cellranger_filtered_manual_ec_div30_core_samples_freeze_ccdifference_seurat_order_v1
```

Slurm execution:

```text
Job ID: 51278701
State: COMPLETED
ExitCode: 0:0
Elapsed: 00:19:15
MaxRSS: 30605012K
Node: gl3076
```

Executed notebook:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook01/executed/01_notebook00_checkpoint_regression_comparison.cellranger_filtered_manual_ec_div30_core_samples_freeze_ccdifference_seurat_order_v1.executed.ipynb
```

Slurm logs:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/execute-notebook01-exec-nb01-regression-51278701.out
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/execute-notebook01-exec-nb01-regression-51278701.err
```

Completion checks:

```text
notebook01_branch_summary.tsv rows: 14
branches: not_regressed=7, regressed_ccdifference=7
regressed_ccdifference h5ad paths: 7/7 non-empty
not_regressed h5ad paths: 0/7, by design
all branches n_hvg_genes=4000
matrix-flow validation lines in Slurm .out: 231
```

Combined saved carry-forward object:

```text
RUN_DIR/h5ad/combined/regressed_ccdifference/analysis_hvg_scaled_umap.h5ad
shape: 97658 cells x 4000 HVGs
.X: scaled CCDifference-regressed residuals for 4000 HVGs
.layers["counts"]: raw counts for the same 4000 HVGs
.obsm["X_pca"]: 97658 x 50
.obsm["X_umap"]: 97658 x 2
.varm["PCs"]: 4000 x 50
```

Final carry-forward `.h5ad` outputs:

```text
RUN_DIR/h5ad/combined/regressed_ccdifference/analysis_hvg_scaled_umap.h5ad
RUN_DIR/h5ad/per_sample/9853-MW-1/regressed_ccdifference/analysis_hvg_scaled_umap.h5ad
RUN_DIR/h5ad/per_sample/9853-MW-2/regressed_ccdifference/analysis_hvg_scaled_umap.h5ad
RUN_DIR/h5ad/per_sample/9853-MW-3/regressed_ccdifference/analysis_hvg_scaled_umap.h5ad
RUN_DIR/h5ad/per_sample/9853-MW-4/regressed_ccdifference/analysis_hvg_scaled_umap.h5ad
RUN_DIR/h5ad/per_sample/9853-MW-5/regressed_ccdifference/analysis_hvg_scaled_umap.h5ad
RUN_DIR/h5ad/per_sample/9853-MW-6/regressed_ccdifference/analysis_hvg_scaled_umap.h5ad
```

### Source Code Updated For Seurat-Aligned Reruns

The source implementation has been updated so Notebook 01 follows the intended
Jeyoon/Seurat ordering more closely:

```text
Default RUN_LABEL:
cellranger_filtered_manual_ec_div30_core_samples_freeze_ccdifference_seurat_order_v1
```

```text
1. Start from Notebook 00 combined normalized/log1p checkpoint
2. Score S_score, G2M_score, phase, and CCDifference
3. For each branch, copy the full normalized/log1p gene matrix
4. For regressed_ccdifference, run sc.pp.regress_out(..., keys=["CCDifference"])
   on the full branch .X matrix
5. Use Scanpy's default regress_out n_jobs behavior; no explicit n_jobs is set
6. Run Seurat-v3 HVG selection after the branch regression step
7. Select 4,000 HVGs from .layers["counts"] with batch_key=None
8. Subset the branch to those HVGs
9. Scale, PCA, neighbors, UMAP, and Leiden from the branch-specific HVG .X
```

### Source Defaults Updated After Completed Run

As of 2026-06-02, the source defaults for the next Notebook 01 run are:

```text
n_pcs=10
n_neighbors=20
neighbors_use_rep=X_pca
leiden_resolution=0.8
```

The completed Slurm job `51278701` above was run before this default update;
treat that output as the completed Seurat-order 4,000-HVG run, but not as a
Leiden-resolution-0.8 rerun unless the notebook is executed again.

The next run uses this explicit neighbors call:

```python
sc.pp.neighbors(
    branch_adata,
    n_neighbors=20,
    n_pcs=10,
    use_rep="X_pca",
    random_state=0,
)
```

and this Leiden call:

```python
sc.tl.leiden(
    branch_adata,
    resolution=0.8,
    key_added="leiden",
    random_state=0,
)
```

The next run also keeps all existing single-branch UMAP/PCA plots and adds
side-by-side branch comparison UMAPs:

```text
plots/combined/comparison/umap_compare_<color>.png
plots/per_sample/<run_sample_id>/comparison/umap_compare_<color>.png
```

These comparison figures place `not_regressed` on the left and
`regressed_ccdifference` on the right. They are side-by-side, not overlaid,
because each branch has its own PCA, neighbor graph, and UMAP coordinate system.

The exact combined-data matrix order and commands are also documented in:

```text
python_notebooks/NOTEBOOK01_COMBINED_MATRIX_ORDER.txt
```

The run emits explicit matrix-flow validation lines to stdout and the executed
notebook. Search the Slurm `.out` log for:

```text
[Notebook01 matrix flow]
```

Those lines confirm:

```text
PCA input:       scaled branch_adata.X
PCA saved in:    branch_adata.obsm["X_pca"], branch_adata.varm["PCs"], branch_adata.uns["pca"]
Neighbors input: branch_adata.obsm["X_pca"] with use_rep="X_pca"
Neighbors saved: branch_adata.uns["neighbors"] plus Scanpy obsp graph matrices
UMAP input:      neighbors graph
UMAP saved in:   branch_adata.obsm["X_umap"]
Leiden input:    neighbors graph
Leiden saved in: branch_adata.obs["leiden"]
```

The same mapping is also recorded in `notebook01_branch_summary.tsv` columns
and in each saved branch object's
`uns["notebook01_branch"]["matrix_flow_validation"]`.

Updated Scanpy HVG command for the next run:

```python
sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=4000,
    layer="counts",
    batch_key=None,
    subset=False,
)
```

Updated CCDifference regression call for the next run:

```python
sc.pp.regress_out(adata, keys=["CCDifference"])
```

Reason for removing explicit `n_jobs`:

```text
The previous explicit n_jobs value was a performance assumption, not a
correctness or reproducibility requirement. Future Notebook 01 code uses the
Scanpy default.
```

Reason for `n_top_genes=4000`:

```text
Jeyoon's Seurat workflow uses FindVariableFeatures(selection.method="vst",
nfeatures=4000). Notebook 01 now defaults to 4,000 Seurat-v3 HVGs.
```

Reason for `batch_key=None`:

```text
batch_key=None most closely matches the pasted Seurat
FindVariableFeatures(..., selection.method="vst") call. A batch key such as
run_sample_id would make Scanpy perform batch-aware HVG selection and is better
treated as a separate integration-feature-selection experiment, not the direct
reproduction of that Seurat step.
```

Seurat-to-Scanpy mapping for the current target implementation:

```text
Seurat merge(SO1..SO6)
  -> Notebook 00 combined checkpoint already contains merged cells.

Seurat NormalizeData(SO1_6)
  -> Notebook 00 normalized/log1p checkpoint is loaded by Notebook 01.

Seurat ScaleData(vars.to.regress="CC.Difference", features=rownames(SO1_6))
  -> Notebook 01 branch copies the full normalized/log1p matrix, then runs
     sc.pp.regress_out(..., keys=["CCDifference"]) on full .X.
  -> Scanpy regression and scaling are separate calls. Notebook 01 scales the
     HVG subset before PCA; this is downstream-equivalent for selected HVGs
     because scaling is per gene, but it is not an all-gene scaled checkpoint.

Seurat saveRDS(...cc_regressed.rds)
  -> Not currently implemented as a separate full all-gene regressed .h5ad.
     Add this explicitly if an intermediate checkpoint is required.

Seurat FindVariableFeatures(selection.method="vst", nfeatures=4000)
  -> sc.pp.highly_variable_genes(..., flavor="seurat_v3",
     n_top_genes=4000, layer="counts", batch_key=None, subset=False)
     after the branch regression step.

Seurat RunPCA(features=VariableFeatures(SO1_6))
  -> sc.pp.pca(...) on the scaled 4,000-HVG branch object.

Seurat IntegrateLayers(method=CCAIntegration, orig.reduction="pca",
new.reduction="integrated.cca")
  -> Not implemented in current Scanpy Notebook 01. Current Notebook 01 runs
     neighbors, UMAP, and Leiden from PCA directly.

Seurat saveRDS(...cc_regressed_integrated.rds)
  -> Not implemented because CCA integration is not yet implemented.
```

Historical-output warning:

```text
The completed ccdifference_v1 output below used the older order:
HVG selection first, 2,000 HVGs, then CCDifference regression on the HVG subset.
Do not treat that output as the Seurat-order reproduction. Use the completed
ccdifference_seurat_order_v1 run above for the 4,000-HVG,
regression-before-HVG outputs.
```

### Historical Completed Run Before Seurat-Ordering Correction

Historical completed run:

```text
NOTEBOOK00_RUN_LABEL=cellranger_filtered_manual_ec_div30_core_samples_freeze
RUN_LABEL=cellranger_filtered_manual_ec_div30_core_samples_freeze_ccdifference_v1
RUN_DIR=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook01/cellranger_filtered_manual_ec_div30_core_samples_freeze_ccdifference_v1
```

Source notebook updated:

```text
/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline/python_notebooks/notebooks/01_notebook00_checkpoint_regression_comparison.ipynb
```

Executed notebook saved:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook01/executed/01_notebook00_checkpoint_regression_comparison.cellranger_filtered_manual_ec_div30_core_samples_freeze_ccdifference_v1.executed.ipynb
```

What was done:

1. Started from the frozen Notebook 00 normalized/log1p checkpoints for
   `cellranger_filtered_manual_ec_div30_core_samples_freeze`.
2. Used the already validated Regev/Tirosh S and G2M cell-cycle gene lists.
3. Scored `S_score`, `G2M_score`, and `phase` on the normalized/log1p object.
4. Created `CCDifference = S_score - G2M_score`.
5. Ran Notebook 01 branch analysis for `combined` plus all six per-sample
   checkpoints.
6. Compared two branches in every scope:
   `not_regressed` and `regressed_ccdifference`.
7. For `regressed_ccdifference`, ran
   `sc.pp.regress_out(..., keys=["CCDifference"])`.
8. Saved final `.h5ad` outputs only for `regressed_ccdifference`.

Exact operation order for the completed Notebook 01 run:

```text
1. Resolve DATA_ROOT and RUN_LABEL
2. Load Notebook 00 normalized/log1p checkpoint
3. Validate .X, .layers["counts"], obs names, and var names
4. Infer run_sample_id values
5. Score cell cycle on the in-memory normalized/log1p object
6. Add S_score, G2M_score, phase, and CCDifference to .obs
7. Save cell-cycle gene and score summary tables
8. Select 2,000 HVGs from .layers["counts"] using Seurat-v3 HVG selection
9. Run not_regressed from normalized/log1p HVG .X
10. Run regressed_ccdifference by regressing .obs["CCDifference"] from HVG .X
11. Scale, PCA, neighbors, UMAP, and Leiden inside each branch
12. Save plots for both branches
13. Save .h5ad only for regressed_ccdifference
14. Write Notebook 01 tables and executed notebook copy
```

Important state detail:

```text
Regression happens only on branch-specific copies after HVG subsetting.
Notebook 00 checkpoints are read-only inputs.
Raw counts remain available in .layers["counts"].
```

Key plots to reopen first:

```text
RUN_DIR/plots/combined/not_regressed/umap_CCDifference.png
RUN_DIR/plots/combined/regressed_ccdifference/umap_CCDifference.png
RUN_DIR/plots/combined/regressed_ccdifference/umap_phase.png
RUN_DIR/plots/combined/regressed_ccdifference/umap_leiden.png
```

Earlier focused CCDifference PCA diagnostic plots remain here:

```text
RUN_DIR/plots/cell_cycle_pca/combined/before_ccdifference_regression/
RUN_DIR/plots/cell_cycle_pca/combined/after_ccdifference_regression/
RUN_DIR/plots/cell_cycle_pca/per_sample/<run_sample_id>/before_ccdifference_regression/
RUN_DIR/plots/cell_cycle_pca/per_sample/<run_sample_id>/after_ccdifference_regression/
```

Main Notebook 01 tables saved:

```text
RUN_DIR/tables/notebook01_run_parameters.tsv
RUN_DIR/tables/notebook01_input_validation.tsv
RUN_DIR/tables/notebook01_analysis_plan.tsv
RUN_DIR/tables/notebook01_hvg_genes.tsv
RUN_DIR/tables/notebook01_hvg_parameters.tsv
RUN_DIR/tables/notebook01_cell_cycle_gene_source.tsv
RUN_DIR/tables/notebook01_cell_cycle_gene_summary.tsv
RUN_DIR/tables/notebook01_cell_cycle_score_summary.tsv
RUN_DIR/tables/notebook01_branch_summary.tsv
RUN_DIR/tables/notebook01_plot_manifest.tsv
RUN_DIR/tables/notebook01_umap_coordinates.tsv
```

Earlier CCDifference diagnostic tables saved:

```text
RUN_DIR/tables/cell_cycle_gene_source.tsv
RUN_DIR/tables/cell_cycle_gene_summary.tsv
RUN_DIR/tables/cell_cycle_genes_present.tsv
RUN_DIR/tables/cell_cycle_genes_missing.tsv
RUN_DIR/tables/cell_cycle_score_summary.tsv
RUN_DIR/tables/cell_cycle_input_validation.tsv
RUN_DIR/tables/cell_cycle_pca_diagnostic_summary.tsv
RUN_DIR/tables/cell_cycle_pca_plot_manifest.tsv
```

Final `.h5ad` outputs saved for carry-forward:

```text
RUN_DIR/h5ad/combined/regressed_ccdifference/analysis_hvg_scaled_umap.h5ad
RUN_DIR/h5ad/per_sample/9853-MW-1/regressed_ccdifference/analysis_hvg_scaled_umap.h5ad
RUN_DIR/h5ad/per_sample/9853-MW-2/regressed_ccdifference/analysis_hvg_scaled_umap.h5ad
RUN_DIR/h5ad/per_sample/9853-MW-3/regressed_ccdifference/analysis_hvg_scaled_umap.h5ad
RUN_DIR/h5ad/per_sample/9853-MW-4/regressed_ccdifference/analysis_hvg_scaled_umap.h5ad
RUN_DIR/h5ad/per_sample/9853-MW-5/regressed_ccdifference/analysis_hvg_scaled_umap.h5ad
RUN_DIR/h5ad/per_sample/9853-MW-6/regressed_ccdifference/analysis_hvg_scaled_umap.h5ad
```

Completion check already verified:

```text
notebook01_branch_summary.tsv has 14 branch rows:
  combined + six samples
  each has not_regressed and regressed_ccdifference

The 7 regressed_ccdifference rows all have non-empty h5ad_path values.
```

Operational note:

```text
The combined h5ad read was killed when attempted outside the Slurm allocation.
The successful run used the active 128 GB Slurm allocation on gl3168.
Do real Notebook 01 reruns through Slurm, not from the login/VS Code side.
```

## Operating Rules For Notebook 01 Execution

Notebook 01 real matrix-scale runs must execute on Great Lakes through Slurm,
not directly on the login node. This follows the same operational rule as
Notebook 00.

Use this conda/Jupyter environment:

```text
/home/elcrespo/miniconda3/envs/mge-organoid-python/bin
```

Use this repo and data root unless the user explicitly says otherwise:

```text
REPO_ROOT=/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
PROJECT_ROOT=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
```

Allowed:

- Load Notebook 00 `.h5ad` checkpoints from
  `PROJECT_ROOT/results/notebook00/<NOTEBOOK00_RUN_LABEL>/h5ad/`.
- Write Notebook 01 outputs under
  `PROJECT_ROOT/results/notebook01/<RUN_LABEL>/`.
- Submit real Notebook 01 execution through a Slurm template.
- Monitor Slurm jobs with `squeue` and `sacct`.
- Inspect Slurm logs under
  `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/`.
- Inspect executed notebooks under
  `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook01/executed/`.

Not allowed unless the user explicitly requests it:

- Do not rerun Notebook 00 loading/filtering from Notebook 01.
- Do not rerun CellBender denoising.
- Do not write large `.h5ad` outputs inside the GitHub checkout.
- Do not overwrite the repo notebook with an executed notebook unless the user
  asks for inline outputs in the repo notebook.
- Do not treat Reddit or issue-thread advice as authoritative without checking
  against the official Scanpy behavior and the actual object states.

## Verified Notebook 00 Inputs

Primary Cell Ranger filtered checkpoint run:

```text
RUN_LABEL=cellranger_filtered_manual_ec_div30_core_samples_freeze
RUN_DIR=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/cellranger_filtered_manual_ec_div30_core_samples_freeze
```

Primary combined input files:

```text
RUN_DIR/h5ad/manual_ec_filtered_counts.h5ad
RUN_DIR/h5ad/manual_ec_filtered_normalized_log1p.h5ad
```

Primary per-sample input files:

```text
RUN_DIR/h5ad/per_sample/<run_sample_id>.manual_ec_filtered_counts.h5ad
RUN_DIR/h5ad/per_sample/<run_sample_id>.manual_ec_filtered_normalized_log1p.h5ad
```

Optional CellBender-denoised comparison checkpoint run:

```text
RUN_LABEL=cellbender_denoised_manual_ec_div30_core_samples_freeze
RUN_DIR=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/cellbender_denoised_manual_ec_div30_core_samples_freeze
```

Notebook 00 validation results:

```text
cellranger_filtered: all 7 checkpoint validation rows passed
cellbender_denoised: all 7 checkpoint validation rows passed
combined h5ad files per run: 2
per-sample h5ad files per run: 12
```

## Input Object Contract

Counts checkpoint:

```text
manual_ec_filtered_counts.h5ad
.X = QC/manual_ec-filtered raw counts
```

Normalized/log1p checkpoint:

```text
manual_ec_filtered_normalized_log1p.h5ad
.X = normalized/log1p expression
.layers["counts"] = QC/manual_ec-filtered raw counts
```

Notebook 01 should preserve this count layer contract in downstream outputs.

## Order Of Operations And Variable Flow

Notebook 01 starts from the normalized/log1p Notebook 00 checkpoint, not the
counts-only checkpoint:

```text
manual_ec_filtered_normalized_log1p.h5ad
```

At load time:

```text
adata.X = normalized/log1p expression
adata.layers["counts"] = raw QC/manual_ec-filtered counts
adata.obs = per-cell metadata and QC covariates
adata.var = per-gene metadata
```

The current run uses the following per-cell `.obs` covariates for QC
regression:

```text
total_counts
pct_counts_mt
```

These covariates are read from `adata.obs`. They are not expression matrices
and they are not stored in `.X`.

The first-pass operation order is:

```text
1. Load Notebook 00 normalized/log1p checkpoint
2. Validate that .layers["counts"] exists and matches .X dimensions
3. Infer run_sample_id values from adata.obs["run_sample_id"]
4. Run Seurat-v3 HVG selection from adata.layers["counts"]
5. For each scope and branch, copy adata[:, hvg_mask]
6. If branch has regress_keys, run sc.pp.regress_out on branch_adata.X
7. Scale branch_adata.X
8. Run PCA from scaled branch_adata.X
9. Build neighbors from PCA
10. Run UMAP from neighbors
11. Run Leiden clustering from neighbors
12. Save tables and plots for that scope/branch
```

Important matrix-state detail:

```text
HVG selection reads .layers["counts"]
Regression modifies branch_adata.X only
Scaling modifies branch_adata.X again
PCA/neighbors/UMAP/Leiden consume the branch-specific state
The original loaded adata object is not used as the final branch object
```

Current branch behavior:

```text
not_regressed
  branch_adata.X starts as normalized/log1p HVG expression
  regression is skipped
  scaling/PCA/neighbors/UMAP/Leiden run from non-regressed .X

regressed_qc
  branch_adata.X starts as normalized/log1p HVG expression
  sc.pp.regress_out(branch_adata, ["total_counts", "pct_counts_mt"]) runs
  branch_adata.X becomes regression residuals
  scaling/PCA/neighbors/UMAP/Leiden run from regressed .X
```

In other words, regression does not delete counts and does not change the
Notebook 00 checkpoint. It changes `.X` inside the temporary branch object used
for PCA/UMAP/clustering. Raw counts remain available in `.layers["counts"]` for
auditability and future methods that need counts.

The code records this sequence in:

```text
python_notebooks/src/mge_organoid_python/notebook01_workflow.py
```

Specifically:

```text
branch_operation_order()
run_regression_embedding_branch()
```

Future runs write the order and `.X` state transitions into:

```text
PROJECT_ROOT/results/notebook01/<RUN_LABEL>/tables/notebook01_branch_summary.tsv
branch_adata.uns["notebook01_branch"]
```

## Notebook 01 Scope

Notebook 01 begins the post-filtering analysis choices:

```text
Seurat-v3 HVG selection from counts
cell-cycle scoring
CCDifference creation
regression comparison
scaling
PCA
neighbors
UMAP
clustering
```

Notebook 01 should also write its own run-scoped outputs under a new Notebook 01
results root, for example:

```text
results/notebook01/<RUN_LABEL>/
results/notebook01/<RUN_LABEL>/tables/
results/notebook01/<RUN_LABEL>/plots/
results/notebook01/<RUN_LABEL>/h5ad/
```

Current driver notebook:

```text
python_notebooks/notebooks/01_notebook00_checkpoint_regression_comparison.ipynb
```

Current Slurm runner:

```text
slurm_templates/14_execute_notebook01_regression.sbatch.template
```

Executed notebooks are written outside the Git checkout:

```text
PROJECT_ROOT/results/notebook01/executed/01_notebook00_checkpoint_regression_comparison.<RUN_LABEL>.executed.ipynb
```

First-pass output tables:

```text
PROJECT_ROOT/results/notebook01/<RUN_LABEL>/tables/notebook01_run_parameters.tsv
PROJECT_ROOT/results/notebook01/<RUN_LABEL>/tables/notebook01_input_validation.tsv
PROJECT_ROOT/results/notebook01/<RUN_LABEL>/tables/notebook01_analysis_plan.tsv
PROJECT_ROOT/results/notebook01/<RUN_LABEL>/tables/notebook01_hvg_genes.tsv
PROJECT_ROOT/results/notebook01/<RUN_LABEL>/tables/notebook01_hvg_parameters.tsv
PROJECT_ROOT/results/notebook01/<RUN_LABEL>/tables/notebook01_branch_summary.tsv
PROJECT_ROOT/results/notebook01/<RUN_LABEL>/tables/notebook01_plot_manifest.tsv
PROJECT_ROOT/results/notebook01/<RUN_LABEL>/tables/notebook01_umap_coordinates.tsv
```

First-pass plots:

```text
PROJECT_ROOT/results/notebook01/<RUN_LABEL>/plots/combined/<branch>/
PROJECT_ROOT/results/notebook01/<RUN_LABEL>/plots/per_sample/<run_sample_id>/<branch>/
```

Branch `.h5ad` outputs are optional and default to off for the first visual
comparison run:

```text
NOTEBOOK01_WRITE_BRANCH_H5AD=0
```

## References For Cell-Cycle Regression Decisions

Use these as guidance, not as hard-coded policy:

- Official Scanpy cell-cycle scoring/regression how-to:
  https://scanpy.readthedocs.io/en/latest/how-to/cell-cycle.html
- theislab/single-cell-tutorial issue on `regress_out` ordering:
  https://github.com/theislab/single-cell-tutorial/issues/35
- Reddit discussion showing common `.X`/layer confusion around Scanpy
  `regress_out`:
  https://www.reddit.com/r/bioinformatics/comments/1ke0uwj/scanpy_regress_out_question/

Practical guidance to encode in Notebook 01:

- Cell-cycle scoring must use explicit S and G2M gene lists stored in tables and
  in `.uns`.
- First regression comparison should not regress cell cycle yet. Use:
  `sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])`.
- Track both first-pass branches:
  `not_regressed` and `regressed_qc`.
- Track the matrix state used by every operation because Scanpy functions
  usually mutate `.X`.
- Treat `CCDifference = S_score - G2M_score` as the named cell-cycle covariate
  that can be added to the regression list later, after cell-cycle scoring is
  implemented.
- Compare downstream PCA/UMAP/clustering from the regressed and non-regressed
  branches rather than silently replacing one with the other.

## Next Goal: Cell-Cycle CCDifference Regression

Use the official Scanpy cell-cycle scoring/regression how-to as the procedural
reference:

```text
https://scanpy.readthedocs.io/en/latest/how-to/cell-cycle.html
```

Important deviation from the Scanpy example:

```text
Scanpy demo regression covariates: S_score, G2M_score
This project's planned regression covariate: CCDifference
```

The scoring setup should follow the Scanpy pattern:

```python
s_genes = [gene for gene in source_s_genes if gene in adata.var_names]
g2m_genes = [gene for gene in source_g2m_genes if gene in adata.var_names]
cell_cycle_genes = [*s_genes, *g2m_genes]

sc.tl.score_genes_cell_cycle(
    adata,
    s_genes=s_genes,
    g2m_genes=g2m_genes,
)

adata.obs["CCDifference"] = adata.obs["S_score"] - adata.obs["G2M_score"]
```

The cell-cycle regression branch should ultimately use:

```python
sc.pp.regress_out(adata, ["CCDifference"])
```

QC plus cell-cycle regression can be tested later as a separate branch:

```python
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt", "CCDifference"])
```

### Standalone Gene-List Validation Step

Before integrating cell-cycle scoring into the main Notebook 01 workflow, use
the importable standalone Python module and one-time Slurm validation runner.

Implemented importable module:

```text
python_notebooks/src/mge_organoid_python/cell_cycle.py
```

That module should own:

- source S-phase and G2M gene lists
- filtering those lists to `adata.var_names`
- creation of `cell_cycle_genes = [*s_genes, *g2m_genes]`
- table generation for present/missing genes
- a helper that runs `sc.tl.score_genes_cell_cycle`
- creation of `adata.obs["CCDifference"]`

Implemented one-time validation script:

```text
python_notebooks/scripts/prepare_notebook01_cell_cycle_genes.py
```

Implemented Slurm runner:

```text
slurm_templates/15_prepare_notebook01_cell_cycle_genes.sbatch.template
```

The validation job should run once against the actual Notebook 00 checkpoint
genes and write tables under a run-specific Notebook 01 output directory:

```text
PROJECT_ROOT/results/notebook01/<RUN_LABEL>/tables/cell_cycle_gene_source.tsv
PROJECT_ROOT/results/notebook01/<RUN_LABEL>/tables/cell_cycle_genes_present.tsv
PROJECT_ROOT/results/notebook01/<RUN_LABEL>/tables/cell_cycle_genes_missing.tsv
PROJECT_ROOT/results/notebook01/<RUN_LABEL>/tables/cell_cycle_gene_summary.tsv
```

Minimum validation fields:

```text
phase
gene
present_in_adata
source
n_source_s_genes
n_source_g2m_genes
n_present_s_genes
n_present_g2m_genes
n_missing_s_genes
n_missing_g2m_genes
```

### PCA Before/After CCDifference Regression

Add a focused diagnostic that proves the cell-cycle regression step changed the
cell-cycle signal in PCA space.

Run this diagnostic separately for:

```text
combined
per_sample
```

For each scope/sample:

```text
1. Load the normalized/log1p checkpoint
2. Score cell cycle using Scanpy's score_genes_cell_cycle
3. Create CCDifference
4. Create a before-regression PCA object using cell_cycle_genes
5. Plot the top three PC pairings before regression
6. Regress CCDifference from .X
7. Create an after-regression PCA object using cell_cycle_genes
8. Plot the top three PC pairings after regression
9. Save tables summarizing phase/score distributions before and after
```

Required PCA pairings:

```text
PC1 vs PC2
PC1 vs PC3
PC2 vs PC3
```

Required colors for before/after PCA diagnostics:

```text
phase
S_score
G2M_score
CCDifference
run_sample_id    # combined only
```

Implemented plot directories:

```text
PROJECT_ROOT/results/notebook01/<RUN_LABEL>/plots/cell_cycle_pca/combined/before_ccdifference_regression/
PROJECT_ROOT/results/notebook01/<RUN_LABEL>/plots/cell_cycle_pca/combined/after_ccdifference_regression/
PROJECT_ROOT/results/notebook01/<RUN_LABEL>/plots/cell_cycle_pca/per_sample/<run_sample_id>/before_ccdifference_regression/
PROJECT_ROOT/results/notebook01/<RUN_LABEL>/plots/cell_cycle_pca/per_sample/<run_sample_id>/after_ccdifference_regression/
```

Implemented diagnostic tables:

```text
PROJECT_ROOT/results/notebook01/<RUN_LABEL>/tables/cell_cycle_score_summary.tsv
PROJECT_ROOT/results/notebook01/<RUN_LABEL>/tables/cell_cycle_pca_diagnostic_summary.tsv
PROJECT_ROOT/results/notebook01/<RUN_LABEL>/tables/cell_cycle_pca_plot_manifest.tsv
```

The before/after PCA diagnostic is a visualization and validation step. It
should not replace the full Notebook 01 branch comparison until the gene-list
validation and scoring behavior are confirmed.

### CCDifference Validation Completed

Cell-cycle gene-list validation completed through Slurm:

```text
Job ID: 51248139
State: COMPLETED
ExitCode: 0:0
Elapsed: 00:00:23
MaxRSS: 3495288K
```

Validation run label:

```text
RUN_LABEL=cellranger_filtered_manual_ec_div30_core_samples_freeze_ccdifference_v1
```

Validated source gene counts:

```text
S phase:  42 present / 43 source genes
G2M:      52 present / 54 source genes
Total:    94 present / 97 source genes
```

Missing genes:

```text
S:   MLF1IP
G2M: FAM64A, HN1
```

Gene validation tables:

```text
PROJECT_ROOT/results/notebook01/cellranger_filtered_manual_ec_div30_core_samples_freeze_ccdifference_v1/tables/cell_cycle_gene_source.tsv
PROJECT_ROOT/results/notebook01/cellranger_filtered_manual_ec_div30_core_samples_freeze_ccdifference_v1/tables/cell_cycle_genes_present.tsv
PROJECT_ROOT/results/notebook01/cellranger_filtered_manual_ec_div30_core_samples_freeze_ccdifference_v1/tables/cell_cycle_genes_missing.tsv
PROJECT_ROOT/results/notebook01/cellranger_filtered_manual_ec_div30_core_samples_freeze_ccdifference_v1/tables/cell_cycle_gene_summary.tsv
```

CCDifference before/after PCA diagnostic completed through Slurm:

```text
Job ID: 51248836
State: COMPLETED
ExitCode: 0:0
Elapsed: 00:02:27
MaxRSS: 14175348K
```

Diagnostic implementation files:

```text
python_notebooks/src/mge_organoid_python/cell_cycle.py
python_notebooks/scripts/prepare_notebook01_cell_cycle_genes.py
python_notebooks/scripts/run_notebook01_ccdifference_pca_diagnostic.py
slurm_templates/15_prepare_notebook01_cell_cycle_genes.sbatch.template
slurm_templates/16_execute_notebook01_ccdifference_pca.sbatch.template
```

Diagnostic outputs:

```text
PROJECT_ROOT/results/notebook01/cellranger_filtered_manual_ec_div30_core_samples_freeze_ccdifference_v1/tables/cell_cycle_input_validation.tsv
PROJECT_ROOT/results/notebook01/cellranger_filtered_manual_ec_div30_core_samples_freeze_ccdifference_v1/tables/cell_cycle_score_summary.tsv
PROJECT_ROOT/results/notebook01/cellranger_filtered_manual_ec_div30_core_samples_freeze_ccdifference_v1/tables/cell_cycle_pca_diagnostic_summary.tsv
PROJECT_ROOT/results/notebook01/cellranger_filtered_manual_ec_div30_core_samples_freeze_ccdifference_v1/tables/cell_cycle_pca_plot_manifest.tsv
PROJECT_ROOT/results/notebook01/cellranger_filtered_manual_ec_div30_core_samples_freeze_ccdifference_v1/plots/cell_cycle_pca/
```

Diagnostic verification:

```text
input validation rows: 7, all passed
score summary rows: 21
PCA summary rows: 14
plot manifest rows: 174
missing plot paths: 0
empty plot paths: 0
```

The diagnostic runs PCA on the present cell-cycle gene subset before and after
regressing `CCDifference`. Regression is applied to the diagnostic
cell-cycle-gene subset only, which is enough to validate the effect on
cell-cycle PCA space. Full Notebook 01 branch correction can later regress
`CCDifference` on the selected downstream feature matrix.

Combined-object PCA variance changed as follows:

```text
before CCDifference regression: PC1=0.336687, PC2=0.101825, PC3=0.025132
after  CCDifference regression: PC1=0.325503, PC2=0.029885, PC3=0.017473
```

Per-sample diagnostics showed the same pattern: PC2 variance dropped after
`CCDifference` regression in all six samples.

## Required Analysis Scopes

Notebook 01 must support two scopes separately:

```text
combined
per_sample
```

The combined scope uses:

```text
results/notebook00/<NOTEBOOK00_RUN_LABEL>/h5ad/manual_ec_filtered_normalized_log1p.h5ad
```

The per-sample scope uses:

```text
results/notebook00/<NOTEBOOK00_RUN_LABEL>/h5ad/per_sample/<run_sample_id>.manual_ec_filtered_normalized_log1p.h5ad
```

The same branch labels should be available in both scopes:

```text
combined.not_regressed
combined.regressed_qc
per_sample.<run_sample_id>.not_regressed
per_sample.<run_sample_id>.regressed_qc
```

Each branch should write its own tables, plots, and optional `.h5ad` output so
we can compare what changes when the selected covariates are regressed.

Initial regression branch:

```python
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
```

Later cell-cycle-aware branch, after scoring creates `CCDifference`:

```python
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt", "CCDifference"])
```

First-pass visible branch plan for the current six-sample DIV30 checkpoint:

```text
scope       run_sample_id  branch          regress_keys
combined                  not_regressed
combined                  regressed_qc    total_counts,pct_counts_mt
per_sample  9853-MW-1      not_regressed
per_sample  9853-MW-1      regressed_qc    total_counts,pct_counts_mt
per_sample  9853-MW-2      not_regressed
per_sample  9853-MW-2      regressed_qc    total_counts,pct_counts_mt
per_sample  9853-MW-3      not_regressed
per_sample  9853-MW-3      regressed_qc    total_counts,pct_counts_mt
per_sample  9853-MW-4      not_regressed
per_sample  9853-MW-4      regressed_qc    total_counts,pct_counts_mt
per_sample  9853-MW-5      not_regressed
per_sample  9853-MW-5      regressed_qc    total_counts,pct_counts_mt
per_sample  9853-MW-6      not_regressed
per_sample  9853-MW-6      regressed_qc    total_counts,pct_counts_mt
```

## Matching Python Modules

Notebook 01 should follow the same pattern as Notebook 00:

```text
python_notebooks/src/mge_organoid_python/notebook01_workflow.py
```

That module should own:

- path classes for Notebook 00 inputs and Notebook 01 outputs
- input validation tables
- branch/scope settings
- regression variant records
- analysis summary tables
- helper functions used by the Notebook 01 driver

The notebook should be a driver. It should not hide major state transitions in
ad hoc notebook-only cells.

## Current First-Pass Implementation

Confirmed current input:

```text
NOTEBOOK01_NOTEBOOK00_RUN_LABEL=cellranger_filtered_manual_ec_div30_core_samples_freeze
```

Confirmed current scopes:

```text
NOTEBOOK01_SCOPES=combined:per_sample
```

Confirmed current branches:

```text
not_regressed
regressed_qc
```

Confirmed current regression covariates:

```text
NOTEBOOK01_REGRESS_KEYS=total_counts:pct_counts_mt
```

Confirmed current HVG settings:

```text
flavor=seurat_v3
n_top_genes=2000
layer=counts
batch_key=None
```

Historical first-pass embedding/clustering settings:

```text
n_pcs=50
n_neighbors=15
leiden_resolution=0.5
random_state=0
```

## Open Decisions After First Run

1. Which Notebook 00 checkpoint should be the primary Notebook 01 input:
   current first pass uses `cellranger_filtered`; later we can repeat against
   `cellbender_denoised`.

2. Whether Notebook 01 should operate on the combined checkpoint only first, or
   also process per-sample checkpoints in parallel. Current first pass runs
   both.

3. Exact Seurat-v3 HVG parameters:
   expected number of HVGs, batch key if any, and whether HVG should be selected
   globally or per sample.

4. Exact cell-cycle gene lists:
   use explicit stored S and G2M gene lists rather than relying on Scanpy to
   expose Seurat lists implicitly.

5. Regression variables:
   first pass is `total_counts` and `pct_counts_mt`. Later add `CCDifference`
   after cell-cycle scoring is implemented.

6. Clustering parameters:
   PCA dimensions, neighbor settings, clustering algorithm, and resolution grid.

## Implementation Guardrails

- Do not modify Notebook 00 as part of Notebook 01 work unless the user asks.
- Do not rerun CellBender.
- Do not write large `.h5ad` outputs inside the Git checkout.
- Use Slurm for real matrix-scale Notebook 01 execution.
- Keep raw counts available in `.layers["counts"]` for downstream auditability.
- Write parameter tables for HVG, cell-cycle scoring, regression, PCA,
  neighbors, UMAP, and clustering.
- Run combined and per-sample analyses separately; do not let per-sample output
  overwrite combined output or vice versa.
- Run and record both regression branches unless the user explicitly narrows the
  comparison.

## Current Execution Verification Step

Run the Notebook 01 driver through Slurm and confirm it:

1. Resolves `DATA_ROOT`.
2. Accepts a Notebook 00 `RUN_LABEL` as input.
3. Loads `results/notebook00/<RUN_LABEL>/h5ad/manual_ec_filtered_normalized_log1p.h5ad`.
4. Confirms `.layers["counts"]` exists and matches expected dimensions.
5. Writes a Notebook 01 input validation table before running downstream
   analysis.
6. Builds a planned analysis table with `combined` and `per_sample` scopes and
   both `not_regressed` and `regressed_qc` variants.
7. Saves UMAP/PCA plots for each branch under the Notebook 01 run directory.
8. Writes an executed notebook with visible plots under
   `PROJECT_ROOT/results/notebook01/executed/`.

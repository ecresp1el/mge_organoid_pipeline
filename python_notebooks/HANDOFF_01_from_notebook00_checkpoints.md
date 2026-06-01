# Handoff: Notebook 01 From Notebook 00 Checkpoints

Branch: `codex/notebook01-from-nb00-checkpoints`

Start point: `main` after the verified Notebook 00 freeze merge.

Notebook 00 is now frozen as a loading, QC, filtering, and checkpoint-creation
notebook. Notebook 01 should start from the `.h5ad` checkpoints written by
Notebook 00 and should not rerun matrix loading, source comparison setup, or
manual_ec filtering.

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

Before integrating cell-cycle scoring into the main Notebook 01 workflow, create
an importable standalone Python module and a one-time Slurm validation runner.

Proposed importable module:

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

Proposed one-time validation script:

```text
python_notebooks/scripts/prepare_notebook01_cell_cycle_genes.py
```

Proposed Slurm runner:

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

Proposed plot directories:

```text
PROJECT_ROOT/results/notebook01/<RUN_LABEL>/plots/cell_cycle_pca/combined/before_ccdifference_regression/
PROJECT_ROOT/results/notebook01/<RUN_LABEL>/plots/cell_cycle_pca/combined/after_ccdifference_regression/
PROJECT_ROOT/results/notebook01/<RUN_LABEL>/plots/cell_cycle_pca/per_sample/<run_sample_id>/before_ccdifference_regression/
PROJECT_ROOT/results/notebook01/<RUN_LABEL>/plots/cell_cycle_pca/per_sample/<run_sample_id>/after_ccdifference_regression/
```

Proposed diagnostic tables:

```text
PROJECT_ROOT/results/notebook01/<RUN_LABEL>/tables/cell_cycle_score_summary.tsv
PROJECT_ROOT/results/notebook01/<RUN_LABEL>/tables/cell_cycle_pca_diagnostic_summary.tsv
PROJECT_ROOT/results/notebook01/<RUN_LABEL>/tables/cell_cycle_pca_plot_manifest.tsv
```

The before/after PCA diagnostic is a visualization and validation step. It
should not replace the full Notebook 01 branch comparison until the gene-list
validation and scoring behavior are confirmed.

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

Confirmed current embedding/clustering settings:

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

# Handoff: Modularize Notebook 00 Data Sources and Analysis Flow

Branch: `codex/modular-00-data-sources`

Primary files discussed so far:

- `python_notebooks/notebooks/00_load_div30_div90_raw_to_anndata.ipynb`
- `scripts/cellbender.sh`
- `python_notebooks/src/mge_organoid_python/data_sources.py`
- `python_notebooks/src/mge_organoid_python/notebook00_workflow.py`
- `python_notebooks/src/mge_organoid_python/notebook00_plots.py`
- `slurm_templates/13_execute_notebook00_source.sbatch.template`

## Operating Rules For Notebook 00 Execution

Notebook 00 real matrix-loading runs must execute on Great Lakes through Slurm,
not directly on the login node. The active execution path is the notebook plus
the Python modules in this repo:

```text
python_notebooks/notebooks/00_load_div30_div90_raw_to_anndata.ipynb
python_notebooks/src/mge_organoid_python/data_sources.py
python_notebooks/src/mge_organoid_python/notebook00_workflow.py
python_notebooks/src/mge_organoid_python/notebook00_plots.py
slurm_templates/13_execute_notebook00_source.sbatch.template
```

Use this conda/Jupyter environment:

```text
/home/elcrespo/miniconda3/envs/mge-organoid-python/bin
```

Use this repo and data root unless the user explicitly says otherwise:

```text
REPO_ROOT=/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
PROJECT_ROOT=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
```

Default Slurm execution command:

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline

sbatch --export=ALL,NOTEBOOK00_LOAD_MATRICES=1,NOTEBOOK00_SAVE_PLOTS=1,NOTEBOOK00_SHOW_PLOTS=0 \
  slurm_templates/13_execute_notebook00_source.sbatch.template
```

Use `NOTEBOOK00_SHOW_PLOTS=0` for normal runs that save PNGs but keep the
executed notebook smaller. Use `NOTEBOOK00_SHOW_PLOTS=1` only when the user
wants the executed notebook, or the repo notebook copied from it, to show
inline plots.

Allowed:

- Submit Notebook 00 runs through `slurm_templates/13_execute_notebook00_source.sbatch.template`.
- Monitor Slurm jobs with `squeue` and `sacct`.
- Inspect Slurm logs under `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/`.
- Check executed notebooks under `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/executed/`.
- Check run outputs under `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/<RUN_LABEL>/`.

Not allowed unless the user explicitly requests it:

- Do not run real matrix-loading Notebook 00 execution directly on the login node.
- Do not rerun CellBender denoising.
- Do not write large data outputs inside the GitHub checkout.
- Do not overwrite the repo notebook with an executed notebook unless the user asks for inline outputs in the repo notebook.
- Do not switch away from the supported Notebook 00 sources without user approval.

Supported Notebook 00 sources:

```text
cellranger_filtered
cellbender_denoised
```

Common source-specific runs:

```bash
# Cell Ranger filtered primary manual_ec run
sbatch --export=ALL,NOTEBOOK00_ACTIVE_SOURCE=cellranger_filtered,NOTEBOOK00_LOAD_MATRICES=1,NOTEBOOK00_SAVE_PLOTS=1,NOTEBOOK00_SHOW_PLOTS=0 \
  slurm_templates/13_execute_notebook00_source.sbatch.template

# Existing CellBender-denoised comparison run
sbatch --export=ALL,NOTEBOOK00_ACTIVE_SOURCE=cellbender_denoised,NOTEBOOK00_LOAD_MATRICES=1,NOTEBOOK00_SAVE_PLOTS=1,NOTEBOOK00_SHOW_PLOTS=0 \
  slurm_templates/13_execute_notebook00_source.sbatch.template
```

If the user asks for the repo notebook to show plots inline, rerun with:

```bash
sbatch --export=ALL,NOTEBOOK00_LOAD_MATRICES=1,NOTEBOOK00_SAVE_PLOTS=1,NOTEBOOK00_SHOW_PLOTS=1 \
  slurm_templates/13_execute_notebook00_source.sbatch.template
```

After that job completes and the executed notebook has no error outputs, copy
the executed notebook back to:

```text
python_notebooks/notebooks/00_load_div30_div90_raw_to_anndata.ipynb
```

Only do this copy-back when the user explicitly wants inline outputs in the repo
notebook.

## Current State For Next Chat

The branch now has a focused Notebook 00 `manual_ec` driver and supporting
Python modules. The forward path no longer builds on the previous custom/MAD
filtering approach, PCA/UMAP/clustering, Harmony, or marker-panel workflow.

The current supported Notebook 00 sources are:

```text
cellranger_filtered: primary manual_ec run from Cell Ranger filtered matrices
cellbender_denoised: comparison manual_ec run from existing CellBender H5 files
```

Important: do not rerun CellBender denoising. The CellBender comparison uses the
existing fixed input files under:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/clean_adata/*_cellbender_denoised.h5
```

Execution should stay on Great Lakes through Slurm, using the conda env through
the notebook template:

```text
/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/jupyter
slurm_templates/13_execute_notebook00_source.sbatch.template
```

The proof scripts are validation scaffolding only. The core notebook runtime is:

```text
python_notebooks/notebooks/00_load_div30_div90_raw_to_anndata.ipynb
python_notebooks/src/mge_organoid_python/data_sources.py
python_notebooks/src/mge_organoid_python/notebook00_workflow.py
python_notebooks/src/mge_organoid_python/notebook00_plots.py
slurm_templates/13_execute_notebook00_source.sbatch.template
```

The current notebook should be continued as the user-facing driver. Most edits
should go into the three Python modules unless the change is specifically about
which source, samples, toggles, or plots a user wants to run.

Completed in this branch:

```text
manual_ec filtered Cell Ranger run for the six DIV30 core samples
manual_ec existing-CellBender-H5 run for the same six DIV30 core samples
manual_ec Cell Ranger filtered vs CellBender-denoised comparison
Notebook 00 source-layer cleanup after manual_ec validation
full post-cleanup Notebook 00 rerun through Slurm
run-specific table and plot directories under results/notebook00/<RUN_LABEL>
executed notebooks under results/notebook00/executed
Slurm-safe comparison-label parsing using ":" or ";" separators
per-sample manual_ec QC plot grids saved as 1xN sample panels
```

Not completed yet:

```text
Notebook 00 freeze as a QC/filter/checkpoint-only notebook
combined and per-sample .h5ad checkpoint writing
collaborator review of manual_ec outputs/plots
raw export manifest cleanup, if no longer needed
physical deletion/archive of remaining historical notes and old MAD/custom filtering scaffolding outside the forward path
```

## Notebook 00 Cleanup Confirmed

The active Notebook 00 code path now supports only the two sources used by the
validated `manual_ec` runs:

```text
cellranger_filtered
cellbender_denoised
```

Cleanup performed after the successful Slurm runs:

```text
removed cellranger_raw from the active Notebook 00 source layer
removed raw Cell Ranger matrix-source resolver support
removed expected_raw_h5ad_path and raw_h5ad reporting from CellBender source tables
removed the legacy load_dataset() tuple wrapper
removed unused DatasetLoadResult.available_samples and DatasetLoadResult.has_skipped_samples properties
deleted python_notebooks/scripts/prove_00_raw_filtered_sources.py
updated remaining Notebook 00 proof scripts so they only refer to active sources
```

The current helper modules required by the active notebook are:

```text
python_notebooks/src/mge_organoid_python/data_sources.py
python_notebooks/src/mge_organoid_python/notebook00_workflow.py
python_notebooks/src/mge_organoid_python/notebook00_plots.py
```

Validation after cleanup used the `mge-organoid-python` conda env:

```text
python -m compileall -q python_notebooks/src/mge_organoid_python python_notebooks/scripts
notebook JSON load passed
source availability passed for cellranger_filtered: 6 available samples
source availability passed for cellbender_denoised: 6 available samples
cellranger_raw is rejected by normalize_data_source()
git diff --check passed
```

No CellBender denoising was run during cleanup.

## Latest Post-Cleanup Slurm Rerun

The full Notebook 00 was rerun after source-layer cleanup. These are pre-freeze
outputs from the currently checked-in notebook. The next confirmed direction is
to remove Notebook 00 HVG/preprocess work and replace it with checkpoint
writing, as described in `Pending Notebook 00 Freeze Contract`.

```text
51191663 nb00-manualec-filtered-postclean    COMPLETED 0:0  00:03:38
51191664 nb00-manualec-cellbender-postclean  COMPLETED 0:0  00:08:03
51191665 nb00-manualec-compare-postclean     COMPLETED 0:0  00:00:16
```

Filtered Cell Ranger post-cleanup run:

```text
RUN_LABEL=cellranger_filtered_manual_ec_div30_core_samples_postcleanup
output=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/cellranger_filtered_manual_ec_div30_core_samples_postcleanup
executed=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/executed/00_load_div30_div90_raw_to_anndata.cellranger_filtered_manual_ec_div30_core_samples_postcleanup.executed.ipynb
manual_ec retained cells=97,658
final genes=17,486
highly variable genes=2,235
```

Existing CellBender-denoised H5 post-cleanup run:

```text
RUN_LABEL=cellbender_denoised_manual_ec_div30_core_samples_postcleanup
output=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/cellbender_denoised_manual_ec_div30_core_samples_postcleanup
executed=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/executed/00_load_div30_div90_raw_to_anndata.cellbender_denoised_manual_ec_div30_core_samples_postcleanup.executed.ipynb
manual_ec retained cells=100,674
final genes=18,549
highly variable genes=2,345
```

Post-cleanup comparison run:

```text
RUN_LABEL=manual_ec_cellranger_filtered_vs_cellbender_denoised_div30_postcleanup
output=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/manual_ec_cellranger_filtered_vs_cellbender_denoised_div30_postcleanup
executed=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/executed/00_load_div30_div90_raw_to_anndata.manual_ec_cellranger_filtered_vs_cellbender_denoised_div30_postcleanup.executed.ipynb
manual_ec_source_comparison_summary.tsv rows=12 data rows
manual_ec_hvg_overlap.tsv Jaccard=0.933305192064162
```

Plots confirmed in both full post-cleanup source runs:

```text
highest_expr_genes_top20.png
manual_ec_qc_violin.png
manual_ec_scatter_total_counts_pct_counts_mt.png
manual_ec_scatter_total_counts_n_genes_by_counts.png
manual_ec_highly_variable_genes.png
```

No CellBender denoising was run during the post-cleanup rerun. The
CellBender-labeled job only read existing `clean_adata/*_cellbender_denoised.h5`
files.

## Current Checked-In Bioinformatics Checkpoint

The checked-in Notebook 00 currently reaches the early single-cell
QC/preprocessing checkpoint:

```text
matrix loading
sample metadata annotation
QC metric calculation
QC visualization
manual_ec cell filtering
total-count normalization
log1p transform
highly variable gene selection
Cell Ranger filtered vs existing CellBender-denoised source comparison
```

This is the state to change next. After the freeze update, Notebook 00 should
stop earlier:

```text
matrix loading
sample metadata annotation
QC metric calculation
QC visualization
manual_ec cell filtering
combined and per-sample .h5ad checkpoint creation
```

Notebook 00 should not include these biological analysis stages:

```text
HVG selection
cell-cycle scoring
CC.Difference creation
scaling/regression
PCA
neighbors
UMAP
clustering
marker analysis
cell type annotation
differential expression
trajectory/pseudotime
integration/batch correction
```

Notebook 01 should start from the Notebook 00 checkpoints with Seurat-v3 HVG
selection from counts, cell-cycle scoring, `CC.Difference`, regression, scaling,
PCA, neighbors, UMAP, and clustering.

Sections below that discuss the original raw/filtered proof, raw export,
pre-cleanup runs, or old MAD/custom filtering are retained as historical context
only. They are not the active Notebook 00 workflow and should not be used as
tomorrow's run instructions.

## Verified Slurm Runs

All notebook execution below used Slurm and the `mge-organoid-python` conda env.

```text
51189568 nb00-manualec-filtered      COMPLETED 0:0
51189668 nb00-manualec-cellbender    COMPLETED 0:0
51189820 nb00-manualec-compare-v2    COMPLETED 0:0
```

Filtered Cell Ranger run:

```text
RUN_LABEL=cellranger_filtered_manual_ec_div30_core_samples
output=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/cellranger_filtered_manual_ec_div30_core_samples
executed=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/executed/00_load_div30_div90_raw_to_anndata.cellranger_filtered_manual_ec_div30_core_samples.executed.ipynb
manual_ec retained cells=97,658
final genes=17,486
highly variable genes=2,235
```

Existing CellBender-denoised H5 run:

```text
RUN_LABEL=cellbender_denoised_manual_ec_div30_core_samples
output=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/cellbender_denoised_manual_ec_div30_core_samples
executed=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/executed/00_load_div30_div90_raw_to_anndata.cellbender_denoised_manual_ec_div30_core_samples.executed.ipynb
manual_ec retained cells=100,674
final genes=18,549
highly variable genes=2,345
```

Comparison run:

```text
RUN_LABEL=manual_ec_cellranger_filtered_vs_cellbender_denoised_div30
output=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/manual_ec_cellranger_filtered_vs_cellbender_denoised_div30
executed=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/executed/00_load_div30_div90_raw_to_anndata.manual_ec_cellranger_filtered_vs_cellbender_denoised_div30.executed.ipynb
manual_ec_source_comparison_summary.tsv rows=12 data rows
manual_ec_hvg_overlap.tsv Jaccard=0.9333051920641621
```

When submitting comparison-only runs through `sbatch --export`, use `:` or `;`
between labels so Slurm does not split the value as separate export tokens:

```text
NOTEBOOK00_COMPARE_RUN_LABELS=cellranger_filtered_manual_ec_div30_core_samples:cellbender_denoised_manual_ec_div30_core_samples
NOTEBOOK00_LOAD_MATRICES=0
```

## Why This Exists

The original Notebook 00 mixed several different responsibilities in one long,
stateful workflow:

1. find the project checkout
2. read sample metadata
3. resolve Cell Ranger raw or filtered matrix paths
4. export raw `.h5ad` files for CellBender
5. run QC metric annotation
6. apply QC filtering
7. concatenate samples
8. preprocess, embed, integrate, cluster, and plot

That was workable for an early teaching notebook, but it was hard to scale
because changing the input data source also changed the mental model of the rest
of the notebook. The current branch addresses this by making the active data
source explicit and making each analysis stage callable on demand.

## Original Data Flow Found During Investigation

The original notebook and Slurm script did not use one single "project folder"
concept.

### Notebook project root

Notebook 00 defines Python `PROJECT_ROOT` by walking upward from the current working directory until it finds:

```text
metadata/div30_div90_sample_id_to_biolabel_map.tsv
```

In local development this usually points at the GitHub checkout, for example:

```text
/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
```

The builder then reads:

```text
PROJECT_ROOT/metadata/div30_div90_sample_id_to_biolabel_map.tsv
```

### Actual Cell Ranger matrix locations

The Cell Ranger matrix paths are not inferred from the GitHub checkout. They are derived from the `per_sample_metrics_csv` column in the metadata TSV.

For each selected sample, the notebook uses:

```python
count_dir = Path(per_sample_metrics_csv).parent / "count"
```

Then it resolves one of these directories:

```text
sample_raw_feature_bc_matrix
raw_feature_bc_matrix
sample_filtered_feature_bc_matrix
filtered_feature_bc_matrix
```

The actual per-cell matrix read is:

```python
sc.read_10x_mtx(matrix_dir, var_names="gene_symbols", make_unique=True)
```

### Original active notebook analysis source

The original builder was instantiated with:

```python
matrix_source="filtered"
```

Therefore this call:

```python
adata_names, adata_list = anndata_builder.per_sample_anndata_list()
```

loaded Cell Ranger filtered 10x matrix directories. The QC plots, QC filtering,
preprocessing, PCA, UMAP, Harmony, Leiden, and comparison plots operated on that
filtered Cell Ranger data.

### Raw export for CellBender

The notebook also has:

```python
save_raw_adatas_from_builder(...)
```

That function ignores the builder's default filtered setting and explicitly loads:

```python
matrix_source="raw"
```

It writes `.h5ad` files to:

```text
$PROJECT_ROOT/raw_adata
```

Important: this is the shell environment variable `PROJECT_ROOT`, not necessarily the notebook's Python `PROJECT_ROOT` variable.

### CellBender script root

`scripts/cellbender.sh` hardcodes:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
```

By default it reads:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/raw_adata/*.h5ad
```

It writes:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/clean_adata/<input_basename>_cellbender_denoised.h5
```

It also writes logs and work files under:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/cellbender_work
```

### Original critical gap

CellBender output is produced by `scripts/cellbender.sh`, but the original
Notebook 00 did not read those CellBender outputs back into the analysis flow.

Original state:

- raw Cell Ranger data can be exported for CellBender
- CellBender can produce denoised `.h5` outputs
- downstream notebook analysis still uses Cell Ranger filtered matrices

## Target Architecture

The next notebook version should make data roots and data sources explicit.

### Separate roots

Use two names instead of one overloaded `PROJECT_ROOT`:

```python
REPO_ROOT = find_repo_root(...)
DATA_ROOT = Path(os.environ.get("MGE_DATA_ROOT", "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"))
```

Expected meanings:

- `REPO_ROOT`: GitHub checkout containing notebooks, scripts, metadata, and source code.
- `DATA_ROOT`: large data directory containing raw `.h5ad`, CellBender outputs, logs, and work files.

### Explicit data sources

Introduce a small data-source registry:

```python
DATA_SOURCES = {
    "cellranger_filtered": {
        "kind": "10x_mtx",
        "matrix_source": "filtered",
    },
    "cellranger_raw": {
        "kind": "10x_mtx",
        "matrix_source": "raw",
    },
    "cellbender_denoised": {
        "kind": "cellbender_h5",
        "base_dir": DATA_ROOT / "clean_adata",
        "suffix": "_cellbender_denoised.h5",
    },
}
```

Then select the active input once:

```python
ACTIVE_SOURCE = "cellranger_filtered"
```

Later this can become:

```python
ACTIVE_SOURCE = "cellbender_denoised"
```

### One loader interface

The notebook should expose one public loader call:

```python
adata_names, adata_list = load_dataset(ACTIVE_SOURCE)
```

Under the hood:

- `cellranger_filtered` and `cellranger_raw` use `sc.read_10x_mtx(...)`
- `cellbender_denoised` uses the CellBender `.h5` output path and an appropriate Scanpy reader
- all loaded objects get the same `.obs` columns: `DIV`, `run_sample_id`, `cell_line`, `matrix_source`, and source path

### Manifest bridge

Raw export should create a manifest that maps notebook sample IDs to files used by CellBender.

Recommended file:

```text
DATA_ROOT/raw_adata/raw_adata_manifest.tsv
```

Recommended columns:

```text
run_sample_id
cell_line
raw_h5ad_path
cellbender_output_h5
source_matrix_dir
n_cells
n_genes
```

This manifest makes the round trip traceable:

```text
metadata TSV -> Cell Ranger raw matrix -> raw .h5ad -> CellBender .h5 -> notebook analysis
```

## Modular Notebook Shape

Notebook 00 should become a driver notebook with sections like this:

1. Configuration
2. Data source selection
3. Load selected dataset
4. QC metric annotation
5. QC filtering, optional
6. Concatenation
7. Preprocessing
8. Embedding
9. Optional Harmony
10. Plot one sample, all samples, or comparisons

The important design rule:

Functions should accept `adata`, `adata_list`, or `adata_names` as inputs instead of assuming hidden global variables.

## Plotting Goal

Plots should be callable individually rather than all firing in one long sequence.

Examples of the intended usage:

```python
plot_qc_sample(qc_runner, run_sample_id="9853-MW-1", metric="pct_counts_mt")
plot_qc_grid(qc_runner, metric="log1p_total_counts")
plot_embedding(adata, basis="umap", color=["run_sample_id", "cell_line"])
plot_marker_panel(adata, markers=["SOX2", "DLX2", "GAD1"])
compare_embeddings(adata_uncorrected, adata_harmony, color="run_sample_id")
```

## Success Criteria

The refactor is successful when all of these are true:

- Notebook 00 clearly distinguishes `REPO_ROOT` from `DATA_ROOT`.
- A reader can see exactly which data source is active from one variable such as `ACTIVE_SOURCE`.
- The notebook can load Cell Ranger filtered data through the new source interface.
- The notebook can load CellBender-denoised data through the source interface for at least the available common samples.
- CellBender output locations are represented in code and missing outputs remain visible in `source_table`.
- The forward `manual_ec` path can run without depending on the historical MAD annotation/filtering path.
- The same `manual_ec` checkpoint can run on `cellranger_filtered` and `cellbender_denoised`.
- The user can compare CellBender-denoised vs Cell Ranger filtered outputs by `run_sample_id`, `cell_line`, cell counts, QC metrics, retained cells, and HVG results.
- The notebook still preserves key `.obs` metadata fields needed for the two-source comparison.
- Historical UMAP, clustering, Harmony, marker panel, MAD filtering, and raw source analysis code is either removed from the forward path or clearly flagged for deletion.

## Implementation Boundaries

The earlier first-pass boundary was to avoid changing scientific thresholds while
modularizing paths and source switching. That boundary is now superseded for the
next analysis pass.

Do not build new work on the historical custom/MAD filtering path. It can remain
in the repository until cleanup, but `manual_ec` should not require:

- `annotate_mad_qc(...)`
- `qc_mad_pass`
- `filter_qc_pass_samples(...)`
- the old `PreprocessSettings()` defaults
- UMAP, clustering, Harmony, or marker plotting

Keep these stable unless explicitly changed:

- `ACTIVE_SOURCE=cellranger_filtered` for `manual_ec`
- `ACTIVE_SOURCE=cellbender_denoised` for the CellBender comparison run
- mitochondrial prefix `MT-`
- target sample IDs
- target DIV selection
- CellBender command-line behavior
- existing CellBender output naming convention

The next pass should focus on reproducing the requested Scanpy QC/preprocess
checkpoint cleanly, with exact parameters and run-specific outputs, then applying
that same checkpoint to CellBender-denoised inputs for comparison.

## Strip/Delete Candidates

Do not delete files until the user approves, but flag these as outside the
forward Notebook 00 path:

```text
historical MAD QC annotation and filtering path
historical basic preprocess/PCA/neighbors/UMAP path
Harmony/integration blocks
Leiden/clustering blocks
marker panel plotting blocks
raw Cell Ranger full-analysis path
generic cross-run comparison scaffolding unrelated to filtered-vs-CellBender
proof scripts once their validation role is replaced by the manual_ec comparison
old run outputs that only document the historical analysis path
```

Keep these pieces because they support the retained comparison:

```text
Cell Ranger filtered source loading
CellBender denoised source loading
source availability and missing-source reporting
sample metadata propagation: DIV, run_sample_id, biological_label, cell_line
run-specific results/notebook00/<RUN_LABEL>/tables and plots directories
CellBender file naming/location conventions
```

## Risks To Watch

- CellBender `.h5` files may need a different reader or postprocessing than Cell Ranger 10x directories.
- Raw `.h5ad` filenames are currently derived from matrix path prefixes, so duplicate-safe naming and `run_sample_id` mapping must remain explicit.
- The notebook currently has both Python `PROJECT_ROOT` and shell `$PROJECT_ROOT`; these must not remain ambiguous.
- Large data should stay under `DATA_ROOT`, not inside the GitHub checkout.
- Any code that rewrites `.obs["run_sample_id"]` from `batch` should be checked carefully so sample identity is preserved.

## Implementation Direction Used

The refactor followed this direction: move reusable workflow code under:

```text
python_notebooks/src/mge_organoid_python/
```

Then keep Notebook 00 as the readable driver. The actual current imports use
the Notebook 00 modules:

```python
from mge_organoid_python.data_sources import Notebook00SourceConfig, load_dataset_result
from mge_organoid_python.notebook00_workflow import (
    PreprocessSettings,
    calculate_qc_metrics,
    concat_samples,
    filter_qc_pass_samples,
    preprocess_basic,
    run_neighbors_umap,
)
from mge_organoid_python.notebook00_plots import (
    PlotConfig,
    plot_embedding,
    plot_marker_panel,
    plot_qc_scatter,
    plot_qc_violin,
    plot_sample_counts,
    plot_source_availability,
)
```

The current branch now has that split:

```text
data_sources.py       path/source resolution, missing-aware loading, source tables
notebook00_workflow.py QC annotation, filtering, concatenation, preprocess, UMAP
notebook00_plots.py    saved/inline plotting helpers
```

This keeps the notebook editable while moving fragile path and analysis
mechanics into testable code.

## Initial Raw vs Filtered Proof

The first implementation pass added:

```text
python_notebooks/src/mge_organoid_python/data_sources.py
python_notebooks/scripts/prove_00_raw_filtered_sources.py
slurm_templates/10_prove_00_raw_filtered_sources.sbatch.template
```

The proof must run on a Great Lakes compute node, not on the login node, because it loads real 10x matrices. Submit it with:

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
sbatch slurm_templates/10_prove_00_raw_filtered_sources.sbatch.template
```

It loads the same sample through both supported sources:

```text
cellranger_raw
cellranger_filtered
```

The proof passes only if:

- raw and filtered matrix directories both exist
- raw and filtered resolve to different directories
- both sources load through Scanpy
- raw has at least as many cell barcodes as filtered
- filtered genes are present in the raw source, even if raw has additional features/genes
- loaded `.obs` contains stable sample metadata such as `run_sample_id` and `cell_line`

Default proof output:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/prove-00-raw-filtered-<job_id>.out
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/prove_00_raw_filtered_sources.tsv
```

## Initial CellBender Location Proof

The CellBender proof builds on the same source/path module and adds:

```text
python_notebooks/scripts/prove_00_cellbender_locations.py
slurm_templates/11_prove_00_cellbender_locations.sbatch.template
```

Submit it with:

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
sbatch slurm_templates/11_prove_00_cellbender_locations.sbatch.template
```

Default sample set:

```text
9853-MW-1
9853-MW-2
9853-MW-3
9853-MW-4
9853-MW-5
9853-MW-6
```

For each sample, it derives:

```text
DATA_ROOT/raw_adata/<run_sample_id>.h5ad
DATA_ROOT/clean_adata/<run_sample_id>_cellbender_denoised.h5
```

This matches the naming behavior in `scripts/cellbender.sh`, where CellBender
outputs are named from the input basename.

The proof reports:

- expected raw `.h5ad` input path
- whether that raw input exists
- expected CellBender denoised `.h5` output path
- whether that CellBender output exists
- output file size for existing outputs
- whether existing outputs can be opened as HDF5
- top-level HDF5 keys and matrix shape when available

This proof does not require every sample to exist. The Slurm template currently
uses `--require-any-existing`, meaning it fails only if no requested CellBender
outputs exist at all or if an existing output cannot be opened as HDF5.

Default proof output:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/prove-00-cellbender-locations-<job_id>.out
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/prove_00_cellbender_locations.tsv
```

Observed proof run:

```text
job_id: 51109820
compute_node: gl3028.arc-ts.umich.edu
log: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/prove-00-cellbender-locations-51109820.out
report: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/prove_00_cellbender_locations.tsv
```

Observed result:

```text
raw_h5ad_found: 6
cellbender_outputs_found: 5
cellbender_outputs_missing: 1
missing_output: 9853-MW-6
```

Existing CellBender outputs for `9853-MW-1` through `9853-MW-5` opened as
HDF5 and had top-level keys:

```text
droplet_latents
global_latents
matrix
metadata
```

Observed matrix shapes:

```text
9853-MW-1: 37143x675447
9853-MW-2: 37143x514516
9853-MW-3: 37143x690463
9853-MW-4: 37143x627418
9853-MW-5: 37143x712267
```

## Current Status Snapshot

Last checked from the login node:

```text
date: 2026-05-29
branch: codex/modular-00-data-sources
login_node: gl-login2.arc-ts.umich.edu
```

What has been done:

- Created the isolated branch `codex/modular-00-data-sources`.
- Added this handoff file to explain the path problem, target architecture, and success criteria.
- Added `python_notebooks/src/mge_organoid_python/data_sources.py` for Notebook 00 source/path handling.
- Exposed the new Notebook 00 source helpers from `python_notebooks/src/mge_organoid_python/__init__.py`.
- Added a compute-node raw-vs-filtered proof script and Slurm template.
- Proved raw vs filtered source switching on Great Lakes with job `51109221`.
- Added a compute-node CellBender location proof script and Slurm template.
- Proved the expected CellBender path model with job `51109820`.
- Confirmed existing CellBender outputs for `9853-MW-1` through `9853-MW-5` are HDF5-readable.
- Added missing-aware source reporting for downstream analysis/plotting.
- Proved that missing-aware reporting identifies available and skipped CellBender samples with job `51111674`.

Raw vs filtered proof result:

```text
job_id: 51109221
compute_node: gl3064.arc-ts.umich.edu
result: PASS
sample: 9853-MW-1
raw_n_obs: 675447
filtered_n_obs: 18047
raw_n_vars: 37143
filtered_n_vars: 18082
filtered_vars_missing_from_raw: 0
```

CellBender status for `9853-MW-6`:

```text
current_job: 51106991_5
state: RUNNING
compute_node: gl1008.arc-ts.umich.edu
input: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/raw_adata/9853-MW-6.h5ad
target_output: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/clean_adata/9853-MW-6_cellbender_denoised.h5
```

The final primary `.h5` output for `9853-MW-6` had not appeared at the last
check. Existing sidecar/partial outputs:

```text
9853-MW-6_cellbender_denoised_posterior.h5
9853-MW-6_cellbender_denoised.pdf
9853-MW-6_cellbender_denoised_cell_barcodes.csv
9853-MW-6_cellbender_denoised.log
```

The current `9853-MW-6` run has passed posterior generation and is computing
denoised counts using the MCKP estimator. The log showed:

```text
Computing denoised counts using mckp estimator
Working on chunk (1/7)
[10.86 mins per chunk]
```

An older `9853-MW-6` CellBender job failed:

```text
job_id: 50803829_5
state: FAILED
failure: ValueError: Workflow hash does not match that of checkpoint.
```

This earlier failure explains why the final output was missing before the
current rerun.

## Missing-Aware Downstream Proof

The missing `9853-MW-6` CellBender output was used as a deliberate test case
for downstream flexibility. The source layer now records availability before
matrix loading, so notebooks can show exactly which samples are usable and which
were skipped.

Added:

```text
python_notebooks/scripts/prove_00_missing_aware_downstream.py
slurm_templates/12_prove_00_missing_aware_downstream.sbatch.template
```

Relevant source-layer fields:

```text
source_path
source_exists
load_status
skip_reason
loaded_in_memory
```

Relevant helper:

```python
result = load_dataset_result(config, load_matrices=False)
source_table = result.source_table
summary_df = result.availability_summary()
```

Observed proof run:

```text
job_id: 51111674
compute_node: gl3180.arc-ts.umich.edu
result: PASS
```

Observed CellBender availability:

```text
available_samples:
  9853-MW-1
  9853-MW-2
  9853-MW-3
  9853-MW-4
  9853-MW-5

skipped_samples:
  9853-MW-6

skip_reason:
  missing_cellbender_output_h5
```

Proof outputs:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/prove-00-missing-aware-51111674.out
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/prove_00_missing_aware_source_status.tsv
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/prove_00_missing_aware_summary.tsv
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/prove_00_missing_aware_source_status.png
```

This is the downstream behavior the notebook should use:

1. Build `source_table`.
2. Display/log available and skipped samples.
3. Plot source availability.
4. Load only available samples when strict mode is disabled.
5. Preserve skipped-sample reporting alongside every downstream analysis result.

## Notebook 00 Driver Refactor

Notebook 00 has been rewritten as a thin driver notebook:

```text
python_notebooks/notebooks/00_load_div30_div90_raw_to_anndata.ipynb
```

The old in-notebook builder class, QC class, raw export utility, Harmony block,
and long plot sequence were removed from the notebook body. Reusable behavior now
lives in:

```text
python_notebooks/src/mge_organoid_python/data_sources.py
python_notebooks/src/mge_organoid_python/notebook00_workflow.py
python_notebooks/src/mge_organoid_python/notebook00_plots.py
```

The notebook now has 19 cells and follows this flow:

1. imports
2. configuration
3. source availability report, without loading matrices
4. load available samples
5. QC annotation
6. QC filtering and concatenation
7. preprocessing and UMAP
8. saved/displayed plots
9. marker panel
10. generated object summary

Default source:

```python
ACTIVE_SOURCE = "cellranger_filtered"
```

The notebook can switch to:

```python
ACTIVE_SOURCE = "cellranger_raw"
ACTIVE_SOURCE = "cellbender_denoised"
```

Missing-source behavior is controlled by:

```python
STRICT_MISSING_SOURCES = False
```

With strict mode disabled, missing samples are logged in `source_table` and
skipped before matrix loading. This is the intended behavior for incomplete
CellBender batches or deliberately omitted samples.

The notebook now only uses objects created in the current run:

```text
source_table
sample_adata_names
sample_adata_list
adata_names
adata_list
combined_adata
manual_ec_adata
```

`sample_adata_names` and `sample_adata_list` are the explicit per-sample
objects. `adata_names` and `adata_list` remain backward-compatible aliases for
older notebook snippets. `combined_adata` is the concatenated loaded object
before manual_ec filtering; `manual_ec_adata` is the filtered and preprocessed
working object.

Plot/table output locations are deterministic and run-specific:

```text
RUN_DIR = DATA_ROOT/results/notebook00/<RUN_LABEL>
TABLE_DIR = RUN_DIR/tables
PLOT_DIR = RUN_DIR/plots
```

With the default source, this resolves to:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/cellranger_filtered_div30_core_samples/tables
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/cellranger_filtered_div30_core_samples/plots
```

Plot helpers save PNGs and display inline through:

```python
plot_config = PlotConfig.from_root(DATA_ROOT, run_label=RUN_LABEL, show=True, save=True)
```

The notebook JSON validates, every code cell parses, and the refactored
notebook has executed successfully on Great Lakes through Slurm with the
`cellranger_filtered` source as a full analysis run.

### Notebook User Inputs

Interactive users can edit the configuration cell directly. Batch users can set
the same inputs with environment variables.

Core inputs:

```text
NOTEBOOK00_ACTIVE_SOURCE
NOTEBOOK00_TARGET_DIVS
NOTEBOOK00_TARGET_RUN_SAMPLE_IDS
NOTEBOOK00_STRICT_MISSING_SOURCES
NOTEBOOK00_LOAD_MATRICES
NOTEBOOK00_RUN_LABEL
```

Analysis toggles:

```text
NOTEBOOK00_APPLY_QC_FILTER
NOTEBOOK00_RUN_PREPROCESS
NOTEBOOK00_RUN_UMAP
```

Plot toggles:

```text
NOTEBOOK00_SHOW_PLOTS
NOTEBOOK00_SAVE_PLOTS
```

Environment/root handling:

```text
REPO_ROOT=/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
PROJECT_ROOT=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
CONDA_ENV_BIN=/home/elcrespo/miniconda3/envs/mge-organoid-python/bin
```

Important naming note:

```text
REPO_ROOT    = GitHub checkout with notebooks, source modules, metadata, Slurm templates
PROJECT_ROOT = large HPC project/data root with raw_adata, clean_adata, logs, results
```

Inside the Python notebook, `PROJECT_ROOT` is still used as the data root
because existing scripts and templates already use that shell variable. The
handoff's conceptual distinction remains: repo files live under `REPO_ROOT`;
large data and results live under `PROJECT_ROOT`.

The Slurm runner executes the notebook with:

```text
kernel: mge-organoid-python
jupyter: /home/elcrespo/miniconda3/envs/mge-organoid-python/bin/jupyter
output: $PROJECT_ROOT/results/notebook00/executed/<notebook>.<RUN_LABEL>.executed.ipynb
```

Reusable Slurm runner:

```text
slurm_templates/13_execute_notebook00_source.sbatch.template
```

Example source-only CellBender run:

```bash
sbatch \
  --job-name=exec-nb00-cellbender \
  --export=ALL,NOTEBOOK00_ACTIVE_SOURCE=cellbender_denoised,NOTEBOOK00_RUN_LABEL=cellbender_denoised_div30_core_samples_source_only,NOTEBOOK00_LOAD_MATRICES=0,NOTEBOOK00_RUN_PREPROCESS=0,NOTEBOOK00_RUN_UMAP=0,NOTEBOOK00_STRICT_MISSING_SOURCES=0 \
  slurm_templates/13_execute_notebook00_source.sbatch.template
```

Known-good Slurm commands from this branch:

```bash
sbatch \
  --job-name=exec-nb00-filtered \
  --export=ALL,NOTEBOOK00_ACTIVE_SOURCE=cellranger_filtered,NOTEBOOK00_RUN_LABEL=cellranger_filtered_div30_core_samples,NOTEBOOK00_LOAD_MATRICES=1,NOTEBOOK00_RUN_PREPROCESS=1,NOTEBOOK00_RUN_UMAP=1,NOTEBOOK00_STRICT_MISSING_SOURCES=0 \
  slurm_templates/13_execute_notebook00_source.sbatch.template

sbatch \
  --job-name=exec-nb00-raw \
  --export=ALL,NOTEBOOK00_ACTIVE_SOURCE=cellranger_raw,NOTEBOOK00_RUN_LABEL=cellranger_raw_div30_core_samples_source_only,NOTEBOOK00_LOAD_MATRICES=0,NOTEBOOK00_RUN_PREPROCESS=0,NOTEBOOK00_RUN_UMAP=0,NOTEBOOK00_STRICT_MISSING_SOURCES=0 \
  slurm_templates/13_execute_notebook00_source.sbatch.template

sbatch \
  --job-name=exec-nb00-cellbender \
  --export=ALL,NOTEBOOK00_ACTIVE_SOURCE=cellbender_denoised,NOTEBOOK00_RUN_LABEL=cellbender_denoised_div30_core_samples_source_only,NOTEBOOK00_LOAD_MATRICES=0,NOTEBOOK00_RUN_PREPROCESS=0,NOTEBOOK00_RUN_UMAP=0,NOTEBOOK00_STRICT_MISSING_SOURCES=0 \
  slurm_templates/13_execute_notebook00_source.sbatch.template
```

Run full matrix work only through Slurm/compute nodes. Source-only reports are
lightweight, but keeping the same Slurm path avoids environment drift between
interactive inspection and batch execution.

### Three-Source Notebook Execution

Three source modes were executed through the refactored notebook:

```text
cellranger_filtered
cellranger_raw
cellbender_denoised
```

Observed jobs:

```text
cellranger_filtered: 51112887, COMPLETED, 0:0, gl3078, 00:05:23
cellranger_raw:      51112888, COMPLETED, 0:0, gl3251, 00:00:16
cellbender_denoised: 51112889, COMPLETED, 0:0, gl3285, 00:00:16
```

The filtered source was run as a full analysis:

```text
NOTEBOOK00_LOAD_MATRICES=1
NOTEBOOK00_RUN_PREPROCESS=1
NOTEBOOK00_RUN_UMAP=1
```

The raw and CellBender sources were run in source-report-only mode:

```text
NOTEBOOK00_LOAD_MATRICES=0
NOTEBOOK00_RUN_PREPROCESS=0
NOTEBOOK00_RUN_UMAP=0
```

Reason: raw and current CellBender outputs are raw-droplet-scale matrices.
Running full concat/preprocess/UMAP on all raw droplets is not a sensible default
for interactive review.

Executed notebooks:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/executed/00_load_div30_div90_raw_to_anndata.cellranger_filtered_div30_core_samples.executed.ipynb
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/executed/00_load_div30_div90_raw_to_anndata.cellranger_raw_div30_core_samples_source_only.executed.ipynb
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/executed/00_load_div30_div90_raw_to_anndata.cellbender_denoised_div30_core_samples_source_only.executed.ipynb
```

Run output folders:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/cellranger_filtered_div30_core_samples
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/cellranger_raw_div30_core_samples_source_only
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/cellbender_denoised_div30_core_samples_source_only
```

Source summaries:

```text
cellranger_filtered: available=6
cellranger_raw:      available=6
cellbender_denoised: available=5, missing_source=1
```

Targeted outputs to inspect:

```text
Filtered full run:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/cellranger_filtered_div30_core_samples/tables/source_table.tsv
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/cellranger_filtered_div30_core_samples/tables/qc_filter_summary.tsv
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/cellranger_filtered_div30_core_samples/tables/embedding_report.tsv
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/cellranger_filtered_div30_core_samples/plots/source_availability.png
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/cellranger_filtered_div30_core_samples/plots/qc_violin_by_sample.png
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/cellranger_filtered_div30_core_samples/plots/qc_scatter_total_counts_pct_counts_mt.png
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/cellranger_filtered_div30_core_samples/plots/umap_sample_cellline_qc.png
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/cellranger_filtered_div30_core_samples/plots/umap_marker_panel.png

Raw source-only run:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/cellranger_raw_div30_core_samples_source_only/tables/source_table.tsv
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/cellranger_raw_div30_core_samples_source_only/tables/source_summary.tsv
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/cellranger_raw_div30_core_samples_source_only/plots/source_availability.png
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/cellranger_raw_div30_core_samples_source_only/plots/loaded_vs_skipped_counts.png

CellBender source-only run:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/cellbender_denoised_div30_core_samples_source_only/tables/source_table.tsv
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/cellbender_denoised_div30_core_samples_source_only/tables/source_summary.tsv
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/cellbender_denoised_div30_core_samples_source_only/plots/source_availability.png
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/cellbender_denoised_div30_core_samples_source_only/plots/loaded_vs_skipped_counts.png
```

The raw and CellBender source-only runs intentionally have empty downstream
analysis tables such as `qc_filter_summary.tsv` and `loaded_sample_shapes.tsv`
because `NOTEBOOK00_LOAD_MATRICES=0`. Their useful outputs are the source
tables, source summaries, and availability plots.

To inspect in VS Code, open the absolute paths above from the Explorer or use
the command palette action "File: Open File..." and paste the path.

Notebook execution logs:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/execute-notebook00-exec-nb00-filtered-51112887.out
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/execute-notebook00-exec-nb00-raw-51112888.out
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/execute-notebook00-exec-nb00-cellbender-51112889.out
```

## Pre-Cleanup Confirmation Runs

Submitted from `gl-login6.arc-ts.umich.edu` on 2026-05-31 before deleting the
historical MAD/custom filtering code.

Important: do not rerun CellBender denoising. The CellBender job below is a
Notebook 00 source-availability check against existing
`clean_adata/*_cellbender_denoised.h5` files only. It does not call
`scripts/cellbender.sh`, does not run CellBender, and does not create new
denoised `.h5` outputs.

Submitted jobs:

```text
51189298  exec-nb00-filtered-precleanup    cellranger_filtered full current notebook run
51189299  exec-nb00-cellbender-precleanup  cellbender_denoised source-only notebook run
```

Submission commands:

```bash
sbatch \
  --job-name=exec-nb00-filtered-precleanup \
  --export=ALL,NOTEBOOK00_ACTIVE_SOURCE=cellranger_filtered,NOTEBOOK00_RUN_LABEL=precleanup_cellranger_filtered_div30_core_samples,NOTEBOOK00_LOAD_MATRICES=1,NOTEBOOK00_RUN_PREPROCESS=1,NOTEBOOK00_RUN_UMAP=1,NOTEBOOK00_APPLY_QC_FILTER=1,NOTEBOOK00_STRICT_MISSING_SOURCES=0 \
  slurm_templates/13_execute_notebook00_source.sbatch.template

sbatch \
  --job-name=exec-nb00-cellbender-precleanup \
  --export=ALL,NOTEBOOK00_ACTIVE_SOURCE=cellbender_denoised,NOTEBOOK00_RUN_LABEL=precleanup_cellbender_denoised_div30_core_samples_source_only,NOTEBOOK00_LOAD_MATRICES=0,NOTEBOOK00_RUN_PREPROCESS=0,NOTEBOOK00_RUN_UMAP=0,NOTEBOOK00_APPLY_QC_FILTER=0,NOTEBOOK00_STRICT_MISSING_SOURCES=0 \
  slurm_templates/13_execute_notebook00_source.sbatch.template
```

Observed final status:

```text
51189298  COMPLETED  0:0  00:06:15  gl3138
51189299  COMPLETED  0:0  00:01:06  gl3065
```

Observed source availability:

```text
cellranger_filtered: available=6
cellbender_denoised: available=6
```

This supersedes the earlier missing-`9853-MW-6` CellBender status in this
handoff. The final CellBender denoised H5 now appears to be present for all six
requested DIV30 samples. Do not rerun CellBender; use these existing H5 outputs.

Pre-cleanup filtered MAD/custom retention summary, for historical confirmation
only:

```text
9853-MW-1: 18047 -> 13998 retained
9853-MW-2:  6060 ->  4086 retained
9853-MW-3: 17928 -> 15903 retained
9853-MW-4: 13047 ->  9280 retained
9853-MW-5: 26408 -> 20590 retained
9853-MW-6: 25530 -> 20399 retained
```

Expected executed notebooks:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/executed/00_load_div30_div90_raw_to_anndata.precleanup_cellranger_filtered_div30_core_samples.executed.ipynb
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/executed/00_load_div30_div90_raw_to_anndata.precleanup_cellbender_denoised_div30_core_samples_source_only.executed.ipynb
```

Expected run directories:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/precleanup_cellranger_filtered_div30_core_samples
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/notebook00/precleanup_cellbender_denoised_div30_core_samples_source_only
```

Key confirmation outputs:

```text
Filtered historical run:
  tables/source_summary.tsv
  tables/source_table.tsv
  tables/qc_filter_summary.tsv
  tables/qc_mad_thresholds.tsv
  tables/preprocess_report.tsv
  tables/embedding_report.tsv
  plots/qc_violin_by_sample.png
  plots/qc_scatter_total_counts_pct_counts_mt.png
  plots/umap_sample_cellline_qc.png
  plots/umap_marker_panel.png

CellBender existing-file source check:
  tables/source_summary.tsv
  tables/source_table.tsv
  plots/source_availability.png
  plots/loaded_vs_skipped_counts.png
```

Expected logs:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/execute-notebook00-exec-nb00-filtered-precleanup-51189298.out
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/execute-notebook00-exec-nb00-filtered-precleanup-51189298.err
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/execute-notebook00-exec-nb00-cellbender-precleanup-51189299.out
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/execute-notebook00-exec-nb00-cellbender-precleanup-51189299.err
```

For future work, CellBender should be treated as a fixed existing input source:

```text
Do not submit scripts/cellbender.sh.
Do not rerun CellBender remove-background.
Do not create new CellBender denoised H5 files unless the user explicitly asks.
Use existing DATA_ROOT/clean_adata/*_cellbender_denoised.h5 outputs for comparison.
```

## Implemented Direction: `manual_ec` Scanpy Path

Status: implemented and Slurm-validated for both Cell Ranger filtered matrices
and existing CellBender-denoised H5 files.

Forget the earlier attempt to make our own Notebook 00 filtering approach. The
historical MAD/custom path was useful for modularization and smoke testing, but
it is not the forward scientific workflow. Future work should keep focusing on
the explicit `manual_ec` Scanpy path for controlled review of Heyoon's
collaborator workflow choices.

Important source rule:

```text
manual_ec primary run uses ACTIVE_SOURCE=cellranger_filtered.
manual_ec comparison run uses ACTIVE_SOURCE=cellbender_denoised.
```

`manual_ec` is not a new physical data source. It is the analysis path applied
after loading matrices through the existing source layer. The first target is
Cell Ranger filtered input; the retained comparison target is CellBender-denoised
input. The data-source layer should remain useful, but `manual_ec` should not
depend on the historical MAD/custom filtering functions.

Implemented run identities:

```text
NOTEBOOK00_ACTIVE_SOURCE=cellranger_filtered
NOTEBOOK00_RUN_LABEL=cellranger_filtered_manual_ec_div30_core_samples

NOTEBOOK00_ACTIVE_SOURCE=cellbender_denoised
NOTEBOOK00_RUN_LABEL=cellbender_denoised_manual_ec_div30_core_samples
```

Do not add more toggles whose purpose is to preserve the older filtering method.
If old code remains temporarily, treat it as historical/cleanup debt rather than
an alternate supported path.

The new path has been stripped of forward dependencies on:

```text
annotate_mad_qc
qc_annotation_summary
qc_mad_pass
filter_qc_pass_samples
PreprocessSettings basic defaults
run_neighbors_umap
```

UMAP, clustering, Harmony, and marker plots can be reintroduced later after the
manual QC/preprocess checkpoint is proven.

### Retained Source Comparison

The notebook keeps a direct comparison between:

```text
cellranger_filtered manual_ec output
cellbender_denoised manual_ec output
```

This comparison uses common `run_sample_id` values only. If any CellBender
output is missing in a future run, the comparison should proceed with the
available common samples and write the missing sample into `source_table`.

Implemented comparison outputs:

```text
tables/manual_ec_source_comparison_summary.tsv
tables/manual_ec_qc_metric_comparison_by_sample.tsv
tables/manual_ec_hvg_overlap.tsv
plots/manual_ec_source_comparison_cell_counts.png
plots/manual_ec_source_comparison_qc_metrics.png
```

Minimum comparison fields:

```text
data_source
run_sample_id
cell_line
starting_n_cells
post_min_gene_cell_filter_n_cells
manual_ec_retained_n_cells
manual_ec_retained_pct
median_n_genes_by_counts
median_total_counts
median_pct_counts_mt
n_highly_variable_genes
```

### Implemented `manual_ec` QC/Plot Checkpoint

After loading and concatenating either selected source, the notebook reaches
this Scanpy checkpoint before any UMAP/cluster work:

```python
sc.pl.highest_expr_genes(adata, n_top=20)

adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True,
)

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
)
sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt")
sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts")
```

The plots are saved under the run-specific Notebook 00 plot directory using
stable names:

```text
highest_expr_genes_top20.png
manual_ec_qc_violin.png
manual_ec_scatter_total_counts_pct_counts_mt.png
manual_ec_scatter_total_counts_n_genes_by_counts.png
```

The notebook now has a separate `manual_ec Per-Sample QC Plot Grids` cell after
the combined manual_ec checkpoint. It loops through
`sample_adata_names`/`sample_adata_list`,
recomputes Scanpy QC metrics on each sample object, and saves per-sample 1xN
Scanpy grids:

```text
manual_ec_per_sample_highest_expr_genes_top20.png
manual_ec_per_sample_qc_violin_n_genes_by_counts.png
manual_ec_per_sample_qc_violin_total_counts.png
manual_ec_per_sample_qc_violin_pct_counts_mt.png
manual_ec_per_sample_scatter_total_counts_pct_counts_mt.png
manual_ec_per_sample_scatter_total_counts_n_genes_by_counts.png
```

### Implemented `manual_ec` Filtering

The requested manual cutoff is implemented as a named, documented filter:

```python
adata = adata[
    (adata.obs["n_genes_by_counts"] > 1500)
    & (adata.obs["n_genes_by_counts"] < 9000)
    & (adata.obs["total_counts"] < 30000)
    & (adata.obs["pct_counts_mt"] < 5),
    :
].copy()
```

The notebook records this output as `manual_ec`, including a table that reports
starting cells, retained cells, removed cells, retained percent, and the exact
thresholds:

```text
tables/manual_ec_filter_summary.tsv
tables/manual_ec_filter_parameters.tsv
```

It also includes these Scanpy tutorial-style count filters before the manual
cutoffs as an explicitly documented preliminary step:

```python
sc.pp.filter_cells(adata, min_genes=20)
sc.pp.filter_genes(adata, min_cells=3)
```

The implementation preserves enough reporting to distinguish:

```text
cellranger_filtered source loading
cellbender_denoised source loading
preliminary min_genes/min_cells filtering
manual_ec biological/QC cutoffs
historical MAD/custom filtering path, if still present in old outputs
```

### Current Checked-In `manual_ec` Preprocessing To Remove

The checked-in notebook currently performs this after `manual_ec` filtering:

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(
    adata,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5,
)
sc.pl.highly_variable_genes(adata)
```

These outputs exist in the pre-freeze implementation but should be removed from
the active Notebook 00 path during the freeze update:

```text
plots/manual_ec_highly_variable_genes.png
tables/manual_ec_preprocess_parameters.tsv
```

### Pending Notebook 00 Freeze Contract

Status: documentation-only planning note. Do not implement until the user gives
explicit confirmation.

The next Notebook 00 edit should freeze the notebook as a checkpoint notebook:

```text
load selected source matrices
annotate sample metadata
calculate QC metrics
save QC plots and QC/filter tables
apply preliminary Scanpy filters
apply manual_ec biological/QC cutoffs
write combined and per-sample .h5ad checkpoints
stop
```

Notebook 00 should not perform:

```text
HVG selection
cell-cycle scoring
CC.Difference creation
regression
scaling
PCA
neighbors
UMAP
clustering
marker analysis
integration
```

The older Notebook 00 HVG artifacts should be removed from the active path:

```text
plots/manual_ec_highly_variable_genes.png
tables/manual_ec_hvg_genes.tsv
tables/manual_ec_hvg_overlap.tsv
HVG-count fields in Notebook 00 comparison tables
per-sample HVG recomputation cells
```

The run-specific output root remains dynamic and must stay based on
`RUN_LABEL`:

```python
RUN_DIR = DATA_ROOT / "results" / "notebook00" / RUN_LABEL
TABLE_DIR = RUN_DIR / "tables"
PLOT_DIR = RUN_DIR / "plots"
H5AD_DIR = RUN_DIR / "h5ad"
PER_SAMPLE_H5AD_DIR = H5AD_DIR / "per_sample"
```

Required combined checkpoints:

```text
h5ad/manual_ec_filtered_counts.h5ad
h5ad/manual_ec_filtered_normalized_log1p.h5ad
```

`manual_ec_filtered_counts.h5ad` should contain the combined QC/manual_ec
filtered AnnData object. In that object:

```text
.X = filtered raw counts
layers = no required new expression layer
```

`manual_ec_filtered_normalized_log1p.h5ad` should be derived from the filtered
counts object. In that object:

```text
.X = normalized/log1p expression
.layers["counts"] = filtered raw counts copied before normalization/log1p
```

Required per-sample checkpoints:

```text
h5ad/per_sample/<run_sample_id>.manual_ec_filtered_counts.h5ad
h5ad/per_sample/<run_sample_id>.manual_ec_filtered_normalized_log1p.h5ad
```

Each per-sample counts file should be a slice of the combined filtered counts
checkpoint for exactly one `run_sample_id`. Each per-sample normalized/log1p
file should be the matching per-sample slice of the combined normalized/log1p
checkpoint. The `.X` and `.layers["counts"]` meaning must match the combined
checkpoint pair.

Notebook 01 should start from the Notebook 00 checkpoints and perform the
post-filtering analysis decisions:

```text
Seurat-v3 HVG selection from counts
cell-cycle scoring
CC.Difference creation
regression
scaling
PCA
neighbors
UMAP
clustering
```

## Recommended Next Steps

1. After user confirmation, implement the Notebook 00 freeze contract above.

   The implementation should remove HVG selection from Notebook 00, write the
   combined and per-sample `.h5ad` checkpoints under
   `results/notebook00/<RUN_LABEL>/h5ad/`, and stop after QC/filtering/checkpoint
   creation.

2. Run the updated Notebook 00 through Slurm and verify the checkpoint outputs.

   Confirm that the manual cutoffs, QC plots, tables, and `.h5ad` files match
   the intended frozen checkpoint:

   ```text
   highest_expr_genes_top20.png
   manual_ec_qc_violin.png
   manual_ec_scatter_total_counts_pct_counts_mt.png
   manual_ec_scatter_total_counts_n_genes_by_counts.png
   manual_ec_per_sample_highest_expr_genes_top20.png
   manual_ec_per_sample_qc_violin_n_genes_by_counts.png
   manual_ec_per_sample_qc_violin_total_counts.png
   manual_ec_per_sample_qc_violin_pct_counts_mt.png
   manual_ec_per_sample_scatter_total_counts_pct_counts_mt.png
   manual_ec_per_sample_scatter_total_counts_n_genes_by_counts.png
   manual_ec_filter_summary.tsv
   manual_ec_filter_parameters.tsv
   manual_ec_filtered_counts.h5ad
   manual_ec_filtered_normalized_log1p.h5ad
   per-sample manual_ec_filtered_counts.h5ad files
   per-sample manual_ec_filtered_normalized_log1p.h5ad files
   ```

3. After Notebook 00 is verified, create a new Notebook 01 planning markdown.

   That planning document should be created after the Notebook 00 branch is
   ready, not mixed into the Notebook 00 freeze implementation.

4. Merge the Notebook 00 branch into `main` only after verification.

   Then create a fresh branch from updated `main` for Notebook 01 work.

5. Run any next notebook update through Slurm, not on the login node.

   Use the same `13_execute_notebook00_source.sbatch.template` path. Do not
   submit `scripts/cellbender.sh`; the CellBender denoised H5 files remain fixed
   inputs.

6. After user confirmation, archive or delete remaining historical workflow
   pieces outside the active Notebook 00 path.

   The forward Notebook 00 path no longer needs historical MAD/custom filtering,
   basic UMAP/clustering, Harmony, marker panels, or raw full-analysis code. Keep
   source loading and missing-aware reporting because they are still needed for
   Cell Ranger filtered vs existing CellBender-H5 comparison.

7. Keep CellBender denoising frozen unless the user explicitly asks otherwise.

   Future runs should reuse the existing `clean_adata/*_cellbender_denoised.h5`
   files. Notebook 00 can compare those fixed files against Cell Ranger filtered
   matrices, but `scripts/cellbender.sh` should not be submitted again for this
   checkpoint.

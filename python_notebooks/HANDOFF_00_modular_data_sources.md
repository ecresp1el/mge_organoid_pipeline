# Handoff: Modularize Notebook 00 Data Sources and Analysis Flow

Branch: `codex/modular-00-data-sources`

Primary files discussed so far:

- `python_notebooks/notebooks/00_load_div30_div90_raw_to_anndata.ipynb`
- `scripts/cellbender.sh`
- `python_notebooks/src/mge_organoid_python/data_sources.py`
- `python_notebooks/src/mge_organoid_python/notebook00_workflow.py`
- `python_notebooks/src/mge_organoid_python/notebook00_plots.py`
- `slurm_templates/13_execute_notebook00_source.sbatch.template`

## Current State For Next Chat

The branch now has a modular Notebook 00 driver and supporting Python modules.
The notebook has been executed on Great Lakes through Slurm for three source
modes:

```text
cellranger_filtered: full matrix load, historical QC/filtering/preprocess/UMAP/plots
cellranger_raw:      source availability only, no raw-droplet-scale analysis
cellbender_denoised: source availability only, missing-aware, no full analysis yet
```

Direction change for next work:

```text
Do not continue building on the previous custom/MAD filtering approach.
Treat it as historical scaffolding from the modularization pass.
The forward analysis target is a standalone manual_ec Scanpy path that starts
from Cell Ranger filtered matrices and reaches QC plots, manual cutoffs,
normalization, log1p, and HVG selection with the requested parameters.
Also retain the ability to run the same manual_ec checkpoint on CellBender
denoised files so CellBender can be compared directly against Cell Ranger
filtered matrices.
Everything outside this scope should be stripped from the forward notebook path
or flagged for deletion.
```

The proof scripts are useful validation scaffolding. The core notebook runtime is:

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
historical filtered full analysis through the modular notebook
raw source availability through the same notebook driver
CellBender source availability through the same notebook driver
missing-source logging/plotting using 9853-MW-6 as the test case
run-specific table and plot directories under results/notebook00/<RUN_LABEL>
```

Not completed yet:

```text
full CellBender matrix loading for all samples
CellBender-backed manual_ec QC/preprocess checkpoint
CellBender-denoised vs Cell Ranger filtered comparison
raw export manifest cleanup
manual_ec filtered-matrix Scanpy checkpoint
strip/delete historical MAD/custom filtering, UMAP/clustering, Harmony, and marker-panel paths from Notebook 00
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
adata_names
adata_list
analysis_names
analysis_list
combined_adata
```

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

## Pending Direction: `manual_ec` Scanpy Path

Status: handoff direction only. Do not implement code until the user approves the
actual change set.

Forget the earlier attempt to make our own Notebook 00 filtering approach. The
historical MAD/custom path was useful for modularization and smoke testing, but
it is not the forward scientific workflow. Future work should focus on an
explicit `manual_ec` Scanpy path for controlled review of Heyoon's collaborator
workflow choices.

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

Proposed run identities:

```text
NOTEBOOK00_ACTIVE_SOURCE=cellranger_filtered
NOTEBOOK00_RUN_LABEL=cellranger_filtered_manual_ec_div30_core_samples

NOTEBOOK00_ACTIVE_SOURCE=cellbender_denoised
NOTEBOOK00_RUN_LABEL=cellbender_denoised_manual_ec_div30_core_samples
```

Do not add more toggles whose purpose is to preserve the older filtering method.
If old code remains temporarily, treat it as historical/cleanup debt rather than
an alternate supported path.

The next implementation should remove dependencies from the new path on:

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

Keep a direct comparison between:

```text
cellranger_filtered manual_ec output
cellbender_denoised manual_ec output
```

This comparison should use common `run_sample_id` values only. If `9853-MW-6` or
any other CellBender output is missing, the comparison should proceed with the
available common samples and write the missing sample into `source_table`.

Recommended comparison outputs:

```text
tables/manual_ec_source_comparison_summary.tsv
tables/manual_ec_qc_metric_comparison_by_sample.tsv
tables/manual_ec_retained_cell_comparison_by_sample.tsv
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

### Proposed `manual_ec` QC/Plot Checkpoint

After loading and concatenating either selected source, the notebook should be
able to get to this Scanpy checkpoint before any UMAP/cluster work:

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

The plots should be saved under the run-specific Notebook 00 plot directory,
using stable names such as:

```text
highest_expr_genes_top20.png
manual_ec_qc_violin.png
manual_ec_scatter_total_counts_pct_counts_mt.png
manual_ec_scatter_total_counts_n_genes_by_counts.png
```

### Proposed `manual_ec` Filtering

The requested manual cutoff should be implemented as a named, documented filter:

```python
adata = adata[
    (adata.obs["n_genes_by_counts"] > 1500)
    & (adata.obs["n_genes_by_counts"] < 9000)
    & (adata.obs["total_counts"] < 30000)
    & (adata.obs["pct_counts_mt"] < 5),
    :
].copy()
```

Record this output as `manual_ec`, including a table that reports starting cells,
retained cells, removed cells, retained percent, and the exact thresholds:

```text
tables/manual_ec_filter_summary.tsv
tables/manual_ec_filter_parameters.tsv
```

Also include these Scanpy tutorial-style count filters, either before the manual
cutoffs or as an explicitly documented preliminary step:

```python
sc.pp.filter_cells(adata, min_genes=20)
sc.pp.filter_genes(adata, min_cells=3)
```

The implementation should preserve enough reporting to distinguish:

```text
cellranger_filtered source loading
cellbender_denoised source loading
preliminary min_genes/min_cells filtering
manual_ec biological/QC cutoffs
historical MAD/custom filtering path, if still present in old outputs
```

### Proposed `manual_ec` Preprocessing

After `manual_ec` filtering, the requested preprocessing checkpoint is:

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

Save the HVG plot and parameters:

```text
plots/manual_ec_highly_variable_genes.png
tables/manual_ec_preprocess_parameters.tsv
```

### Seurat-Equivalent Cell-Cycle Requirement

This is also part of the controlled conversion request, but can be implemented
after the `manual_ec` QC/preprocess checkpoint above is working.

The conversion should preserve Heyoon's biological parameters and use Rich et
al. only to choose the closest Seurat-like Scanpy implementation. For cell-cycle
scoring:

- use Seurat-equivalent S and G2M gene lists
- do not rely on Scanpy implicitly exposing Seurat's `cc.genes`
- explicitly store the exact S and G2M gene lists used for scoring in the output

Proposed output locations:

```text
tables/manual_ec_cell_cycle_s_genes.tsv
tables/manual_ec_cell_cycle_g2m_genes.tsv
tables/manual_ec_cell_cycle_score_parameters.tsv
```

The corresponding AnnData object should also store these lists in `.uns`, for
example:

```python
adata.uns["manual_ec_cell_cycle_s_genes"] = s_genes
adata.uns["manual_ec_cell_cycle_g2m_genes"] = g2m_genes
adata.uns["manual_ec_cell_cycle_note"] = "Seurat-equivalent gene lists used for Scanpy score_genes_cell_cycle."
```

## Recommended Next Steps

1. Get user approval for the `manual_ec` direction above.

   The next implementation should not preserve or extend the historical
   custom/MAD filtering path. It should make `manual_ec` the focused path from
   Cell Ranger filtered matrices to QC plots, manual cutoffs, normalization,
   log1p, and HVG selection, while retaining a direct comparison against
   CellBender-denoised inputs.

2. Strip or bypass historical analysis dependencies in the forward notebook path.

   Flag old MAD/custom filtering, basic PCA/UMAP, clustering, Harmony, marker
   panels, raw full-analysis, and generic comparison scaffolding for deletion.
   Do not remove source loading or missing-source reporting because those are
   still needed for the filtered-vs-CellBender comparison.

3. Implement and validate the `manual_ec` Cell Ranger filtered checkpoint.

   The first successful run should produce, at minimum:

   ```text
   plots/highest_expr_genes_top20.png
   plots/manual_ec_qc_violin.png
   plots/manual_ec_scatter_total_counts_pct_counts_mt.png
   plots/manual_ec_scatter_total_counts_n_genes_by_counts.png
   plots/manual_ec_highly_variable_genes.png
   tables/manual_ec_filter_summary.tsv
   tables/manual_ec_filter_parameters.tsv
   tables/manual_ec_preprocess_parameters.tsv
   ```

4. Apply the same `manual_ec` checkpoint to CellBender-denoised inputs.

   Start with available common samples. If one CellBender output is still
   missing, record it in `source_table` and compare the samples that exist in
   both sources.

5. Add the retained comparison outputs.

   The comparison should summarize CellBender-denoised vs Cell Ranger filtered
   by sample/cell line, using cell counts, QC metrics, retained cells, and HVG
   overlap.

6. Keep the data-source layer and missing-aware reporting behavior.

   Missing files, failed jobs, or intentionally omitted samples should remain
   visible in `source_table`. This is source bookkeeping, not the old filtering
   approach.

7. Recheck the running `9853-MW-6` CellBender task before an all-sample
   CellBender comparison run.

   A partial CellBender comparison can use available common samples first.

   Last checked during this handoff update:

   ```text
   job_id: 51106991_5
   state: RUNNING
   elapsed: 02:05:06
   node: gl1008
   final_h5_present: no
   posterior_h5_present: yes
   ```

   Check again with:

   ```bash
   squeue -j 51106991_5 -o '%.18i %.9T %.20M %.10D %.30R'
   ls -lh /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/clean_adata/9853-MW-6_cellbender_denoised.h5
   sacct -j 51106991_5 --format=JobID,JobName%20,State,ExitCode,Elapsed,NodeList%30 -P
   ```

8. Once `9853-MW-6` finishes, rerun the CellBender location proof.

   ```bash
   sbatch slurm_templates/11_prove_00_cellbender_locations.sbatch.template
   ```

   Success should change from:

   ```text
   cellbender_outputs_found: 5
   cellbender_outputs_missing: 1
   ```

   to:

   ```text
   cellbender_outputs_found: 6
   cellbender_outputs_missing: 0
   ```

9. Rerun the missing-aware proof without expecting `9853-MW-6` to be missing.

   ```bash
   EXPECT_MISSING_SAMPLE= sbatch slurm_templates/12_prove_00_missing_aware_downstream.sbatch.template
   ```

   This proves the same reporting path works when every requested sample is
   available.

10. Prove actual CellBender matrix loading on a compute node before enabling a
    full all-sample CellBender comparison run.

   `data_sources.py` has a CellBender reader path, but the safe proven mode so
   far is source availability only. The next code/proof milestone should be:

   ```python
   config = Notebook00SourceConfig(active_source="cellbender_denoised")
   result = load_dataset_result(config, load_matrices=True)
   ```

   Start with one sample or a small target sample set before loading all six
   raw-droplet-scale outputs.

11. Decide whether the proof scripts should stay as active validation utilities
   or move under an archive/validation area after the notebook is stable.

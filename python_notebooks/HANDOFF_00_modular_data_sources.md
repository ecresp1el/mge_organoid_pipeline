# Handoff: Modularize Notebook 00 Data Sources and Analysis Flow

Branch: `codex/modular-00-data-sources`

Primary files discussed so far:

- `python_notebooks/notebooks/00_load_div30_div90_raw_to_anndata.ipynb`
- `scripts/cellbender.sh`

## Why This Exists

Notebook 00 currently mixes several different responsibilities in one long, stateful workflow:

1. find the project checkout
2. read sample metadata
3. resolve Cell Ranger raw or filtered matrix paths
4. export raw `.h5ad` files for CellBender
5. run QC metric annotation
6. apply QC filtering
7. concatenate samples
8. preprocess, embed, integrate, cluster, and plot

That is workable for an early teaching notebook, but it is hard to scale because changing the input data source also changes the mental model of the rest of the notebook. The next version should make the active data source explicit and make each analysis stage callable on demand.

## Current Data Flow

The notebook and the Slurm script do not currently use one single "project folder" concept.

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

### Current active notebook analysis source

The current builder is instantiated with:

```python
matrix_source="filtered"
```

Therefore this call:

```python
adata_names, adata_list = anndata_builder.per_sample_anndata_list()
```

loads Cell Ranger filtered 10x matrix directories. The QC plots, QC filtering, preprocessing, PCA, UMAP, Harmony, Leiden, and comparison plots currently operate on that filtered Cell Ranger data.

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

### Critical current gap

CellBender output is produced by `scripts/cellbender.sh`, but Notebook 00 does not currently read those CellBender outputs back into the analysis flow.

Current state:

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
- The notebook can load Cell Ranger raw data through the same source interface.
- CellBender output locations are represented in code, even if the files are not present on the local machine.
- Raw export writes or previews a manifest that maps `run_sample_id` to raw `.h5ad` and expected CellBender `.h5` output.
- QC annotation remains separate from QC filtering.
- Running a QC plot does not require rerunning preprocessing, UMAP, Harmony, or every other plot.
- Running an embedding plot does not require rerunning raw export or CellBender-related code.
- The notebook still preserves the key `.obs` metadata fields across all supported sources.
- The user can compare at least two loaded/processed sources by `run_sample_id`, `cell_line`, or QC metrics without rewriting path logic.

## Implementation Boundaries

First implementation pass should avoid changing scientific thresholds unless explicitly requested.

Keep these stable unless there is a clear reason:

- MAD threshold rules
- mitochondrial prefix
- target sample IDs
- target DIV selection
- CellBender command-line behavior
- existing CellBender output naming convention

The first pass should focus on structure, naming, path clarity, and source switching.

## Risks To Watch

- CellBender `.h5` files may need a different reader or postprocessing than Cell Ranger 10x directories.
- Raw `.h5ad` filenames are currently derived from matrix path prefixes, so duplicate-safe naming and `run_sample_id` mapping must remain explicit.
- The notebook currently has both Python `PROJECT_ROOT` and shell `$PROJECT_ROOT`; these must not remain ambiguous.
- Large data should stay under `DATA_ROOT`, not inside the GitHub checkout.
- Any code that rewrites `.obs["run_sample_id"]` from `batch` should be checked carefully so sample identity is preserved.

## Proposed Next Step

Before changing the notebook cells, create or identify a small Python module for reusable workflow code, likely under:

```text
python_notebooks/src/mge_organoid_python/
```

Then keep Notebook 00 as the readable driver:

```python
from mge_organoid_python.data_sources import load_dataset, export_raw_for_cellbender
from mge_organoid_python.qc import PerSampleMADQC, filter_qc_pass
from mge_organoid_python.plots import plot_qc_sample, plot_embedding, compare_embeddings
```

This keeps the notebook editable while moving fragile path and analysis mechanics into testable code.

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

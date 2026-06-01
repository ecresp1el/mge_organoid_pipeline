# Handoff: Notebook 01 From Notebook 00 Checkpoints

Branch: `codex/notebook01-from-nb00-checkpoints`

Start point: `main` after the verified Notebook 00 freeze merge.

Notebook 00 is now frozen as a loading, QC, filtering, and checkpoint-creation
notebook. Notebook 01 should start from the `.h5ad` checkpoints written by
Notebook 00 and should not rerun matrix loading, source comparison setup, or
manual_ec filtering.

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

## Notebook 01 Scope

Notebook 01 begins the post-filtering analysis choices:

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

Notebook 01 should also write its own run-scoped outputs under a new Notebook 01
results root, for example:

```text
results/notebook01/<RUN_LABEL>/
results/notebook01/<RUN_LABEL>/tables/
results/notebook01/<RUN_LABEL>/plots/
results/notebook01/<RUN_LABEL>/h5ad/
```

## Decisions To Confirm Before Implementation

1. Which Notebook 00 checkpoint should be the primary Notebook 01 input:
   `cellranger_filtered` or `cellbender_denoised`.

2. Whether Notebook 01 should operate on the combined checkpoint only first, or
   also process per-sample checkpoints in parallel.

3. Exact Seurat-v3 HVG parameters:
   expected number of HVGs, batch key if any, and whether HVG should be selected
   globally or per sample.

4. Exact cell-cycle gene lists:
   use explicit stored S and G2M gene lists rather than relying on Scanpy to
   expose Seurat lists implicitly.

5. Regression variables:
   at minimum confirm whether to regress `CC.Difference`, total counts, percent
   mitochondrial counts, or any batch/source covariates.

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

## First Implementation Step

Create or update the Notebook 01 driver so it:

1. Resolves `DATA_ROOT`.
2. Accepts a Notebook 00 `RUN_LABEL` as input.
3. Loads `results/notebook00/<RUN_LABEL>/h5ad/manual_ec_filtered_normalized_log1p.h5ad`.
4. Confirms `.layers["counts"]` exists and matches expected dimensions.
5. Writes a Notebook 01 input validation table before running downstream
   analysis.

# PCDH19 scRNA-seq primary-processing handoff

## Purpose

Build a reusable, object-oriented Python pipeline for primary scRNA-seq
processing of the PCDH19 dataset. This workflow proceeds one major step at a
time and is independent of the earlier reference-mapping, probe, and genotype-
classification workflows.

## Approval rule

For every major step:

1. Implement only the current approved step.
2. Run and validate that step.
3. Save its code, configuration, outputs, logs, and checkpoint.
4. Report inputs, parameters, changes, dimensions, warnings, and review plots
   or tables.
5. Stop and wait for explicit approval.

A successful computation is `IN_REVIEW`, not `APPROVED`. Only the exact run
explicitly approved by the user may become the input to the next major step.
Requested diagnostics, parameter changes, and reruns remain part of the same
step until approval.

## Project and execution conventions

- Canonical code: `paper3_pcdh19/scripts/primary_processing/`.
- Thin shell entry points and Great Lakes SLURM wrappers contain no scientific
  analysis logic.
- Existing Great Lakes environment, logging, versioned-run, and atomic-
  publication conventions are reused.
- Exact Python, shell, SLURM, configuration, handoff, and metadata files are
  copied into each run before submission.
- SLURM executes those frozen copies, not mutable repository files.
- Submitted scheduler scripts: `PAPER3_ROOT/jobs/`.
- Scheduler logs: `PAPER3_ROOT/logs/primary_processing/`.
- Scientific results:
  `PAPER3_ROOT/results/primary_processing/<STEP>/<RUN_ID>/`.
- A new versioned run is the default. Guarded `--replace-run RUN_ID` is limited
  to one named inactive run within the same step.
- Raw counts and processing provenance are preserved.

## Input distinction

For every sample, `sample_filtered_feature_bc_matrix.h5` is the Cell Ranger
called-cell matrix and the primary expression input. Only `Gene Expression`
features enter the canonical AnnData.

Each `sample_raw_feature_bc_matrix.h5` contains a larger barcode universe. Its
path and structural provenance are registered for possible Step 04 ambient-RNA
assessment, but its uncalled droplets do not enter the Step 00 AnnData.
Probe-level matrices are not transcriptome inputs.

The authoritative biological metadata are in `config/sample_key.csv`.
`config/sample_manifest_draft.tsv` supplies registered technical identifiers
and vendor metrics. Fields absent from these registered tables remain absent
or explicitly unavailable; sample names are not used to infer embryo, litter,
or batch structure.

## Major steps

### Step 00 — input validation and canonical AnnData

Create one canonical AnnData from the 12 filtered matrices and registered
metadata. Preserve sparse, unnormalized integer counts; use globally unique
cell IDs; validate features, barcodes, counts, metadata joins, dimensions, and
uniqueness; write an H5AD; and reopen it for validation.

Do not filter cells or genes, calculate biological QC metrics, normalize,
scale, reduce dimensions, construct neighbors, cluster, or annotate. Stop for
review after the computation.

### Step 01 — QC metrics

Calculate and visualize cell-level QC metrics without removing cells. Do not
select final thresholds. Stop for review.

### Step 02 — QC filtering

Use approved Step 01 evidence to propose and then apply reviewed filters.
Preserve pass/fail state and exclusion reasons. Stop for review.

### Step 03 — doublet assessment

Generate reproducible scores and candidate calls. Do not automatically remove
tool-called cells. Stop for review.

### Step 04 — ambient RNA and contamination assessment

First determine whether meaningful ambient contamination exists. Do not apply
CellBender or another correction automatically. Review evidence, method, and
consequences before correction. Stop for review.

### Step 05 — broad biological contaminant assessment

Assess clearly unwanted broad populations without treating this as fine cell-
type annotation. Any biological removal requires review. Stop for review.

### Step 06 — technical, sample, and batch diagnostics

Use diagnostic normalization and reductions when needed to determine whether
technical structure matters. Do not apply Harmony, scVI, Seurat integration,
or another correction without first demonstrating need and obtaining review.
Stop for review.

### Step 07 — freeze validated HiCAT input

Publish the exact representation required by HiCAT with raw-count provenance,
cell/gene metadata, QC and doublet evidence, exclusion history, configuration,
and final dimensions. Confirm HiCAT requirements before finalization. Do not
begin HiCAT without explicit instruction.

## Approval ledger

Maintain a human-readable ledger with step, run ID, status, input/output
checkpoint, timestamps, code/configuration identity, dimensions, parameters,
findings, requested additions, reruns, outstanding questions, approval
decision, and approved run ID.

The only valid states for a computed checkpoint are `IN_REVIEW` and
`APPROVED`; submission and failure states may also be recorded operationally.
Only `APPROVED` can be consumed by the next major step.

## Current authorization

Step 00 has been computed and is **IN REVIEW**. Do not begin Step 01 or any
later processing until the user explicitly approves this exact run.

### Step 00 run awaiting review

- Run ID:
  `00_input_validation_and_canonical_anndata_20260830_113749_d8b6bf7`
- Successful Great Lakes job: `59279775` (`COMPLETED`, exit `0:0`, 2 minutes
  12 seconds, approximately 28.97 GB maximum resident memory).
- Frozen executable commit: `7fe06d5e2d7034deda5fd3fcc9dd63c544cca750`.
- Canonical checkpoint:
  `objects/pcdh19_step00_canonical_raw_counts.h5ad`.
- Checkpoint size: 6,601,426,576 bytes.
- Dimensions: 450,788 cells × 19,071 genes.
- Sparse nonzero entries: 1,295,361,777.
- Validation: 221 PASS, 0 FAIL.
- Documentation audit: 66 module/class/function definitions PASS, 0 FAIL.
- H5AD SHA-256:
  `c420d501cc0acca6567014d8f6c7962a6b7ba1c80036c5a7c39fe593033848a6`.

The raw matrices contain 19,404 Gene Expression features and therefore 333
genes not present in the filtered matrices. All 19,071 canonical filtered
genes occur in the raw matrices with their feature definitions and order
preserved. The raw matrices remain path-registered only; their extra genes and
uncalled droplets do not enter Step 00.

The first scheduler attempt, job `59279713`, failed before scientific code ran
because a space-containing Cell Ranger path was not shell-escaped in the
generated environment file. Commit `7fe06d5` fixed that execution defect. The
authorized guarded replacement retained the failed job ID and logs, froze the
corrected code, and produced successful job `59279775` in the same run.

### Step 00 implementation map

- `scripts/primary_processing/models.py`: immutable settings, path resolution,
  source-structure records, and validation-ledger records.
- `scripts/primary_processing/loaders.py`: biological/technical registry join,
  non-loading 10x HDF5 inspection, exact barcode/metrics checks, and sparse
  raw-count AnnData loading.
- `scripts/primary_processing/validation.py`: structural, matrix-state,
  forbidden-analysis, and Python-docstring audits.
- `scripts/primary_processing/publishing.py`: atomic staging/publication,
  software and output manifests, review report, and approval-ledger updates.
- `scripts/primary_processing/workflow.py`: Step 00 orchestration only.
- `scripts/primary_processing/cli.py`: explicit frozen command-line interface.
- `bin/submit_primary_processing_step_00.sh`: validation, run freezing,
  versioning/guarded replacement, provenance capture, and SLURM submission.
- `slurm/primary_processing_00_input_validation_and_canonical_anndata.sbatch`:
  thin Great Lakes executor with runtime/failure markers; no scientific logic.

Every frozen Python module, class, and function has a docstring. The run-level
`tables/documentation_audit.tsv` enforces that contract mechanically.

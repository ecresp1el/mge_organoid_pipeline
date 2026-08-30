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

Implement and run Step 00 only. Present its validation results and stop. Do
not begin Step 01 or any later processing until the user explicitly approves
the Step 00 run.

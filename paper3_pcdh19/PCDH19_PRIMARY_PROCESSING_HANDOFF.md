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

Step 00 was explicitly **APPROVED** by the user on 2026-08-30. Step 01 QC
metrics have been computed and remain **IN_REVIEW**. Before approval, the user
authorized a separate Step 01a per-sample MAD sensitivity amendment. Do not
begin Step 02 or any later processing until the user explicitly approves the
combined evidence from the exact Step 01 and Step 01a runs.

### Approved Step 00 run

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

### Step 01 computed run — IN_REVIEW

- Run ID: `01_qc_metrics_20260830_115715_2b57907`.
- Successful Great Lakes job: `59281063` (`COMPLETED`, exit `0:0`, 5 minutes
  21 seconds, 27,382,536 KB maximum resident memory).
- Frozen executable commit:
  `2b57907b4a31e0f99e7d235715ccf6084750de11`.
- Checkpoint: `objects/pcdh19_step01_qc_metrics.h5ad`.
- Checkpoint size: 6,617,737,368 bytes.
- Dimensions: 450,788 cells x 19,071 genes; 1,295,361,777 sparse nonzeros.
- Validation: 40 PASS, 0 FAIL, including exact logical sparse-matrix
  fingerprint preservation and a successful H5AD round trip.
- Output H5AD SHA-256:
  `5bf137f0b19b9a734fe4085d47fb28b073deca36dbd00c9b2c2cba6e1ebeaa0e`.
- Figures: 15 PNG plus matching 15 PDF files: one pooled diagnostic, one for
  each of the 12 samples, and sample/design-group comparisons.
- Filtering, thresholds, exclusions, normalization, reductions, graphs,
  clustering, and annotation: none.

Pooled medians are 4,014 total counts, 2,601 detected genes, and 0.648%
mitochondrial counts. The pooled 5th-95th percentile intervals are
1,647-12,560 counts, 1,281-5,358 genes, and 0.145%-1.819% mitochondrial
counts. Across samples, medians range from 2,978 to 4,591 counts, 2,042 to
2,867 genes, and 0.285% to 0.980% mitochondrial counts. These are descriptive
review observations only and do not define high/low quality or a threshold.

The run-local report, all plots, complete summary tables, software versions,
validation ledger, frozen executed code/configuration, and output manifest
are inside the run directory. Request additional Step 01 diagnostics or
explicitly approve this exact run before any Step 02 work.

### Step 01a authorized scope

Step `01a_qc_mad_sensitivity` consumes the exact in-review Step 01 H5AD in
backed read-only mode. It calculates candidate boundaries independently for
each of the 12 technical samples; design group, genotype, and sex cannot
contribute to boundary definition. The evaluated stringencies are 3, 4, and
5 scaled MAD, where scaled MAD is `1.4826 * median(abs(x - median(x)))`.

Low total counts and low detected genes are evaluated on `log1p` values and
back-transformed to the original metric scale. High mitochondrial percentage
is evaluated on the observed percentage scale. Outputs include every
boundary, per-cell candidate flags, individual and joint counts/percentages,
exact overlap patterns, per-sample distribution plots with all boundaries,
and cross-sample boundary/flag/overlap visualizations in PNG and PDF.

Step 01a must not calculate an upper total-count or upper detected-gene
candidate rule. High-complexity cells are reserved for the later Scrublet
step. It must not remove cells or genes, write a replacement H5AD, define a QC
pass/fail decision, start Step 02, or treat a candidate flag as an exclusion.
Its successful result remains `IN_REVIEW` with Step 01.

### Step 01 implemented scope

Step 01 must load only the approved Step 00 H5AD and use
`scanpy.pp.calculate_qc_metrics()` to add descriptive QC metadata. It must
calculate total counts, detected genes, and mitochondrial counts/fraction;
preserve every cell, gene, and raw count; create pooled, each-sample, and
design-group summaries/plots; and stop `IN_REVIEW` without thresholds or
filtering.

Mitochondrial genes are the 13 delivered symbols beginning `mt-`. The targeted
Flex panel has no `Rpl`/`Rps` genes, so Step 01 records ribosomal fraction as
unavailable instead of fabricating an all-zero metric. `percent_top=None`
avoids an expensive, unrequested top-gene calculation and does not affect the
required totals, detected-gene, or mitochondrial metrics.

Step 01 implementation files are `step01_models.py`, `step01_metrics.py`,
`step01_plots.py`, `step01_validation.py`, `step01_publishing.py`,
`step01_workflow.py`, and `step01_cli.py`. The thin submitter is
`bin/submit_primary_processing_step_01.sh`; the thin Great Lakes wrapper is
`slurm/primary_processing_01_qc_metrics.sbatch`. Every frozen Python
definition is docstring-audited in the run.

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

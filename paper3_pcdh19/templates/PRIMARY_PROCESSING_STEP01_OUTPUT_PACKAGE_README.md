# PCDH19 primary processing: Step 01 QC-metrics run

This versioned package calculates Scanpy QC metrics from one explicitly
approved Step 00 checkpoint. It does not filter cells or genes and does not
establish thresholds.

- `code/`: exact frozen object-oriented Python package, submitter, and SLURM
  wrapper executed by the job;
- `config/`: environment/settings, approved Step 00 status/manifest, approval-
  ledger snapshot, requirements, and authoritative handoff;
- `objects/`: Step 01 H5AD with unchanged raw counts and added QC metadata;
- `tables/`: pooled/sample/design summaries, feature-set audit, metric
  dictionary, validation checks, documentation audit, plot manifest, software
  versions, object state, and output manifest;
- `figures/`: pooled, each-sample, sample-comparison, and design-group QC
  diagnostics;
- `logs/`: links to scheduler logs under
  `PAPER3_ROOT/logs/primary_processing/`;
- `provenance/`: repository, submission, runtime, and scheduler identities.

`COMPUTATION_SUCCESS.txt` means only that computation and validation passed.
The result remains `IN_REVIEW` until the user explicitly approves this exact
run. It cannot feed Step 02 before that approval.

# PCDH19 primary processing: Step 00 run package

This versioned package contains the canonical raw-count AnnData construction
for one `00_input_validation_and_canonical_anndata` run.

Directory roles:

- `code/`: exact Python package, submitter, and SLURM wrapper frozen before
  submission and executed by the job;
- `config/`: exact path/resource settings, package pins, biological sample key,
  technical manifest, and authoritative handoff used by this run;
- `objects/`: canonical H5AD, published only after all validations pass;
- `tables/`: input summaries, metadata dictionary, documentation audit,
  validation checks, environment versions, object summary, and output manifest;
- `logs/`: links to scheduler logs stored under
  `PAPER3_ROOT/logs/primary_processing/`;
- `provenance/`: submission command, repository state, job identity, and runtime
  environment.

`COMPUTATION_SUCCESS.txt` means that the frozen code ran and passed its
structural checks. It does not mean the step is scientifically approved.
`STEP_STATUS.tsv` and the workflow-level `APPROVAL_LEDGER.tsv` remain
`IN_REVIEW` until the user explicitly approves this exact run.

Step 00 performs no QC filtering, normalization, scaling, reduction, graph
construction, clustering, or annotation. Raw Cell Ranger counts are stored in
the canonical object's `.X`; they are not redundantly copied to another layer.

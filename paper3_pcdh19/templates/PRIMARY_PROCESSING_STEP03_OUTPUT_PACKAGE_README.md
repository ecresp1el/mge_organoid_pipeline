# Primary processing Step 03 run package

This directory is a frozen, versioned execution package for approved
scDblFinder detection. The job executes the Python, R, shell, SLURM,
configuration, and metadata copies under this run directory—not mutable
repository files.

The project-scoped R library is reproducible with the frozen copy of
`install_primary_processing_step03_r_environment.sh`; the run also freezes the
expected key-package versions and the complete installed-package inventory.

## Scope

- Input: exact approved Step 02 H5AD with sparse raw integer counts.
- Capture: one independently processed `GEX_1` capture; 12 Probe Barcode
  technical samples remain reporting metadata.
- Detector: native R scDblFinder with `clusters=TRUE`, `dbr.sd=1`, no supplied
  `dbr`, and otherwise model defaults.
- Reproducibility: one second-seed run with identical model arguments.
- Filtering: none. Calls and scores are retained for review.

The primary output-only `returnType="full"` retains scDblFinder's exact
internal PCA for the requested expression-space diagnostic. Artificial
doublets are excluded from the saved PCA and checkpoint. This return choice
does not change classifier training or thresholding.

## Outputs

- `objects/pcdh19_step03_scdblfinder.h5ad`: unchanged raw counts plus primary
  and repeat scores/calls and the diagnostic internal PCA.
- `tables/scdblfinder_per_cell_details.tsv.gz`: complete native-R real-cell
  result fields.
- `tables/`: call fractions, sample/design/cluster composition,
  reproducibility, versions, validation, method contract, and dictionaries.
- `figures/`: score, called-fraction, reproducibility, and internal-PCA plots.
- `STEP03_SCDBLFINDER_REPORT.md`: the review report.
- `STEP_STATUS.tsv`: remains `IN_REVIEW` until explicit user approval.

`config/`, `code/`, `intermediate/`, and `provenance/` preserve the executable
run and cross-language bridge. A deliberate rerun normally creates a new run;
guarded replacement is limited to one named inactive Step 03 run.

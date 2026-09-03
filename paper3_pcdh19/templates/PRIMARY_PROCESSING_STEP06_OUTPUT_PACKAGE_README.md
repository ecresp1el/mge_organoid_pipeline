# Primary processing Step 06 run package

This immutable, versioned package executes the unintegrated technical,
biological-sample, and batch diagnostics specified in
`STEP06_DIAGNOSTIC_SPECIFICATION.txt`.

## Lineage

- Input: exact approved Step 02 raw-count H5AD.
- Step 03: rejected; no scDblFinder field is consumed.
- Steps 04 and 05: explicitly skipped; no checkpoint exists.
- Step 06: diagnostic only and unintegrated.

The scientific runtime is Python-only and operates directly on the approved
AnnData/H5AD with `anndata`, `scanpy`, NumPy, pandas, SciPy, scikit-learn, and
Matplotlib. It invokes no R runtime or R-based single-cell/plotting package.

## Outputs

- `objects/pcdh19_step06_unintegrated_diagnostics.h5ad`: unchanged raw
  counts plus validated descriptive annotations, PCA, one global UMAP, and
  the unintegrated neighbor graph.
- `figures/step06_primary_diagnostic_report.pdf`: coherent A-L review report.
- `tables/`: full numerical companion tables, parameters, versions,
  validation ledger, output manifest, and rendering-cell IDs.
- `STEP06_DIAGNOSTIC_REPORT.md`: concise methods and provisional outcome.
- `STEP_STATUS.tsv`: remains `IN_REVIEW`.

All numerical analyses use every cell. The deterministic balanced subset in
`tables/rendering_cell_ids.tsv.gz` affects rendering only. This run performs
no integration, batch correction, final annotation, or cell/gene removal.

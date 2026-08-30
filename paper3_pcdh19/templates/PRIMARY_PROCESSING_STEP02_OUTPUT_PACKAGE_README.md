# Primary processing Step 02 run package

This directory is a frozen, versioned execution package for approved QC
filtering. The job executes the copies under `code/`, not mutable repository
files.

## Scope

- Input: exact approved Step 01 raw-count H5AD.
- Decision input: exact approved Step 01a per-cell candidate flags.
- Filter: per-technical-sample 5-MAD low counts OR low genes OR high mt%.
- Genes removed: none.
- Normalization, reduction, clustering, integration, scDblFinder: none.

## Outputs

- `objects/pcdh19_step02_qc_filtered.h5ad`: retained raw integer counts.
- `tables/step02_per_cell_disposition.tsv.gz`: every original cell and every
  applicable exclusion reason.
- `tables/step02_before_after_by_sample.tsv`: exact sample counts.
- `tables/step02_before_after_by_design_group.tsv`: reporting-only design counts.
- `STEP02_QC_FILTERING_REPORT.md`: filtering review report.
- `STEP_STATUS.tsv`: remains `IN_REVIEW` until explicit user approval.

`config/`, `code/`, and `provenance/` preserve the complete executable run
package. A deliberate rerun normally creates a new run ID; guarded replacement
is limited to one named inactive Step 02 run.

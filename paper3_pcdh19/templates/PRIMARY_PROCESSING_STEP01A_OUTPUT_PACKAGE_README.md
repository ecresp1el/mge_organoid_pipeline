# Primary processing Step 01a output package

This versioned package is a non-filtering technical-outlier sensitivity
amendment to the exact Step 01 checkpoint. It calculates 3, 4, and 5 scaled-
MAD candidate boundaries independently within each technical sample.

The source H5AD is opened read-only. No cell or gene is removed, no output
H5AD is created, and no design/genotype group contributes to a boundary.
Upper total-count and detected-gene thresholds are out of scope because
high-complexity cells will be assessed separately with Scrublet.

`code/` and `config/` contain exact submitted copies. `tables/` contains
boundaries, per-cell flags, individual/joint summaries, validation, software,
and the output manifest. `figures/` contains every boundary and candidate-
percentage visualization in PNG and PDF. `STEP_STATUS.tsv` remains
`IN_REVIEW` until explicit user approval.


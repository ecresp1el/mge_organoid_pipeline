# GSE94641 source cache

This directory is the immutable local cache for the processed GSE94641 mouse
embryonic MGE reference. `source_files/` contains the checksum-locked GEO
archives and metadata; `processed_expression/` contains the 225 safely
extracted per-cell processed-expression files. The large downloaded assets are
intentionally ignored by Git and are fully identified by
`config/gse94641_reference_validation.lock.json` and the validation package's
`source_manifest.tsv`.

Run `../../bin/run_step_00_gse94641_reference_validation.sh` from any working
directory. This component validates the reference and stops before query
loading or label transfer.

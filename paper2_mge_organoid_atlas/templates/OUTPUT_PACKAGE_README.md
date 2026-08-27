# Paper 2 output package

This is an immutable, run-scoped output package for the Paper 2 cross-study
MGE organoid atlas.

The exact submitted code is under `code/`, submitted configuration under
`config/`, tabular outputs under `tables/`, logs under `logs/`, and run/Git/
environment records under `provenance/`.

For `00_input_audit`, all registered source objects are read-only. The audit
inventories their file identities, object structures, expression
representations, features, metadata, and reductions. It does not repeat QC,
alter cells, choose an integration method, or create the frozen master input.

Run completion requires `SUCCESS.txt`. If `FAILED.txt` exists, review it and
the SLURM logs before creating a new run; never overwrite this package.


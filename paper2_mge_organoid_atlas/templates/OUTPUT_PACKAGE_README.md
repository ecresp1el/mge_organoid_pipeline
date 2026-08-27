# Paper 2 output package

This is a run-scoped output package for the Paper 2 cross-study MGE organoid
atlas. `config/resolved.env` records whether it was created in immutable
versioned mode or explicit working replacement mode.

The exact submitted code is under `code/`, submitted configuration under
`config/`, tabular outputs under `tables/`, logs under `logs/`, and run/Git/
environment records under `provenance/`.

For `00_input_audit`, all registered source objects are read-only. The audit
inventories their file identities, object structures, expression
representations, features, metadata, and reductions. It does not repeat QC,
alter cells, choose an integration method, or create the frozen master input.

The audit also regenerates cluster UMAP plots from the saved embeddings. These
are all-cell input inventories, not newly computed embeddings. DIV90 uses the
established plotting-only vertical orientation but retains current clusters
6/7. Tables preserve original coordinates, raw cluster IDs, displayed cluster
names, mapping sources, and raw versus displayed cluster counts.

PNG and PDF UMAP outputs are mandatory and fixed at 300 dpi. SVG is optional;
when enabled its rasterized point layers are also fixed at 300 dpi and text
remains editable.

Run completion requires `SUCCESS.txt`. If `FAILED.txt` exists, review it and
the SLURM logs before retrying. Frozen/versioned packages must not be
overwritten; a package may be replaced only through the explicit guarded
`--replace-run` working mode.

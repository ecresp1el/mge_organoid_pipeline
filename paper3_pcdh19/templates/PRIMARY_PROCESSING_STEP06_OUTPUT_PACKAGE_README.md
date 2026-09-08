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
- `provenance/step06_progress_events.jsonl`: append-only live events with the
  exact inputs, parameters, outputs, elapsed time, and peak memory for every
  major function.
- `provenance/step06_progress_latest.json`: atomically refreshed current event.

## Reusable AnnData contents

The H5AD contains raw sparse integer counts in `.X`, original sample/QC
metadata, per-cell mean-log marker-program scores, descriptive Leiden labels,
and provisional cluster-state labels in `.obs`. `.var` retains all genes and
adds `highly_variable_step06`. `.obsm['X_pca']` and `.obsm['X_umap']` contain
the all-cell coordinates; `.obsp` contains `distances` and `connectivities`.

`.layers` is empty and `.raw` is absent: normalized/log-transformed/scaled
matrices are not saved, and raw counts are in `.X`. PCA gene loadings are not
saved in `.varm`. Saved coordinates can be redrawn without refitting, but
are not a complete fitted PCA model for projecting new cells.

`.uns` contains three dictionaries:

| Key | Contents |
| --- | --- |
| `neighbors` | Graph-key names and neighbor settings: metric, neighborhood size, PC count, method, seed and unintegrated status. |
| `primary_processing` | Inherited Step 00–02 input policies, QC calculation history and QC filtering counts/rules. |
| `step06_diagnostics` | Upstream run ID/checksum, selected normalization/HVG/PCA/graph/UMAP/Leiden settings, seed, bypass/no-correction records and the automatic provisional outcome code/label. |

The inherited `primary_processing.latest_step` remains `02_qc_filtering`;
`step` and `forbidden_operations` retain Step 00's original scope. They are
historical fields, not accurate current-stage or completed-operation
indicators. `step06_diagnostics.normalization` describes the working matrix
used to compute diagnostics, not the saved `.X`. Its input checksum refers
to Step 02; the output H5AD checksum is in `tables/output_manifest.tsv`.

The complete parameter set, PCA variance ratios, full HVG statistics,
cluster-level standardized marker scores, mixing metrics, code/software
identity and review status are in companion tables, frozen code/config and
`STEP_STATUS.tsv`; they are not all embedded in `.uns`. Keep the full run
package alongside the H5AD. The source repository's `STEP06_SAVED_ANNDATA.md`
provides the field-by-field inventory and exact September 3 metadata snapshot.

Watch a running package, including full function inputs and outputs, with:

```bash
PYTHONPATH="RUN_DIR/code" /home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python \
  -m primary_processing.step06_progress_cli \
  --run-dir "RUN_DIR" --follow --details
```

All numerical analyses use every cell. The deterministic balanced subset in
`tables/rendering_cell_ids.tsv.gz` affects rendering only. This run performs
no integration, batch correction, final annotation, or cell/gene removal.

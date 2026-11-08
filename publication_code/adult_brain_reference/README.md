# Adult brain reference mapping

- `reference_transfer.py`: shared-gene normalization, feature selection, reference projection, and label transfer.
- `plot_reference_rivers.py`: regional and fine-subtype river panels.
- `reference_labels.example.tsv`: reference-label hierarchy format.
- `figures/`: combined river-plot PNG.

Dependencies: Python 3.10 or newer and the packages in `requirements.txt`.

## Inputs

The transfer function accepts sparse raw counts in cells-by-genes orientation,
row-aligned cell metadata, and a common ordered gene list. Reference metadata
contains `unified_leaf_subtype`; query metadata contains its cell identifiers
and `div90_broad_class` annotations.

Reference atlas: Siletti et al. (2023).

The reference comprises adult MGE interneurons, LAMP5-LHX6/chandelier neurons,
and the selected CHAT-like Splatter population. CGE, medium spiny, and eccentric
medium spiny populations are excluded before transfer. The CHAT-like subset is
cluster 400, subclusters 1634–1638 and 1640–1642, annotated `NT-CHOL NT-VGLUT3`.
The selected CHAT-like reference cells also express GABAergic markers.

The river plotter accepts a cell-level prediction table with `predicted.id` and
`div90_broad_class`, plus a reference hierarchy with three columns:
`unified_leaf_subtype`, `unified_major_subtype_roi`, and
`unified_pallial_subpallial_bin`.

## Usage

Run from this folder. With matrices, metadata, and shared gene names loaded:

```python
from reference_transfer import transfer_labels

result = transfer_labels(reference_counts, query_counts,
                         reference_metadata, query_metadata, genes)
result["predictions"].to_csv("predictions.tsv", sep="\t", index=False)
result["label_scores"].to_csv("label_scores.tsv", sep="\t", index=False)
result["selected_features"].to_csv("selected_features.tsv", sep="\t", index=False)
```

```sh
python3 plot_reference_rivers.py --predictions predictions.tsv \
  --label-map reference_labels.tsv --outdir mapping_outputs
```

## Analysis and display

Counts are normalized per cell to 10,000 across the supplied shared genes and
log1p-transformed. The 3,000 highest-variance genes are selected jointly across
reference and query. A 50-component truncated SVD is fitted to reference cells;
query cells are projected into that basis. L2-normalized coordinates are used
for cosine nearest-neighbor search. Each query cell receives the label with the
greatest similarity-weighted support among its 50 nearest reference cells.
Vote scores measure relative neighbor support. All eligible reference cells
are used; the transfer function does not subsample cells.

Regional and major-subtype assignments are derived from the winning fine label
through the reference hierarchy. `Pallial/cortical` groups reference cells from
cerebral cortex and hippocampus; `Subpallial` groups the other selected sampled
regions. These labels describe reference sampling regions.

Query classes use one shared order, ranked by descending Pallial/cortical
assignment fraction before ribbon filtering. Ribbons with fewer than 10 query
cells are omitted. Node heights and ribbon widths are recalculated from the
displayed edges at each resolution, so displayed totals can differ between
panels. Full and displayed edge counts are exported separately. Pallial/cortical assignments are purple and subpallial assignments are pink.
The combined
figure is written as a 300 dpi PNG, PDF, and SVG.

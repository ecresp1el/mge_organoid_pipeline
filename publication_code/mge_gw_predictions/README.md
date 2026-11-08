# Cell-type and gestational-week prediction scores

- `plot_prediction_scores.py`: both prediction-score UMAP grids.
- `figures/`: cell-type and gestational-week prediction-score PNGs.

Dependencies: Python 3.10 or newer and the packages in `requirements.txt`.

## Inputs

Prepared cell-level TSV or TSV.gz tables containing:

| Columns | Content |
| --- | --- |
| `study_id`, `cell_id` | Study and cell identifiers |
| `sample`, `cluster` | Sample and cluster annotations |
| `umap_1`, `umap_2` | Cell coordinates |
| `shi_seurat_full_prediction_score_<label>` | Individual cell-type scores |
| `shi_seurat_full_week_prediction_score_GW<week>` | Individual gestational-week scores |

The script reads prediction scores and embeddings already present in these
tables. Scores are displayed on a fixed 0–1 grey-to-blue scale. DIV90 clusters
6 and 7 are excluded; Samarasinghe samples containing `Ctrl` are included.
DIV90 UMAP2 is reflected for display. Cells with nonfinite coordinates are
omitted from the relevant study row. All other eligible cells are plotted.

## Usage

Run from this folder:

```sh
python3 plot_prediction_scores.py --input predictions.tsv.gz \
  --outdir prediction_outputs
```

Multiple input tables can follow `--input`. The default output is two 150 dpi
PNGs. Use `--dpi` to adjust resolution or `--formats png pdf svg` for vector
exports. Per-panel cell counts and display-filter counts are written alongside
the generated plots.

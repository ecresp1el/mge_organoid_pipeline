# Step 06 saved AnnData: asset and metadata inventory

This inventory was checked directly against the saved H5AD, including decoding its `.uns` with AnnData's HDF5 reader. It describes the **existing September 3 checkpoint**, not an intended future schema. Documentation does not change the object or its review state.

## Main reusable asset

- Run: `06_technical_sample_batch_diagnostics_20260903_114110_6087eeb`.
- Job: `59983804`.
- File: `objects/pcdh19_step06_unintegrated_diagnostics.h5ad`.
- Size: **7,001,220,104 bytes** (7.00 GB; approximately 6.52 GiB).
- Shape: **446,349 cells × 19,071 genes**.
- Status: **IN_REVIEW**, recorded in the run's `STEP_STATUS.tsv` and the workflow approval ledger, not as a dedicated status field in this H5AD's `.uns`.
- Full path:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/primary_processing/06_technical_sample_batch_diagnostics/06_technical_sample_batch_diagnostics_20260903_114110_6087eeb/objects/pcdh19_step06_unintegrated_diagnostics.h5ad
```

The key output is a raw-count AnnData with reusable diagnostic coordinates, graphs and annotations. You can reload it to inspect cells/clusters or redraw the saved PCA/UMAP. It is not a complete persisted representation of every intermediate used in the analysis.

## What is stored in each AnnData slot

| Slot | Actual stored contents | Interpretation |
| --- | --- | --- |
| `.X` | Sparse CSR `int32` raw counts, 446,349 × 19,071 | Approved Step 02 counts, with no further cell/gene removal or count correction. |
| `.obs` | Registered sample/genotype/sex/design fields; source metadata and raw-count QC; Step 02 pass/exclusion metadata; ten `program_*` columns; `step06_leiden`; `step06_provisional_state` | Program columns contain per-cell **unstandardized mean log1p expression**. Provisional state is the cluster's winning standardized marker program, not a validated cell type. Only retained cells are in this object; complete Step 02 exclusion history remains in the upstream disposition table. |
| `.var` | Gene IDs, symbols, feature type, genome and inherited QC; `highly_variable_step06` | All 19,071 genes remain; the boolean flag selects the 3,000 diagnostic HVGs. |
| `.obsm['X_pca']` | 446,349 × 50 coordinate array | Saved joint, unintegrated PC scores for every cell. |
| `.obsm['X_umap']` | 446,349 × 2 coordinate array | Saved joint, unintegrated UMAP for every cell. |
| `.obsp['distances']` | Sparse cell-by-cell neighbor-distance graph | Stored edges supply adjacency for sample-entropy calculations. |
| `.obsp['connectivities']` | Sparse cell-by-cell connectivity graph | Shared unintegrated graph used in the diagnostic analysis. |
| `.uns` | Exactly three top-level dictionaries: `neighbors`, `primary_processing`, `step06_diagnostics` | Metadata; detailed inventory below. |
| `.layers` | Empty | No normalized, log-transformed, scaled or corrected expression layer. |
| `.raw` | Absent | Raw counts are in `.X`; absence of `.raw` does not mean counts were lost. |
| `.varm`, `.varp` | Empty | In particular, no PCA gene loadings are saved in `.varm`. |

## Read `.uns` as three separate records

`.uns` is AnnData's unstructured metadata dictionary. It can contain nested dictionaries, arrays and scalar values. It is not itself an expression layer, and it does not automatically provide a complete analysis history.

| Dictionary | Origin | Meaning |
| --- | --- | --- |
| `neighbors` | Added by Step 06 output assembly | Names of graph arrays in `.obsp` and selected graph-construction parameters. |
| `primary_processing` | Inherited unchanged from Step 02 | Original input policies and the earlier Step 00–02 processing history. Some inherited fields are stale if read as describing the entire Step 06 object. |
| `step06_diagnostics` | Added by Step 06 output assembly | Exact upstream input identity, selected diagnostic parameters, bypass/no-correction records, and the automatic outcome label. |

### Important inherited-field ambiguities

- `primary_processing.latest_step` is **still `02_qc_filtering`**. It was not advanced when Step 06 diagnostic fields were attached. It describes the inherited record and is not a reliable current-stage indicator for this object.
- `primary_processing.step` is **still `00_input_validation_and_canonical_anndata`**. It identifies original construction, not the latest computation.
- `primary_processing.forbidden_operations` is the **original Step 00 scope restriction**. Its entries include QC filtering, PCA, UMAP and clustering, although later steps did perform QC and diagnostics. It must not be read as proof those operations never occurred.
- `primary_processing.raw_matrix_policy` mentions registration for Step 04. This is historical registration of raw droplet matrices, not evidence Step 04 ran; Step 04 was subsequently skipped.
- `primary_processing.matrix_state` correctly describes the **saved `.X`** as raw integer counts. `step06_diagnostics.normalization` describes the **temporary working expression** used for diagnostics. Both can be true: the workflow reloads raw counts before saving its output.
- `step06_diagnostics.cells_removed = 0` means **Step 06** removed none. The nested Step 02 history separately records the earlier **4,439** QC removals.

For the existing object, use the run identity, `STEP_STATUS.tsv`, and `step06_diagnostics` together. The current approval state is external to `.uns`. The historical object has not been rewritten to repair these schema ambiguities.

### `.uns['neighbors']`: all fields

| Field | Saved value | Meaning |
| --- | --- | --- |
| `connectivities_key` | `connectivities` | Look up `.obsp['connectivities']`; this is a key name, not the graph itself. |
| `distances_key` | `distances` | Look up `.obsp['distances']`. |
| `params.integrated` | `False` | Graph was constructed without integration/batch correction. |
| `params.method` | `umap` | Neighbor-connectivity construction method; this alone is not a UMAP parameter record. |
| `params.metric` | `cosine` | Cell-to-cell distance metric in PCA space. |
| `params.n_neighbors` | `30` | Requested neighborhood size; stored nonself edge count can differ. |
| `params.n_pcs` | `50` | PCs supplied to the neighbor calculation. |
| `params.random_state` | `20260903` | Graph-construction random seed. |

### `.uns['primary_processing']`: all fields

| Field | Meaning / saved information |
| --- | --- |
| `canonical_input` | Original per-sample filtered Cell Ranger H5 matrices. |
| `cell_id_policy` | Global cell IDs join technical sample ID and Cell Ranger barcode. |
| `feature_policy` | Only `Gene Expression` features entered the canonical object. |
| `forbidden_operations` | Historical Step 00 restrictions; see ambiguity note above. |
| `genome` | `GRCm39`. |
| `latest_step` | Inherited `02_qc_filtering`; not updated to Step 06. |
| `matrix_state` | Saved `.X` remains sparse, unnormalized integer Cell Ranger counts. |
| `raw_matrix_policy` | Historical registration-only policy for raw droplet matrices. |
| `step` | Original Step 00 construction identifier. |
| `var_names_policy` | Ensembl gene IDs index genes; symbols are in `.var['gene_symbol']`. |
| `step01_qc_metrics.calculator` | `scanpy.pp.calculate_qc_metrics`. |
| `step01_qc_metrics.filtering` | `none` during that QC-metric step. |
| `step01_qc_metrics.mitochondrial_prefix` | `mt-`. |
| `step01_qc_metrics.percent_top` | Python `None` (JSON `null` below): no top-gene-percentage metrics requested. |
| `step01_qc_metrics.qc_vars` | `['mt']`. |
| `step01_qc_metrics.ribosomal_metrics_created` | `False`; no ribosomal fraction was manufactured for absent panel genes. |
| `step01_qc_metrics.ribosomal_prefixes` | `['Rpl', 'Rps']`, the prefix definitions. |
| `step01_qc_metrics.thresholds` | `none` during Step 01. |
| `step02_qc_filtering.cells_before` | `450788`. |
| `step02_qc_filtering.cells_after` | `446349`. |
| `step02_qc_filtering.cells_removed` | `4439`. |
| `step02_qc_filtering.combination` | Union of low counts, low detected genes, or high mitochondrial fraction. |
| `step02_qc_filtering.doublet_filtering` | `none`. |
| `step02_qc_filtering.flag_columns` | Three exact `*_5mad` source columns; listed in the snapshot below. |
| `step02_qc_filtering.genes_removed` | `0`. |
| `step02_qc_filtering.normalization` | `none` during Step 02. |
| `step02_qc_filtering.source` | Approved Step 01a per-sample MAD candidate flags. |
| `step02_qc_filtering.stringency_mad` | `5`. Exact thresholds and per-cell disposition are upstream assets. |

### `.uns['step06_diagnostics']`: all fields

| Field | Saved value | Meaning |
| --- | --- | --- |
| `input_step02_run_id` | `02_qc_filtering_20260830_124611_97e1bb5` | Approved upstream run. This is not the Step 06 output run ID. |
| `input_sha256` | Full digest in the snapshot below | Hash of the **input Step 02 file**, not this output H5AD. |
| `steps_03_05_bypassed` | `True` | Rejected Step 03 and skipped Steps 04–05 were not consumed. |
| `normalization` | `normalize_total target_sum=10000; log1p` | Temporary working expression used for diagnostics; `.X` is still raw. |
| `hvg_flavor` | `seurat_v3` | Feature selection on raw counts. |
| `n_top_genes` | `3000` | Selected diagnostic HVGs. |
| `pca_components` | `50` | Number of saved PCs. |
| `n_neighbors` | `30` | Requested graph neighborhood size. |
| `neighbor_metric` | `cosine` | Cell distance metric. |
| `umap_min_dist` | `0.3` | UMAP layout setting. |
| `umap_spread` | `1.0` | UMAP layout setting. |
| `leiden_resolution` | `1.0` | Clustering resolution; backend/iterations are not included here. |
| `random_seed` | `20260903` | Shared diagnostic seed. |
| `batch_correction` | `none` | No correction applied. |
| `integration` | `none` | No integrated representation fitted. |
| `cells_removed` | `0` | No additional removals in Step 06. |
| `provisional_outcome_code` | `1` | Automatic heuristic branch, not approval or statistical validation. |
| `provisional_outcome_label` | `Minimal sample-associated structure` | Unreviewed headline; its numerical evidence/rule must be read separately. |

## What `.uns` does not contain, and where to find it

| Missing from `.uns` / H5AD | Where it is available or consequence |
| --- | --- |
| Normalized/log/scaled expression matrices | Not persisted. Recompute from raw `.X` using the frozen parameters if needed. |
| PCA gene loadings and complete fitted PCA model | Not saved; `.varm` is empty. Coordinates support redraw/inspection, but are insufficient alone to project new cells with the identical fitted model. |
| Standard Scanpy `uns['pca']`, `uns['umap']`, `uns['leiden']`, `uns['log1p']` records | Absent. This object was rebuilt from raw input and selected diagnostics were attached; it is not the entire working Scanpy object. |
| PCA variance ratios | `tables/pca_variance_ratio.tsv`. |
| Full HVG statistics | `tables/highly_variable_genes.tsv.gz`; `.var` retains only the added boolean `highly_variable_step06` flag. |
| Marker lists, scoring/winner formula, scaling/solver settings, Leiden backend/iterations, headline thresholds | Frozen `code/primary_processing/step06_analysis.py` and `step06_metrics.py`; explanations in the [code and tuning guide](STEP06_CODE_AND_TUNING_GUIDE.md). |
| Cluster-level standardized marker means | `tables/marker_program_scores_by_cluster.tsv`; saved `.obs['program_*']` values are unstandardized per-cell means. |
| Entropy values, sample-center distances, state proportions, per-sample review flags | Numerical companion tables under `tables/`; the headline is not the evidence. |
| All configurable settings, including centroid PCs, rendering cap and DPI | `tables/analysis_parameters.tsv`, `config/submitted_step06.env`, `config/resolved.env`. |
| Output run ID as a dedicated Step 06 `.uns` field, review state, completion time, code identity, package versions | Run directory/name, `STEP_STATUS.tsv`, shared approval ledger, `provenance/`, `tables/software_versions.tsv`, frozen code/config. |
| Output H5AD SHA-256 | `tables/output_manifest.tsv`, whose H5AD row describes this output. |

**Keep the full run package with the H5AD.** The H5AD is the main reusable analysis object; the tables, frozen code/config, status and manifests supply provenance and intermediate summaries not embedded in it.

## Read just the metadata without loading the count matrix

With the run's pinned AnnData environment:

```python
from pathlib import Path
from pprint import pprint

import h5py
from anndata.io import read_elem

run_dir = Path('/absolute/path/to/the/Step06/run')
with h5py.File(
    run_dir / 'objects/pcdh19_step06_unintegrated_diagnostics.h5ad', 'r'
) as handle:
    # Decode only /uns using AnnData's reader; do not read the large /X matrix.
    metadata = read_elem(handle['uns'])
pprint(metadata)
```

For full analysis, use `anndata.read_h5ad(...)`; for read-only inspection of the large object, `anndata.read_h5ad(..., backed='r')` keeps `.X` on disk, although metadata and embeddings may still occupy memory. Close its backed file handle with `adata.file.close()` when finished. These inspection operations do not update the historical H5AD.

## Exact decoded `.uns` snapshot

The following is the directly observed `.uns` represented as JSON: NumPy scalar/array values have been converted to ordinary JSON values and Python `None` becomes `null`. The dictionaries, field names and values are preserved; this is not a proposed replacement schema.

```json
{
  "neighbors": {
    "connectivities_key": "connectivities",
    "distances_key": "distances",
    "params": {
      "integrated": false,
      "method": "umap",
      "metric": "cosine",
      "n_neighbors": 30,
      "n_pcs": 50,
      "random_state": 20260903
    }
  },
  "primary_processing": {
    "canonical_input": "per-sample sample_filtered_feature_bc_matrix.h5",
    "cell_id_policy": "<technical_sample_id>_<cellranger_barcode>",
    "feature_policy": "Gene Expression",
    "forbidden_operations": [
      "QC filtering",
      "normalization",
      "scaling",
      "PCA",
      "neighbors",
      "UMAP",
      "clustering",
      "annotation"
    ],
    "genome": "GRCm39",
    "latest_step": "02_qc_filtering",
    "matrix_state": "sparse_unnormalized_integer_cellranger_counts",
    "raw_matrix_policy": "registered_for_step04_not_loaded_into_canonical_object",
    "step": "00_input_validation_and_canonical_anndata",
    "step01_qc_metrics": {
      "calculator": "scanpy.pp.calculate_qc_metrics",
      "filtering": "none",
      "mitochondrial_prefix": "mt-",
      "percent_top": null,
      "qc_vars": [
        "mt"
      ],
      "ribosomal_metrics_created": false,
      "ribosomal_prefixes": [
        "Rpl",
        "Rps"
      ],
      "thresholds": "none"
    },
    "step02_qc_filtering": {
      "cells_after": 446349,
      "cells_before": 450788,
      "cells_removed": 4439,
      "combination": "low_total_counts OR low_n_genes_by_counts OR high_pct_counts_mt",
      "doublet_filtering": "none",
      "flag_columns": [
        "low_total_counts_5mad",
        "low_n_genes_by_counts_5mad",
        "high_pct_counts_mt_5mad"
      ],
      "genes_removed": 0,
      "normalization": "none",
      "source": "approved Step 01a per-sample MAD candidate flags",
      "stringency_mad": 5
    },
    "var_names_policy": "Ensembl gene ID; gene symbol in var['gene_symbol']"
  },
  "step06_diagnostics": {
    "batch_correction": "none",
    "cells_removed": 0,
    "hvg_flavor": "seurat_v3",
    "input_sha256": "fadba4a25a7b6b7320219b21c189b6325687493519ba0fe1bd27efd79606b103",
    "input_step02_run_id": "02_qc_filtering_20260830_124611_97e1bb5",
    "integration": "none",
    "leiden_resolution": 1.0,
    "n_neighbors": 30,
    "n_top_genes": 3000,
    "neighbor_metric": "cosine",
    "normalization": "normalize_total target_sum=10000; log1p",
    "pca_components": 50,
    "provisional_outcome_code": 1,
    "provisional_outcome_label": "Minimal sample-associated structure",
    "random_seed": 20260903,
    "steps_03_05_bypassed": true,
    "umap_min_dist": 0.3,
    "umap_spread": 1.0
  }
}
```

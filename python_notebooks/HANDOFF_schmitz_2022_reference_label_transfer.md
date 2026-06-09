# Handoff: Schmitz 2022 Reference Metadata And Label Transfer

Standalone workflow seed for using Schmitz et al. 2022 as a new reference in
cross-study comparisons.

This is separate from the Shi-reference workflow. Shi remains the existing
human fetal brain reference. Schmitz will be evaluated as a new inhibitory
neuron reference, with its own notebook and data-staging logic.

Study:

```text
Schmitz et al. 2022
The development and evolution of inhibitory neurons in primate cerebrum
Nature
DOI: 10.1038/s41586-022-04510-w
GitHub: https://github.com/mtvector/dev-and-evo-of-primate-inhibitory-neurons
UCSC Cell Browser: https://dev-inhibitory-neurons.cells.ucsc.edu/dev-inhibitory-neurons
```

## Current Run Target

Default run label:

```text
schmitz_2022_reference_metadata_umap_v1
```

Run directory:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/schmitz_2022_reference_metadata_umap/schmitz_2022_reference_metadata_umap_v1
```

## Code Added

Reusable module:

```text
python_notebooks/src/mge_organoid_python/schmitz_reference.py
```

CLI wrapper:

```text
python_notebooks/scripts/fetch_schmitz_2022_metadata_umap.py
```

Notebook:

```text
python_notebooks/notebooks/schmitz_2022_reference_metadata_umap.ipynb
```

## Data Staged

Scope:

```text
Metadata, Cell Browser configs, and decoded Cell Browser UMAP coordinates only.
Expression matrices were intentionally not downloaded.
```

Reason:

```text
The public UCSC `exprMatrix.tsv.gz` files are Float32 processed expression
matrices with negative centered values, not raw UMI count matrices. Treat raw
count/Seurat-reference use as a separate decision.
```

Staged datasets:

| dataset_id | cells | label field | coordinates |
| --- | ---: | --- | --- |
| `schmitz_macaque_dev` | 109,111 | `class` | `coords_0` / UMAP |
| `schmitz_mouse_dev` | 76,804 | `class` | `coords_0`, `coords_1`, `coords_2` |
| `schmitz_mouse_adult` | 141,069 | `class` | `coords_0`, `coords_1`, `coords_2` |

Main staged files:

```text
metadata/schmitz_macaque_dev.metadata_with_umap.tsv.gz
metadata/schmitz_mouse_dev.metadata_with_umap.tsv.gz
metadata/schmitz_mouse_adult.metadata_with_umap.tsv.gz

coords/{dataset_id}/{coord_name}.coords.bin
coords/{dataset_id}/{coord_name}.coords.tsv.gz

config/{dataset_id}.dataset.json
config/{dataset_id}.desc.json
```

Coordinate note:

```text
UMAP coordinates were decoded from UCSC Cell Browser UInt16 `coords.bin` files.
The output tables include raw browser UInt16 coordinates and normalized 0..1
browser-scale coordinates. These are good for reproducing the browser UMAP
geometry, but they are not the original Scanpy UMAP coordinate units.
```

## Metadata And Annotation Columns

Primary reference annotation:

```text
class
```

This corresponds to the paper/code concept named `supervised_name` in the
Schmitz GitHub analysis scripts.

Primary sample columns:

```text
schmitz_macaque_dev: batch_name, file_name
schmitz_mouse_dev:   batch_name
schmitz_mouse_adult: dataset_name, batch_name
```

Primary biology columns:

```text
timepoint
region
class
leiden
hires_leiden, where present
phase, where present
latent_time
Expressed Genes
UMI Count
latent_cell_probability
```

The Schmitz `class` labels are more granular than the Shi major-type labels.
They encode developmental source/trajectory plus marker anchors, e.g.
`MGE_LHX6/MAF`, `MGE_LHX6/NPY`, `LGE_FOXP1/ISL1`, `CGE_NR2F2/PROX1`,
`Transition`, and cell-cycle classes.

## Reference-Relevant MGE Labels

Current broad MGE-like rule used only for inventory:

```text
class starts with MGE_
or class starts with Ctx_LHX6
or class starts with Ctx_SST
```

Inventory result:

| dataset_id | MGE-like cells | percent |
| --- | ---: | ---: |
| `schmitz_macaque_dev` | 34,088 | 31.24% |
| `schmitz_mouse_dev` | 22,250 | 28.97% |
| `schmitz_mouse_adult` | 32,295 | 22.89% |

Major MGE-like classes:

```text
schmitz_macaque_dev:
  MGE_LHX6/MAF 26,967
  MGE_LHX6/NPY 4,077
  MGE_CRABP1/MAF 1,554
  MGE_CRABP1/TAC3 1,490

schmitz_mouse_dev:
  MGE_LHX6/MAF 19,564
  MGE_LHX6/NPY 1,945
  MGE_CRABP1/MAF 741

schmitz_mouse_adult:
  MGE_LHX6/MAF 16,391
  Ctx_LHX6/SST 6,822
  Ctx_LHX6/PVALB 4,589
  MGE_LHX6/NPY 1,968
  Ctx_LHX6/LAMP5 924
  Ctx_SST/NDNF 894
  MGE_CRABP1/MAF 707
```

## Output Tables

Key inventory tables:

```text
tables/schmitz_2022_datasets.tsv
tables/schmitz_2022_metadata_fields.tsv
tables/schmitz_2022_metadata_value_counts.tsv
tables/schmitz_2022_class_composition.tsv
tables/schmitz_2022_reference_label_summary.tsv
tables/schmitz_2022_umap_files.tsv
tables/schmitz_2022_reference_metadata_umap_complete.tsv
tables/previews/*.tsv
```

Completion table:

```text
dataset_cells:schmitz_macaque_dev            109111
metadata_with_umap_rows:schmitz_macaque_dev  109111
dataset_cells:schmitz_mouse_dev               76804
metadata_with_umap_rows:schmitz_mouse_dev     76804
dataset_cells:schmitz_mouse_adult            141069
metadata_with_umap_rows:schmitz_mouse_adult  141069
```

## Planned Schmitz-Reference Comparisons

New comparison direction:

```text
Use Schmitz as the reference and map project/cross-study target cells onto
Schmitz inhibitory-neuron annotations.
```

Initial target:

```text
Metadata/annotation comparison and marker-expression sanity checks before any
anchor-based label transfer.
```

Candidate reference choices:

```text
1. `schmitz_macaque_dev` as the main developmental primate inhibitory-neuron
   reference.
2. `schmitz_mouse_dev` as a developmental cross-species comparator, not a direct
   human/primate transfer reference unless ortholog handling is explicit.
3. `schmitz_mouse_adult` for adult interneuron endpoint context, not the first
   developmental transfer reference.
```

Potential target datasets:

```text
Varela DIV30
Varela DIV90
Walsh
Bershteyn 2025
Bershteyn 2023
Samarasinghe 2021
Siebert 2026
```

Do not fold this into the existing Shi prediction plots yet. The Schmitz
reference requires a separate output namespace and different label vocabulary.

## Next Implementation Step

Build a dedicated Schmitz-reference notebook/module stage:

```text
1. Load staged Schmitz metadata-with-UMAP tables.
2. Plot Schmitz UMAPs colored by `class`, `region`, `timepoint`, and MGE-like
   inventory labels.
3. Decide whether Schmitz reference labels should be:
   - full `class`
   - MGE/LGE/CGE coarse class
   - restricted MGE-like classes only
4. Decide expression source:
   - accept UCSC processed Float32 expression for marker/reference plotting only
   - reconstruct/download raw counts for Seurat-style label transfer
5. Only after expression-source validation, implement Schmitz-as-reference
   transfer with its own output root.
```

Suggested future output root:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/schmitz_2022_reference_label_transfer/
```

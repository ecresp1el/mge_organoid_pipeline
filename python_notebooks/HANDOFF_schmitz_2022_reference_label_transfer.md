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

Slurm template:

```text
slurm_templates/29_fetch_schmitz_2022_metadata_umap.sbatch.template
```

## Runtime Environment

Run this workflow from the repo checkout:

```text
/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
```

Large runtime data and outputs live under:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
```

Use the project conda environment:

```text
CONDA_ENV_BIN=/home/elcrespo/miniconda3/envs/mge-organoid-python/bin
python=/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python
env=mge-organoid-python
```

Required shell exports for direct command-line runs:

```bash
export PROJECT_ROOT=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
export PYTHONPATH=/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline/python_notebooks/src:${PYTHONPATH:-}
export PATH=/home/elcrespo/miniconda3/envs/mge-organoid-python/bin:${PATH}
```

Do not rely on the login shell's generic `python`; on this system it may be
missing or too old. Use the explicit conda-env Python path above.

## Slurm Management

Use Slurm for normal Schmitz metadata/UMAP reruns so the environment, logs, and
outputs follow the rest of the project.

Default Slurm template:

```text
slurm_templates/29_fetch_schmitz_2022_metadata_umap.sbatch.template
```

Prepared job file:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/29_fetch_schmitz_2022_metadata_umap.sbatch
```

Prepare and submit:

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
mkdir -p /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs
mkdir -p /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs

cp slurm_templates/29_fetch_schmitz_2022_metadata_umap.sbatch.template \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/29_fetch_schmitz_2022_metadata_umap.sbatch

sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/29_fetch_schmitz_2022_metadata_umap.sbatch
```

Template defaults:

```text
account: parent0
partition: standard
cpus-per-task: 4
memory: 32G
time: 04:00:00
job name: schmitz-ref-umap
```

Environment variables accepted by the template:

```text
REPO_ROOT=/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
PROJECT_ROOT=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
CONDA_ENV_BIN=/home/elcrespo/miniconda3/envs/mge-organoid-python/bin
SCHMITZ_2022_RUN_LABEL=schmitz_2022_reference_metadata_umap_v1
SCHMITZ_2022_OVERWRITE=0
```

Use `SCHMITZ_2022_OVERWRITE=1` only when intentionally redownloading staged
metadata/config/coordinate files:

```bash
sbatch --export=ALL,SCHMITZ_2022_OVERWRITE=1 \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/29_fetch_schmitz_2022_metadata_umap.sbatch
```

Expected logs:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/schmitz-ref-umap-schmitz-ref-umap-<jobid>.out
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/schmitz-ref-umap-schmitz-ref-umap-<jobid>.err
```

Monitor:

```bash
squeue -j <jobid> -o '%.18i %.9P %.28j %.8u %.2t %.10M %.6D %R'
sacct -j <jobid> --format=JobID,JobName%28,State,ExitCode,Elapsed,MaxRSS,ReqMem -P
tail -n 80 /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/schmitz-ref-umap-schmitz-ref-umap-<jobid>.out
tail -n 80 /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/schmitz-ref-umap-schmitz-ref-umap-<jobid>.err
```

The Slurm job runs:

```bash
/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python \
  python_notebooks/scripts/fetch_schmitz_2022_metadata_umap.py \
  --project-root "$PROJECT_ROOT" \
  --run-label "$SCHMITZ_2022_RUN_LABEL"
```

It writes the standardized table registry, UMAP plot manifest, and completion
table to stdout at the end of the job.

Current execution status:

```text
The current `schmitz_2022_reference_metadata_umap_v1` outputs were generated
with the same conda environment by direct CLI execution while developing the
workflow. The Slurm template and prepared job file are in place for future
reruns, but no duplicate Slurm rerun has been submitted yet.
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

tables/per_dataset/schmitz_macaque_dev_schmitz_reference_umap_cells.tsv.gz
tables/per_dataset/schmitz_mouse_dev_schmitz_reference_umap_cells.tsv.gz
tables/per_dataset/schmitz_mouse_adult_schmitz_reference_umap_cells.tsv.gz
tables/schmitz_2022_reference_umap_cells.tsv.gz

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

## Standardized Cell Tables

The staging script now writes cross-study-style standardized cell tables for
each Schmitz dataset and one combined table.

Standardized table columns:

```text
cell_id
study_id
dataset_id
species
study_label
sample
batch_name
dataset_name
region
timepoint
class
cluster
leiden
hires_leiden
phase
latent_time
expressed_genes
umi_count
coord_name
umap_1
umap_2
umap_1_u16
umap_2_u16
```

Important Schmitz-specific rule:

```text
Keep `species` and `dataset_id` attached to every cell. Unlike most other
cross-study objects, Schmitz is not one single-species UMAP. The combined table
contains macaque and mouse rows, but the UMAP coordinates are dataset-local and
should be plotted as facets/per-dataset panels, not overlaid as one shared
embedding.
```

Standardized table registry:

```text
tables/schmitz_2022_standardized_cell_tables.tsv
```

Current standardized table row counts:

```text
schmitz_macaque_dev   109,111
schmitz_mouse_dev      76,804
schmitz_mouse_adult   141,069
combined              326,984
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
tables/schmitz_2022_standardized_cell_tables.tsv
tables/schmitz_2022_umap_plot_manifest.tsv
tables/schmitz_2022_reference_metadata_umap_complete.tsv
tables/previews/*.tsv
```

## UMAP Plots

UMAP output directory:

```text
plots/umap/
```

Current plot manifest:

```text
tables/schmitz_2022_umap_plot_manifest.tsv
```

Plots generated:

```text
plots/umap/schmitz_2022_umap_facets_by_species_dataset.png

plots/umap/schmitz_macaque_dev/schmitz_macaque_dev_umap_by_class.png
plots/umap/schmitz_macaque_dev/schmitz_macaque_dev_umap_by_region.png
plots/umap/schmitz_macaque_dev/schmitz_macaque_dev_umap_by_timepoint.png

plots/umap/schmitz_mouse_dev/schmitz_mouse_dev_umap_by_class.png
plots/umap/schmitz_mouse_dev/schmitz_mouse_dev_umap_by_region.png
plots/umap/schmitz_mouse_dev/schmitz_mouse_dev_umap_by_timepoint.png

plots/umap/schmitz_mouse_adult/schmitz_mouse_adult_umap_by_class.png
plots/umap/schmitz_mouse_adult/schmitz_mouse_adult_umap_by_region.png
plots/umap/schmitz_mouse_adult/schmitz_mouse_adult_umap_by_timepoint.png
```

Plot interpretation:

```text
The species/dataset overview is faceted by Schmitz dataset. Do not interpret it
as a joint integrated mouse-macaque embedding.
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

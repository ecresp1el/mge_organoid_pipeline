# DIV30 Paper Cluster Annotations and First URD Handoff

## Purpose

This handoff records the DIV30 cell-state labels that are mapped from the paper/manual Seurat cluster annotations and documents the first URD lineage run setup.

The immediate lineage question is whether a radial-glia-rooted ordering is biologically coherent before comparing alternative roots:

1. Jia RGC1 root
2. Jia RGC1 + RGC2 root
3. Shi high-RGC-score root

The first URD pass uses the paper/manual Radial glia annotation as the root because the user clarified that the published annotations should be mapped directly onto the DIV30 Seurat cluster IDs.

## Source Data

Primary DIV30 Seurat object:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/varela_this_paper/varela_this_paper_seurat.rds
```

Jia DIV30 score table used as the metadata sidecar:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/jia_program_div30_scoring/jia_program_div30_scoring_v1/tables/div30_jia_program_scores_obs.tsv
```

UMAP table used for plotting the mapped labels:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/varela_div30_2f0j5mwk/umap.tsv
```

Shi Seurat label-transfer metadata available for later root comparison:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_reference_div30_seurat_label_transfer/shi_reference_div30_seurat_label_transfer_v1/tables/div30_shi_seurat_label_transfer_obs.tsv.gz
```

## Cluster Columns

The Jia score table contains both `seurat_clusters` and `RNA_snn_res.0.2`.
They were checked cell-by-cell and match for all 90,631 DIV30 cells, so the paper/manual annotation can be mapped without ambiguity.

Counts by cluster:

| Seurat cluster | Cells | Paper/manual annotation |
|---:|---:|---|
| 0 | 19,148 | Radial glia |
| 1 | 18,907 | SST+ cIN |
| 2 | 17,287 | MGE subpallial neurons |
| 3 | 15,932 | Radial glia |
| 4 | 13,750 | PV neuron precursor |
| 6 | 3,696 | Inhibitory progenitors |
| 7 | 1,911 | Radial glia |

Collapsed paper/manual counts:

| Paper/manual annotation | Cells |
|---|---:|
| Radial glia | 36,991 |
| Inhibitory progenitors | 3,696 |
| PV neuron precursor | 13,750 |
| SST+ cIN | 18,907 |
| MGE subpallial neurons | 17,287 |

## Mapping Step

Script:

```text
python_notebooks/scripts/map_div30_paper_cluster_annotations.py
```

Output root:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_paper_cluster_annotations/div30_paper_cluster_annotations_v1
```

Key outputs:

```text
tables/div30_paper_cluster_annotation_mapping.tsv
tables/div30_paper_cluster_annotation_counts.tsv
tables/div30_paper_cluster_annotation_counts_by_sample.tsv
tables/div30_jia_scores_with_paper_cluster_annotations.tsv.gz
plots/div30_umap_paper_cluster_annotations_and_seurat_clusters.png
```

The script is intentionally a sidecar. It does not rewrite the Seurat object or Jia score table.

## First URD Run Definition

First root definition:

```text
paper_cluster_annotation == "Radial glia"
```

This corresponds to Seurat clusters:

```text
0, 3, 7
```

The first URD pass is intended to produce a root-specific diffusion/flood pseudotime before full branching/tip-walk reconstruction. The expected outputs are:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_v1/
```

with pseudotime and summary tables/plots under that run directory.

The URD runner now writes the imposed parameter table before doing heavy work:

```text
tables/div30_first_urd_parameters.tsv
```

This first pass is intentionally a stratified pilot by default:

```text
MAX_CELLS=30000
```

It samples across the mapped paper/manual annotations, records the exact selected cells in the input metadata, and can be promoted to a full-cell run by setting:

```text
MAX_CELLS=0
```

before submitting the Slurm job.

Workflow files:

```text
python_notebooks/scripts/export_div30_first_urd_inputs.py
scripts/14_div30_first_urd.R
slurm_templates/31_div30_first_urd.sbatch.template
```

The exporter reads the existing Matrix Market count export from the AnnData conversion rather than reopening the large Seurat object. The URD runner reads only plain Matrix Market counts plus cell metadata, which avoids mixing Seurat-dependent R packages with the older URD package.

## URD Processing Logic and Imposed Parameters

The workflow is split into two explicit stages.

### Stage 1: Python input export

Script:

```text
python_notebooks/scripts/export_div30_first_urd_inputs.py
```

Inputs:

```text
results/python_anndata/varela_div30_2f0j5mwk/matrix_counts.mtx
results/python_anndata/varela_div30_2f0j5mwk/barcodes.tsv
results/python_anndata/varela_div30_2f0j5mwk/features.tsv
results/python_anndata/varela_div30_2f0j5mwk/umap.tsv
results/div30_paper_cluster_annotations/div30_paper_cluster_annotations_v1/tables/div30_jia_scores_with_paper_cluster_annotations.tsv.gz
results/shi_reference_div30_seurat_label_transfer/shi_reference_div30_seurat_label_transfer_v1/tables/div30_shi_seurat_label_transfer_obs.tsv.gz
```

Logic imposed by the exporter:

```text
Join key: cell_id
Matrix orientation: features x cells
Default selection: stratified pilot capped at MAX_CELLS=30000
Stratification column: paper_cluster_annotation
Sampling seed: SEED=7
Root candidate flag: urd_root_candidate = paper_cluster_annotation == "Radial glia"
```

The exporter does not normalize expression, choose variable genes, compute PCA, compute diffusion maps, or run pseudotime. It only writes a plain input bundle and a manifest.

Stage 1 outputs:

```text
inputs/div30_first_urd_counts.mtx
inputs/div30_first_urd_features.tsv
inputs/div30_first_urd_barcodes.tsv
inputs/div30_first_urd_cell_metadata.tsv
inputs/div30_first_urd_input_manifest.tsv
```

### Stage 2: R/URD lineage run

Script:

```text
scripts/14_div30_first_urd.R
```

The R script has been refactored into named processing functions so the lineage logic is isolated:

```text
build_config()
parameter_table()
read_urd_input_bundle()
create_filtered_urd()
select_variable_genes()
run_urd_geometry()
run_flood_pseudotime()
build_pseudotime_table()
write_tables()
write_plots()
run_pipeline()
```

Current/default URD parameters:

| Stage | Parameter | Value | Meaning |
|---|---:|---:|---|
| root | root_label | Radial glia | Root candidates are cells flagged by the exporter from paper/manual annotations. |
| filter | min_genes | 500 | `createURD()` keeps cells with at least this many detected genes. |
| filter | min_cells | 3 | `createURD()` keeps genes detected in at least this many cells. |
| filter | min_counts | 10 | `createURD()` keeps genes with at least this many total counts. |
| variable genes | num_variable_genes | 3000 | Top variance genes from `urd@logupx.data` stored in `urd@var.genes`. |
| PCA | pca_mp_factor | 2 | `calcPCA(mp.factor = 2)`. |
| diffusion map | knn | 100 | `calcDM(knn = 100)`. |
| diffusion map | sigma | local | `calcDM(sigma = "local")`. |
| flood pseudotime | n_floods | 20 | `floodPseudotime(n = 20)`. |
| flood pseudotime | minimum_cells_flooded | 2 | `floodPseudotime(minimum.cells.flooded = 2)`. |
| flood processing | max_frac_NA | 0.4 | `floodPseudotimeProcess(max.frac.NA = 0.4)`. |
| flood processing | pseudotime_fun | mean | Flood replicate pseudotimes are summarized with `mean`. |
| flood processing | stability_div | min(10, n_floods) | Default is 10 for the current 20-flood run. |
| seed | seed | 7 | Seed set before PCA/diffusion/flood steps. |

Stage 2 outputs:

```text
tables/div30_first_urd_parameters.tsv
tables/div30_first_urd_pseudotime.tsv
tables/div30_first_urd_summary.tsv
tables/div30_first_urd_marker_correlations.tsv
tables/div30_first_urd_pseudotime_by_paper_cluster.tsv
plots/div30_first_urd_umap_pseudotime.png
plots/div30_first_urd_pseudotime_by_paper_cluster.png
div30_first_urd_object.rds
```

Note: job `51670320` was already running when the code was refactored to write `div30_first_urd_parameters.tsv`, so that specific job may not emit the new parameter table unless rerun. The parameters above document the effective defaults used for the submitted job.

## URD Environment Status

URD was not installed in the conda R environment:

```text
/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/Rscript
```

The first install attempt there failed because the conda environment did not expose the compiler wrapper executables needed for source packages at that time:

```text
x86_64-conda-linux-gnu-cc
x86_64-conda-linux-gnu-c++
x86_64-conda-linux-gnu-gfortran
```

A later check showed those wrapper paths are now visible inside `mge-organoid-python`, so the conda R route may be viable after resolving `destiny`:

```text
/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/x86_64-conda-linux-gnu-cc
/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/x86_64-conda-linux-gnu-c++
/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/x86_64-conda-linux-gnu-gfortran
```

The Great Lakes module R path does expose working compilers:

```text
module load R/4.4.3
/sw/pkgs/arc/stacks/gcc/13.2.0/R/4.4.3/bin/Rscript
/sw/pkgs/arc/gcc/13.2.0/bin/gcc
/sw/pkgs/arc/gcc/13.2.0/bin/g++
/sw/pkgs/arc/gcc/13.2.0/bin/gfortran
```

The current URD installation attempt uses the official GitHub package:

```text
farrellja/URD@v1.1.1
```

Resolved package blocker:

```text
module load R/4.4.3
URD 1.1.1
destiny 2.14.0
Matrix 1.7.2
ggplot2 4.0.3
```

`destiny` was installed from the Bioconductor 3.9 source archive and patched to remove the hard `VIM::hotdeck` import. That keeps this first URD run aligned with URD's older dependency era and avoids the modern `VIM`/`xgboost` dependency chain. The input matrix is expected to contain no missing values, so the removed imputation path should not be used for this workflow.

The first submission attempt failed before job creation because `standard` nodes advertise about 184 GB memory and the initial template requested `--mem=220G`. The template now requests `--mem=160G`, matching the memory profile used by the other large DIV30/Seurat jobs in this repo.

## First URD Submission

Submitted job:

```text
51670320
```

Submission command:

```text
cp slurm_templates/31_div30_first_urd.sbatch.template /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/31_div30_first_urd.sbatch
sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/31_div30_first_urd.sbatch
```

Immediate queue state after submission:

```text
Initially PENDING on standard, reason: Priority.
Then RUNNING on standard node gl3327.
```

Expected log path once the job starts:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/31_div30_first_urd_51670320.log
```

Observed launch progress:

```text
Export completed with selected_cells=30000.
Selected paper/manual annotation counts:
  Radial glia               12246
  SST+ cIN                   6258
  MGE subpallial neurons     5722
  PV neuron precursor        4551
  Inhibitory progenitors     1223

The R package gate passed for Matrix, URD, and ggplot2.
The URD script reached count loading and createURD().
createURD() completed.
Root cells retained after filtering: 12246.
Variable genes stored in URD object: 3000.
PCA completed; 24 PCs exceeded 2x the Marchenko-Pastur null upper bound and 48 PCs were stored.
The job then entered calcDM(knn = 100, sigma = "local").
```

## Interpretation Checks for First URD

After the first URD pseudotime exists, evaluate whether pseudotime increases from root/progenitor toward neuronal states by comparing:

```text
jia_score_RGC1
jia_score_RGC2
jia_score_IPC
DLX2 expression
ASCL1 expression
DCX expression
neuronal maturation marker expression
paper_cluster_annotation
shi_seurat_full_developmental_class
```

The desired pattern for a coherent developmental ordering is high root/RGC signal at low pseudotime, increasing IPC/progenitor transition signal after the root, and increasing neuronal/maturation markers later in pseudotime.

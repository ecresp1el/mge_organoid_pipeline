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

## URD Environment Status

URD was not installed in the conda R environment:

```text
/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/Rscript
```

Installing there failed because the conda environment does not have the compiler wrapper executables needed for source packages:

```text
x86_64-conda-linux-gnu-cc
x86_64-conda-linux-gnu-c++
x86_64-conda-linux-gnu-gfortran
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

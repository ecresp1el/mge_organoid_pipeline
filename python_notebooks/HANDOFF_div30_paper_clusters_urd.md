# DIV30 Paper Cluster Annotations and First URD Handoff

## Purpose

This handoff records the DIV30 cell-state labels that are mapped from the paper/manual Seurat cluster annotations and documents the first URD lineage run setup.

The immediate lineage question is whether a radial-glia-rooted ordering is biologically coherent before comparing alternative roots:

1. Jia RGC1 root
2. Jia RGC1 + RGC2 root
3. Shi high-RGC-score root

The first URD pass uses the paper/manual Radial glia annotation as the root because the user clarified that the published annotations should be mapped directly onto the DIV30 Seurat cluster IDs.

DIV90 companion workflow:

```text
python_notebooks/HANDOFF_div90_jia_lineage_urd_plan.md
```

The DIV90 workflow now has a corrected v2 Jia-lineage smoke run for later DIV30-vs-DIV90 comparison. It uses cluster 12 as the Jia RootScore root, excludes glial/OPC/stress clusters, uses clusters `0+5+8`, `1`, and `2` as Jia endpoint tips, and retains clusters `3`/`11` as unassigned candidates for post-tree marker projection rather than tips.

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

The `python_anndata/varela_div30_2f0j5mwk/` directory is the Matrix Market export directory created during the Seurat-to-AnnData conversion job, not an independent data source. Its count matrix traces back to the same DIV30 Seurat object:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/varela_this_paper/varela_this_paper_seurat.rds
```

That project-local Seurat object is documented in `WORKFLOW.md` as a true copy of the original legacy DIV30 object:

```text
/nfs/turbo/umms-parent/mgeo_scRNAseq/day30_old/Day30.rds
```

Conversion log:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/08_convert_python_anndata_50452847.log
```

The log records this exact conversion command:

```text
Rscript python_notebooks/scripts/seurat_export_for_anndata.R \
  --seurat /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/varela_this_paper/varela_this_paper_seurat.rds \
  --outdir /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/varela_div30_2f0j5mwk \
  --assay RNA \
  --reduction umap \
  --expression_layer data
```

The export manifest beside the matrix says:

```text
assay                  RNA
expression_layer_requested data
layer_used_data        data
layer_used_counts      counts
reduction              umap
n_features             18082
n_cells                90631
```

Therefore `matrix_counts.mtx` is the Seurat `RNA` assay `counts` layer from the DIV30 Seurat object, exported as features x cells. The parallel `matrix_data.mtx` and legacy `matrix.mtx` files are the Seurat `RNA` assay `data` layer, but the URD input exporter intentionally uses `matrix_counts.mtx`.

Important integration clarification:

```text
The UMAP coordinates are taken as-is from the Seurat object's existing `umap`
reduction. That UMAP may reflect the prior Seurat processing/integration/batch
handling used to build the object.

The URD expression input is different: it is the Seurat `RNA` assay `counts`
layer exported from that same object after SeuratObject::JoinLayers().
These are raw count values stored in the post-QC/post-processing Seurat object,
not an integrated or batch-corrected count matrix. Seurat integration changes
the low-dimensional/reduction space or integrated assay/data representation; it
does not create "integrated raw counts" in the RNA counts layer.
```

So this first URD run uses:

```text
UMAP for plotting/context: existing Seurat umap reduction
Expression for URD: RNA counts layer
Metadata/root labels: sidecar annotations joined by cell_id
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

### Current Runtime Diagnosis

The long-running step is not input export, URD object creation, normalization, variable-gene selection, or PCA. Those completed on the 30k run.

The bottleneck is:

```text
URD::calcDM(urd, knn = 100, sigma = "local")
```

In URD 1.1.1, `calcDM()` does not use the PCA coordinates produced by `calcPCA()`. Its source code does:

```r
data.use <- t(object@logupx.data[genes.use, cells.use])
data.use <- as.matrix(data.use)
dm <- DiffusionMap(data.use, sigma = sigma.use, k = knn, ...)
```

Because `genes.use` defaults to `object@var.genes`, the current run computes the diffusion map on:

```text
30,000 cells x 3,000 variable genes
```

rather than on:

```text
30,000 cells x 48 stored PCs
```

This explains the observed behavior:

```text
calcPCA completed.
calcDM is still running.
The R process uses about one CPU core.
Memory use is modest relative to the 160G request.
```

Main issue: the current URD diffusion-map implementation densifies the selected log-normalized expression matrix and calls `destiny::DiffusionMap()` on 3,000 genes with `k=100`. The slowdown is algorithmic/single-threaded, not a Slurm memory problem and not something fixed by simply requesting more CPUs.

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

## Parallel Smoke URD Submission

Submitted a separate smoke run while leaving job `51670320` running.

Smoke job:

```text
51674431
```

Submission command:

```text
sbatch \
  --job-name=div30_urd_sm5k \
  --export=ALL,RUN_LABEL=div30_first_urd_paper_radial_glia_smoke5k_knn100_v1,MAX_CELLS=5000,SEED=7,URD_SEED=7,ROOT_LABEL="Radial glia",URD_KNN=100 \
  slurm_templates/31_div30_first_urd.sbatch.template
```

What stays biologically the same:

```text
Same source DIV30 RNA counts export.
Same paper/manual Radial glia root definition.
Same paper/Jia/Shi metadata sidecars.
Same URD_KNN=100.
Same default URD filters, PCA parameter, flood count, and flood processing parameters.
```

What changes only to make it a smoke test:

```text
RUN_LABEL=div30_first_urd_paper_radial_glia_smoke5k_knn100_v1
MAX_CELLS=5000
```

Immediate queue state after submission:

```text
Initially PENDING on standard, reason: Priority.
Then RUNNING on standard node gl3398.
Completed successfully with ExitCode 0:0 after 00:11:55.
```

Expected smoke log path:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/31_div30_first_urd_51674431.log
```

Expected smoke output root:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/
```

Smoke outputs generated:

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

Smoke result summary:

```text
n_cells_urd       5000
n_root_cells_urd  2040
n_variable_genes  3000
knn               100
sigma             local
n_floods          20
finite pseudotime 5000/5000
pseudotime range  0.0000 to 0.5440
```

Median smoke pseudotime by paper/manual annotation:

| Paper/manual annotation | Median URD pseudotime | Cells |
|---|---:|---:|
| Radial glia | 0.000 | 2040 |
| Inhibitory progenitors | 0.159 | 204 |
| SST+ cIN | 0.294 | 1043 |
| MGE subpallial neurons | 0.335 | 954 |
| PV neuron precursor | 0.354 | 759 |

Smoke marker/program correlations with URD pseudotime:

| Feature | Spearman rho |
|---|---:|
| `jia_score_RGC1` | -0.677 |
| `jia_score_RGC2` | -0.797 |
| `jia_score_IPC` | 0.438 |
| `logupx_DLX2` | 0.268 |
| `logupx_ASCL1` | -0.204 |
| `logupx_DCX` | 0.423 |
| `neuronal_maturation_score` | 0.607 |
| `shi_seurat_full_prediction_score_MGE` | 0.860 |
| `shi_seurat_full_prediction_score_progenitor` | -0.887 |
| `shi_seurat_full_expected_shi_week_numeric` | -0.185 |

Initial interpretation: the 5k smoke run did not fail. The paper/manual Radial glia root has pseudotime 0, RGC scores decrease with pseudotime, IPC/DCX/neuronal maturation increase with pseudotime, and later neuronal paper/manual annotations have higher median pseudotime than Radial glia. This supports the workflow mechanics and suggests the root is biologically plausible in this smoke run, but the full 30k run is still needed for the larger pilot.

### Smoke Lineage Decision-Tree Report

Report script:

```text
scripts/15_div30_urd_lineage_decision_report.R
```

Command run on the completed smoke URD object:

```text
Rscript scripts/15_div30_urd_lineage_decision_report.R \
  --urd-rds /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/div30_first_urd_object.rds \
  --outdir /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/lineage_decision_report \
  --annotation-col paper_cluster_annotation \
  --top-n 50 \
  --correlation-genes variable
```

Report output root:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/lineage_decision_report/
```

Report outputs:

```text
lineage_decision_tree_report.md
tables/root_annotation_composition.tsv
tables/pseudotime_ordering_by_annotation.tsv
tables/top_negative_pseudotime_genes.tsv
tables/top_positive_pseudotime_genes.tsv
tables/branch_structure_status.tsv
tables/decision_genes_between_branches.tsv
tables/flood_stability_summary.tsv
tables/gene_cascade_heatmap_matrix.tsv
plots/umap_pseudotime.png
plots/diffusion_map_pseudotime.png
plots/diffusion_map_annotation.png
plots/flood_stability.png
plots/tree_visualization.png
plots/gene_cascade_heatmap.png
plots/lineage_decision_tree.png
```

Branch status:

```text
The URD object has a `tree` slot, but it is empty (`tree_length = 0`).
This object contains flood pseudotime, not a reconstructed URD branch tree.
Therefore branch-specific decision genes are not available yet.
```

Top root-associated genes by negative Spearman correlation with pseudotime:

```text
NOTCH1, SPARC, ADGRV1, TTYH1, FKBP10, FGFR2, TEAD2, HES5
```

Top terminal-associated genes by positive Spearman correlation with pseudotime:

```text
LRRC7, CELF3, GAD1, L1CAM, AKAP6, FNBP1L, GPC2, KALRN
```

These are the current smoke-run transition-driving genes because they are ranked directly by Spearman correlation between expression and URD pseudotime across the variable genes used in the diffusion map. Negative genes are root-associated; positive genes are terminal/late-associated.

Flood stability summary: pseudotime estimates converge smoothly as walks per cell increase. Spearman correlation to the final pseudotime rises from 0.917 at 2 walks/cell to 1.000 at 20 walks/cell; median absolute delta falls from 0.0194 to 0.

Requested validation layers now generated:

```text
1. UMAP pseudotime: plots/umap_pseudotime.png
2. Pseudotime by annotation: tables/pseudotime_ordering_by_annotation.tsv and plots/lineage_decision_tree.png
3. Diffusion map plots: plots/diffusion_map_pseudotime.png and plots/diffusion_map_annotation.png
4. Flood stability: tables/flood_stability_summary.tsv and plots/flood_stability.png
5. Tree visualization: plots/tree_visualization.png
6. Branch-specific genes: tables/decision_genes_between_branches.tsv documents not available because no populated tree exists
7. Gene cascade heatmap: plots/gene_cascade_heatmap.png and tables/gene_cascade_heatmap_matrix.tsv
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

## Actual URD Tree Smoke Run

Goal for this step: stop adding pseudotime-only visualizations and ask whether the three neuronal paper/manual populations form distinct URD branches after the radial glia-rooted pseudotime has been established.

Tree build scripts:

```text
scripts/16_div30_urd_build_lineage_tree.R
scripts/17_div30_urd_finalize_lineage_tree_report.R
```

Why two scripts: `buildTree()` and random-walk processing can leave a large in-memory URD object. The build script now stops cleanly after saving the tree RDS and writes a `finalize_lineage_tree_report.command.txt`. The finalizer then starts a fresh R process to make tables, branch-gene summaries, and lightweight tree plots from the saved tree object. This avoids carrying the random-walk/buildTree memory footprint into plotting and branch-gene reporting.

Tree Slurm template:

```text
slurm_templates/32_div30_urd_lineage_tree.sbatch.template
```

The first URD Slurm template now also runs the smoke panel automatically after `div30_first_urd_object.rds` is saved:

```text
slurm_templates/31_div30_first_urd.sbatch.template
```

Every future first-URD run should therefore emit:

```text
lineage_decision_report/plots/umap_pseudotime.png
lineage_decision_report/plots/diffusion_map_pseudotime.png
lineage_decision_report/plots/diffusion_map_annotation.png
lineage_decision_report/plots/flood_stability.png
lineage_decision_report/plots/gene_cascade_heatmap.png
lineage_decision_report/tables/top_positive_pseudotime_genes.tsv
lineage_decision_report/tables/top_negative_pseudotime_genes.tsv
```

Tree root and tips used:

| Role | Annotation | URD tree ID |
|---|---|---:|
| Root | Radial glia | root cells, not a tip |
| Tip | SST+ cIN | 1 |
| Tip | PV neuron precursor | 2 |
| Tip | MGE subpallial neurons | 3 |

Numeric tip IDs are intentional. URD's tree internals convert segment names to numeric values in `buildTree()`, so the biological labels are stored in `tables/tree_tip_mapping.tsv` while the actual tree uses IDs `1`, `2`, and `3`.

Tree command run on the completed 5k smoke URD object:

```text
Rscript scripts/16_div30_urd_build_lineage_tree.R \
  --urd-rds /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/div30_first_urd_object.rds \
  --outdir /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/lineage_tree_radial_glia_tips_v1 \
  --annotation-col paper_cluster_annotation \
  --pseudotime-name paper_radial_glia_root \
  --root-label "Radial glia" \
  --tip-labels "SST+ cIN,PV neuron precursor,MGE subpallial neurons" \
  --n-per-tip 5000 \
  --seed 7
```

Finalizer command run from the saved tree object:

```text
Rscript scripts/17_div30_urd_finalize_lineage_tree_report.R \
  --tree-rds /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/lineage_tree_radial_glia_tips_v1/div30_urd_lineage_tree_object.rds \
  --outdir /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/lineage_tree_radial_glia_tips_v1 \
  --annotation-col paper_cluster_annotation \
  --pseudotime-name paper_radial_glia_root \
  --tip-mapping /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/lineage_tree_radial_glia_tips_v1/tables/tree_tip_mapping.tsv \
  --top-n-branch-genes 50
```

Tree output root:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/lineage_tree_radial_glia_tips_v1/
```

Tree outputs now generated:

```text
div30_urd_lineage_tree_object.rds
div30_urd_tip_random_walks.rds
finalize_lineage_tree_report.command.txt
urd_lineage_tree_report.md
tables/lineage_tree_parameters.tsv
tables/pseudotime_logistic_parameters.tsv
tables/tree_tip_mapping.tsv
tables/tree_status.tsv
tables/tree_segment_joins.tsv
tables/tree_tip_composition.tsv
tables/branch_specific_genes.tsv
plots/pseudotime_logistic.png
plots/breakpoint_decisions-1.pdf
plots/breakpoint_decisions-3.pdf
plots/urd_tree_annotation.png
plots/urd_tree_pseudotime.png
```

Core URD calls completed:

```text
pseudotimeDetermineLogistic()
pseudotimeWeightTransitionMatrix()
simulateRandomWalksFromTips()
processRandomWalksFromTips()
buildTree()
```

Important tree parameters imposed:

| Parameter | Value |
|---|---:|
| `n.per.tip` | 5000 |
| `optimal.cells.forward` | 100 |
| `max.cells.back` | 40 |
| `pseudotime.direction` | `<` |
| `minimum.visits` | 10 |
| `visit.threshold` | 0.7 |
| `cells.per.pseudotime.bin` | 80 |
| `bins.per.pseudotime.window` | 5 |
| `divergence.method` | `ks` |

Tree result:

```text
tree_slot_length      18
n_requested_tips      3
n_segment_joins       4
n_segments            5
has_distinct_branching TRUE
```

Segment joins:

| Parent | Child | Pseudotime |
|---:|---:|---:|
| 4 | 2 | 0.1626288 |
| 4 | 3 | 0.1626288 |
| 5 | 1 | 0.0000000 |
| 5 | 4 | 0.0000000 |

Interpretation of the smoke tree: PV neuron precursor (`2`) and MGE subpallial neurons (`3`) form a shared branch at pseudotime ~0.163. SST+ cIN (`1`) joins at pseudotime 0 in this smoke tree. So the first URD tree does support some separation of neuronal populations, but the cleanest split is PV precursor plus MGE subpallial versus SST+ cIN rather than three independent terminal branches.

URD warnings to remember:

```text
No obvious breakpoint between 1 and 3.
1144 cells were not visited by a branch that exists at their pseudotime and were not assigned.
```

These warnings do not mean the tree failed. They mean that at least one branch comparison is not sharply separable in this 5k smoke tree, and a substantial subset of cells did not receive a branch assignment under the current visit threshold.

Tip composition after correcting NA-safe tip counting:

| Tip ID | Annotation | Cells | Median pseudotime |
|---:|---|---:|---:|
| 1 | SST+ cIN | 1043 | 0.294 |
| 2 | PV neuron precursor | 759 | 0.354 |
| 3 | MGE subpallial neurons | 954 | 0.335 |

Top branch-specific genes from the smoke tree finalizer:

| Tip | Annotation | Top genes |
|---:|---|---|
| 1 | SST+ cIN | `PLS3`, `GALNT14`, `ERO1B`, `ASCL1`, `NFIB` |
| 2 | PV neuron precursor | `MEGF10`, `GBX1`, `GLCE`, `SLC2A13`, `ASIC1` |
| 3 | MGE subpallial neurons | `CDH23`, `SPOCK1`, `CDH13`, `ZFHX3`, `VAT1L` |

Branch-specific genes are computed after the tree exists by contrasting each requested tip against the other requested tips across URD variable genes. The finalizer reports mean logUPX, expression fraction, log2 fold change, expression-fraction delta, and a branch-specificity score.

## Current Full 30k URD Status

The original 30k job is still running as of the last check:

```text
JOBID     51670320
STATE     RUNNING
ELAPSED   5:39:25
NODE      gl3327
STAGE     calcDM(knn = 100, sigma = "local")
```

Do not stop it unless explicitly deciding to abandon the full diffusion-map run. The 5k smoke object is currently the completed object used for the tree proof of workflow.

## Jia VZ-RGC-Like Alternative Root

Goal: use the Jia hierarchy to identify the earliest VZ-like progenitor compartment inside paper/manual Radial glia, recompute URD flood pseudotime from that stricter root, keep the same neuronal tips, and compare tree topology against the broad Radial glia-root tree.

Script added:

```text
scripts/18_div30_urd_jia_vz_root_pseudotime.R
```

This script does not recompute PCA or the diffusion map. It reads an existing URD object, selects the Jia VZ-RGC-like root from metadata and expression, then reruns only `floodPseudotime()` and `floodPseudotimeProcess()` from the new root.

Root definition used on the 5k smoke object:

| Step | Rule | Threshold |
|---|---|---:|
| Start | `paper_cluster_annotation == "Radial glia"` | 2040 cells |
| High RGC1 | `jia_score_RGC1 >=` Radial glia 75th percentile | 0.2815764 |
| Low RGC2 | `jia_score_RGC2 <=` Radial glia 50th percentile | 0.3066419 |
| Low IPC | `jia_score_IPC <=` Radial glia 50th percentile | 0.2369850 |
| Jia hierarchy candidates | all three Jia rules above | 137 cells |
| Proliferative refinement | top half by proliferation marker score among Jia candidates | threshold 0.338125 |
| Final VZ-RGC-like root | Jia hierarchy candidate plus proliferative refinement | 69 cells |

Proliferation marker genes found and used:

```text
MKI67, TOP2A, PCNA, MCM2, MCM3, MCM4, MCM5, MCM6, MCM7,
STMN1, HMGB2, UBE2C, AURKB, CDK1, CCNB1, CCNB2
```

Alternative-root output root:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_vz_root_v1/
```

Alternative-root outputs:

```text
div30_urd_jia_vz_root_object.rds
tables/jia_vz_root_summary.tsv
tables/jia_vz_root_thresholds.tsv
tables/jia_vz_root_cells.tsv
plots/jia_vz_root_umap.png
lineage_decision_report/
lineage_tree_sst_pv_mge_tips_v1/
```

Command run:

```text
Rscript scripts/18_div30_urd_jia_vz_root_pseudotime.R \
  --urd-rds /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/div30_first_urd_object.rds \
  --outdir /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_vz_root_v1 \
  --annotation-col paper_cluster_annotation \
  --radial-glia-label "Radial glia" \
  --pseudotime-name jia_vz_rgc_root \
  --rgc1-high-quantile 0.75 \
  --rgc2-low-quantile 0.50 \
  --ipc-low-quantile 0.50 \
  --seed 7 \
  --n-floods 20
```

The tree builder was updated to support an explicit logical root column:

```text
--root-col jia_vz_rgc_root_candidate
```

Alternative-root tree command:

```text
Rscript scripts/16_div30_urd_build_lineage_tree.R \
  --urd-rds /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_vz_root_v1/div30_urd_jia_vz_root_object.rds \
  --outdir /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_vz_root_v1/lineage_tree_sst_pv_mge_tips_v1 \
  --annotation-col paper_cluster_annotation \
  --pseudotime-name jia_vz_rgc_root \
  --root-label "Jia VZ-RGC-like" \
  --root-col jia_vz_rgc_root_candidate \
  --tip-labels "SST+ cIN,PV neuron precursor,MGE subpallial neurons" \
  --n-per-tip 5000 \
  --seed 7
```

Alternative-root tree output root:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_vz_root_v1/lineage_tree_sst_pv_mge_tips_v1/
```

Tree topology comparison:

| Root | First join | First join pseudotime | Second join | Second join pseudotime | Unassigned warning |
|---|---|---:|---|---:|---:|
| Broad Radial glia | PV neuron precursor + MGE subpallial neurons | 0.1626288 | SST+ cIN joins PV/MGE branch | 0.0000000 | 1144 cells |
| Jia VZ-RGC-like | PV neuron precursor + MGE subpallial neurons | 0.3072650 | SST+ cIN joins PV/MGE branch | 0.2228481 | 875 cells |

Interpretation: the Jia VZ-RGC-like root does not change the major topology in the 5k smoke tree. Both roots group PV neuron precursor with MGE subpallial neurons first, with SST+ cIN joining as the other arm. The Jia VZ-root is cleaner in two practical ways: branchpoints are no longer collapsed to pseudotime 0, and fewer cells are unassigned by branch visitation (`875` vs `1144`). This suggests the stricter Jia/proliferative VZ-RGC-like root may be a stronger root definition for URD tree reconstruction than broad Radial glia, while preserving the same biological branch relationship among the tips.

Alternative-root pseudotime ordering by paper/manual annotation:

| Annotation | Median pseudotime |
|---|---:|
| Radial glia | 0.291 |
| Inhibitory progenitors | 0.346 |
| SST+ cIN | 0.455 |
| MGE subpallial neurons | 0.467 |
| PV neuron precursor | 0.491 |

Top VZ-root pseudotime genes:

| Direction | Genes |
|---|---|
| Negative/root-associated | `FZD5`, `SPARC`, `CELSR1`, `FKBP10`, `GPC3`, `LRP2`, `NOTCH1`, `SPON1`, `TEAD2`, `ADGRV1` |
| Positive/late-associated | `LRRC7`, `GAD1`, `CELF3`, `L1CAM`, `FNBP1L`, `RUNX1T1`, `AKAP6`, `NRXN3`, `GPC2`, `NSG2` |

## Jia RootScore Root, No Cluster-Based Selection

Goal: identify a Jia-like proliferative VZ-RGC root population without using Seurat or paper/manual clusters for selection.

Scripts added:

```text
scripts/19_div30_jia_rootscore_candidates.R
scripts/20_div30_urd_reflood_from_root_column.R
```

RootScore formula:

```text
RootScore = z(jia_score_RGC1) + z(proliferation_score) - z(jia_score_IPC)
```

The score was computed for every cell. Cluster annotations were not used to compute or select candidates. The only ingredients were `jia_score_RGC1`, `jia_score_IPC`, proliferation gene expression, UMAP coordinates for plotting, and existing URD pseudotime for positional diagnostics.

Proliferation genes used:

```text
MKI67, TOP2A, CENPF, UBE2C, PCNA, MCM2, MCM3, MCM4, MCM5,
MCM6, MCM7, STMN1, HMGB2, AURKB, CDK1, CCNB1, CCNB2
```

Diagnostic marker genes summarized in candidate roots:

```text
HES1, FGFR2, NOTCH1, NOTCH2, VIM, NES, DACH1, DLX1, DLX2, ASCL1
```

RootScore output root:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1/
```

Command run:

```text
Rscript scripts/19_div30_jia_rootscore_candidates.R \
  --urd-rds /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/div30_first_urd_object.rds \
  --outdir /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1 \
  --pseudotime-name paper_radial_glia_root \
  --candidate-pcts 1,2,5,10 \
  --selected-pct 2
```

RootScore outputs:

```text
div30_urd_jia_rootscore_object.rds
jia_rootscore_candidate_report.md
tables/root_score_all_cells.tsv
tables/root_score_distribution.tsv
tables/root_score_candidate_counts.tsv
tables/root_score_candidate_cells_top1pct.tsv
tables/root_score_candidate_cells_top2pct.tsv
tables/root_score_candidate_cells_top5pct.tsv
tables/root_score_candidate_cells_top10pct.tsv
tables/root_score_marker_expression_by_candidate_set.tsv
tables/root_score_pseudotime_by_candidate_set.tsv
tables/root_score_proliferation_genes_used.tsv
tables/root_score_marker_genes_used.tsv
plots/root_score_distribution.png
plots/root_score_umap.png
plots/root_score_candidate_umap.png
plots/root_score_by_pseudotime.png
```

Candidate root set summary:

| Candidate set | Cells | Min RootScore | Median RootScore | Median RGC1 | Median IPC | Median proliferation | Median broad-root URD pseudotime |
|---|---:|---:|---:|---:|---:|---:|---:|
| Top 1% | 50 | 4.1997 | 4.3529 | 0.2918 | 0.1769 | 1.4538 | 0 |
| Top 2% | 100 | 3.9951 | 4.1961 | 0.2882 | 0.1769 | 1.4100 | 0 |
| Top 5% | 250 | 3.6296 | 3.8936 | 0.2819 | 0.1653 | 1.1974 | 0 |
| Top 10% | 500 | 3.2034 | 3.6282 | 0.2756 | 0.1682 | 1.0629 | 0 |

RootScore distribution:

```text
min     -5.7560
median  -0.6383
mean     0.0000
sd       2.2720
95%      3.6270
98%      3.9806
99%      4.1926
max      5.0343
```

Marker behavior in candidate roots:

| Candidate set | VZ/RGC markers | IPC/neurogenic markers |
|---|---|---|
| Top 1% | High `FGFR2`, `NOTCH1`, `VIM`, `NES`, `HES1`; `NOTCH1` expressed in 100%, `FGFR2`/`NES` in 98% | Low `DLX1`, `DLX2`, `ASCL1`; medians are 0 |
| Top 2% | High `FGFR2`, `NOTCH1`, `VIM`, `NES`, `HES1`; `NOTCH1` expressed in 100%, `NES` in 99% | Low `DLX1`, `DLX2`, `ASCL1`; medians are 0 |
| Top 5% | Still VZ/RGC-high, but proliferation score is lower than top 1-2% | Low IPC/neurogenic markers persist |
| Top 10% | Broader and less stringent; still VZ/RGC-like but less enriched | Low IPC/neurogenic markers persist |

Selected root: top 2% by RootScore.

Rationale: top 1% is maximally stringent but only 50 cells. Top 2% keeps the same biological profile, remains strongly proliferative and low IPC-like, and gives 100 root cells, which is a more stable URD root size. Top 5% and 10% are useful sensitivity sets but broaden the compartment.

Selected root metadata columns written into `div30_urd_jia_rootscore_object.rds`:

```text
jia_rootscore
jia_rootscore_z_RGC1
jia_rootscore_z_IPC
jia_rootscore_z_proliferation
jia_rootscore_proliferation_score
jia_rootscore_top1pct_candidate
jia_rootscore_top2pct_candidate
jia_rootscore_top5pct_candidate
jia_rootscore_top10pct_candidate
jia_rootscore_selected_root
jia_rootscore_selected_top_percent
```

The selected top 2% RootScore root was then used as an actual alternative URD root by rerunning flood pseudotime from the metadata column:

```text
Rscript scripts/20_div30_urd_reflood_from_root_column.R \
  --urd-rds /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1/div30_urd_jia_rootscore_object.rds \
  --outdir /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1/reflood_top2pct_v1 \
  --root-col jia_rootscore_selected_root \
  --pseudotime-name jia_rootscore_top2pct_root \
  --seed 7 \
  --n-floods 20
```

Reflooded output root:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1/reflood_top2pct_v1/
```

Reflood summary:

```text
root_col        jia_rootscore_selected_root
pseudotime      jia_rootscore_top2pct_root
n_root_cells    100
n_floods        20
```

Top 2% RootScore-root pseudotime ordering:

| Annotation | Median pseudotime |
|---|---:|
| Radial glia | 0.236 |
| Inhibitory progenitors | 0.299 |
| SST+ cIN | 0.444 |
| MGE subpallial neurons | 0.474 |
| PV neuron precursor | 0.495 |

Root-associated genes under the top 2% RootScore root:

```text
SPARC, NOTCH1, FKBP10, ADGRV1, TTYH1, FGFR2, LRP2, TEAD2, CRB2, CELSR1, NOTCH3
```

Terminal-associated genes under the top 2% RootScore root:

```text
LRRC7, CELF3, GAD1, L1CAM, FNBP1L, AKAP6, GPC2, KALRN, RUNX1T1, NRXN3, NSG2
```

## Jia Program Trends Over Jia RootScore-Rooted Pseudotime

Goal: determine whether the Jia-rooted URD pseudotime supports a trunk ordering of `RGC1 -> RGC2 -> IPC`, `RGC1/RGC2 overlap -> IPC`, or another structure.

Script added:

```text
scripts/21_div30_jia_program_pseudotime_trends.R
```

Input object:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1/reflood_top2pct_v1/div30_urd_reflood_object.rds
```

Pseudotime used:

```text
jia_rootscore_top2pct_root
```

Command run:

```text
Rscript scripts/21_div30_jia_program_pseudotime_trends.R \
  --urd-rds /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1/reflood_top2pct_v1/div30_urd_reflood_object.rds \
  --outdir /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1/reflood_top2pct_v1/jia_program_trends \
  --pseudotime-name jia_rootscore_top2pct_root \
  --n-bins 30
```

Outputs:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1/reflood_top2pct_v1/jia_program_trends/
```

Files:

```text
jia_program_pseudotime_trends_report.md
plots/jia_program_scores_vs_pseudotime.png
plots/jia_program_binned_trends.png
tables/jia_program_scores_by_cell.tsv
tables/jia_program_binned_means.tsv
tables/jia_program_peak_pseudotimes.tsv
tables/jia_program_ordering_interpretation.tsv
```

Peak pseudotimes from binned mean scores:

| Program | Peak bin | Peak pseudotime | Peak mean score | Cells in peak bin |
|---|---:|---:|---:|---:|
| RGC1 | 1 | 0.000 | 0.2962 | 100 |
| RGC2 | 1 | 0.000 | 0.3571 | 100 |
| IPC | 20 | 0.4065 | 0.4516 | 221 |

Interpretation:

```text
RGC1/RGC2 overlap -> IPC
```

The Jia RootScore-rooted URD pseudotime does not support a clean sequential `RGC1 -> RGC2 -> IPC` progression in this 5k smoke object. Instead, RGC1 and RGC2 both peak in the root bin at pseudotime 0, then both decline as IPC rises and peaks later around pseudotime 0.406. This suggests the selected VZ-like root is RGC1/RGC2-high and IPC-low, with IPC emerging downstream along the developmental trunk.

## Jia RootScore Top-2% URD Tree Decision

Goal: keep the Jia RootScore top-2% root provisional until a full URD tree is reconstructed with the same tips used in the broad Radial glia-root tree.

Input object:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1/reflood_top2pct_v1/div30_urd_reflood_object.rds
```

Root:

```text
jia_rootscore_selected_root
```

Pseudotime:

```text
jia_rootscore_top2pct_root
```

Tips unchanged:

| Tip ID | Annotation |
|---:|---|
| 1 | SST+ cIN |
| 2 | PV neuron precursor |
| 3 | MGE subpallial neurons |

Tree command run:

```text
Rscript scripts/16_div30_urd_build_lineage_tree.R \
  --urd-rds /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1/reflood_top2pct_v1/div30_urd_reflood_object.rds \
  --outdir /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1/reflood_top2pct_v1/lineage_tree_sst_pv_mge_tips_v1 \
  --annotation-col paper_cluster_annotation \
  --pseudotime-name jia_rootscore_top2pct_root \
  --root-label "Jia RootScore top2pct" \
  --root-col jia_rootscore_selected_root \
  --tip-labels "SST+ cIN,PV neuron precursor,MGE subpallial neurons" \
  --n-per-tip 5000 \
  --seed 7
```

Finalizer command run:

```text
Rscript scripts/17_div30_urd_finalize_lineage_tree_report.R \
  --tree-rds /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1/reflood_top2pct_v1/lineage_tree_sst_pv_mge_tips_v1/div30_urd_lineage_tree_object.rds \
  --outdir /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1/reflood_top2pct_v1/lineage_tree_sst_pv_mge_tips_v1 \
  --annotation-col paper_cluster_annotation \
  --pseudotime-name jia_rootscore_top2pct_root \
  --tip-mapping /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1/reflood_top2pct_v1/lineage_tree_sst_pv_mge_tips_v1/tables/tree_tip_mapping.tsv \
  --top-n-branch-genes 50
```

Tree outputs:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1/reflood_top2pct_v1/lineage_tree_sst_pv_mge_tips_v1/
```

Files:

```text
div30_urd_lineage_tree_object.rds
div30_urd_tip_random_walks.rds
urd_lineage_tree_report.md
tables/tree_status.tsv
tables/tree_segment_joins.tsv
tables/tree_tip_composition.tsv
tables/branch_specific_genes.tsv
plots/pseudotime_logistic.png
plots/breakpoint_decisions-1.pdf
plots/breakpoint_decisions-3.pdf
plots/urd_tree_annotation.png
plots/urd_tree_pseudotime.png
```

Core URD tree calls completed:

```text
pseudotimeDetermineLogistic()
pseudotimeWeightTransitionMatrix()
simulateRandomWalksFromTips()
processRandomWalksFromTips()
buildTree()
```

Topology comparison:

| Root | First join | First join pseudotime | Second join | Second join pseudotime | Unassigned warning |
|---|---|---:|---|---:|---:|
| Broad Radial glia | PV neuron precursor + MGE subpallial neurons | 0.1626288 | SST+ cIN joins PV/MGE branch | 0.0000000 | 1144 cells |
| Jia hierarchy/proliferative VZ-RGC | PV neuron precursor + MGE subpallial neurons | 0.3072650 | SST+ cIN joins PV/MGE branch | 0.2228481 | 875 cells |
| Jia RootScore top 2% | PV neuron precursor + MGE subpallial neurons | 0.3714297 | SST+ cIN joins PV/MGE branch | 0.1716992 | 618 cells |

Branch-specific genes:

| Tip | Annotation | Top genes |
|---:|---|---|
| 1 | SST+ cIN | `PLS3`, `GALNT14`, `ERO1B`, `ASCL1`, `NFIB` |
| 2 | PV neuron precursor | `MEGF10`, `GBX1`, `GLCE`, `SLC2A13`, `ASIC1` |
| 3 | MGE subpallial neurons | `CDH23`, `SPOCK1`, `CDH13`, `ZFHX3`, `VAT1L` |

The branch-specific genes match the previous tip-level branch-gene summaries because the current finalizer computes each tip versus the other requested tips over variable genes. The important tree-specific comparison is the preserved segment topology plus improved branchpoint behavior.

Decision: finalize the Jia RootScore top-2% population as the production root for this DIV30 URD lineage workflow and stop further root optimization on the 5k smoke object.

Rationale: the selected RootScore top-2% cells have the expected proliferative VZ-RGC-like marker profile, reflood pseudotime is biologically ordered, Jia program trends show an RGC1/RGC2-high root followed by IPC emergence, and the full URD tree preserves the same major branch topology as the broad Radial glia-root tree while producing nonzero branchpoints and fewer unassigned cells.

## Full 30k Run Note

Slurm job `51670320` completed successfully according to `sacct`:

```text
State     COMPLETED
ExitCode  0:0
Elapsed   06:17:26
MaxRSS    55459644K
```

However, as of the last direct filesystem check, the expected full-run outputs were not present under:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_v1/
```

The log tail ends at:

```text
[2026-06-11 18:39:43] Processing flood pseudotime
```

There is no final `Saved URD object` line in the checked log. Treat this as a full-run output anomaly to investigate before applying the finalized production root to the 30k object.

## Branchpoint Decision Genes And Jia Program Localization

Goal: after finalizing the Jia RootScore top-2% root on the 5k smoke workflow, localize Jia progenitor programs on the tree before interpreting branch-decision genes biologically.

Important conceptual separation:

```text
Decision genes answer:
What molecular programs define branch divergence?

Jia scores answer:
Which progenitor programs occupy specific parts of the lineage tree?
```

Script added:

```text
scripts/22_div30_urd_branchpoint_decision_and_jia_localization.R
```

Input tree:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1/reflood_top2pct_v1/lineage_tree_sst_pv_mge_tips_v1/div30_urd_lineage_tree_object.rds
```

Command run:

```text
Rscript scripts/22_div30_urd_branchpoint_decision_and_jia_localization.R \
  --tree-rds /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1/reflood_top2pct_v1/lineage_tree_sst_pv_mge_tips_v1/div30_urd_lineage_tree_object.rds \
  --outdir /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1/reflood_top2pct_v1/branchpoint_decision_jia_localization_v1 \
  --pseudotime-name jia_rootscore_top2pct_root \
  --window 0.06 \
  --min-cells 40 \
  --top-n 50 \
  --genes variable
```

Output root:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1/reflood_top2pct_v1/branchpoint_decision_jia_localization_v1/
```

Outputs:

```text
branchpoint_decision_and_jia_localization_report.md
tables/jia_program_by_tree_segment.tsv
tables/tree_segment_cell_assignments.tsv
tables/branchpoint_cell_selection_summary.tsv
tables/branchpoint1_upstream_cells.tsv
tables/branchpoint1_downstream_cells.tsv
tables/branchpoint1_decision_genes.tsv
tables/branchpoint2_upstream_cells.tsv
tables/branchpoint2_downstream_cells.tsv
tables/branchpoint2_decision_genes.tsv
plots/jia_program_by_tree_segment_heatmap.png
plots/branchpoint1_decision_gene_heatmap.png
plots/branchpoint2_decision_gene_heatmap.png
```

### Jia Program Localization By URD Segment

Segment meaning:

| Segment | Meaning |
|---:|---|
| 5 | root/trunk before SST vs PV/MGE split |
| 4 | shared PV/MGE trunk before PV vs MGE split |
| 1 | SST branch |
| 2 | PV branch |
| 3 | MGE branch |

Jia mean scores by segment:

| Segment | Cells | Median pseudotime | Mean RGC1 | Mean RGC2 | Mean IPC |
|---:|---:|---:|---:|---:|---:|
| 5 | 330 | 0.137 | 0.263 | 0.326 | 0.262 |
| 4 | 678 | 0.249 | 0.220 | 0.256 | 0.256 |
| 1 | 1775 | 0.365 | 0.119 | 0.130 | 0.417 |
| 2 | 720 | 0.497 | 0.097 | 0.037 | 0.345 |
| 3 | 879 | 0.479 | 0.124 | 0.054 | 0.364 |

Interpretation: the trunk is not a pure RGC compartment. Segment 5 is the strongest RGC-like segment, especially RGC2-high. Segment 4 is a mixed progenitor/trunk segment where RGC1/RGC2 decline and IPC begins to rise. The SST branch has the strongest IPC signal among terminal segments. PV and MGE terminal segments retain lower RGC1/RGC2 and moderate IPC. This supports the caveat that the current progenitor trunk is a mixed Jia-program compartment.

### Branchpoint Cell Selection

The script selected cells in a pseudotime window of 0.06 immediately around each split. It expanded only if needed, but the initial 0.06 window met the minimum cell requirement for every region.

| Branchpoint | Region | Split pseudotime | Cells | Median selected pseudotime |
|---|---|---:|---:|---:|
| SST vs PV/MGE | upstream segment 5 | 0.1717 | 223 | 0.1515 |
| SST vs PV/MGE | downstream SST segment 1 | 0.1717 | 340 | 0.2081 |
| SST vs PV/MGE | downstream PV/MGE segment 4 | 0.1717 | 254 | 0.2096 |
| PV vs MGE | upstream segment 4 | 0.3714 | 134 | 0.3501 |
| PV vs MGE | downstream PV segment 2 | 0.3714 | 83 | 0.4012 |
| PV vs MGE | downstream MGE segment 3 | 0.3714 | 207 | 0.4036 |

### Branchpoint 1: SST vs PV/MGE

DE comparison:

```text
SST_downstream vs PV_MGE_downstream
```

Output:

```text
tables/branchpoint1_decision_genes.tsv
plots/branchpoint1_decision_gene_heatmap.png
```

Top genes by adjusted p-value:

| Gene | logFC | adjusted p-value | Direction |
|---|---:|---:|---|
| HNRNPA0 | 0.738 | 5.25e-12 | SST higher |
| KPNA2 | 0.719 | 9.41e-12 | SST higher |
| MARCKSL1 | 0.766 | 1.37e-11 | SST higher |
| CKB | 0.616 | 2.54e-11 | SST higher |
| MARCKS | 0.687 | 1.68e-10 | SST higher |
| VIM | 0.973 | 1.80e-10 | SST higher |
| TUBB2B | 0.749 | 1.80e-10 | SST higher |
| NAALAD2 | -0.668 | 1.80e-10 | PV/MGE higher |
| H1FX | 0.734 | 2.89e-10 | SST higher |
| ASCL1 | 0.670 | 2.89e-10 | SST higher |
| PRTG | -0.874 | 1.89e-09 | PV/MGE higher |
| SOX2 | 0.665 | 1.17e-08 | SST higher |
| NNAT | 0.985 | 1.29e-08 | SST higher |
| CCND2 | 0.668 | 1.88e-08 | SST higher |
| HMGB2 | 0.659 | 2.85e-08 | SST higher |

Interpretation: the SST-side early branch is enriched for a proliferative/progenitor-like and neurogenic transition signature (`VIM`, `SOX2`, `ASCL1`, `CCND2`, `HMGB2`, `NNAT`) relative to the PV/MGE shared branch. This reinforces the working hypothesis that the first split may reflect a broader lineage-state divergence rather than a simple interneuron subtype decision.

### Branchpoint 2: PV vs MGE

DE comparison:

```text
PV_downstream vs MGE_downstream
```

Output:

```text
tables/branchpoint2_decision_genes.tsv
plots/branchpoint2_decision_gene_heatmap.png
```

Top genes by adjusted p-value:

| Gene | logFC | adjusted p-value | Direction |
|---|---:|---:|---|
| HDAC9 | -1.564 | 2.21e-19 | MGE higher |
| ISLR2 | -1.507 | 7.20e-17 | MGE higher |
| GAD1 | 1.761 | 2.42e-14 | PV higher |
| STMN2 | -1.168 | 2.42e-14 | MGE higher |
| GAP43 | -1.122 | 2.56e-14 | MGE higher |
| PCDH17 | -1.352 | 8.40e-14 | MGE higher |
| GLCE | 1.202 | 8.40e-14 | PV higher |
| NRXN1 | -1.180 | 2.34e-13 | MGE higher |
| ZFHX3 | -1.347 | 4.02e-13 | MGE higher |
| ASIC1 | 1.141 | 7.67e-12 | PV higher |
| DNER | -1.014 | 8.23e-12 | MGE higher |
| SPOCK1 | -1.259 | 1.04e-11 | MGE higher |
| NRXN3 | 1.162 | 1.42e-11 | PV higher |
| TLE4 | -0.877 | 1.44e-11 | MGE higher |
| RGMB | -1.311 | 2.51e-11 | MGE higher |

Interpretation: the PV vs MGE split has a sharper neuronal/axon-guidance and subtype-separation signal. PV-side cells are higher for `GAD1`, `GLCE`, `ASIC1`, and `NRXN3`; MGE-side cells are higher for `HDAC9`, `ISLR2`, `STMN2`, `GAP43`, `PCDH17`, `ZFHX3`, `SPOCK1`, and `TLE4`.

Working model after this analysis: Jia progenitor programs occupy the trunk in a mixed state, with RGC-like programs highest before the first split and IPC strongest along the SST branch. The SST vs PV/MGE split may reflect a broader progenitor-state or pallial/subpallial divergence, while the PV vs MGE split shows clearer downstream branch-specific neuronal programs.

## Jia Lineage-Focused Branchpoint Figure And Interpretation

Date: 2026-06-11

Script:

```text
scripts/23_div30_urd_jia_lineage_branchpoint_figure.R
```

Purpose: generate a publication-oriented branchpoint figure and reinterpret the URD branches using Jia lineage programs. This step intentionally does not center the biology on whether the tree recovers a clean RGC1 -> RGC2 -> IPC sequence. The main question is which Jia inhibitory neuron lineage programs are supported by the genes enriched at each URD branchpoint.

Input tree:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1/reflood_top2pct_v1/lineage_tree_sst_pv_mge_tips_v1/div30_urd_lineage_tree_object.rds
```

Output directory:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1/reflood_top2pct_v1/jia_lineage_branchpoint_interpretation_v1/
```

Run command:

```bash
module load R/4.4.3
Rscript scripts/23_div30_urd_jia_lineage_branchpoint_figure.R \
  --tree-rds /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1/reflood_top2pct_v1/lineage_tree_sst_pv_mge_tips_v1/div30_urd_lineage_tree_object.rds \
  --outdir /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1/reflood_top2pct_v1/jia_lineage_branchpoint_interpretation_v1 \
  --pseudotime-name jia_rootscore_top2pct_root \
  --window 0.06 \
  --min-cells 40 \
  --top-ns 50,100,250 \
  --genes variable
```

### Programmatic logic

1. Load the saved URD tree object.
2. Use URD `cells.in.segment` and `segment.joins` to define branchpoint-local windows.
3. For branchpoint 1, compare cells immediately downstream of the split:
   - SST branch: segment 1
   - PV/MGE branch: segment 4
4. For branchpoint 2, compare cells immediately downstream of the split:
   - PV branch: segment 2
   - MGE branch: segment 3
5. Run Wilcoxon DE over the 3,000 URD variable genes in each branchpoint comparison.
6. Export all DE genes and branch-enriched top 50, top 100, and top 250 lists for each daughter branch.
7. Compute hypergeometric overlap between each branch-enriched gene list and Jia lineage marker/TF programs.
8. Compute per-cell Jia lineage program scores from mean logUPX expression of available lineage genes.
9. Generate the publication panel and supporting lineage heatmaps/tree score plots.

### Main outputs

Publication figure:

```text
plots/urd_branchpoint_jia_publication_panel.png
plots/urd_branchpoint_jia_publication_panel.pdf
```

Individual panels:

```text
plots/panel_a_urd_tree_branchpoints.png
plots/panel_b_branchpoint1_top20_heatmap.png
plots/panel_c_branchpoint2_top20_heatmap.png
plots/panel_d_jia_rgc_ipc_segment_heatmap.png
```

Jia lineage support plots:

```text
plots/jia_lineage_program_scores_on_tree.png
plots/jia_lineage_marker_segment_heatmap.png
```

Decision gene and interpretation tables:

```text
tables/branchpoint1_all_de_sst_vs_pv_mge.tsv
tables/branchpoint2_all_de_pv_vs_mge.tsv
tables/sst_vs_pv_mge_sst_enriched_top50.tsv
tables/sst_vs_pv_mge_sst_enriched_top100.tsv
tables/sst_vs_pv_mge_sst_enriched_top250.tsv
tables/sst_vs_pv_mge_pv_mge_enriched_top50.tsv
tables/sst_vs_pv_mge_pv_mge_enriched_top100.tsv
tables/sst_vs_pv_mge_pv_mge_enriched_top250.tsv
tables/pv_vs_mge_pv_enriched_top50.tsv
tables/pv_vs_mge_pv_enriched_top100.tsv
tables/pv_vs_mge_pv_enriched_top250.tsv
tables/pv_vs_mge_mge_enriched_top50.tsv
tables/pv_vs_mge_mge_enriched_top100.tsv
tables/pv_vs_mge_mge_enriched_top250.tsv
tables/jia_lineage_overlap_enrichment.tsv
tables/branchpoint_jia_lineage_interpretation.tsv
tables/jia_lineage_marker_inventory.tsv
tables/jia_lineage_marker_expression_by_segment.tsv
tables/jia_lineage_marker_expression_by_branchpoint_region.tsv
```

### Branchpoint windows used

| Branchpoint | Region | Split pseudotime | n cells | Median pseudotime |
|---|---|---:|---:|---:|
| SST vs PV/MGE | upstream segment 5 | 0.1717 | 223 | 0.1515 |
| SST vs PV/MGE | downstream SST segment 1 | 0.1717 | 340 | 0.2081 |
| SST vs PV/MGE | downstream PV/MGE segment 4 | 0.1717 | 254 | 0.2096 |
| PV vs MGE | upstream segment 4 | 0.3714 | 134 | 0.3501 |
| PV vs MGE | downstream PV segment 2 | 0.3714 | 83 | 0.4012 |
| PV vs MGE | downstream MGE segment 3 | 0.3714 | 207 | 0.4036 |

### Jia RGC/RGC2/IPC localization

| Segment | Interpretation | n cells | mean RGC1 | mean RGC2 | mean IPC |
|---|---|---:|---:|---:|---:|
| 5 | root/trunk | 330 | 0.263 | 0.326 | 0.262 |
| 4 | PV/MGE trunk | 678 | 0.220 | 0.256 | 0.256 |
| 1 | SST | 1775 | 0.119 | 0.130 | 0.417 |
| 2 | PV | 720 | 0.097 | 0.037 | 0.345 |
| 3 | MGE | 879 | 0.124 | 0.054 | 0.364 |

Interpretation: RGC-like scores are strongest in the root/trunk and PV/MGE trunk, while IPC score is strongest on the SST branch. This supports treating the trunk as a mixed progenitor compartment rather than a pure RGC state.

### Jia lineage overlap interpretation

The corrected interpretation table reports no forced lineage call when a branch has zero Jia marker overlap.

| Branchpoint | Branch | Top n | Best Jia lineage | Supporting genes | adjusted p-value |
|---|---|---:|---|---|---:|
| SST vs PV/MGE | SST | 50 | EPHA5/MEF2C | ASCL1, DLX2 | 0.250 |
| SST vs PV/MGE | SST | 100 | EPHA5/MEF2C | ASCL1, DLX2, DLX1 | 0.176 |
| SST vs PV/MGE | SST | 250 | EPHA5/MEF2C | ASCL1, DLX2, DLX1, ARX, HMGA2, PBX3 | 0.009 |
| SST vs PV/MGE | PV/MGE | 50/100/250 | no Jia marker overlap | none | NA |
| PV vs MGE | PV | 50 | EPHA5/MEF2C | BCL11A, EPHA5 | 0.250 |
| PV vs MGE | PV | 100 | LHX6/NFIA | LHX6, TCF4 | 0.333 |
| PV vs MGE | PV | 250 | LHX6/NFIA | LHX6, TCF4, ARX, NFIB | 0.151 |
| PV vs MGE | MGE | 50 | CRABP1/ANGPT2 | SPOCK1 | 0.536 |
| PV vs MGE | MGE | 100/250 | LHX8/ISL1 | ID4 | 0.536 to 0.799 |

Working interpretation: the strongest Jia-lineage support in this smoke-scale tree is the SST-side enrichment for the EPHA5/MEF2C SVZ-derived cortex-bound inhibitory neuron program, especially in the top-250 branch-enriched list. PV has suggestive LHX6/NFIA support at broader top-n thresholds. MGE shows only weak single-gene support from SPOCK1 or ID4. The PV/MGE side at branchpoint 1 has no Jia marker overlap in the exported top lists, so it should not be assigned a Jia lineage from this result alone.

Important caveat: several Jia markers are present in the expression matrix but absent from the 3,000-variable-gene DE universe. See `tables/jia_lineage_marker_inventory.tsv` before overinterpreting absent overlaps.

## Refined Figure Sets: Jia Localization, Decision Genes, Marker Validation

Date: 2026-06-11

Script:

```text
scripts/24_div30_urd_refined_branchpoint_figures.R
```

Purpose: refine the figure set around the current biological model rather than doing more root optimization. The root-optimization phase is considered sufficient for the smoke-scale tree because the Jia RootScore top-2% root and broad Radial glia root produced nearly identical topologies. This stage prioritizes figure refinement and interpretation:

```text
RGC1/RGC2-enriched progenitor trunk
to IPC-enriched transition/branches
to neuronal lineage divergence
```

Output directory:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1/reflood_top2pct_v1/refined_branchpoint_figure_sets_v1/
```

Run command:

```bash
module load R/4.4.3
Rscript scripts/24_div30_urd_refined_branchpoint_figures.R \
  --tree-rds /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1/reflood_top2pct_v1/lineage_tree_sst_pv_mge_tips_v1/div30_urd_lineage_tree_object.rds \
  --outdir /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1/reflood_top2pct_v1/refined_branchpoint_figure_sets_v1 \
  --pseudotime-name jia_rootscore_top2pct_root \
  --window 0.06 \
  --min-cells 40 \
  --genes variable
```

### Figure Set 1: Jia lineage localization by segment

This replaces continuous mean-logUPX tree overlays as the primary Jia lineage localization display. For each URD segment, the script calculates mean Jia lineage scores for:

```text
LHX8_ISL1
NR2F1_NR2F2
EPHA5_MEF2C
LHX6_NFIA
CRABP1_ANGPT2
```

Outputs:

```text
plots/figset1_segment_lineage_score_heatmap.png
plots/figset1_segment_lineage_zscore_heatmap.png
plots/figset1_segment_lineage_enrichment_summary.png
tables/jia_lineage_scores_by_cell.tsv
tables/jia_lineage_scores_by_segment.tsv
tables/jia_lineage_segment_enrichment_summary.tsv
```

Top positive segment-lineage enrichments:

| Segment | Lineage program | Delta mean | adjusted p-value |
|---|---|---:|---:|
| SST | LHX6/NFIA | 0.346 | 2.97e-245 |
| SST | EPHA5/MEF2C | 0.471 | 1.54e-187 |
| MGE | LHX8/ISL1 | 0.278 | 6.97e-110 |
| SST | NR2F1/NR2F2 | 0.215 | 3.16e-90 |
| MGE | CRABP1/ANGPT2 | 0.175 | 1.60e-54 |
| root/trunk | NR2F1/NR2F2 | 0.282 | 1.71e-45 |
| PV/MGE trunk | NR2F1/NR2F2 | 0.174 | 8.83e-37 |

Interpretation: segment-level whole-program scoring shows that the SST branch carries strong EPHA5/MEF2C signal but also high LHX6/NFIA signal. MGE is the clearest segment-level match for LHX8/ISL1 and CRABP1/ANGPT2. NR2F1/NR2F2 is strongest in the root/trunk and PV/MGE trunk, consistent with progenitor/trunk localization. Because full lineage programs contain broad TFs, these segment-level scores should be read as localization patterns, not direct branch-decision genes.

### Figure Set 2: branchpoint decision genes

Decision genes are kept separate from Jia lineage scores. The script reruns branchpoint-local DE using the same branchpoint windows:

```text
Branchpoint 1: SST vs PV/MGE
Branchpoint 2: PV vs MGE
```

Outputs:

```text
plots/figset2_branchpoint1_volcano.png
plots/figset2_branchpoint2_volcano.png
plots/figset2_branchpoint1_top20_heatmap.png
plots/figset2_branchpoint1_top50_heatmap.png
plots/figset2_branchpoint2_top20_heatmap.png
plots/figset2_branchpoint2_top50_heatmap.png
tables/branchpoint1_sst_vs_pv_mge_ranked_marker_table.tsv
tables/branchpoint2_pv_vs_mge_ranked_marker_table.tsv
tables/branchpoint1_top20_decision_genes.tsv
tables/branchpoint1_top50_decision_genes.tsv
tables/branchpoint2_top20_decision_genes.tsv
tables/branchpoint2_top50_decision_genes.tsv
tables/branchpoint1_sst_enriched_top20.tsv
tables/branchpoint1_pv_mge_enriched_top20.tsv
tables/branchpoint2_pv_enriched_top20.tsv
tables/branchpoint2_mge_enriched_top20.tsv
```

Ranking note: the ranked marker tables are ordered by adjusted p-value first and effect size second. A combined `rank_score = abs(logFC) * -log10(adjusted_p_value)` is included for later prioritization if the analysis needs effect-size-weighted ranking.

### Figure Set 3: Jia marker validation

Rather than plotting entire lineage programs on the tree, this figure set validates branch identity using Jia-derived markers that appeared in branch DE or marker-validation logic.

Candidate marker groups:

```text
PV: LHX6, ARX, NFIB, TCF4
SST: EPHA5, ASCL1, DLX1, DLX2, HMGA2
MGE: SPOCK1, LHX8, ISL1, ID4
```

Outputs:

```text
plots/figset3_jia_marker_tree_overlays.png
plots/figset3_jia_marker_dotplot_by_segment.png
plots/figset3_jia_marker_heatmap_by_segment.png
tables/jia_marker_validation_by_segment.tsv
tables/jia_marker_validation_inventory.tsv
```

Marker-validation readout:

- SST branch shows strong `ASCL1`, `DLX1`, `DLX2`, and `EPHA5` localization, supporting the EPHA5/MEF2C-like interpretation.
- PV branch shows `LHX6` and `EPHA5` signal, but several broad TFs are not PV-specific at the segment level.
- MGE branch shows strong `SPOCK1`, `ISL1`, and `LHX8`, supporting weaker but visible CRABP1/ANGPT2 and LHX8/ISL1-related identity.
- `ID4` is strongest in the root/trunk and PV/MGE trunk, so it should be treated as progenitor/trunk-associated rather than a clean MGE terminal marker in this tree.

### Summary panel

Combined refined summary:

```text
plots/refined_urd_branchpoint_summary_panel.png
plots/refined_urd_branchpoint_summary_panel.pdf
```

Current prioritized conclusion:

```text
RGC1/RGC2-enriched progenitor trunk
to IPC-enriched transition/branches
to neuronal lineage divergence
```

Likely branch assignments, with current evidence:

| Branch | Interpretation | Evidence strength |
|---|---|---|
| SST | EPHA5/MEF2C-like lineage | strongest marker/DE support |
| PV | LHX6/NFIA-like decision signal | suggestive, strongest from branchpoint DE overlap rather than whole-program segment score |
| MGE | CRABP1/ANGPT2 and LHX8/ISL1-associated | weaker, supported by `SPOCK1`, `LHX8`, `ISL1`, and segment-level lineage scores |

Next natural step: scale this exact workflow from the 5k smoke object to the 30k/full dataset once the full URD tree object is available.

## Jia Fig. S11-Style Marker Validation

Date: 2026-06-11

Script:

```text
scripts/25_div30_urd_jia_fig_s11_marker_validation.R
```

SLURM workflow template:

```text
slurm_templates/33_div30_urd_jia_fig_s11_marker_validation.sbatch.template
```

Workflow placement: run after the URD lineage tree has been built and finalized. This is a downstream validation/figure-generation step in the DIV30 URD workflow, not an exploratory root-optimization or branchpoint-DE step.

Purpose: reproduce the structure of Jia Fig. S11 as closely as possible on the organoid URD tree. This is a marker-expression validation figure only.

Important constraints followed:

- No branchpoint DE genes.
- No lineage program averages.
- No marker substitutions.
- No z-score scaling across genes.
- Same URD tree coordinates for every marker panel.
- Exact marker order requested.

Input tree:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1/reflood_top2pct_v1/lineage_tree_sst_pv_mge_tips_v1/div30_urd_lineage_tree_object.rds
```

Output directory:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1/reflood_top2pct_v1/jia_fig_s11_style_marker_validation_v1/
```

Run command:

```bash
module load R/4.4.3
Rscript scripts/25_div30_urd_jia_fig_s11_marker_validation.R \
  --tree-rds /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1/reflood_top2pct_v1/lineage_tree_sst_pv_mge_tips_v1/div30_urd_lineage_tree_object.rds \
  --outdir /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_first_urd/div30_first_urd_paper_radial_glia_smoke5k_knn100_v1/jia_rootscore_root_v1/reflood_top2pct_v1/jia_fig_s11_style_marker_validation_v1 \
  --point-size 0.28
```

SLURM submit pattern:

```bash
cp slurm_templates/33_div30_urd_jia_fig_s11_marker_validation.sbatch.template \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/33_div30_urd_jia_fig_s11_marker_validation.sbatch
sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/33_div30_urd_jia_fig_s11_marker_validation.sbatch
```

Expression source: URD `logupx.data`. In the saved tree object, `logupx.data` had missing dimnames after serialization, so the script repairs row/column names from `count.data` when dimensions match. The plotted values are still `logupx.data` expression values, not z-scored values.

### Panel A: developmental progression markers

Exact order:

```text
HES1 | CACNA1E | DLX2 | DCX
```

Interpretation order:

```text
VZ RGC -> SVZ RGC -> IPC -> newborn neuron
```

Output:

```text
plots/jia_fig_s11_panel_a_developmental_markers.png
plots/jia_fig_s11_panel_a_developmental_markers.pdf
```

### Panel B: Jia inhibitory neuron lineage markers

Exact order:

```text
LHX8 | NR2F1 | EPHA5 | MEF2C | CRABP1
```

Panel B lineage labels now follow the Jia lineage interpretation scheme:

```text
VZ-RGC
to SVZ-RGC
to Neurogenesis
to Jia inhibitory neuron lineage branches
```

Displayed marker and lineage labels:

| Plotted marker | Lineage label | Terminal identity |
|---|---|---|
| LHX8 | LHX8/ISL1 | Subpallial cholinergic interneuron |
| NR2F1 | NR2F1/NR2F2 | Subpallial GABAergic interneuron |
| EPHA5 | EPHA5/MEF2C | Cortical GABAergic interneuron |
| MEF2C | LHX6/NFIA | Cortical interneuron lineage |
| CRABP1 | CRABP1/ANGPT2 | Subpallial GABAergic interneuron |

Output:

```text
plots/jia_fig_s11_panel_b_lineage_markers.png
plots/jia_fig_s11_panel_b_lineage_markers.pdf
```

### Combined output

```text
plots/jia_fig_s11_style_urd_marker_validation.png
plots/jia_fig_s11_style_urd_marker_validation.pdf
```

Tables:

```text
tables/jia_fig_s11_marker_order.tsv
tables/jia_fig_s11_marker_expression_summary.tsv
tables/jia_fig_s11_marker_expression_by_cell.tsv.gz
```

All requested markers were present.

Marker-expression summary:

| Gene | Category | pct expressed | mean logUPX | max logUPX |
|---|---|---:|---:|---:|
| HES1 | VZ RGC | 0.176 | 0.289 | 4.62 |
| CACNA1E | SVZ RGC | 0.385 | 0.649 | 4.38 |
| DLX2 | IPC | 0.663 | 1.619 | 6.13 |
| DCX | Newborn neuron | 0.906 | 2.751 | 6.39 |
| LHX8 | Jia lineage marker | 0.645 | 1.860 | 5.47 |
| NR2F1 | Jia lineage marker | 0.578 | 1.417 | 6.13 |
| EPHA5 | Jia lineage marker | 0.319 | 0.698 | 5.42 |
| MEF2C | Jia lineage marker | 0.164 | 0.208 | 3.46 |
| CRABP1 | Jia lineage marker | 0.456 | 0.945 | 5.51 |

Interpretation note: this figure should be interpreted as a Jia Fig. S11-style marker validation overlay on the URD tree. Branchpoint DE remains a separate analysis and should not be mixed into this figure.

## Production Status And Required Figure Checklist

Date: 2026-06-11

Current status: the DIV30 work above is a mature smoke-scale workflow, not the final production lineage run. The current tree is based on the 5k DIV30 smoke object with the Jia RootScore top-2% root. It is useful for establishing the workflow, figure set, and biological interpretation, but it should not be treated as the final manuscript-scale analysis until the same workflow is run on the larger/final DIV30 object and then mirrored on DIV90.

Required interpretation model to preserve across future runs:

```text
RGC1/RGC2-enriched progenitor trunk
to IPC-enriched transition/branches
to neuronal lineage divergence
```

The goal is not to force a clean RGC1 -> RGC2 -> IPC sequence. The priority is to test whether progenitor/trunk state, IPC transition, and Jia inhibitory lineage marker patterns are coherent on each URD tree.

### Must-have figures for every URD lineage run

These outputs should be generated for DIV30 full-scale and DIV90 smoke/full runs before a handoff is considered complete.

1. URD smoke/quality panel
   - UMAP pseudotime
   - Diffusion map pseudotime
   - Diffusion map annotation
   - Flood stability
   - Top pseudotime gene cascade
   - Top positive pseudotime genes table
   - Top negative pseudotime genes table

2. URD tree topology panel
   - URD tree with branchpoints labeled
   - Segment joins and branch pseudotimes table
   - Tip/root composition table
   - Tree colored by annotation
   - Tree colored by pseudotime

3. Branchpoint decision-gene panel
   - Branchpoint 1 volcano plot
   - Branchpoint 1 top 20 heatmap
   - Branchpoint 1 top 50 heatmap
   - Branchpoint 1 ranked marker table
   - Branchpoint 2 volcano plot
   - Branchpoint 2 top 20 heatmap
   - Branchpoint 2 top 50 heatmap
   - Branchpoint 2 ranked marker table
   - Top 20/50 branch-enriched tables per daughter branch

4. Jia lineage localization by segment
   - Segment x Jia lineage mean-score heatmap
   - Segment x Jia lineage z-score heatmap
   - Segment x Jia lineage enrichment summary
   - Segment-level table for lineage scores

5. Jia marker validation by segment
   - Jia marker tree overlays
   - Jia marker dotplot by segment
   - Jia marker heatmap by segment
   - Marker inventory table

6. Must-have Jia Fig. S11-style marker-expression validation
   - This is now a required validation figure for every production/smoke URD run.
   - Exact Panel A order: `HES1 | CACNA1E | DLX2 | DCX`
   - Exact Panel B order: `LHX8 | NR2F1 | EPHA5 | MEF2C | CRABP1`
   - Panel B labels must include:
     - `LHX8/ISL1`: Subpallial cholinergic interneuron
     - `NR2F1/NR2F2`: Subpallial GABAergic interneuron
     - `EPHA5/MEF2C`: Cortical GABAergic interneuron
     - `LHX6/NFIA`: Cortical interneuron lineage
     - `CRABP1/ANGPT2`: Subpallial GABAergic interneuron
   - Expression only; no z-score scaling across genes.
   - No marker substitutions.
   - No branchpoint DE genes.
   - No lineage program averages.

### Workflow order to preserve

1. Confirm source object, UMAP, metadata labels, raw/log expression slots.
2. Export/create URD input object with annotations carried into `@meta`.
3. Run root/tip smoke test.
4. Validate pseudotime and diffusion geometry.
5. Build URD tree.
6. Generate branchpoint decision-gene figures.
7. Generate Jia lineage segment-localization figures.
8. Generate Jia marker-validation figures.
9. Generate Jia Fig. S11-style marker-expression figure.
10. Update handoff with exact source paths, parameter values, root/tip definitions, and every figure/table output.

## DIV90 Next Workflow Target

DIV90 has not yet been run through this URD lineage workflow. Before shipping the major DIV30/DIV90 handoffs, DIV90 needs its own smoke-test pass, root/tip selection, tree generation, and required figure set.

### DIV90 source objects found

Likely source Seurat object with real labels:

```text
/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/umap_props_output/clustered_day90_with_cluster_names_2.rds
```

Cached AnnData:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/varela_div90.h5ad
```

Existing Seurat/AnnData inventory:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/seurat_anndata_umap_inventory/seurat_anndata_umap_inventory_v1/tables/seurat_object_inventory/varela_div90/
```

Existing DIV90 UMAP plots:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/seurat_anndata_umap_inventory/seurat_anndata_umap_inventory_v1/plots/umaps/varela_div90/
```

Most important existing UMAP for label inspection:

```text
umap_X_umap_seurat_cluster_number_name.png
```

### DIV90 UMAP cluster-label audit

Date: 2026-06-11

Script:

```text
scripts/26_div90_umap_cluster_label_audit.R
```

Purpose: upstream metadata audit before DIV90 URD. This is the figure to use when deciding how the Seurat cluster numbers map to exact metadata names and paper-level biology groups. It is not a URD output and should be kept separate from lineage reconstruction.

Inputs:

```text
obs:  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/varela_div90_srznf5gg/obs.tsv
umap: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/varela_div90_srznf5gg/umap.tsv
```

Output directory:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_umap_cluster_label_audit/div90_umap_cluster_label_audit_v1/
```

Run command:

```bash
module load R/4.4.3
Rscript scripts/26_div90_umap_cluster_label_audit.R \
  --obs /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/varela_div90_srznf5gg/obs.tsv \
  --umap /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/varela_div90_srznf5gg/umap.tsv \
  --outdir /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_umap_cluster_label_audit/div90_umap_cluster_label_audit_v1 \
  --label-col cluster_number_name \
  --cluster-col cluster_id
```

Outputs:

```text
plots/div90_umap_cluster_number_name_labeled.png
plots/div90_umap_cluster_number_name_labeled.pdf
plots/div90_umap_cluster_numbers_only.png
plots/div90_umap_cluster_numbers_only.pdf
plots/div90_umap_cluster_overlay_grid.png
plots/div90_umap_cluster_overlay_grid.pdf
tables/div90_cluster_number_name_to_biology_mapping.tsv
div90_umap_cluster_label_audit_report.md
```

This plot is the current best answer to: "show the DIV90 UMAP with the cluster number and exact name in metadata."

The overlay-grid plot is the current best answer to: "show each individual DIV90 cluster overlaid on the grey UMAP one by one in the same image."

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_umap_cluster_label_audit/div90_umap_cluster_label_audit_v1/plots/div90_umap_cluster_overlay_grid.png
```

Important reconciliation note: the metadata has 13 Seurat cluster IDs, but these collapse naturally toward 10 paper-level biology groups because three biology categories are split across two Seurat clusters:

- `4` and `10` are both Pre-Astrocytes/Astrocytes states.
- `5` and `8` are both LHX8+ vMGE GABAergic Striatal/GP-fated states.
- `6` and `7` are both Stressed Cells.

Therefore:

```text
13 Seurat cluster IDs
minus duplicated/split biology states for astrocytes, LHX8/vMGE, and stressed cells
= 10 broader biology groups
```

Exact DIV90 metadata-to-biology mapping for URD planning:

| Cluster | Exact metadata name | Biology name | URD smoke role guess | Paper-level reconciliation |
|---:|---|---|---|---|
| 0 | `0 - MGE Striatal/GP Fated` | MGE Striatal/GP Fated | candidate MGE/LHX8/subpallial tip | one Seurat cluster label |
| 1 | `1 - SST+, NPY +, Cortical Fated` | SST+, NPY +, Cortical Fated | candidate SST/cortical interneuron tip | one Seurat cluster label |
| 2 | `2 - CRABP1+/PV Precursors` | CRABP1+/PV Precursors | candidate CRABP1/PV-associated branch | one Seurat cluster label |
| 3 | `3 - PV precursors/Migrating cells/Cortical-fated` | PV precursors/Migrating cells/Cortical-fated | candidate PV lineage tip | one Seurat cluster label |
| 4 | `4 - Pre-Astrocytes/Astrocytes 1` | Pre-Astrocytes/Astrocytes 1 | exclude from neuronal URD: astrocyte lineage | collapse with cluster 10 for paper-level biology |
| 5 | `5 - LHX8+ vMGE GABergic Striatal/GP fated 1` | LHX8+ vMGE GABergic Striatal/GP fated 1 | candidate MGE/LHX8/subpallial tip | collapse with cluster 8 for paper-level biology |
| 6 | `6 - Stressed Cells` | Stressed Cells | exclude/track separately: stressed | collapse with cluster 7 for paper-level biology |
| 7 | `7 - Stressed Cells` | Stressed Cells | exclude/track separately: stressed | collapse with cluster 6 for paper-level biology |
| 8 | `8 - LHX8+ vMGE GABergic Striatal/GP fated 2` | LHX8+ vMGE GABergic Striatal/GP fated 2 | candidate MGE/LHX8/subpallial tip | collapse with cluster 5 for paper-level biology |
| 9 | `9 - Pre-OPCs/OPCs` | Pre-OPCs/OPCs | exclude from neuronal URD: OPC lineage | one Seurat cluster label |
| 10 | `10 - Pre-Astrocytes/Astrocytes 2` | Pre-Astrocytes/Astrocytes 2 | exclude from neuronal URD: astrocyte lineage | collapse with cluster 4 for paper-level biology |
| 11 | `11 - PV Precursors` | PV Precursors | candidate PV lineage tip | one Seurat cluster label |
| 12 | `12 - Dividing cells` | Dividing cells | candidate root/proliferative progenitor | one Seurat cluster label |

### Stiletti comparison placeholder

There is now an explicit upstream placeholder for a separate Stiletti/comparative UMAP process. That future process should use the same DIV90 UMAP/metadata audit outputs above and the matching DIV30 paper-cluster UMAP outputs as inputs. It should not be mixed into the URD root/tip smoke-test work.

Future Stiletti comparison inputs to collect in the separate process:

```text
DIV30 UMAP coordinates
DIV30 Seurat cluster IDs
DIV30 paper/manual biology labels
DIV90 UMAP coordinates
DIV90 cluster_id
DIV90 cluster_number_name
paper-level collapsed biology groups for DIV30 and DIV90
```

### DIV90 object inventory

From the existing inventory:

| Field | Value |
|---|---|
| Cells | 22,338 |
| Features | 18,082 |
| Source Seurat reductions | `pca`, `integrated.cca`, `umap` |
| Cached AnnData UMAP | `X_umap` |
| Default assay | `RNA` |
| Graphs | `RNA_nn`, `RNA_snn` |
| Real label column | `cluster_number_name` |
| Cluster ID columns | `seurat_clusters`, `cluster_id`, `RNA_snn_res.0.5` |
| Sample column | `orig.ident` |

DIV90 metadata columns:

```text
orig.ident
nCount_RNA
nFeature_RNA
percent.mt
S.Score
G2M.Score
Phase
old.ident
CC.Difference
RNA_snn_res.0.5
seurat_clusters
cluster_id
cluster_number_name
```

DIV90 `cluster_number_name` labels:

| Cluster | Label | n cells |
|---:|---|---:|
| 0 | MGE Striatal/GP Fated | 3,601 |
| 1 | SST+, NPY +, Cortical Fated | 3,548 |
| 2 | CRABP1+/PV Precursors | 3,503 |
| 3 | PV precursors/Migrating cells/Cortical-fated | 2,283 |
| 4 | Pre-Astrocytes/Astrocytes 1 | 1,987 |
| 5 | LHX8+ vMGE GABergic Striatal/GP fated 1 | 1,924 |
| 6 | Stressed Cells | 1,226 |
| 7 | Stressed Cells | 1,063 |
| 8 | LHX8+ vMGE GABergic Striatal/GP fated 2 | 922 |
| 9 | Pre-OPCs/OPCs | 915 |
| 10 | Pre-Astrocytes/Astrocytes 2 | 583 |
| 11 | PV Precursors | 425 |
| 12 | Dividing cells | 358 |

Samples:

| Sample | n cells |
|---|---:|
| 10496-MW-4 | 6,314 |
| 10496-MW-6 | 4,850 |
| 10496-MW-2 | 3,533 |
| 10496-MW-1 | 3,095 |
| 10496-MW-3 | 2,714 |
| 10496-MW-5 | 1,832 |

### DIV90 Jia-lineage-driven smoke-test plan

Status: planned, not yet run.

This plan supersedes the earlier broad SST/PV/MGE tip strategy. DIV30 used SST, PV, and MGE labels as practical smoke-test tips because the DIV30 paper-cluster framework was broad. DIV90 has more informative terminal labels, so the DIV90 URD tree should be built around Jia fetal MGE lineage endpoint logic instead of the paper annotation hierarchy.

Main DIV90 question:

```text
Which Jia fetal MGE lineage programs are represented by the terminal DIV90 states?
```

Not the primary DIV90 question:

```text
How do progenitors generate SST, PV, and MGE neurons?
```

Machine-readable planning files:

```text
metadata/div90_jia_lineage_urd_plan.tsv
metadata/div90_jia_rootscore_markers.tsv
metadata/jia_lineage_modules.tsv
```

Dedicated DIV90 handoff:

```text
python_notebooks/HANDOFF_div90_jia_lineage_urd_plan.md
```

First DIV90 Jia-lineage smoke run completed:

```text
Slurm job: 51685516
Template: slurm_templates/34_div90_jia_lineage_urd_smoke.sbatch.template
Run label: div90_urd_jia_lineage_smoke5k_knn100_v1
Output root: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_urd_jia_lineage_smoke5k_knn100_v1/
Log: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/34_div90_jia_lineage_urd_smoke_51685516.log
State: COMPLETED
Exit code: 0:0
Elapsed: 00:20:25
```

The run emitted the same smoke-panel outputs as DIV30 under `lineage_decision_report/`, built the Jia endpoint-tip tree under `lineage_tree_jia_endpoint_tips_v1/`, and ran the Jia Fig. S11-style marker validation under `jia_fig_s11_style_marker_validation_v1/`.

First-pass topology: `tip_epha5_mef2c` and `tip_lhx6_nfia` fuse first; the combined branch then splits from `tip_lhx8_isl1` at pseudotime ~0.267. This smoke result supports a stronger LHX8/ISL1-like branch versus a combined EPHA5/MEF2C/LHX6/NFIA-like branch, not a clean separation among all three requested Jia endpoint tips.

#### DIV90 root

Use cluster `12 - Dividing cells` as the root candidate pool, then derive a Jia-style progenitor root inside that pool.

RootScore:

```text
z(RGC1_score) + z(RGC2_score) + z(HES1) + z(VIM) + z(NES) - z(IPC_score) - z(DLX1) - z(DLX2) - z(ASCL1)
```

Root cells for first smoke tree:

```text
top 1-2% RootScore cells within cluster 12
```

Interpretation of the root: proliferative, VZ/SVZ-like progenitor cells that are high for RGC programs and core progenitor markers, but low for IPC/neurogenic markers.

Important dependency: if `RGC1_score`, `RGC2_score`, and `IPC_score` are not already present in the DIV90 object, compute them with the same Jia scoring logic used for DIV30 before root selection. Do not substitute Seurat cluster identity for the root score.

#### DIV90 tips

Do not define DIV90 tips as SST/PV/MGE broad labels. Define tips as candidate Jia lineage endpoints:

| Tip ID | Jia lineage target | DIV90 clusters | Exact metadata names | Validation genes | Rationale |
|---|---|---|---|---|---|
| `tip_lhx8_isl1` | LHX8/ISL1-like lineage | `0`, `5`, `8` | `0 - MGE Striatal/GP Fated`; `5 - LHX8+ vMGE GABergic Striatal/GP fated 1`; `8 - LHX8+ vMGE GABergic Striatal/GP fated 2` | `LHX8`, `ISL1`, `GBX2`, `TAC1` | strongest candidate for Jia LHX8/ISL1 lineage |
| `tip_epha5_mef2c` | EPHA5/MEF2C-like lineage | `2`, `3`, `11` | `2 - CRABP1+/PV Precursors`; `3 - PV precursors/Migrating cells/Cortical-fated`; `11 - PV Precursors` | `MEF2C`, `MAFB`, `ETV1`, `ERBB4` | strongest candidate for Jia EPHA5/MEF2C lineage |
| `tip_lhx6_nfia` | LHX6/NFIA-like lineage | `1` | `1 - SST+, NPY +, Cortical Fated` | `LHX6`, `SST`, `NPY`, `ERBB4`, `CXCR4`, `ARX` | strongest candidate for Jia LHX6/NFIA lineage |

For the first neuronal DIV90 URD smoke test, exclude or track separately:

| Cluster(s) | Exact metadata names | Reason |
|---|---|---|
| `4`, `10` | `Pre-Astrocytes/Astrocytes 1`; `Pre-Astrocytes/Astrocytes 2` | glial lineage, not a neuronal terminal lineage endpoint |
| `6`, `7` | `Stressed Cells`; `Stressed Cells` | stress-associated states can distort graph geometry |
| `9` | `Pre-OPCs/OPCs` | OPC lineage, not part of the first inhibitory neuron lineage test |

#### DIV90 post-tree required analyses

After the DIV90 URD tree exists:

1. Score the actual Jia lineage modules on the tree:
   - `LHX8/ISL1`
   - `NR2F1/NR2F2`
   - `EPHA5/MEF2C`
   - `LHX6/NFIA`
   - `CRABP1/ANGPT2`
2. Determine whether the chosen tips actually correspond to the expected Jia lineages.
3. Generate segment x lineage score and z-score heatmaps.
4. Generate branchpoint decision-gene tables and figures separately from Jia lineage score plots.
5. Generate the Jia Fig. S11-style marker-expression validation figure using the exact marker order:
   - Panel A: `HES1 | CACNA1E | DLX2 | DCX`
   - Panel B: `LHX8 | NR2F1 | EPHA5 | MEF2C | CRABP1`
6. Keep marker-expression validation, lineage module localization, and branchpoint DE as separate figure classes.

#### DIV90 open checks before first URD tree

- Confirm that the DIV90 object has, or can receive, Jia `RGC1_score`, `RGC2_score`, and `IPC_score` columns using the same scoring method as DIV30.
- Confirm that cluster 12 is connected to the neuronal manifold in the source UMAP/diffusion space.
- Confirm whether cluster 2 (`CRABP1+/PV Precursors`) behaves as EPHA5/MEF2C-like, CRABP1/ANGPT2-like, or transitional after lineage module scoring.
- Confirm whether NR2F1/NR2F2 and CRABP1/ANGPT2 need explicit terminal tips after the first Jia-lineage smoke tree, or whether they localize as internal/branch-associated programs.

### DIV90 handoff requirement

Before shipping the DIV90 handoff, record:

- exact source Seurat/H5AD path,
- annotation column used,
- clusters retained/excluded,
- root definition,
- tip definitions,
- URD parameters,
- run labels,
- Slurm job IDs,
- every required figure/table path,
- interpretation and caveats.

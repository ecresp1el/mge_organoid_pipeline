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
ELAPSED   4:22:09
NODE      gl3327
STAGE     calcDM(knn = 100, sigma = "local")
```

Do not stop it unless explicitly deciding to abandon the full diffusion-map run. The 5k smoke object is currently the completed object used for the tree proof of workflow.

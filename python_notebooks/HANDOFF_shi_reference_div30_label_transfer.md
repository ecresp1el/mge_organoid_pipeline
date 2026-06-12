# Handoff: Shi Reference DIV30 Label Transfer

Standalone workflow for mapping DIV30 AnnData cells to Shi et al. reference
cell-level labels.

This is not Notebook 02. It is a separate reference-mapping workflow.

## Active Add-On: Wang et al. 2025 hnbMO Cell Ranger Rerun

This handoff also tracks the active Wang/Liu 2025 hnbMO rerun because it is
currently running on Great Lakes.

Study:

```text
Wang et al. 2025 / GSE286235 / PRJNA1206345
Generation of human nucleus basalis organoids with functional nbM-cortical
cholinergic projections in transplanted assembloids
```

Goal:

```text
Rerun Cell Ranger from healthy SRA reads so downstream Seurat can use
Cell Ranger filtered_feature_bc_matrix outputs instead of the raw-H5
nFeature_RNA >= 1360 proxy filter.
```

Healthy samples:

```text
BF_H9_D36      H9 day 36
BF_H9_D63      H9 day 63
BFCO_IMR_D63   IMR90-4 day 63
```

Code/config added:

```text
config/liu_2025_hnbmo_healthy_sra_runs.tsv
scripts/01g_gse286235_sra_to_10x_fastqs.sh
scripts/05m_liu_2025_hnbmo_cellranger_count.sh
slurm_templates/01g_gse286235_sra_to_10x_fastqs_array.sbatch.template
slurm_templates/05m_liu_2025_hnbmo_cellranger_count_array.sbatch.template
```

Active Slurm jobs as of 2026-06-08 status check:

```text
FASTQ reconstruction array: 51484325
  51484325_1 BF_H9_D36      COMPLETED, ExitCode 0:0, elapsed 01:19:54
  51484325_2 BF_H9_D63      COMPLETED, ExitCode 0:0, elapsed 01:35:38
  51484325_3 BFCO_IMR_D63   COMPLETED, ExitCode 0:0, elapsed 01:55:22

Original Cell Ranger count array: 51484326
  51484326_[1-3] FAILED immediately, ExitCode 1:0.
  Cause: Great Lakes cellranger/6.1.2 module sourceme references unset
  _RUN10X while Slurm template had `set -u` active.

Corrected Cell Ranger count array: 51485738
  Submitted after patching slurm_templates/05m_liu_2025_hnbmo_cellranger_count_array.sbatch.template
  to relax nounset only around module load.
  51485738_1 BF_H9_D36      COMPLETED, ExitCode 0:0, elapsed 01:12:55
  51485738_2 BF_H9_D63      COMPLETED, ExitCode 0:0, elapsed 01:27:43
  51485738_3 BFCO_IMR_D63   COMPLETED, ExitCode 0:0, elapsed 02:08:48
```

Check status later:

```bash
squeue -j 51484325,51484326,51485738 -o '%.18i %.9P %.28j %.8u %.2t %.10M %.6D %R'
sacct -j 51484325,51484326,51485738 --format=JobID,JobName%28,State,ExitCode,Elapsed,MaxRSS,ReqMem -P
tail -n 60 /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/01g_gse286235_sra_fastq_51484325_3.log
tail -n 60 /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/05m_liu_2025_cellranger_count_51485738_1.log
```

Input/output locations:

```text
SRA + reconstructed FASTQs:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/data/raw/liu_2025_hnbmo_sra/

FASTQs:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/data/raw/liu_2025_hnbmo_sra/fastqs/{BF_H9_D36,BF_H9_D63,BFCO_IMR_D63}/

Cell Ranger outputs, once dependency releases:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/liu_2025_hnbmo_cellranger_counts/{BF_H9_D36,BF_H9_D63,BFCO_IMR_D63}/outs/filtered_feature_bc_matrix.h5
```

Cell Ranger 6.1.2 called-cell counts:

```text
BF_H9_D36       6,648 cells
BF_H9_D63       4,115 cells
BFCO_IMR_D63    8,438 cells
Total          19,201 cells
```

This is higher than the paper's reported 14,245 healthy cells, so the rerun is
successful technically but not numerically identical to the authors' Cell
Ranger v3.1/reference/filtering state.

FASTQ reconstruction assumptions:

```text
Use healthy SRA runs only.
Treat each SRA run as one lane for its sample.
Use fasterq-dump --split-files --include-technical.
Rename to Cell Ranger-style names: <sample>_S1_L00x_R{1,2}_001.fastq.gz.
```

Cell Ranger rerun assumptions:

```text
Great Lakes module: cellranger/6.1.2
Authors reported: Cell Ranger v3.1
Chemistry: --chemistry SC3Pv3
Reference: /nfs/turbo/umms-parent/Manny_human_ref/refdata-gex-GRCh38-2020-A
Storage saver: --no-bam
Cell calling: Cell Ranger default; no --force-cells and no nFeature proxy filter
```

Important caveat:

```text
This is closer to the authors' upstream workflow than filtering GEO raw H5s,
but it is not exact-author because Great Lakes does not expose Cell Ranger 3.1
and the paper did not report the exact GRCh38 reference bundle.
```

### Wang 2025 hnbMO Exploratory Seurat Processing

Status as of 2026-06-08:

```text
Exploratory only. Do not include this Wang object in the broader study/cross-study
analysis yet.
```

New script/template:

```text
scripts/05n_wang_2025_hnbmo_cellranger_seurat.R
slurm_templates/05n_wang_2025_hnbmo_cellranger_seurat.sbatch.template
```

Jobs:

```text
51508301 FAILED during integration after merged UMAP completed.
  Cause: duplicated 10x barcodes across samples confused split/integration path.

51508462 FAILED during integration.
  Cause: Seurat/future globals max size limit during FindIntegrationAnchors.

51508627 COMPLETED, ExitCode 0:0, elapsed 00:07:16, MaxRSS ~11.9G.
```

Fixes applied before successful job:

```text
Prefix cell names by sample at object creation with RenameCells(add.cell.id=sample_id).
Integrate from original per-sample objects rather than SplitObject(merged).
Set options(future.globals.maxSize = 8 * 1024^3).
```

Output root:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/wang_2025_hnbmo_cellranger_seurat_exploratory
```

Main objects:

```text
wang_2025_hnbmo_cellranger_merged_seurat.rds
wang_2025_hnbmo_cellranger_integrated_seurat.rds
```

Main plots:

```text
plots/merged_vs_integrated_umap_by_sample.{png,pdf}
plots/merged_umap_by_sample.{png,pdf}
plots/integrated_umap_by_sample.{png,pdf}
plots/merged_umap_by_cluster.{png,pdf}
plots/integrated_umap_by_cluster.{png,pdf}
plots/qc_cells_kept_filtered_by_sample.{png,pdf}
plots/qc_metric_violins_kept_vs_filtered.{png,pdf}
```

QC tracking:

```text
Cell Ranger filtered cells:        19,201
After additional Seurat QC:        18,705
Additional Seurat QC removed:         496

Additional Seurat QC defaults:
  nFeature_RNA >= 200
  nFeature_RNA < Inf
  percent.mt < 20
  gene min.cells = 3
```

Per-sample QC:

```text
BF_H9_D36      Cell Ranger 6,648 -> Seurat QC 6,534; removed 114
BF_H9_D63      Cell Ranger 4,115 -> Seurat QC 3,735; removed 380
BFCO_IMR_D63   Cell Ranger 8,438 -> Seurat QC 8,436; removed   2
```

Useful tables:

```text
tables/seurat_qc_filtering_by_sample.tsv
tables/seurat_qc_per_cell.tsv
tables/cell_count_totals.tsv
tables/analysis_assumptions.tsv
tables/merged_cluster_counts_by_sample.tsv
tables/integrated_cluster_counts_by_sample.tsv
```

UMAP interpretation:

```text
Merged/no-integration UMAP separates strongly by sample/cell line/timepoint.
Seurat integrated UMAP mixes samples onto a shared manifold and is closer in
spirit to the paper's Figure 1C presentation.
```

Requested follow-up as of 2026-06-08:

```text
Paper-reported healthy Figure 1C cell total: 14,245 cells.

New exploratory filter requested:
  nFeature_RNA >= 200
  nFeature_RNA < 6000
  percent.mt < 20

New output branch:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/wang_2025_hnbmo_cellranger_seurat_exploratory_maxFeature6000

New plots added by script:
  merged_umap_nFeature_RNA_after_qc.{png,pdf}
  merged_umap_nCount_RNA_after_qc.{png,pdf}
  integrated_umap_nFeature_RNA_after_qc.{png,pdf}
  integrated_umap_nCount_RNA_after_qc.{png,pdf}
  merged_vs_integrated_umap_qc_metrics_after_qc.{png,pdf}
  post_qc_violin_nFeature_nCount_by_sample.{png,pdf}

Slurm job:
  51509995 COMPLETED, ExitCode 0:0, elapsed 00:08:11, MaxRSS ~12.2G.

Result:
  Cell Ranger filtered cells: 19,201
  After added Seurat QC with nFeature_RNA < 6000: 18,638
  Total removed by added Seurat QC: 563

Per-sample result:
  BF_H9_D36      6,648 -> 6,506; removed 142
  BF_H9_D63      4,115 -> 3,712; removed 403
  BFCO_IMR_D63   8,438 -> 8,420; removed 18

User-requested violin plot:
  post-QC kept cells only, after nFeature_RNA < 6000 branch.
  Panels: nFeature_RNA and nCount_RNA by H9 D36, H9 D63, IMR90-4 D63.
  PNG:
    /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/wang_2025_hnbmo_cellranger_seurat_exploratory_maxFeature6000/plots/post_qc_violin_nFeature_nCount_by_sample.png
  PDF:
    /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/wang_2025_hnbmo_cellranger_seurat_exploratory_maxFeature6000/plots/post_qc_violin_nFeature_nCount_by_sample.pdf

New percent.mt < 10 branch requested on 2026-06-08:
  Keep the nFeature_RNA cap from the prior branch.
  Filters:
    nFeature_RNA >= 200
    nFeature_RNA < 6000
    percent.mt < 10

  Output branch:
    /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/wang_2025_hnbmo_cellranger_seurat_exploratory_maxFeature6000_percentMT10

  Slurm:
    51519761 CANCELLED before start; original request was 128G and sat pending.
    51519875 FAILED after 00:00:42, ExitCode 1:0.
      Cause: script generated the new post-QC violin before `sample_colors` was defined.
      Partial QC tables were written before failure.
    51520086 COMPLETED, ExitCode 0:0, elapsed 00:07:53, MaxRSS ~11.1G.
      Submitted after patching scripts/05n_wang_2025_hnbmo_cellranger_seurat.R
      to define `sample_colors` immediately after `sample_info`.

  Result:
    Cell Ranger filtered cells: 19,201
    After added Seurat QC with nFeature_RNA < 6000 and percent.mt < 10: 17,860
    Total removed by added Seurat QC: 1,341

  Per-sample result:
    BF_H9_D36      6,648 -> 6,149; removed 499
    BF_H9_D63      4,115 -> 3,312; removed 803
    BFCO_IMR_D63   8,438 -> 8,399; removed 39

  Median post-QC metrics:
    BF_H9_D36      median nFeature_RNA 2,642;   median nCount_RNA 7,133;   median percent.mt 5.28
    BF_H9_D63      median nFeature_RNA 2,164.5; median nCount_RNA 5,777.5; median percent.mt 4.51
    BFCO_IMR_D63   median nFeature_RNA 1,660;   median nCount_RNA 3,871;   median percent.mt 0.76

  Expected key plots when complete:
    plots/post_qc_violin_nFeature_nCount_by_sample.{png,pdf}
    plots/merged_vs_integrated_umap_by_sample.{png,pdf}
    plots/merged_vs_integrated_umap_qc_metrics_after_qc.{png,pdf}

  User-requested violin plot:
    /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/wang_2025_hnbmo_cellranger_seurat_exploratory_maxFeature6000_percentMT10/plots/post_qc_violin_nFeature_nCount_by_sample.png

  Check commands:
    squeue -j 51520086 -o '%.18i %.9P %.28j %.8u %.2t %.10M %.6D %R'
    sacct -j 51519761,51519875,51520086 --format=JobID,JobName%28,State,ExitCode,Elapsed,MaxRSS,ReqMem -P
    tail -n 80 /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/05n_wang_2025_hnbmo_cellranger_seurat_maxFeature6000_mt10_51520086.log
```

## Core Rule

DIV30 `seurat_clusters` are not used for prediction.

Predictions are made per DIV30 cell using expression similarity to Shi reference
cells. Existing DIV30 `seurat_clusters` are used only afterward for summaries,
stacked bars, heatmaps, and interpretation.

## Current Run Target

Default run label:

```text
shi_reference_div30_label_transfer_v2
```

Default run directory:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_reference_div30_label_transfer/shi_reference_div30_label_transfer_v2
```

Executed notebook path:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_reference_div30_label_transfer/executed/shi_reference_div30_label_transfer.shi_reference_div30_label_transfer_v2.executed.ipynb
```

## Code Added

Reusable module:

```text
python_notebooks/src/mge_organoid_python/shi_label_transfer.py
```

Notebook:

```text
python_notebooks/notebooks/shi_reference_div30_label_transfer.ipynb
```

Slurm template:

```text
slurm_templates/20_execute_shi_reference_div30_label_transfer.sbatch.template
```

## Runtime Environment And Slurm Management

Use the same project environment conventions as Notebook 00, the Seurat/AnnData
inventory, marker-expression plotting, and the Schmitz staging workflow.

Repo checkout:

```text
REPO_ROOT=/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
```

Large data/results root:

```text
PROJECT_ROOT=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
```

Python/conda environment:

```text
CONDA_ENV_BIN=/home/elcrespo/miniconda3/envs/mge-organoid-python/bin
env=mge-organoid-python
python=/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python
jupyter=/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/jupyter
kernel=mge-organoid-python
```

Required exports for direct CLI/debug runs:

```bash
export PROJECT_ROOT=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
export PATH=/home/elcrespo/miniconda3/envs/mge-organoid-python/bin:${PATH}
export PYTHONPATH=/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline/python_notebooks/src:${PYTHONPATH:-}
```

Do not rely on the login shell's generic `python`; use the explicit
`mge-organoid-python` interpreter or submit through Slurm templates.

### Standalone DIV30 Shi kNN/Notebook Workflow

Template:

```text
slurm_templates/20_execute_shi_reference_div30_label_transfer.sbatch.template
```

Prepared job file pattern:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/20_execute_shi_reference_div30_label_transfer.sbatch
```

Submit:

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
mkdir -p /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs
mkdir -p /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs

cp slurm_templates/20_execute_shi_reference_div30_label_transfer.sbatch.template \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/20_execute_shi_reference_div30_label_transfer.sbatch

sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/20_execute_shi_reference_div30_label_transfer.sbatch
```

Template defaults:

```text
account: parent0
partition: standard
cpus-per-task: 8
memory: 160G
time: 16:00:00
job name: shi-div30-xfer
```

Main environment variables accepted by the template:

```text
SHI_LABEL_TRANSFER_RESULTS_DIRNAME=shi_reference_div30_label_transfer
SHI_LABEL_TRANSFER_RUN_LABEL=shi_reference_div30_label_transfer_v2
SHI_LABEL_TRANSFER_QUERY_H5AD=$PROJECT_ROOT/results/python_anndata/varela_div30.h5ad
SHI_LABEL_TRANSFER_REFERENCE_H5AD=$PROJECT_ROOT/results/python_anndata/shi_2019_paper_qc.h5ad
SHI_LABEL_TRANSFER_TABLE_S2_XLSX=$PROJECT_ROOT/reference/shi_2021_tables_s2_to_s9/science.abj6641_table_s2.xlsx
SHI_LABEL_TRANSFER_N_TOP_VARIABLE_GENES=3000
SHI_LABEL_TRANSFER_N_PCS=50
SHI_LABEL_TRANSFER_N_NEIGHBORS=31
SHI_LABEL_TRANSFER_RESTRICTED_LABELS=MGE,LGE,CGE
SHI_LABEL_TRANSFER_WRITE_H5AD=1
SHI_LABEL_TRANSFER_SAVE_PLOTS=1
```

Expected executed notebook:

```text
$PROJECT_ROOT/results/shi_reference_div30_label_transfer/executed/shi_reference_div30_label_transfer.<RUN_LABEL>.executed.ipynb
```

Expected logs:

```text
$PROJECT_ROOT/logs/shi-label-transfer-shi-div30-xfer-<jobid>.out
$PROJECT_ROOT/logs/shi-label-transfer-shi-div30-xfer-<jobid>.err
```

### Cross-Study Shi Seurat Transfer Workflow

Templates:

```text
slurm_templates/26_cross_study_shi_seurat_label_transfer.sbatch.template
slurm_templates/27_cross_study_shi_seurat_label_transfer_array.sbatch.template
slurm_templates/28_finalize_cross_study_shi_prediction_plots.sbatch.template
slurm_templates/28_cross_study_shi_prediction_plots_plot_only.sbatch.template
```

These Seurat-transfer jobs load the Great Lakes Seurat module:

```text
module load Bioinformatics
SEURAT5_MODULE=r-seurat/5.1.0-R-4.4.1-c3m7yfq
module load "$SEURAT5_MODULE"
```

Single-job submit pattern:

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
cp slurm_templates/26_cross_study_shi_seurat_label_transfer.sbatch.template \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/26_cross_study_shi_seurat_label_transfer.sbatch
sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/26_cross_study_shi_seurat_label_transfer.sbatch
```

Array/finalizer pattern:

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
cp slurm_templates/27_cross_study_shi_seurat_label_transfer_array.sbatch.template \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/27_cross_study_shi_seurat_label_transfer_array.sbatch
cp slurm_templates/28_finalize_cross_study_shi_prediction_plots.sbatch.template \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/28_finalize_cross_study_shi_prediction_plots.sbatch

array_job=$(sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/27_cross_study_shi_seurat_label_transfer_array.sbatch | awk '{print $4}')
sbatch --dependency=afterok:${array_job} \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/28_finalize_cross_study_shi_prediction_plots.sbatch
```

Common cross-study environment variables:

```text
CROSS_STUDY_SHI_RUN_LABEL=cross_study_shi_seurat_label_transfer_v1
CROSS_STUDY_SHI_STUDIES=varela_div30,varela_div90,siebert_2026,walsh,bershteyn_2025,bershteyn_2023,samarasinghe_2021
CROSS_STUDY_SHI_REUSE_EXISTING=true
CROSS_STUDY_SHI_FORCE_RERUN=false
CROSS_STUDY_SHI_DIMS=50
CROSS_STUDY_SHI_MIN_SHARED_FEATURES=500
SHI_REFERENCE_RDS=$PROJECT_ROOT/results/shi_2019_paper_qc/shi_2019_seurat.rds
SHI_REFERENCE_LABELS_TSV=$PROJECT_ROOT/results/shi_reference_div30_seurat_label_transfer/shi_reference_div30_seurat_label_transfer_v1/seurat/shi_reference_labels_for_seurat.tsv
```

Expected cross-study logs:

```text
$PROJECT_ROOT/logs/cross-shi-transfer-cross-shi-transfer-<jobid>.out
$PROJECT_ROOT/logs/cross-shi-transfer-cross-shi-transfer-<jobid>.err
$PROJECT_ROOT/logs/cross-shi-xfer-array-<array_jobid>_<taskid>.out
$PROJECT_ROOT/logs/cross-shi-xfer-array-<array_jobid>_<taskid>.err
$PROJECT_ROOT/logs/cross-shi-final-<jobid>.out
$PROJECT_ROOT/logs/cross-shi-final-<jobid>.err
```

### Monitoring Commands

Use these checks for Shi-related Slurm jobs:

```bash
squeue -j <jobid> -o '%.20i %.9P %.30j %.8u %.2t %.10M %.6D %R'
sacct -j <jobid> --format=JobID,JobName%28,State,ExitCode,Elapsed,MaxRSS,ReqMem -P
tail -n 80 /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/<logfile>.out
tail -n 80 /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/<logfile>.err
```

## Inputs

DIV30 query AnnData:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/varela_div30.h5ad
```

Shi reference AnnData:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/shi_2019_paper_qc.h5ad
```

Shi Table S2 cell labels:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/reference/shi_2021_tables_s2_to_s9/science.abj6641_table_s2.xlsx
```

Table S2 provides cell-level `Major types` labels. The workflow joins those
labels to the cached Shi AnnData by normalized barcode.

Observed Table S2 labels during setup:

```text
MGE
LGE
CGE
progenitor
Excitatory IPC
Excitatory neuron
Thalamic neurons
Microglia
OPC
Endothelial
```

## Comparisons

The workflow runs two independent per-cell label-transfer comparisons.

### A. Full Relevant Developmental Reference

As of `shi_reference_div30_label_transfer_v2`, the default full reference
includes all 10 observed Shi Table S2 major-type labels:

```text
MGE
LGE
CGE
progenitor
Excitatory IPC
Excitatory neuron
Thalamic neurons
Microglia
OPC
Endothelial
```

### B. Restricted MGE/LGE/CGE Reference

Default included labels:

```text
MGE, LGE, CGE
```

This is intentionally compared against the full-reference result to identify
cells whose MGE/LGE/CGE label is potentially being forced by removing other
developmental states.

## Prediction Score Source of Truth

Use only the Seurat `TransferData` score contract for current Shi prediction
scores.

Authoritative implementation paths:

```text
scripts/13_run_cross_study_shi_seurat_label_transfer.R
python_notebooks/scripts/seurat_shi_label_transfer_export.R
python_notebooks/scripts/run_shi_seurat_label_transfer_smoke.py
python_notebooks/src/mge_organoid_python/cross_study_shi_prediction_plots.py
python_notebooks/src/mge_organoid_python/shi_prediction_schema.py
```

Exact Seurat-side score syntax:

```r
predictions <- as.data.frame(Seurat::TransferData(
  anchorset = anchors,
  refdata = reference@meta.data[[label_col]],
  dims = dims_use,
  k.weight = transfer_k_weight,
  verbose = TRUE
), stringsAsFactors = FALSE)

shi_seurat_full_prediction_score <- as.numeric(predictions$prediction.score.max)
shi_seurat_full_uncertainty_score <- 1 - shi_seurat_full_prediction_score

label_score_cols_raw <- setdiff(
  grep("^prediction\\.score\\.", colnames(predictions), value = TRUE),
  "prediction.score.max"
)
```

The per-label score columns are exported from Seurat's
`prediction.score.<label>` columns and renamed to
`shi_seurat_full_prediction_score_<label_token>`. The winning-label score is
`prediction.score.max`, i.e. Seurat's max anchor-weighted support score across
the returned label score columns.

Gestational-week scores use the same `TransferData` contract with
`refdata = <Shi week label>`. Downstream plotting may collapse duplicate Shi GW
score columns by canonical GW label and sum their support, then recompute the
canonical max GW label/score. That is a week-label granularity correction, not a
different prediction-score formula.

Do not use `GetTransferPredictions()` for this workflow unless an explicit
confidence policy is added. Current outputs intentionally do not apply its
default `score.filter = 0.75` unassignment behavior.

Non-authoritative legacy score path:

```text
python_notebooks/src/mge_organoid_python/shi_label_transfer.py::run_knn_label_transfer
python_notebooks/notebooks/shi_reference_div30_label_transfer.ipynb
```

That legacy path computes `prediction_score` as the fraction of k nearest Shi
neighbors voting for the winning label after sparse scaling, SVD projection, and
cosine kNN. It is retained only to explain the completed
`shi_reference_div30_label_transfer_v2` run and now requires an explicit
`allow_legacy_knn_scoring=True` opt-in in code.

## Original Shi 2019 Source/Reference Contents

Keep this distinction explicit for methods and discussion:

```text
The original Shi 2019 GEO download is not a deposited Seurat object with
prediction scores. It is a raw-count matrix with gestational-week information
encoded in cell/barcode suffixes, plus GEO series metadata:

data/raw/shi_2019_geo_files/suppl/GSE135827_GE_mat_raw_count_with_week_info.txt.gz
data/raw/shi_2019_geo_files/suppl/GSE135827_RAW.tar
data/raw/shi_2019_geo_files/matrix/GSE135827_series_matrix.txt.gz
data/raw/shi_2019_geo_files/miniml/GSE135827_family.xml.tgz
data/raw/shi_2019_geo_files/soft/GSE135827_family.soft.gz
```

Standalone Shi GEO object built from that matrix:

```text
results/shi_2019/shi_2019_seurat.rds

n_cells:         56,412
n_features:      21,289
reductions:      pca, umap
clusters:        25
weeks detected:  GW09, GW12_01, GW12_02, GW13, GW16, GW18_01

metadata columns:
  orig.ident
  nCount_RNA
  nFeature_RNA
  raw_cell_id
  has_duplicated_raw_cell_id
  study_id
  week_label
  week_numeric
  barcode
  percent.mt
  RNA_snn_res.0.8
  seurat_clusters
```

Paper-QC Shi GEO object:

```text
results/shi_2019_paper_qc/shi_2019_seurat.rds

cells before QC:        56,412
cells after paper QC:   56,136
features after QC:      21,191
reductions:             pca, umap
clusters:               25
weeks detected:         GW09, GW12_01, GW12_02, GW13, GW16, GW18_01

additional metadata compared with the default object:
  percent.hb
```

Shi major labels used as Seurat TransferData reference labels:

```text
The major cell-type labels were attached from the Shi supplementary annotation
table, not inferred from a downloaded prediction-score matrix.

Join summary:
  reference cells:              56,136
  table S2 label rows:          56,412
  matched reference cells:      55,704
  unmatched reference cells:       432
  matched reference fraction:   99.23%

Full matched Table S2 label counts:
  MGE                15,321
  progenitor         12,041
  LGE                11,767
  CGE                 7,128
  Thalamic neurons    4,212
  Excitatory neuron   3,736
  Excitatory IPC      1,065
  Microglia             232
  OPC                   163
  Endothelial            39
```

Cross-study v2 Shi reference actually used by Seurat TransferData:

```text
For cross-study v2, the reference was restricted to cells with both a major
Shi label and a gestational-week label. This gives 38,831 reference cells:

  MGE                10,446
  progenitor          8,278
  LGE                 6,879
  CGE                 4,322
  Thalamic neurons    4,212
  Excitatory neuron   3,388
  Excitatory IPC        937
  Microglia             196
  OPC                   136
  Endothelial            37

Week labels in the cross-study v2 reference:
  GW09      5,251
  GW12_01   6,704
  GW12_02  10,454
  GW13      5,426
  GW16      6,179
  GW18_01   4,817
```

GE-only age reference used in cross-study v2:

```text
The GE-only week classifier reruns FindTransferAnchors/TransferData using only
MGE/LGE/CGE Shi reference cells. Counts by label/week:

  MGE GW09:     1,564   LGE GW09:     1,025   CGE GW09:       86
  MGE GW12_01:  1,310   LGE GW12_01:  1,425   CGE GW12_01:   477
  MGE GW12_02:  1,513   LGE GW12_02:  1,484   CGE GW12_02:   471
  MGE GW13:     1,947   LGE GW13:     1,063   CGE GW13:      532
  MGE GW16:     2,110   LGE GW16:     1,359   CGE GW16:    1,317
  MGE GW18_01:  2,002   LGE GW18_01:    523   CGE GW18_01: 1,439
```

Methods implication:

```text
Shi supplies the biological reference labels and week labels used as refdata.
Our `shi_seurat_*prediction_score*` columns are newly computed Seurat
TransferData support scores from those reference labels. They are not original
Shi author score columns downloaded from GEO.
```

2026-06-09 implementation update:

```text
Codex refactored the active Shi scoring paths so the canonical score contract is
now enforced close to export/plot entry points.

Files changed:
  scripts/13_run_cross_study_shi_seurat_label_transfer.R
  python_notebooks/scripts/seurat_shi_label_transfer_export.R
  python_notebooks/scripts/run_shi_seurat_label_transfer_smoke.py
  python_notebooks/src/mge_organoid_python/cross_study_shi_prediction_plots.py
  python_notebooks/src/mge_organoid_python/shi_label_transfer.py
  python_notebooks/src/mge_organoid_python/shi_prediction_schema.py

New/strengthened invariants:
  - `predicted.id` and `prediction.score.max` must be present after TransferData.
  - At least one Seurat `prediction.score.<label>` column must be present.
  - `prediction.score.max` must be finite, in [0,1], and equal to the row-wise
    maximum of the per-label support scores.
  - `shi_seurat_full_uncertainty_score` must equal
    `1 - shi_seurat_full_prediction_score`.
  - Canonical exported score columns must use
    `shi_seurat_full_prediction_score_<label_token>`.
  - Python plotting/smoke code now uses the same label-token sanitizer via
    `shi_prediction_schema.py`.
  - Reused/old tables are validated against the canonical identities before
    plotting, so stale kNN-style or raw-Seurat-named score tables may now fail
    early instead of silently flowing into plots.

QC-only additions:
  - `shi_seurat_full_score_delta`
  - `shi_seurat_full_low_confidence_flag`
  - `n_anchors`
  - `n_unique_query_anchor_cells`
  - `k_weight_used`
  - `k_weight_reason`

Diagnostics now also log:
  - `k_weight_requested`
  - `k_weight_used`
  - `k_weight_reason`
  - `median_prediction_score_max`
  - `fraction_prediction_score_max_ge_0_75`
  - winner-take-all composition by top transferred Shi label

The 0.75 field is a confidence/QC metric only. It is not used as a default
composition or plotting filter.
```

## Legacy v2 Method

This describes the completed `shi_reference_div30_label_transfer_v2` run only.
Do not use this Scanpy-side kNN vote-fraction method for the next scoring
target.

1. Load DIV30 query AnnData and Shi reference AnnData.
2. Read Shi Table S2 `Major types` labels.
3. Join Table S2 labels to Shi reference cells by barcode.
4. Harmonize query/reference genes with exact and unambiguous case-insensitive
   gene matching.
5. For each comparison:
   - subset Shi reference cells to the relevant label set
   - select top variable genes from the Shi reference among shared genes
   - fit sparse StandardScaler and TruncatedSVD on the Shi reference
   - project DIV30 cells into the same SVD space
   - run kNN from each DIV30 cell to Shi reference cells
   - predict labels by neighbor-vote at the single-cell level
   - predict Shi gestational week from the same neighbor set
6. Add predictions to DIV30 `adata.obs`.
7. Summarize predictions by DIV30 `seurat_clusters`.
8. Summarize Shi gestational-week predictions by DIV30 sample.

Prediction score is the fraction of k nearest Shi reference neighbors voting for
the winning label. Uncertainty is `1 - prediction_score`. Entropy is normalized
neighbor-label entropy. The gestational-week score uses the same vote-fraction
logic, but over Shi neighbor `week_label`; mean/median neighbor week use the
numeric Shi gestational week.

## Next Target: Seurat Prediction Scores

Default next run label:

```text
shi_reference_div30_seurat_label_transfer_v1
```

This next target should replace Scanpy-side prediction scoring with Seurat
anchor-based label transfer. Scanpy/AnnData should remain the plotting and
summary environment after Seurat predictions are exported and joined back onto
DIV30 cells.

Core rule for the next target:

```text
Do not compute prediction_score from Scanpy kNN neighbor-vote fractions.
Use Seurat's transferred prediction scores as the prediction-score source of
truth.
```

Recommended first runnable target:

1. Run one Seurat label-transfer comparison using the full Shi reference label
   set.
2. Export per-cell Seurat predictions and the full label-score matrix.
3. Import those tables into the existing DIV30 AnnData.
4. Generate a small smoke-test plot set:
   - UMAP colored by predicted Shi label
   - UMAP colored by maximum Seurat prediction score
   - UMAP panels for MGE, LGE, and CGE score columns
   - ridge/density-style score distributions by DIV30 sample
   - stacked bar of predicted labels by DIV30 `seurat_clusters`
5. Only after that runs, add the MGE/LGE/CGE-restricted comparison and
   gestational-week transfer.

Expected new output roots:

```text
RUN_DIR/seurat/
RUN_DIR/tables/
RUN_DIR/plots/
RUN_DIR/h5ad/
```

Expected Seurat-side files:

```text
seurat/div30_shi_seurat_full_predictions.tsv.gz
seurat/div30_shi_seurat_full_prediction_scores.tsv.gz
seurat/div30_shi_seurat_full_transfer_diagnostics.tsv
```

Expected AnnData-side files:

```text
tables/div30_shi_seurat_label_transfer_obs.tsv.gz
tables/div30_shi_seurat_label_transfer_label_scores_long.tsv.gz
h5ad/div30_shi_seurat_label_transfer_predictions.h5ad
```

Required new `adata.obs` columns for the first target:

```text
shi_seurat_full_predicted_shi_label
shi_seurat_full_prediction_score
shi_seurat_full_uncertainty_score
shi_seurat_full_broad_region_class
shi_seurat_full_developmental_class
```

Required per-label score columns in `adata.obs` for plotting:

```text
shi_seurat_full_prediction_score_MGE
shi_seurat_full_prediction_score_LGE
shi_seurat_full_prediction_score_CGE
shi_seurat_full_prediction_score_progenitor
shi_seurat_full_prediction_score_Excitatory_IPC
shi_seurat_full_prediction_score_Excitatory_neuron
shi_seurat_full_prediction_score_Thalamic_neurons
shi_seurat_full_prediction_score_Microglia
shi_seurat_full_prediction_score_OPC
shi_seurat_full_prediction_score_Endothelial
```

Column-name sanitization rule: preserve the original Shi label in a companion
long table, but replace spaces and punctuation with underscores for wide
`adata.obs` score columns.

Implementation added for the first smoke target:

```text
python_notebooks/scripts/seurat_shi_label_transfer_export.R
python_notebooks/scripts/run_shi_seurat_label_transfer_smoke.py
slurm_templates/21_run_shi_seurat_label_transfer_smoke.sbatch.template
```

Submitted smoke run:

```text
Slurm job: 51405243
Job script:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/21_run_shi_seurat_label_transfer_smoke.sbatch
Expected run directory:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_reference_div30_seurat_label_transfer/shi_reference_div30_seurat_label_transfer_v1
Logs:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/shi-seurat-label-transfer-shi-seurat-xfer-51405243.out
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/shi-seurat-label-transfer-shi-seurat-xfer-51405243.err
```

Completed Seurat transfer/plot runs:

```text
DIV30 refreshed plot run:
  Slurm job: 51410249
  State: COMPLETED
  ExitCode: 0:0
  Output root:
    /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_reference_div30_seurat_label_transfer/shi_reference_div30_seurat_label_transfer_v1
  Completion:
    n_query_cells = 90631
    n_score_columns = 10
    n_plots = 8
  Notes:
    - reused existing Seurat prediction tables
    - regenerated plots with non-squished cluster UMAP panel layout
    - subtype overlay panels now keep all Shi labels, including zero-count CGE/Endothelial

DIV30 unexpected-label marker summary:
  Slurm job: 51410250
  State: COMPLETED
  ExitCode: 0:0
  Output table:
    tables/markers/shi_seurat_unexpected_label_marker_summary.tsv
  Result summary:
    Thalamic neurons: 556 cells, mean max score 0.604
    Excitatory neuron: 54 cells, mean max score 0.534
    Excitatory IPC: 1 cell, mean max score 0.473
  Interpretation note:
    These are low-confidence Seurat edge calls. The "Excitatory neuron" cells
    are high for SST/GAD1, so do not interpret that label literally without
    checking markers and score thresholds.

DIV90 Seurat transfer:
  Slurm job: 51410263
  State: COMPLETED
  ExitCode: 0:0
  Elapsed: 00:06:08
  MaxRSS: 22621688K
  Output root:
    /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_reference_div90_seurat_label_transfer/shi_reference_div90_seurat_label_transfer_v1
  Query inputs:
    AnnData: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/varela_div90.h5ad
    Seurat: /nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/umap_props_output/clustered_day90_with_cluster_names_2.rds
  Completion:
    n_query_cells = 22338
    n_score_columns = 10
    n_plots = 8
  Diagnostics:
    reference_cells_labelled = 55704
    shared_features_used = 2542
    anchors = 3249
```

Completed Seurat plot inventory for each timepoint:

```text
<slug>_umap_shi_seurat_full_predicted_shi_label.png
<slug>_umap_shi_seurat_full_predicted_shi_label_with_subtype_overlays.png
<slug>_umap_shi_seurat_full_prediction_score.png
<slug>_umap_shi_seurat_full_mge_lge_cge_score_panel.png
<slug>_shi_seurat_full_prediction_score_density_by_sample.png
<slug>_shi_seurat_full_expected_shi_gw_ridge_by_sample.png
<slug>_shi_seurat_full_shi_label_stacked_bar_by_seurat_clusters.png
<slug>_shi_seurat_full_shi_label_stacked_bar_by_seurat_clusters_with_cluster_umaps.png
```

Rsync completed Seurat plot folders from a Mac:

```bash
mkdir -p /Users/ecrespo/Downloads/shi_reference_div30_seurat_label_transfer_v1/plots
mkdir -p /Users/ecrespo/Downloads/shi_reference_div90_seurat_label_transfer_v1/plots

rsync -avh --progress \
  elcrespo@greatlakes.arc-ts.umich.edu:/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_reference_div30_seurat_label_transfer/shi_reference_div30_seurat_label_transfer_v1/plots/ \
  /Users/ecrespo/Downloads/shi_reference_div30_seurat_label_transfer_v1/plots/

rsync -avh --progress \
  elcrespo@greatlakes.arc-ts.umich.edu:/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_reference_div90_seurat_label_transfer/shi_reference_div90_seurat_label_transfer_v1/plots/ \
  /Users/ecrespo/Downloads/shi_reference_div90_seurat_label_transfer_v1/plots/
```

Completed DIV30/DIV90 shared-axis comparison plots:

```text
Script:
  python_notebooks/scripts/plot_shi_seurat_timepoint_comparisons.py
Output root:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_reference_div30_div90_seurat_label_transfer_comparison/shi_reference_div30_div90_seurat_comparison_v1
Inputs:
  DIV30 tables/div30_shi_seurat_label_transfer_obs.tsv.gz
  DIV90 tables/div90_shi_seurat_label_transfer_obs.tsv.gz
Completion:
  n_cells = 112969
  n_plots = 8
  numeric_gw_range = [9, 18]
  even_gw_range = [0, 4]
  kde_bw_method = 0.08
  kde_min_density = 0.01
Plots:
  plots/div30_div90_shi_seurat_full_expected_shi_gw_ridge_by_sample_shared_axis.png
  plots/div30_div90_shi_seurat_full_expected_shi_gw_density_by_timepoint_shared_axis.png
  plots/div30_div90_shi_seurat_full_expected_shi_gw_ridge_by_sample_shared_axis_even_gw_spacing.png
  plots/div30_div90_shi_seurat_full_expected_shi_gw_density_by_timepoint_shared_axis_even_gw_spacing.png
  plots/div30_div90_shi_seurat_full_expected_shi_gw_ridge_by_sample_shared_axis_predicted_mge.png
  plots/div30_div90_shi_seurat_full_expected_shi_gw_density_by_timepoint_shared_axis_predicted_mge.png
  plots/div30_div90_shi_seurat_full_expected_shi_gw_ridge_by_sample_shared_axis_even_gw_spacing_predicted_mge.png
  plots/div30_div90_shi_seurat_full_expected_shi_gw_density_by_timepoint_shared_axis_even_gw_spacing_predicted_mge.png
Tables:
  tables/div30_div90_expected_shi_gw_summary_by_sample.tsv
  tables/div30_div90_expected_shi_gw_summary_by_sample_predicted_mge.tsv
  tables/div30_div90_shi_label_counts_by_timepoint.tsv
  tables/div30_div90_shi_label_counts_by_timepoint_predicted_mge.tsv
```

Rsync the shared-axis comparison plots from a Mac:

```bash
mkdir -p /Users/ecrespo/Downloads/shi_reference_div30_div90_seurat_comparison_v1/plots

rsync -avh --progress \
  elcrespo@greatlakes.arc-ts.umich.edu:/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_reference_div30_div90_seurat_label_transfer_comparison/shi_reference_div30_div90_seurat_comparison_v1/plots/ \
  /Users/ecrespo/Downloads/shi_reference_div30_div90_seurat_comparison_v1/plots/
```

## Cross-Study Prediction Inventory

This inventory was run before bringing the cross-study UMAP objects into the
notebook and before extending Shi-style Seurat prediction scores beyond the
currently plotted Varela DIV30/DIV90 objects.

Current inclusion decision:

- **Include as cross-study UMAP / target datasets:** Varela DIV30, Varela DIV90,
  Walsh, Bershteyn 2025, Bershteyn 2023, Samarasinghe 2021, and Siebert 2026.
- **Do not include yet:** Xiang et al. The object exists, but it is not ready for
  this comparison set because sample metadata/biology require additional
  reconstruction and validation.
- **Reference only:** Shi et al. 2019 is the label-transfer reference. Do not use
  Shi as a cross-study target UMAP/comparison dataset.

## Cross-Study Shi Prediction Plot Workflow

The prediction-score UMAP grids are handled separately from marker expression:

```text
Python module:
  python_notebooks/src/mge_organoid_python/cross_study_shi_prediction_plots.py
Notebook:
  python_notebooks/notebooks/cross_study_shi_prediction_scores.ipynb
CLI wrapper:
  python_notebooks/scripts/run_cross_study_shi_prediction_plots.py
Full finalizer Slurm template:
  slurm_templates/28_finalize_cross_study_shi_prediction_plots.sbatch.template
UMAP-only plot Slurm template:
  slurm_templates/28_cross_study_shi_prediction_plots_plot_only.sbatch.template
```

Current prediction run label:

```text
cross_study_shi_seurat_label_transfer_v1
```

Prediction UMAP aesthetic update on 2026-06-08:

- Continuous Shi prediction-score UMAP grids now use the same broad aesthetic as
  the marker-expression v12 grids: studies as rows, score classes as columns,
  tight row/column spacing, grey background cells, fixed-size whiteBlue score
  overlays, and aligned 60%-width per-column colorbars.
- Prediction scores keep their natural score scale of `0` to `1`; unlike
  marker expression, they are not q99-scaled per feature.
- The plotting code no longer sorts cells by score before drawing. Continuous
  score overlays use table/cell order. Categorical predicted-label UMAPs also
  draw cells in table/cell order rather than category order.
- Plot-only Slurm job `51530906` wrote the refreshed UMAP grids and was then
  canceled while it was spending extra time on summary plots. The UMAP outputs
  were complete before cancellation. The plot-only Slurm template has since been
  updated to call the new `plot-umap` command so future aesthetic rerenders only
  regenerate UMAP grids.
- Prediction UMAP rerender job `51536925` completed after adding
  `bershteyn_2025` back into the prediction UMAP grids. `samarasinghe_2021`
  remains filtered internally to control samples only. The filter audit table is:
  `tables/cross_study_shi_umap_internal_plot_filter_summary.tsv`.
- Plot-filtered predicted-GW composition outputs were added for the paper-style
  GW18 comparison across studies. The denominator is each study's cells passing
  the internal plot filters, not a score-thresholded subset:
  `tables/cross_study_shi_predicted_gw_label_counts_by_study_plot_filtered.tsv`,
  `tables/cross_study_shi_predicted_gw_label_fractions_by_study_plot_filtered.tsv`,
  and `tables/cross_study_shi_predicted_gw18_fraction_by_study_plot_filtered.tsv`.
  Plots are:
  `plots/summary/cross_study_shi_predicted_gw_label_fractions_by_study_plot_filtered.png`
  and `plots/summary/cross_study_shi_predicted_gw18_fraction_by_study_plot_filtered.png`.
- A second denominator-matched version was added for cells whose top transferred
  Shi major label is `MGE`. This is still not score-thresholded. Tables are:
  `tables/cross_study_shi_predicted_gw_label_counts_by_study_mge_predicted_plot_filtered.tsv`,
  `tables/cross_study_shi_predicted_gw_label_fractions_by_study_mge_predicted_plot_filtered.tsv`,
  and `tables/cross_study_shi_predicted_gw18_fraction_by_study_mge_predicted_plot_filtered.tsv`.
  Plots are:
  `plots/summary/cross_study_shi_predicted_gw_label_fractions_by_study_mge_predicted_plot_filtered.png`
  and `plots/summary/cross_study_shi_predicted_gw18_fraction_by_study_mge_predicted_plot_filtered.png`.
- `Max label score` and `Max GW score` are no longer included as UMAP grid
  columns. There is no maximum-prediction-score threshold/filter for including
  cells in prediction UMAPs.

Prediction-score policy after 2026-06-09 refactor:

- The exported Shi scores come from Seurat `TransferData` output: `predicted.id`,
  `prediction.score.max`, and the per-class `prediction.score.<label>` columns.
  This is Seurat's anchor-weighted label-support score.
- The pipeline does **not** call `GetTransferPredictions()` and therefore does
  **not** apply its default `score.filter = 0.75` behavior that would return
  `"Unassigned"` for lower-confidence cells.
- Current UMAP plots therefore show all cells that pass the plot-level filters
  using the raw Seurat prediction/support scores. There is no score cutoff.
- Winner-take-all `TransferData` labels/scores are the default composition
  logic. Confidence-filtered outputs such as `score >= 0.75` must remain
  explicitly labeled QC/sensitivity outputs and must not replace the default
  published-style composition unless requested.
- For gestational-week score panels, duplicate Shi week classes are collapsed
  after Seurat transfer by summing default per-class support scores that map to
  the same canonical GW label, e.g. `GW12`, `GW12_01`, and `GW12_02` become
  one plotted `GW12` score column; `GW18` and `GW18_01` become one plotted
  `GW18` score column. The individual GW score grid now uses the ordered unique
  GW columns `GW09`, `GW12`, `GW13`, `GW16`, `GW18`. This is a
  label-granularity correction, not a new anchor-weighting method.

## Cross-Study Marker Expression Plot Workflow

This is the Python-side workflow for the first major cross-study marker-expression
figure. It is **not** a Shi prediction-score plot. It plots marker-gene expression
directly on each study's own UMAP coordinates.

Code entrypoints:

```text
Python module:
  python_notebooks/src/mge_organoid_python/cross_study_marker_expression.py
Notebook:
  python_notebooks/notebooks/cross_study_marker_expression.ipynb
CLI wrapper:
  python_notebooks/scripts/run_cross_study_marker_expression.py
Seurat exporter:
  scripts/12_export_cross_study_marker_expression_tables.R
Slurm template:
  slurm_templates/25_cross_study_marker_expression.sbatch.template
```

Current run label:

```text
cross_study_marker_expression_v12
```

Data provenance and expression scale:

- Varela DIV30 and Varela DIV90 are read from existing AnnData/H5AD caches:
  `results/python_anndata/varela_div30.h5ad` and
  `results/python_anndata/varela_div90.h5ad`.
- Walsh, Bershteyn 2025, Bershteyn 2023, Samarasinghe 2021, and Siebert 2026
  are exported from the registered Seurat objects using
  `scripts/12_export_cross_study_marker_expression_tables.R`.
- The standardized per-cell tables are written under:
  `results/cross_study_marker_expression/<run_label>/tables/per_study/`.
- Each table contains `cell_id`, `study_id`, `study_label`, `sample`, `cluster`,
  `umap_1`, `umap_2`, and one column per marker gene.
- Expression values are log-normalized expression used for plotting. Varela,
  Walsh, Bershteyn 2025, Bershteyn 2023, and Samarasinghe use the configured
  Seurat/AnnData RNA `data` values. Siebert 2026 is exported from RNA `counts`
  and transformed during export to `log1p(CP10K)`.
- The figure label is:
  `log1p(CP10K); each gene column has its own scale`.

Current marker panels:

- Core ON-target genes:
  `DCX`, `GAD2`, `DLX5`, `LHX6`, `MAF`, `SST`, `PVALB`, `ERBB4`, `MEF2C`,
  `MAFB`, `LHX8`, `NKX2-1`.
- Core OFF-target genes:
  `SP8`, `PAX6`, `NEUROD2`, `ISL1`, `ACHE`, `NKX6-2`, `MKI67`.
- PV Precursors genes:
  `MAFB`, `MEF2C`, `ERBB4`, `ETV1`, `CRABP1`, `TAC1`, `ST18`, `PVALB`.
- PV Precursors paired OFF-target genes:
  `SP8`, `EBF1`, `NKX2-2`, `RAX`, `HMX3`, `DBH`.
- The marker panels are defined in
  `python_notebooks/src/mge_organoid_python/cross_study_marker_expression.py`
  as named `MarkerGenePanel` objects. They are also exported to
  `tables/cross_study_marker_expression_gene_panels.tsv`, so the exact panel
  membership used by each run is auditable.

Layout logic:

- Studies are rows.
- Genes are columns.
- Study rows are ordered as:
  `This Study, DIV 30`, `This Study, DIV 90`, `Siebert et al. 2026`,
  `Walsh et al. 2025`, `Bershteyn et al. 2025`, `Bershteyn et al. 2023`,
  `Samarasinghe et al. 2021`.
- Row labels include the number of cells plotted for that study.
- All cells for each study are drawn in grey as the UMAP background.
- Cells with expression greater than zero are overlaid using a Seurat-like
  `whiteBlue` expression colormap, and each per-gene colorbar uses that same
  white-to-blue scale.
- Marker-expression v7 uses 3x larger dots than v6 for both the grey background
  layer and the colored expression overlay.
- Marker-expression v8 adds an internal UMAP plotting filter for
  `samarasinghe_2021`: only samples whose `sample` value contains `Ctrl` are
  visualized. The row label is marked `Controls only`, and the cell count shown
  on the figure is the post-filter control-cell count.
- The Samarasinghe control-only filter is deliberately internal to
  `cross_study_marker_expression.py`, because it defines the biological subset
  for this figure rather than acting as a general user-facing plot option. The
  filter validates that only control samples remain and writes an audit table:
  `tables/cross_study_marker_expression_internal_plot_filter_summary.tsv`.
- Each combined panel draws a vertical divider between that panel's ON/PV
  precursor genes and its paired OFF-target genes.
- Each gene column has its own colorbar. That color scale is shared across all
  study rows for that gene, so comparisons are fair within a gene across
  studies. Do **not** interpret color intensity as directly comparable between
  different gene columns, because each gene has its own scale.
- Marker-expression v9 passes raw expression values directly to the color layer.
  For each gene, the colorbar spans raw expression units from `0` to the
  99th percentile of positive plotted expression values for that gene. Values
  above that 99th percentile are drawn at the top color only, so outliers do not
  stretch the colorbar. The expression table itself is not modified.
- Marker-expression v10 uses `viridis` for the expression overlay/colorbars and
  writes an audit table computed from the exact same plotted data and colorbar
  rules:
  `tables/cross_study_marker_expression_distribution_audit.tsv`. This table
  summarizes expression distributions/ranges by study and gene, including
  positive-cell counts, quantiles, raw maxima, per-gene colorbar maxima, and the
  number/percentage of values above the plotted colorbar maximum.
- Marker-expression v11 keeps the original core ON/OFF plots and adds a second
  named panel for `PV Precursors` plus its paired OFF-target genes. New plot
  tokens are `pv_precursors_on_target`, `pv_precursors_off_target`, and
  `pv_precursors_on_off_target`. The extraction gene table now contains the
  unique union of all genes needed by both panels.
- Marker-expression v12 tightens row and column spacing, switches the
  expression overlay/colorbars to `whiteBlue`, removes `ZFHX3` from the PV
  paired OFF-target panel, and sets the PV panel order to:
  `MAFB`, `MEF2C`, `ERBB4`, `ETV1`, `CRABP1`, `TAC1`, `ST18`, `PVALB` /
  `SP8`, `EBF1`, `NKX2-2`, `RAX`, `HMX3`, `DBH`.
- A v12 plot-only rerender keeps the same marker tables/data but shortens the
  bottom per-gene colorbars to 60% of each column width while keeping them
  aligned in one row. The expression overlay uses the table/cell order directly;
  it does not sort cells by expression before drawing. The plot-only Slurm
  template is:
  `slurm_templates/25_cross_study_marker_expression_plot_only.sbatch.template`.

Plot-study toggle:

- Study exclusion is a plotting filter, not a data-export deletion.
- The Slurm/notebook environment variable is:
  `CROSS_STUDY_MARKER_EXCLUDE_STUDIES`.
- Default for the current iteration:
  `CROSS_STUDY_MARKER_EXCLUDE_STUDIES=bershteyn_2025`.
- To include all studies in a future run, set:
  `CROSS_STUDY_MARKER_EXCLUDE_STUDIES=none`.
- To exclude multiple studies, use comma-, semicolon-, or space-separated study
  IDs, for example:
  `CROSS_STUDY_MARKER_EXCLUDE_STUDIES="bershteyn_2025 walsh"`.

Script:

```text
scripts/07_cross_study_prediction_inventory.R
```

Slurm template:

```text
slurm_templates/22_cross_study_prediction_inventory.sbatch.template
```

Submitted job:

```text
Slurm job: 51414106
State: COMPLETED
ExitCode: 0:0
Elapsed: 00:03:39
MaxRSS: 37875628K
Node: gl3420
Job script:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/22_cross_study_prediction_inventory.sbatch
Logs:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/cross-study-inventory-cross-study-inv-51414106.out
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/cross-study-inventory-cross-study-inv-51414106.err
```

Output table directory:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_prediction_inventory/cross_study_prediction_inventory_v1/tables
```

Inventory tables:

```text
cross_study_prediction_inventory_complete.tsv
cross_study_prediction_object_summary.tsv
cross_study_prediction_readiness.tsv
cross_study_prediction_metadata_columns.tsv
cross_study_prediction_metadata_value_counts.tsv
cross_study_prediction_primary_sample_counts.tsv
cross_study_prediction_primary_cluster_counts.tsv
div30_div90_existing_shi_prediction_table_summary.tsv
div30_div90_existing_shi_prediction_sample_counts.tsv
div30_div90_existing_shi_prediction_cluster_counts.tsv
div30_div90_existing_shi_prediction_label_counts.tsv
```

Objects inspected:

| study_id | object | cells | RNA features | UMAP for plotting | shared Shi genes |
| --- | --- | ---: | ---: | --- | ---: |
| `varela_div30` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/varela_this_paper/varela_this_paper_seurat.rds` | 90,631 | 18,082 | `umap` | 14,354 |
| `varela_div90` | `/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/umap_props_output/clustered_day90_with_cluster_names_2.rds` | 22,338 | 18,082 | `umap` | 14,354 |
| `walsh` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/walsh_day75/walsh_day75_final_annotated.rds` | 4,519 | 20,194 | `umap_sel` | 14,945 |
| `bershteyn_2025` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/bershteyn_2025/bershteyn_2025_seurat.rds` | 124,583 | 45,068 | `umap` | 20,021 |
| `bershteyn_2023` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/bershteyn_2023/bershteyn_2023_seurat.rds` | 98,042 | 45,068 | `umap` | 20,021 |
| `xiang_2018` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/xiang_2018/xiang_2018_seurat.rds` | 58,950 | 23,287 | `umap` | **exclude for now; not ready for cross-study inclusion** |
| `samarasinghe_2021` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/samarasinghe_2021_zenodo_processed_object/samarasinghe_2021_zenodo_seurat.rds` | 49,942 | 27,379 | `umap` | Shi overlap not re-audited after Zenodo registration |
| `siebert_2026` | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siebert_2026/siebert_2026_seurat.rds` | 64,676 | 32,131 | `umap` | 15,305 |

Readiness conclusion:

```text
Seven target cross-study objects are currently usable: Varela DIV30, Varela
DIV90, Walsh, Bershteyn 2025, Bershteyn 2023, Samarasinghe 2021, and Siebert
2026. Samarasinghe was updated to use the official Zenodo processed `datExpr`
object instead of rebuilding from GEO counts. Re-run the cross-study
prediction inventory/transfer step before interpreting Samarasinghe Shi
label-transfer scores, because the older handoff audit was performed before
this Zenodo object was registered.

Xiang is **not ready for inclusion**. The existing `xiang_2018` Seurat object
has only `orig.ident = Xiang2018`, and the useful sample biology must be
reconstructed from GSE98201 Cell Ranger aggregation suffixes `-1` through `-8`.
The raw Xiang rownames are also Ensembl IDs, so marker/label-transfer work
requires the GSE98201 feature map:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/data/raw/xiang_2018_geo_files/suppl/GSE98201_genes.tsv.gz

Until the Xiang metadata and biology mapping are validated, exclude Xiang from
cross-study UMAP comparisons and Shi-style target label-transfer runs.

Xiang 2017 tSNE reproduction run:

```text
Slurm job: 51483745
Script: scripts/10_xiang_2017_tsne_reproduction.R
Slurm template: slurm_templates/24_xiang_2017_tsne_reproduction.sbatch.template
Output root:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/xiang_2017_tsne_reproduction/xiang_2017_tsne_reproduction_v1
```

This is a modern approximation of the Xiang et al. Seurat v1.4 workflow unless
an old Seurat v1.4 runtime is added later. The run uses paper-style filters
(`min.cells = 4`, `min.features = 201`), per-cell log2 normalization, variable
genes with dispersion > 0, regression of `sample` and `nCount_RNA`, 20 PCs, and
tSNE using PC1-PC5 only.

Recovered Xiang barcode suffix metadata:

| barcode_suffix | sample | condition | timepoint | replicate | GEO accession | raw barcodes |
| --- | --- | --- | --- | --- | --- | ---: |
| `-1` | `hMGEO_d30_rep1` | `hMGEO` | `d30` | `rep1` | `GSM2589129` | 6,480 |
| `-2` | `hCO_d30_rep1` | `hCO` | `d30` | `rep1` | `GSM2589130` | 2,969 |
| `-3` | `hMGEO_d72_rep1` | `hMGEO` | `d72` | `rep1` | `GSM2589131` | 9,722 |
| `-4` | `hCO_d72_rep1` | `hCO` | `d72` | `rep1` | `GSM2589132` | 10,258 |
| `-5` | `hMGEO_d30_rep2` | `hMGEO` | `d30` | `rep2` | `GSM2684867` | 5,438 |
| `-6` | `hCO_d30_rep2` | `hCO` | `d30` | `rep2` | `GSM2684868` | 6,355 |
| `-7` | `hMGEO_d79_rep2` | `hMGEO` | `d79` | `rep2` | `GSM2684869` | 8,821 |
| `-8` | `hCO_d79_rep2` | `hCO` | `d79` | `rep2` | `GSM2684870` | 9,193 |

Varela DIV30 and DIV90 already have exported Shi Seurat prediction obs tables.
Walsh, Bershteyn, Samarasinghe, and Siebert do not yet have Shi prediction
metadata columns, so they need new Seurat TransferData runs before Scanpy-side
plotting. Xiang is excluded until it is made ready. Shi itself remains the
reference and should not be treated as a cross-study target.
```

Primary metadata columns to use:

| study_id | sample column | cluster column | sample columns present | cluster columns present |
| --- | --- | --- | --- | --- |
| `varela_div30` | `orig.ident` | `seurat_clusters` | `orig.ident` | `seurat_clusters`, `RNA_snn_res.0.2` |
| `varela_div90` | `orig.ident` | `seurat_clusters` | `orig.ident` | `seurat_clusters`, `cluster_id`, `cluster_number_name`, `RNA_snn_res.0.5` |
| `walsh` | `sample_id` | `seurat_clusters` | `sample_id`, `orig.ident` | `seurat_clusters`, `RNA_snn_res.2` |
| `bershteyn_2025` | `sample` | `seurat_clusters` | `sample` | `seurat_clusters`, `predicted.GEcluster`, `predicted.GEtype`, `predicted.GEgws` |
| `bershteyn_2023` | `orig.ident` | `seurat_clusters` | `orig.ident`, `samples` | `seurat_clusters`, `celltype`, `process` |
| `xiang_2018` | **exclude for now** | **exclude for now** | `orig.ident` only in current object | `seurat_clusters`, `RNA_snn_res.0.5`; sample biology requires suffix reconstruction |
| `samarasinghe_2021` | `orig.ident` | `seurat_clusters` | `orig.ident`, `Time`, `Genotype` | `RNA_snn_res.0.3`, `seurat_clusters`, `new.cluster.ids` |
| `siebert_2026` | `orig.ident` | `seurat_clusters` | `orig.ident`, `sample` | `seurat_clusters`, `SCT_snn_res.1`, `SCT_snn_res.0.8` |

Primary sample inventory:

| study_id | samples from primary metadata column |
| --- | --- |
| `varela_div30` | `9583-MW-6` 22,528; `9583-MW-5` 22,230; `9583-MW-3` 15,773; `9583-MW-1` 14,564; `9583-MW-4` 11,004; `9583-MW-2` 4,532 |
| `varela_div90` | `10496-MW-4` 6,314; `10496-MW-6` 4,850; `10496-MW-2` 3,533; `10496-MW-1` 3,095; `10496-MW-3` 2,714; `10496-MW-5` 1,832 |
| `walsh` | `GSM7979671` / `MEL1_dFB_d75` 2,273; `GSM7979672` / `MEL1_vFB_d75` 2,246 |
| `bershteyn_2025` | `010720S` 11,834; `200520S2` 11,123; `010519S1` 10,778; `280120S` 10,776; `010519S2` 10,661; `070120S` 10,537; `220720S1` 10,345; `100620S` 9,885; `220720S2` 8,802; `150120S` 8,294; `251219S` 6,929; `200520S1` 6,756; `111219S` 5,722; `200319S` 2,141 |
| `bershteyn_2023` | `MB279` 10,009; `MS35r41` 9,208; `r41v2ym` 8,722; `MS35mock` 8,403; `mockv2ym` 8,208; `mockv2dw` 8,206; `D0` 8,118; `MB528` 7,127; `r41v2dw` 6,656; `MB460` 6,447; `MB527` 4,933; `D14` 4,851; `MB280` 4,205; `MB461` 2,949 |
| `xiang_2018` | **excluded for now**; current object has only `orig.ident = Xiang2018` 58,950, while usable sample biology requires validated suffix reconstruction |
| `samarasinghe_2021` | `D70_Ctrl_docked_1_seurat` 10,931; `D56_Ctrl_unfused_1_seurat` 9,306; `D56_Rett_unfused_1_seurat` 9,186; `D100_Rett_docked_2_seurat` 7,561; `D100_Ctrl_docked_2_seurat` 6,698; `D70_Rett_docked_1_seurat` 6,260 |
| `siebert_2026` | `Old_1` 16,606; `Young_2` 16,377; `Old_2` 16,073; `Young_1` 15,620 |

Sample alias and metadata notes:

- Walsh `sample_id` values are GEO accessions. Use `MEL1_dFB_d75` as the
  readable alias for `GSM7979671` and `MEL1_vFB_d75` as the readable alias for
  `GSM7979672`. Both are day-75/DIV75 units; the object also carries
  `orig.ident` values `dFB` and `vFB`, plus `domain` values `dFB_domain` and
  `vFB_domain`.
- Siebert 2026 has four primary sample values in both `orig.ident` and
  `sample`: `Old_1`, `Young_2`, `Old_2`, and `Young_1`. Additional informative
  metadata columns include `cellLine`, `donor_batch`, and author annotation
  columns `predictions` and `predictions_class`. Use `orig.ident`/`sample` for
  primary sample-level plots; use `cellLine` or `donor_batch` only when the
  analysis needs line- or donor/batch-level stratification.

Siebert 2026 metadata value inventory was extracted by Slurm job `51479753`
using:

```text
scripts/08_extract_siebert_metadata_values.R
slurm_templates/23_extract_siebert_metadata_values.sbatch.template
```

Output tables:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siebert_2026/metadata_value_inventory
```

Key Siebert 2026 metadata values:

- `cellLine`: `CellLine1` 28,080; `CellLine3` 20,719; `CellLine2` 14,993;
  `CellLine4` 884.
- `donor_batch`: `Old_1_CellLine1` 7,554; `Old_2_CellLine1` 7,429;
  `Young_2_CellLine1` 6,792; `Young_1_CellLine1` 6,305;
  `Young_2_CellLine3` 5,707; `Young_1_CellLine3` 5,484;
  `Old_1_CellLine3` 4,880; `Old_2_CellLine3` 4,648;
  `Old_1_CellLine2` 4,172; `Old_2_CellLine2` 3,996;
  `Young_1_CellLine2` 3,413; `Young_2_CellLine2` 3,412;
  `Young_2_CellLine4` 466; `Young_1_CellLine4` 418.
- `predictions`: 16 author labels. Top labels are `progenitors` 44,680,
  `MGE_IN` 10,200, `CGE_NR2F2/PROX1` 3,775, and `MGE_SST` 1,673.
- `predictions_class`: 21 author labels. Top labels are `IPC_late` 21,578,
  `MGE_IN` 10,899, `Radial_glia` 10,386, `Neuroblast` 7,185, and
  `CGE_NR2F2/PROX1` 3,850.

Siebert 2026 organoid-age interpretation:

- The paper describes single-cell RNA-seq of MGE-enriched GE organoids at
  `DIV51` and `DIV164`, derived from four independent cell lines. Treat
  `Young_*` as the younger organoid-age group and `Old_*` as the older
  organoid-age group. Based on the methods, the likely mapping is
  `Young = DIV51` and `Old = DIV164`.
- The `_1` and `_2` suffixes are best interpreted as sequencing/library
  lane-style batches within each age group, because the methods state that
  cell lines from the same age were pooled and run in two lanes.
- The methods list cell lines `H28126`, `H1`, and `WTCII` at both `DIV51` and
  `DIV164`, and `13234` at `DIV51` only. In the current Seurat object, only
  `CellLine4` is absent from the `Old_*` groups, so `CellLine4` is consistent
  with the DIV51-only line (`13234`). The object does not encode a direct
  mapping from `CellLine1`-`CellLine3` to `H28126`/`H1`/`WTCII`, so keep those
  labels anonymized unless an external demultiplexing key is found.

Primary cluster inventory:

| study_id | cluster summary |
| --- | --- |
| `varela_div30` | `seurat_clusters`, 7 clusters: 0, 1, 2, 3, 4, 6, 7 |
| `varela_div90` | `seurat_clusters`, 13 clusters: 0 through 12 |
| `walsh` | `seurat_clusters`, 24 clusters: 0 through 23 |
| `bershteyn_2025` | `seurat_clusters`, 9 clusters: 0 through 8 |
| `bershteyn_2023` | `seurat_clusters`, 6 clusters: 0 through 5 |
| `xiang_2018` | **excluded for now**; `seurat_clusters`, 28 clusters: 0 through 27 exist, but sample biology/metadata are not ready |
| `samarasinghe_2021` | `seurat_clusters`, 14 clusters: 0 through 13; manual `new.cluster.ids`, 10 labels |
| `siebert_2026` | `seurat_clusters`, 28 clusters: 0 through 27 |

Existing Varela DIV30/DIV90 Shi Seurat prediction tables:

| timepoint | cells | samples | clusters | label score columns | week score columns | predicted label column | predicted week column |
| --- | ---: | ---: | ---: | ---: | ---: | --- | --- |
| DIV30 | 90,631 | 6 | 7 | 11 | 6 | `shi_seurat_full_predicted_shi_label` | `shi_seurat_full_predicted_shi_week_label` |
| DIV90 | 22,338 | 6 | 13 | 11 | 6 | `shi_seurat_full_predicted_shi_label` | `shi_seurat_full_predicted_shi_week_label` |

Current Varela predicted-label inventory:

```text
DIV30:
  MGE 49,251
  progenitor 40,476
  Thalamic neurons 556
  OPC 266
  Excitatory neuron 54
  LGE 25
  Microglia 2
  Excitatory IPC 1

DIV90:
  MGE 17,525
  progenitor 3,247
  LGE 590
  CGE 490
  OPC 269
  Thalamic neurons 180
  Endothelial 20
  Microglia 13
  Excitatory neuron 4
```

What is needed to extend prediction scores across studies:

```text
1. Keep using the Shi reference Seurat object as the reference:
   /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_2019_paper_qc/shi_2019_seurat.rds

2. For each target object, run Seurat FindTransferAnchors/TransferData using
   RNA features shared with Shi. Do not use sample IDs or clusters as prediction
   inputs.

   Do not run Xiang as a target yet. If Xiang is revisited later, first
   validate sample metadata from barcode suffixes and harmonize Ensembl rownames
   to gene symbols using:
   /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/data/raw/xiang_2018_geo_files/suppl/GSE98201_genes.tsv.gz

   Do not include Shi itself as a target object. Shi remains the reference-only
   dataset for these label-transfer comparisons.

   For Samarasinghe, use the official Zenodo processed object:
   /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/samarasinghe_2021_zenodo_processed_object/samarasinghe_2021_zenodo_seurat.rds

3. Export per-cell predicted Shi label, max label score, full per-label score
   matrix, predicted Shi gestational-week label, max week score, and full
   per-week score matrix.

4. Bring the exported tables back into AnnData/Scanpy only for plotting and
   summaries, using the primary sample and cluster columns listed above.
```

Rsync inventory tables from a Mac:

```bash
mkdir -p /Users/ecrespo/Downloads/cross_study_prediction_inventory_v1/tables

rsync -avh --progress \
  elcrespo@greatlakes.arc-ts.umich.edu:/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_prediction_inventory/cross_study_prediction_inventory_v1/tables/ \
  /Users/ecrespo/Downloads/cross_study_prediction_inventory_v1/tables/
```

## HTML Reference For This Target

Local inspiration file from the Mac:

```text
file:///Users/ecrespo/Downloads/ms2024.html
```

Public Dataverse source:

```text
https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/ZQSPKW
```

Dataverse API metadata endpoint:

```text
https://dataverse.harvard.edu/api/datasets/:persistentId/?persistentId=doi:10.7910/DVN/ZQSPKW
```

The API lists one unrestricted HTML file:

```text
filename: ms2024.html
dataFile id: 11713377
contentType: text/html
filesize: 20,925,576 bytes
license: CC0 1.0
```

Download command for the analysis environment:

```bash
curl -L \
  "https://dataverse.harvard.edu/api/access/datafile/11713377" \
  -o /tmp/ms2024.html
```

The rendered Dataverse landing page may require JavaScript/browser
verification, but the API and direct file endpoint are readable from the Linux
workspace.

Use this HTML specifically for the Seurat prediction-score plotting pattern.
Do not copy the biological interpretation directly, because the context differs.
The useful transferable pattern is:

1. Treat Seurat `prediction.score.*` assays as score matrices.
2. Plot selected score rows on UMAP.
3. Melt full score matrices to long tables for density/ridge-style plots.
4. Use metadata `predicted.*` and `predicted.*.score` columns for winner-label
   summaries and confidence distributions.

```text
Observed relevant pattern in ms2024.html:
  prediction.score.adSupercluster
  prediction.score.adROIGroup
  prediction.score.mtg_type
  FeaturePlot(..., features = rownames(prediction.score.* assay), ...)
  reshape2::melt(t(GetAssayData(prediction.score.* assay)))
  FetchData(..., vars = c("predicted.*.score", "predicted.*", "mpt"))
```

The Seurat route should preserve the same core rule: DIV30 `seurat_clusters`
must not be used as prediction input. They are only used after prediction for
cluster summaries, stacked bars, heatmaps, and interpretation.

Candidate Seurat-side workflow for this target:

```r
anchors <- FindTransferAnchors(
  reference = shi_reference_seurat,
  query = div30_query_seurat,
  normalization.method = "LogNormalize",
  reference.assay = "RNA",
  query.assay = "RNA",
  dims = 1:50
)

label_predictions <- TransferData(
  anchorset = anchors,
  refdata = shi_reference_seurat$shi_label,
  dims = 1:50
)

div30_query_seurat <- AddMetaData(div30_query_seurat, label_predictions)
```

Expected Seurat prediction fields to export:

```text
cell_id
predicted.id or predicted.<transfer_name>
prediction.score.max or predicted.<transfer_name>.score
prediction.score.<label> columns or prediction.score.<transfer_name> assay rows
```

Suggested AnnData column names after import:

```text
shi_seurat_full_predicted_shi_label
shi_seurat_full_prediction_score
shi_seurat_full_uncertainty_score
shi_seurat_full_prediction_score_MGE
shi_seurat_full_prediction_score_LGE
shi_seurat_full_prediction_score_CGE
...
```

For gestational-week transfer, run a second `TransferData` call with
`refdata = shi_reference_seurat$week_label` and import analogous columns:

```text
shi_seurat_full_predicted_shi_week_label
shi_seurat_full_week_prediction_score
shi_seurat_full_week_uncertainty_score
```

After import, Scanpy plotting should remain unchanged in spirit: color the DIV30
Seurat UMAP by `shi_seurat_full_predicted_shi_label`,
`shi_seurat_full_prediction_score`, individual
`shi_seurat_full_prediction_score_*` label-score columns, broad/developmental
classes derived from the predicted label, and later week prediction columns.
Cluster-level summaries should use DIV30 `seurat_clusters` only after these
per-cell predictions have been joined.

## Legacy v2 Output Contract

This output contract belongs to the completed Scanpy-side
`shi_reference_div30_label_transfer_v2` run. For the next Seurat score target,
use the expected output roots and files listed in
`Next Target: Seurat Prediction Scores`.

Run outputs:

```text
RUN_DIR/tables/
RUN_DIR/plots/
RUN_DIR/h5ad/
```

For `shi_reference_div30_label_transfer_v2`, the completed run contains:

```text
plots/*.png: 14 files
tables/*: 34 files
```

Key tables:

```text
tables/shi_reference_table_s2_label_join_summary.tsv
tables/shi_reference_table_s2_label_counts.tsv
tables/shi_reference_comparison_label_sets.tsv
tables/shi_label_transfer_gene_harmonization_summary.tsv
tables/shi_label_transfer_gene_harmonization_detail.tsv
tables/div30_shi_label_transfer_obs.tsv.gz
tables/div30_shi_label_transfer_predictions_long.tsv.gz
tables/div30_shi_label_transfer_cluster_summaries.tsv
tables/div30_shi_label_transfer_label_fractions_by_cluster.tsv
tables/div30_shi_full_vs_mge_lge_cge_predictions.tsv.gz
tables/div30_shi_full_vs_mge_lge_cge_summary_by_seurat_clusters.tsv
tables/div30_shi_week_prediction_counts_by_sample.tsv
tables/div30_shi_week_prediction_summary_by_sample.tsv
tables/shi_reference_week_counts.tsv
tables/shi_reference_week_metadata.tsv
tables/shi_label_transfer_output_manifest.tsv
tables/shi_label_transfer_complete.tsv
```

Annotated AnnData:

```text
h5ad/div30_shi_label_transfer_predictions.h5ad
```

Prediction columns added to `adata.obs`:

```text
shi_full_predicted_shi_label
shi_full_prediction_score
shi_full_uncertainty_score
shi_full_prediction_entropy
shi_full_broad_region_class
shi_full_developmental_class
shi_full_predicted_shi_week_label
shi_full_week_prediction_score
shi_full_week_uncertainty_score
shi_full_week_prediction_entropy
shi_full_mean_neighbor_week_numeric
shi_full_median_neighbor_week_numeric
shi_full_std_neighbor_week_numeric

shi_mge_lge_cge_predicted_shi_label
shi_mge_lge_cge_prediction_score
shi_mge_lge_cge_uncertainty_score
shi_mge_lge_cge_prediction_entropy
shi_mge_lge_cge_broad_region_class
shi_mge_lge_cge_developmental_class
shi_mge_lge_cge_predicted_shi_week_label
shi_mge_lge_cge_week_prediction_score
shi_mge_lge_cge_week_uncertainty_score
shi_mge_lge_cge_week_prediction_entropy
shi_mge_lge_cge_mean_neighbor_week_numeric
shi_mge_lge_cge_median_neighbor_week_numeric
shi_mge_lge_cge_std_neighbor_week_numeric
```

## Legacy v2 Plot Inventory

`shi_reference_div30_label_transfer_v2` generated exactly 14 PNG plots:
7 plots for `full_relevant` and 7 plots for `mge_lge_cge_only`.

Run plot directory:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_reference_div30_label_transfer/shi_reference_div30_label_transfer_v2/plots
```

Shared plot inputs:

```text
DIV30 query AnnData:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/varela_div30.h5ad
  UMAP coordinates from .obsm["X_umap_seurat"]
  sample groups from .obs["orig.ident"]
  post hoc cluster groups from .obs["seurat_clusters"]

Shi reference AnnData:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/shi_2019_paper_qc.h5ad

Shi Table S2:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/reference/shi_2021_tables_s2_to_s9/science.abj6641_table_s2.xlsx
  cell-level Major types and gestational-week-derived labels
```

Plot helper functions are defined in:

```text
python_notebooks/src/mge_organoid_python/shi_label_transfer.py
```

Full-reference plots:

```text
1. div30_umap_shi_full_predicted_shi_label.png
   Input: query .obsm["X_umap_seurat"] + query.obs["shi_full_predicted_shi_label"]
   Function: plot_umap_categorical

2. div30_umap_shi_full_prediction_score.png
   Input: query .obsm["X_umap_seurat"] + query.obs["shi_full_prediction_score"]
   Function: plot_umap_continuous

3. div30_umap_shi_full_broad_region_class.png
   Input: query .obsm["X_umap_seurat"] + query.obs["shi_full_broad_region_class"]
   Function: plot_umap_categorical

4. div30_shi_full_shi_label_stacked_bar_by_seurat_clusters.png
   Input table: tables/div30_shi_full_label_counts.tsv
   Source columns: seurat_clusters, shi_full_predicted_shi_label, fraction
   Function: plot_stacked_bar

5. div30_shi_full_shi_label_heatmap_by_seurat_clusters.png
   Input table: tables/div30_shi_full_label_counts.tsv
   Source columns: seurat_clusters, shi_full_predicted_shi_label, fraction
   Function: plot_heatmap

6. div30_shi_full_mean_neighbor_gw_age_density_by_sample.png
   Input: query.obs["orig.ident"] + query.obs["shi_full_mean_neighbor_week_numeric"]
   Function: plot_overlaid_density_by_group

7. div30_shi_full_gw_prediction_score_density_by_sample.png
   Input: query.obs["orig.ident"] + query.obs["shi_full_week_prediction_score"]
   Function: plot_overlaid_density_by_group
```

MGE/LGE/CGE-restricted plots:

```text
8. div30_umap_shi_mge_lge_cge_predicted_shi_label.png
   Input: query .obsm["X_umap_seurat"] + query.obs["shi_mge_lge_cge_predicted_shi_label"]
   Function: plot_umap_categorical

9. div30_umap_shi_mge_lge_cge_prediction_score.png
   Input: query .obsm["X_umap_seurat"] + query.obs["shi_mge_lge_cge_prediction_score"]
   Function: plot_umap_continuous

10. div30_umap_shi_mge_lge_cge_broad_region_class.png
    Input: query .obsm["X_umap_seurat"] + query.obs["shi_mge_lge_cge_broad_region_class"]
    Function: plot_umap_categorical

11. div30_shi_mge_lge_cge_shi_label_stacked_bar_by_seurat_clusters.png
    Input table: tables/div30_shi_mge_lge_cge_label_counts.tsv
    Source columns: seurat_clusters, shi_mge_lge_cge_predicted_shi_label, fraction
    Function: plot_stacked_bar

12. div30_shi_mge_lge_cge_shi_label_heatmap_by_seurat_clusters.png
    Input table: tables/div30_shi_mge_lge_cge_label_counts.tsv
    Source columns: seurat_clusters, shi_mge_lge_cge_predicted_shi_label, fraction
    Function: plot_heatmap

13. div30_shi_mge_lge_cge_mean_neighbor_gw_age_density_by_sample.png
    Input: query.obs["orig.ident"] + query.obs["shi_mge_lge_cge_mean_neighbor_week_numeric"]
    Function: plot_overlaid_density_by_group

14. div30_shi_mge_lge_cge_gw_prediction_score_density_by_sample.png
    Input: query.obs["orig.ident"] + query.obs["shi_mge_lge_cge_week_prediction_score"]
    Function: plot_overlaid_density_by_group
```

Notes for tomorrow:

```text
tables/shi_reference_comparison_label_sets.tsv confirms full_relevant has all
10 Shi labels, and mge_lge_cge_only has MGE, LGE, and CGE.

tables/div30_shi_full_label_counts.tsv keeps zero-count label rows, so CGE,
Microglia, OPC, and Endothelial do not silently disappear when their predicted
fraction is 0.

The overlaid GW-age plots are not based on DIV30 collection age directly. They
are inferred from the Shi reference neighbors for each DIV30 cell and then
overlaid by DIV30 sample (`orig.ident`).
```

## Slurm Command

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
cp slurm_templates/20_execute_shi_reference_div30_label_transfer.sbatch.template \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/20_execute_shi_reference_div30_label_transfer.sbatch
sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/20_execute_shi_reference_div30_label_transfer.sbatch
```

## Rsync From Mac

From `/Users/ecrespo/Downloads`:

```bash
mkdir -p ./shi_reference_div30_label_transfer_v2

rsync -avh --progress \
  elcrespo@greatlakes.arc-ts.umich.edu:/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_reference_div30_label_transfer/shi_reference_div30_label_transfer_v2/plots/ \
  ./shi_reference_div30_label_transfer_v2/plots/

rsync -avh --progress \
  elcrespo@greatlakes.arc-ts.umich.edu:/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_reference_div30_label_transfer/shi_reference_div30_label_transfer_v2/tables/ \
  ./shi_reference_div30_label_transfer_v2/tables/
```

## Slurm Status

Latest corrected run:

```text
51359356:
  Run label: shi_reference_div30_label_transfer_v2
  State: COMPLETED
  ExitCode: 0:0
  Elapsed: 00:06:22
  MaxRSS: 24756096K
  Node: gl3076
  Change from v1:
    - full_relevant now includes all 10 Shi Table S2 major-type labels
    - zero-count labels are retained in label-count tables and plot legends
    - added Shi gestational-week prediction columns and sample-overlaid density plots
  Output log: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/shi-label-transfer-shi-div30-xfer-51359356.out
  Error log: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/shi-label-transfer-shi-div30-xfer-51359356.err
```

Setup/debug attempt:

```text
51358066:
  State: FAILED
  Reason: Table S2 join initially matched only 9.8% of Shi reference cells.
  Fix: updated label join to infer Table S2 numeric sample suffixes and use
       per-cell unique barcode keys before suffix-level fallback.
```

Corrected submitted run:

```text
51358138:
  State: COMPLETED
  ExitCode: 0:0
  Elapsed: 00:06:14
  MaxRSS: 24623288K
  Node: gl3118
  CPUs: 8
  Memory: 160G
  Output log: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/shi-label-transfer-shi-div30-xfer-51358138.out
  Error log: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/shi-label-transfer-shi-div30-xfer-51358138.err
  Pre-submit reference-label join check: 55,704 / 56,136 matched (99.23%)
  Completion marker: present
```

## 2026-06-07 Cross-Study Shi Seurat Transfer Status

New workflow target:

```text
Run label:
  cross_study_shi_seurat_label_transfer_v1

Output root:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_shi_seurat_label_transfer/cross_study_shi_seurat_label_transfer_v1
```

Code added/updated for this cross-study workflow:

```text
scripts/13_run_cross_study_shi_seurat_label_transfer.R
python_notebooks/src/mge_organoid_python/cross_study_shi_prediction_plots.py
python_notebooks/scripts/run_cross_study_shi_prediction_plots.py
python_notebooks/notebooks/cross_study_shi_prediction_scores.ipynb
slurm_templates/26_cross_study_shi_seurat_label_transfer.sbatch.template
slurm_templates/27_cross_study_shi_seurat_label_transfer_array.sbatch.template
slurm_templates/28_finalize_cross_study_shi_prediction_plots.sbatch.template
```

Important implementation details:

```text
- Shi is reference only.
- Xiang is excluded.
- Varela DIV30 and DIV90 reuse existing Seurat TransferData score exports.
- Reused Varela tables do not carry UMAP columns, so Python augments UMAP/sample/cluster
  metadata from the existing Varela H5AD caches during combine/plot.
- Non-Varela targets run Seurat FindTransferAnchors + TransferData from the Shi RDS.
- The Shi RDS has week_label but not a major shi_label column, so the R workflow attaches
  major labels from the validated existing TSV:
    /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_reference_div30_seurat_label_transfer/shi_reference_div30_seurat_label_transfer_v1/seurat/shi_reference_labels_for_seurat.tsv
- Prediction inputs remain RNA expression only. Samples, clusters, author labels, and UMAP
  coordinates are not used for prediction.
- One Seurat parameter required explicit handling during cross-study transfer:
  `TransferData(k.weight)`. The initial Siebert 2026 serial run produced too
  few anchors for Seurat's default `k.weight = 50`, so the workflow now uses an
  adaptive guard:
    n_anchors >= 100 -> k.weight = 50
    n_anchors < 100  -> k.weight = min(4, n_anchors - 1)
- This `k.weight` change controls the number of anchors used for Seurat's
  weighting step. It is not a prediction-score cutoff and does not invoke
  `GetTransferPredictions(score.filter = 0.75)`.
```

Serial job attempts:

```text
51486655:
  State: FAILED
  Partition: largemem
  Reason:
    siebert_2026 produced only 25 anchors.
    TransferData failed because default k.weight was larger than the anchor count:
      "Please set k.weight to be smaller than the number of anchors (25)."

51486702:
  State: FAILED
  Partition: largemem
  Change:
    set k.weight = min(50, n_anchors - 1).
  Reason:
    siebert_2026 still failed.
    Seurat reported:
      "Number of anchor cells is less than k.weight. Consider lowering k.weight to less than 5..."

51486807:
  State: CANCELLED
  Reason:
    cancelled manually after deciding to switch from one serial all-study job to a per-study array.
```

Parallel array jobs:

```text
Array job:
  51486835

Final combine/plot job:
  51486836
  Dependency:
    afterok:51486835

Array mapping:
  51486835_1  siebert_2026
  51486835_2  walsh
  51486835_3  bershteyn_2025
  51486835_4  bershteyn_2023
  51486835_5  samarasinghe_2021
```

Status as of 2026-06-07 22:09 EDT:

```text
Completed per-study outputs:
  varela_div30  reused existing Seurat predictions
  varela_div90  reused existing Seurat predictions
  walsh         completed new Seurat transfer successfully

Running:
  51486835_3  bershteyn_2025
  51486835_4  bershteyn_2023

Pending:
  51486835_5  samarasinghe_2021
  51486836    final combine/plot job, currently held by dependency

Failed:
  51486835_1  siebert_2026
```

Current per-study obs tables present:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_shi_seurat_label_transfer/cross_study_shi_seurat_label_transfer_v1/tables/per_study/varela_div30_shi_seurat_label_transfer_obs.tsv.gz
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_shi_seurat_label_transfer/cross_study_shi_seurat_label_transfer_v1/tables/per_study/varela_div90_shi_seurat_label_transfer_obs.tsv.gz
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_shi_seurat_label_transfer/cross_study_shi_seurat_label_transfer_v1/tables/per_study/walsh_shi_seurat_label_transfer_obs.tsv.gz
```

Walsh completed details:

```text
Shared features: 14,945
Anchors: 3,183
TransferData k.weight: 50
Warning:
  zero predicted cells for Excitatory IPC, Microglia, Endothelial
Interpretation:
  This is not a runtime failure. Keep zero-count labels in summaries/legends.
```

Siebert failure details:

```text
Shared features: 15,305
Anchors: 25
TransferData k.weight after patch: 4
Failure:
  Error in data.frame(predicted.id = prediction.ids, prediction.score = as.matrix(prediction.scores),  :
    arguments imply differing number of rows: 0, 64676

Interpretation:
  Siebert is not just a default-k.weight problem. The low/poor anchor set is causing
  Seurat TransferData to return no predicted IDs while retaining a score matrix with
  64,676 query rows. Treat siebert_2026 as a special failed target for now.
```

Recommended next actions:

```text
1. Let 51486835_3, 51486835_4, and 51486835_5 finish.
2. Because 51486835_1 failed, 51486836 will not run automatically under afterok.
3. After the successful array tasks finish, manually run the Python final combine/plot
   for available successful studies, excluding siebert_2026 unless it is debugged.
4. Debug siebert_2026 separately. Possible directions:
     - inspect RNA assay/layers and feature naming in the Siebert object
     - try fewer dims or a different anchor reduction strategy
     - verify whether its RNA data layer is sparse/empty/odd after JoinLayers
     - consider excluding Siebert from this first cross-study Seurat-transfer figure
       if anchors remain too weak.
```

Status/check commands:

```bash
squeue -j 51486835,51486836 -o '%.20i %.9P %.30j %.8u %.2t %.10M %.6D %R'

sacct -j 51486835,51486836 \
  --format=JobID,JobName%30,State,ExitCode,Elapsed,MaxRSS,ReqMem -P

for i in 1 2 3 4 5; do
  echo "==== task ${i} ===="
  tail -n 80 /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/cross-shi-xfer-array-51486835_${i}.err
done

find /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_shi_seurat_label_transfer/cross_study_shi_seurat_label_transfer_v1/tables/per_study \
  -maxdepth 1 -type f -name '*_shi_seurat_label_transfer_obs.tsv.gz' -printf '%f\t%TY-%Tm-%Td %TH:%TM\t%s\n' | sort
```

## 2026-06-08 Siebert Cross-Study Shi Transfer Fix

Supersedes the 2026-06-07 note that treated `siebert_2026` as failed/excluded.

Root cause:

```text
A. Low anchors:
   Siebert RNA counts were populated, but RNA data/logcounts were empty:
     counts: 32131 x 64676, nnzero 126244258
     data:   0 x 0,       nnzero 0

   The transfer helper accepted the empty data layer as usable, so NormalizeData()
   was skipped. After treating 0 x 0 layers as missing and normalizing counts,
   Siebert anchors increased from 25 to 4089.

B. 0 predicted IDs:
   The 25-anchor run had valid anchor labels, so this was not solved by only
   naming refdata with a "_reference" suffix. The degenerate low-anchor run
   produced invalid/NA TransferData weights/scores, causing Seurat to construct
   zero predicted IDs for 64676 query rows.
```

Code patch:

```text
scripts/13_run_cross_study_shi_seurat_label_transfer.R

get_data_layer() now ignores empty or shape-mismatched data/logcounts/counts
layers. A layer must match assay dimensions and have finite positive nnzero
before it is considered usable.
```

Successful Siebert rerun:

```text
Slurm job: 51515010
State: COMPLETED
Elapsed: 00:10:09
MaxRSS: 48001428K

Log:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/siebert-xfer-fixed-51515010.err

Key log lines:
  Ignoring empty or shape-mismatched data layer for assay RNA (layer dim=0x0, assay dim=32131x64676, nnzero=0)
  Ignoring empty or shape-mismatched logcounts layer for assay RNA (layer dim=0x0, assay dim=32131x64676, nnzero=0)
  No data/logcounts layer found for assay RNA; running NormalizeData
  Found 4089 anchors
  TransferData k.weight for siebert_2026: 50 (anchors=4089)
```

Siebert output now present:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_shi_seurat_label_transfer/cross_study_shi_seurat_label_transfer_v1/tables/per_study/siebert_2026_shi_seurat_label_transfer_obs.tsv.gz
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_shi_seurat_label_transfer/cross_study_shi_seurat_label_transfer_v1/seurat/per_study/siebert_2026_shi_seurat_full_transfer_diagnostics.tsv
```

Siebert diagnostics:

```text
n_query_cells: 64676
n_reference_cells: 38831
n_shared_features: 15305
n_anchors: 4089
transfer_k_weight: 50
n_missing_prediction: 0
Warning: zero predicted cells for Shi labels Microglia and Endothelial.
```

Combined cross-study table was regenerated explicitly for all seven target
studies:

```text
varela_div30
varela_div90
siebert_2026
walsh
bershteyn_2025
bershteyn_2023
samarasinghe_2021

Combined rows: 454731
```

All five UMAP grid plot files were regenerated and verified to include all seven
studies:

```text
plots/umap_grids/cross_study_umap_shi_seurat_full_predicted_shi_label_grid.png
plots/umap_grids/cross_study_umap_shi_seurat_full_label_score_grid.png
plots/umap_grids/cross_study_umap_shi_seurat_full_all_label_scores_grid.png
plots/umap_grids/cross_study_umap_shi_seurat_full_gw_prediction_grid.png
plots/umap_grids/cross_study_umap_shi_seurat_full_individual_gw_scores_grid.png
```

Important plotting note:

```text
The all-in-one Python plotting command refreshed the UMAP grids, but the
sample-level summary-plot phase was cancelled after hanging. The existing
summary plots and cross_study_shi_seurat_label_transfer_complete.tsv should be
treated as stale until summary plotting is split or optimized. The UMAP grid
manifest is current and includes all seven studies.
```

## 2026-06-09 Shi TransferData Score-Schema Refactor

Reason:

```text
The active workflows needed one enforced Shi prediction-score schema. The only
active prediction-score source should be Seurat TransferData anchor-weighted
support scores, not the legacy Scanpy/kNN vote fraction. This may affect any
downstream code that expects raw `prediction.score.<label>` score-matrix column
names or stale tables where max score/uncertainty do not match the canonical
per-label score matrix.
```

Code changes:

```text
scripts/13_run_cross_study_shi_seurat_label_transfer.R
  - validates TransferData output immediately after major-label and GW transfer
  - exports canonical `shi_seurat_full_prediction_score_<label_token>` columns
  - validates max-score and uncertainty identities before writing/plotting
  - logs k.weight request/use/reason, anchors, unique query anchor cells, median
    max score, fraction with max score >= 0.75, and winner-take-all composition
  - adds per-cell QC columns: score delta, low-confidence flag, anchor/k.weight QC

python_notebooks/scripts/seurat_shi_label_transfer_export.R
  - uses the same canonical score-matrix naming and identity checks for the
    standalone/smoke Seurat export path

python_notebooks/src/mge_organoid_python/shi_prediction_schema.py
  - centralizes Python-side canonical column names, label-token sanitization,
    and score identity validation

python_notebooks/src/mge_organoid_python/cross_study_shi_prediction_plots.py
python_notebooks/scripts/run_shi_seurat_label_transfer_smoke.py
python_notebooks/src/mge_organoid_python/shi_label_transfer.py
  - now reuse the shared token/validation helper
  - legacy `run_knn_label_transfer()` still raises unless
    `allow_legacy_knn_scoring=True`
```

Local verification before Slurm submission:

```text
/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/Rscript parse check: OK
/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python -m py_compile: OK
git diff --check: OK
Legacy kNN guard check: OK
Repo search:
  no active GetTransferPredictions
  no score.filter
  no allow_legacy_knn_scoring=True outside the guard message
```

Submission plan:

```text
Submit the cross-study array workflow rather than the single serial template.
The array template force-reruns the five non-Varela target transfers:
  siebert_2026
  walsh
  bershteyn_2025
  bershteyn_2023
  samarasinghe_2021

Then submit the finalizer with afterok dependency on the array job. The reused
Varela DIV30/DIV90 tables are still validated during final combine/plot, so this
run should expose both fresh-transfer and reused-table schema implications.
```

Submitted Slurm jobs:

```text
Array job:
  51559750

Finalizer job:
  51559751
  Dependency:
    afterok:51559750

Immediate status after submission:
  51559750_[1-5]  PENDING on largemem
  51559751        PENDING on standard, Dependency

Prepared job files:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/27_cross_study_shi_seurat_label_transfer_array.sbatch
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/28_finalize_cross_study_shi_prediction_plots.sbatch
```

Monitor:

```bash
squeue -j 51559750,51559751 -o '%.20i %.9P %.30j %.8u %.2t %.10M %.6D %R'

sacct -j 51559750,51559751 \
  --format=JobID,JobName%30,State,ExitCode,Elapsed,MaxRSS,ReqMem -P

for i in 1 2 3 4 5; do
  tail -n 80 /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/cross-shi-xfer-array-51559750_${i}.err
done

tail -n 80 /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/cross-shi-final-51559751.err
```

Follow-up status:

```text
Transfer array 51559750:
  All five non-Varela transfer tasks completed successfully.

Finalizer 51559751:
  FAILED quickly, ExitCode 1:0.
  Cause:
    New Python score validation had a NumPy broadcasting bug when checking
    one max-score vector against the per-label score matrix:
      operands could not be broadcast together with shapes (454731,) (454731,10)

Patch:
  python_notebooks/src/mge_organoid_python/shi_prediction_schema.py
  now checks max-score range and per-label score-matrix range separately.

Important execution note:
  A direct local `combine` command was run once after the failed finalizer to
  validate the patch against the full table and completed successfully:
    combined_obs_rows 454731
  For the official handoff-style run, continue using Slurm jobs.

Official Slurm finalizer rerun after schema patch:
  51638316
  Job file:
    /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/28_finalize_cross_study_shi_prediction_plots.schema_fix_rerun.sbatch
  Immediate status:
    PENDING on standard
```

Monitor rerun:

```bash
squeue -j 51638316 -o '%.20i %.9P %.35j %.8u %.2t %.10M %.6D %R'
sacct -j 51638316 --format=JobID,JobName%35,State,ExitCode,Elapsed,MaxRSS,ReqMem -P
tail -n 120 /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/cross-shi-final-51638316.err
tail -n 120 /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/cross-shi-final-51638316.out
```

## Shi Threshold/Stage Composition Plot Reruns

Dedicated Slurm template:

```text
slurm_templates/30_shi_threshold_tradeoff_plots.sbatch.template
```

This regenerates `python_notebooks/scripts/plot_shi_mge_threshold_tradeoff.py`
through Slurm, including:

```text
winner-take-all versus score cutoff summaries
MGE-only, MGE/LGE/CGE, and all Shi major-label score scopes
predicted-stage sample-composition stacked bars
combined all-study predicted-stage sample-composition stacked bar
```

Important sample-label/order note:

```text
The transfer obs tables carry raw sample IDs for Bershteyn 2023. The plotting
script now relabels and reorders those samples from GSE208672 metadata:
  DIV0 hESC
  DIV14 MGE progenitor
  DIV42 EOP lots 1..3 as paired U/S samples
  additional GSE DIV42 EOP batches 4..9 kept visible as U samples

GSE208672 gives explicit D0/D14 annotations and marks the remaining samples as
"End of Process"; study-design notes indicate these EOP cells are week 6, which
is plotted as DIV42. S/U denotes sorted/unsorted. They are ordered after D14 by
lot/batch rather than by cell-count abundance.

Bershteyn 2025 samples from GSE283775 are also "End of Process" sorted batches.
No explicit DIV was found in that GEO metadata, so those plot labels are marked
"DIV unknown EOP" and ordered batch 1..14 rather than assigned an invented DIV.
```

The predicted-stage stacked bars use the fixed Shi reference-panel palette:

```text
GW09 #231611
GW12 #3F1C6A
GW13 #A02E6B
GW16 #EB5840
GW18 #FCC031
```

Submit:

```bash
cp slurm_templates/30_shi_threshold_tradeoff_plots.sbatch.template \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/30_shi_threshold_tradeoff_plots.sbatch
sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/30_shi_threshold_tradeoff_plots.sbatch
```

Completed Slurm reruns after adding Bershteyn sample relabeling/order and the
combined all-study predicted-stage plot:

```text
51641793  COMPLETED  first generated combined all-study plot
51641821  COMPLETED  corrected combined-plot title/legend/study-label spacing
51642033  COMPLETED  converted Bershteyn 2023 week-6 labels to DIV42 and S/U
51643511  COMPLETED  repeated predicted-age plots for all Shi, MGE, MGE/LGE/CGE
```

Current display/order choices:

```text
Varela samples are labeled as "Varela et al. / this paper / DIV30" and
"Varela et al. / this paper / DIV90".

The combined predicted-stage composition plot leaves an explicit x-axis gap
between studies and draws the study name/year above each group.

Bershteyn 2023 is plotted before Bershteyn 2025.

Predicted-age composition plots are now repeated for three winner-take-all
major-label cell sets:
  all_shi_major_labels
  mge
  mge_lge_cge

Those outputs live in subfolders under:
  plots/summary/predicted_age_sample_composition/
```

## GE-Only Age Classifier Addition

Major classifier logic change:

```text
Run label:
  cross_study_shi_seurat_label_transfer_v2_ge_only_age

Whole-Shi transfer outputs remain in the existing namespace:
  shi_seurat_full_predicted_shi_label
  shi_seurat_full_prediction_score*
  shi_seurat_full_predicted_shi_week_label
  shi_seurat_full_week_prediction_score*

New GE-only age classifier outputs use a separate namespace:
  shi_seurat_ge_only_predicted_shi_week_label
  shi_seurat_ge_only_week_prediction_score
  shi_seurat_ge_only_week_uncertainty_score
  shi_seurat_ge_only_week_prediction_score_*
  shi_seurat_ge_only_expected_shi_gw_numeric
  shi_seurat_ge_only_expected_shi_gw_even
```

Classifier inventory for this run:

```text
1. Whole-Shi major-label classifier
   reference classes: all Shi major cell-type labels
   target columns: shi_seurat_full_predicted_shi_label,
                   shi_seurat_full_prediction_score*

2. Whole-Shi week/age classifier
   reference classes: Shi GW09, GW12, GW13, GW16, GW18 across the full
   reference object
   target columns: shi_seurat_full_predicted_shi_week_label,
                   shi_seurat_full_week_prediction_score*

3. GE-only week/age classifier
   reference classes: Shi GW09, GW12, GW13, GW16, GW18 after subsetting the
   reference to major labels MGE, LGE, and CGE only
   target columns: shi_seurat_ge_only_predicted_shi_week_label,
                   shi_seurat_ge_only_week_prediction_score*
```

Computation:

```text
The R transfer script now creates a Shi reference subset containing only
reference cells whose major label is MGE, LGE, or CGE. It reruns
FindTransferAnchors/TransferData for Shi gestational-week labels against this
GE-only reference. This is a distinct classifier, not a post-hoc subset of the
whole-Shi week scores.

The original major-label classifier is still whole-Shi and still provides the
denominator filters:
  all_shi_major_labels
  mge
  mge_lge_cge

Predicted-age sample-composition plotting now uses:
  all_shi_major_labels -> whole-Shi week classifier
  mge                  -> GE-only week classifier
  mge_lge_cge          -> GE-only week classifier

The 0.90 cutoff remains a separate major-label support analysis. It is not used
to threshold week/age calls in the stacked age-composition plots.
```

Important PCA/anchor caveat:

```text
The current v2 GE-only implementation recomputes PCA after subsetting the Shi
reference to MGE/LGE/CGE cells. That was done to make the GE-only week
classifier a fully separate Seurat reference-space classifier:

  subset reference to MGE/LGE/CGE cells
  recompute ScaleData/RunPCA on that subset
  recompute FindTransferAnchors from GE-only reference to each query
  TransferData only Shi GW/week labels from those GE-only anchors

This is not a minor downstream filter. Recomputing PCA on the GE subset changes
the reference coordinate system used to find anchors. Axes that separated
whole-Shi broad classes are removed/reweighted, and GE-internal variation can
dominate the anchor space. That can improve a GE-focused age classifier in
principle, but it can also make results worse if the whole-Shi PCA provided a
more stable or better-calibrated cross-study alignment.

Current interpretation of v2:
  shi_seurat_ge_only_* columns answer:
    "What week does this query cell resemble when the reference and anchors are
     rebuilt using only Shi MGE/LGE/CGE cells?"

They do not answer:
    "Among cells called MGE/LGE/CGE by the old classifier, what is the old
     whole-Shi week score?"

If the GE-only PCA makes the biology look worse, the next candidate logic should
be documented as a separate v3 mode, for example:
  - keep the whole-Shi anchor/PCA space and only change which cells are included
    in age summaries; or
  - test a GE-cell reference subset without recomputing a GE-specific PCA, if
    Seurat supports the intended reduction cleanly; or
  - use whole-Shi week predictions for all age plots and treat MGE/MGE-LGE-CGE
    only as denominator filters.
```

Current function/file map:

```text
Main computation file:
  scripts/13_run_cross_study_shi_seurat_label_transfer.R

Reference setup:
  main()
    - loads Shi reference RDS
    - attaches/chooses major-label metadata column
    - chooses GW/week metadata column
    - filters reference cells with both major label and week
    - writes:
        diagnostics/shi_reference_labels_used_by_seurat.tsv
        diagnostics/shi_reference_weeks_used_by_seurat.tsv
    - creates ge_reference by:
        ge_reference_cells <- cells with major label in MGE/LGE/CGE
        ge_reference <- subset(reference, cells = ge_reference_cells)
    - writes:
        diagnostics/shi_ge_only_reference_labels_by_week_used_by_seurat.tsv

Per-study transfer:
  run_transfer_one(study, reference, ge_reference, ...)
    - loads target Seurat object
    - joins Assay5 RNA layers if needed
    - ensures reference/query RNA data layer exists

Whole-Shi anchor path:
  run_transfer_one()
    - shared_features <- intersect(full reference genes, query genes)
    - ScaleData(reference, features = shared_features)
    - RunPCA(reference, features = shared_features)
    - FindTransferAnchors(reference = reference, query = query, ...)
    - TransferData(anchorset = anchors, refdata = major labels)
    - TransferData(anchorset = anchors, refdata = GW/week labels)
    - writes shi_seurat_full_* columns

GE-only anchor path:
  run_transfer_one()
    - ge_shared_features <- intersect(GE-only reference genes, query genes)
    - ScaleData(ge_reference, features = ge_shared_features)
    - RunPCA(ge_reference, features = ge_shared_features)
    - FindTransferAnchors(reference = ge_reference, query = query, ...)
    - TransferData(anchorset = ge_anchors, refdata = GW/week labels)
    - writes shi_seurat_ge_only_* columns

Final table/UMAP plotting:
  python_notebooks/src/mge_organoid_python/cross_study_shi_prediction_plots.py
  python_notebooks/scripts/run_cross_study_shi_prediction_plots.py

Threshold and age-composition plotting:
  python_notebooks/scripts/plot_shi_mge_threshold_tradeoff.py
```

Workflow map:

```text
Step 1: Per-study Seurat label transfer
  code:
    scripts/13_run_cross_study_shi_seurat_label_transfer.R
  Slurm:
    slurm_templates/27_cross_study_shi_seurat_label_transfer_array.sbatch.template
    slurm_templates/26_cross_study_shi_seurat_label_transfer.sbatch.template
  what runs:
    - whole-Shi major-label FindTransferAnchors/TransferData
    - whole-Shi GW/week TransferData on the same whole-Shi anchors
    - GE-only FindTransferAnchors/TransferData after subsetting the Shi
      reference to MGE/LGE/CGE cells
  outputs:
    results/cross_study_shi_seurat_label_transfer/
      cross_study_shi_seurat_label_transfer_v2_ge_only_age/
        tables/per_study/*_shi_seurat_label_transfer_obs.tsv.gz
        seurat/per_study/*_shi_seurat_full_predictions.tsv.gz
        seurat/per_study/*_shi_seurat_full_prediction_scores.tsv.gz
        seurat/per_study/*_shi_seurat_full_week_predictions.tsv.gz
        seurat/per_study/*_shi_seurat_full_week_prediction_scores.tsv.gz
        seurat/per_study/*_shi_seurat_ge_only_week_predictions.tsv.gz
        seurat/per_study/*_shi_seurat_ge_only_week_prediction_scores.tsv.gz
        diagnostics/shi_ge_only_reference_labels_by_week_used_by_seurat.tsv

Step 2: Combine/validate tables and draw cross-study UMAP grids
  code:
    python_notebooks/src/mge_organoid_python/cross_study_shi_prediction_plots.py
    python_notebooks/scripts/run_cross_study_shi_prediction_plots.py
  Slurm:
    slurm_templates/28_finalize_cross_study_shi_prediction_plots.sbatch.template
  what runs:
    - reads all per-study obs tables
    - validates canonical max score and uncertainty for whole-Shi major labels,
      whole-Shi weeks, and GE-only weeks
    - writes combined/long score tables
    - renders cross-study UMAP grids
  outputs:
    results/cross_study_shi_seurat_label_transfer/
      cross_study_shi_seurat_label_transfer_v2_ge_only_age/
        tables/cross_study_shi_seurat_label_transfer_obs.tsv.gz
        tables/cross_study_shi_seurat_label_scores_long.tsv.gz
        tables/cross_study_shi_seurat_week_scores_long.tsv.gz
        tables/cross_study_shi_seurat_ge_only_week_scores_long.tsv.gz
        tables/cross_study_shi_transfer_diagnostics_summary.tsv
        plots/umap_grids/

Step 3: Winner/90% tradeoff and predicted-age sample composition plots
  code:
    python_notebooks/scripts/plot_shi_mge_threshold_tradeoff.py
  Slurm:
    slurm_templates/30_shi_threshold_tradeoff_plots.sbatch.template
  what runs:
    - winner-take-all versus 0.90 major-label support summaries
    - predicted-age stacked bars by sample
    - all_shi_major_labels uses shi_seurat_full_predicted_shi_week_label
    - mge and mge_lge_cge use shi_seurat_ge_only_predicted_shi_week_label
  outputs:
    results/cross_study_shi_seurat_label_transfer/
      cross_study_shi_seurat_label_transfer_v2_ge_only_age/
        plots/summary/threshold_tradeoff/
        plots/summary/predicted_age_sample_composition/all_shi_major_labels/
        plots/summary/predicted_age_sample_composition/mge/
        plots/summary/predicted_age_sample_composition/mge_lge_cge/
        tables/threshold_tradeoff/
        tables/predicted_age_sample_composition/
```

Slurm templates updated for this addition:

```text
slurm_templates/27_cross_study_shi_seurat_label_transfer_array.sbatch.template
  now runs all seven studies, including Varela DIV30/DIV90, under the v2 run
  label with --reuse-existing false --force-rerun true. The array walltime is
  3 hours per task on the standard partition with 160G. The previous
  largemem/240G request was unnecessarily conservative after observed peak RSS
  from the v2 run.

slurm_templates/28_finalize_cross_study_shi_prediction_plots.sbatch.template
  defaults to the v2 run label and validates both whole-Shi and GE-only week
  score schemas. The finalizer walltime is 4 hours.

slurm_templates/28_cross_study_shi_prediction_plots_plot_only.sbatch.template
  defaults to the v2 run label for any plot-only rerender from existing tables.

slurm_templates/30_shi_threshold_tradeoff_plots.sbatch.template
  defaults to the v2 run label. Age-composition outputs for MGE and MGE/LGE/CGE
  use the GE-only week classifier.
```

Submit the v2 run:

```bash
cp slurm_templates/27_cross_study_shi_seurat_label_transfer_array.sbatch.template \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/27_cross_study_shi_seurat_label_transfer_array.v2_ge_only_age.sbatch
cp slurm_templates/28_finalize_cross_study_shi_prediction_plots.sbatch.template \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/28_finalize_cross_study_shi_prediction_plots.v2_ge_only_age.sbatch
cp slurm_templates/30_shi_threshold_tradeoff_plots.sbatch.template \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/30_shi_threshold_tradeoff_plots.v2_ge_only_age.sbatch

array_job=$(sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/27_cross_study_shi_seurat_label_transfer_array.v2_ge_only_age.sbatch | awk '{print $4}')
final_job=$(sbatch --dependency=afterok:${array_job} \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/28_finalize_cross_study_shi_prediction_plots.v2_ge_only_age.sbatch | awk '{print $4}')
plot_job=$(sbatch --dependency=afterok:${final_job} \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/30_shi_threshold_tradeoff_plots.v2_ge_only_age.sbatch | awk '{print $4}')
```

Monitor:

```bash
squeue -j ${array_job},${final_job},${plot_job} -o '%.20i %.9P %.35j %.8u %.2t %.10M %.6D %R'
sacct -j ${array_job},${final_job},${plot_job} --format=JobID,JobName%35,State,ExitCode,Elapsed,MaxRSS,ReqMem -P
```

Submitted v2 GE-only age classifier Slurm chain:

```text
array_job 51644156
  cross-shi-xfer-array
  seven-study array: Varela DIV30, Varela DIV90, Siebert 2026, Walsh 2025,
  Bershteyn 2023, Bershteyn 2025, Samarasinghe 2021

final_job 51644157
  dependency: afterok:51644156

plot_job 51644158
  dependency: afterok:51644157
```

Repair note:

```text
Initial v2 array 51644156 completed five studies but failed on the two
Bershteyn tasks after Seurat TransferData completed:

  51644156_5 bershteyn_2023 failed: missing configured sample column `sample`
  51644156_6 bershteyn_2025 failed: missing configured sample column `orig.ident`

The default study metadata config had the Bershteyn sample columns swapped.
The transfer script now uses:

  bershteyn_2023 sample_col = orig.ident|samples
  bershteyn_2025 sample_col = sample

Targeted Slurm repair chain submitted, rerunning only array tasks 5-6:

  repair_job 51648985
    submitted with --array=5-6%2

  final_job 51648986
    dependency: afterok:51648985

  plot_job 51648987
    dependency: afterok:51648986

Live repair resource correction:

  repair_job 51648985 was updated from largemem/240G/6h to
  standard/128G/2h after the failed Bershteyn attempts showed peak RSS of
  ~65G for 2023 and ~89G for 2025.

Plot-job dependency correction:

  plot_job 51648987 was released from its finalizer dependency after all seven
  per-study obs tables existed. The threshold/age plot script reads
  tables/per_study/*_shi_seurat_label_transfer_obs.tsv.gz directly, so it does
  not need to wait for the UMAP-grid finalizer.
```

## New Adjacent Reference Direction: Schmitz 2022

Schmitz et al. 2022 will be developed as a separate reference workflow where
Schmitz is the reference, not a target dataset for Shi prediction plots.

Dedicated handoff:

```text
python_notebooks/HANDOFF_schmitz_2022_reference_label_transfer.md
```

Dedicated metadata/UMAP staging code:

```text
python_notebooks/src/mge_organoid_python/schmitz_reference.py
python_notebooks/scripts/fetch_schmitz_2022_metadata_umap.py
python_notebooks/notebooks/schmitz_2022_reference_metadata_umap.ipynb
slurm_templates/29_fetch_schmitz_2022_metadata_umap.sbatch.template
```

Prepared Slurm job file:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/29_fetch_schmitz_2022_metadata_umap.sbatch
```

Completed staging run:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/schmitz_2022_reference_metadata_umap/schmitz_2022_reference_metadata_umap_v1
```

The run downloaded Cell Browser configs, cell-level metadata, and decoded UMAP
coordinates for:

```text
schmitz_macaque_dev   109,111 cells
schmitz_mouse_dev      76,804 cells
schmitz_mouse_adult   141,069 cells
```

The Schmitz staging run now also writes cross-study-style standardized cell
tables and UMAP plots:

```text
tables/per_dataset/*_schmitz_reference_umap_cells.tsv.gz
tables/schmitz_2022_reference_umap_cells.tsv.gz
tables/schmitz_2022_standardized_cell_tables.tsv
tables/schmitz_2022_umap_plot_manifest.tsv
plots/umap/
```

Schmitz-specific plotting rule:

```text
Track `species` and `dataset_id` explicitly. The combined Schmitz table contains
macaque and mouse rows, but the Cell Browser UMAP coordinates are dataset-local.
Use per-dataset/faceted UMAP panels rather than overlaying macaque and mouse as
one shared embedding.
```

Primary Schmitz annotation column:

```text
class
```

This is the Cell Browser/public metadata version of the author-script
`supervised_name` labels. Use `class` as the initial Schmitz reference label
column; keep `leiden` and `hires_leiden` as unsupervised cluster IDs.

Important caveat:

```text
The public UCSC `exprMatrix.tsv.gz` files are Float32 processed expression
matrices with negative centered values. They were not downloaded in the staging
run. Do not run Seurat-style Schmitz label transfer until the expression source
is explicitly chosen and validated.
```

Do not add Schmitz to `cross_study_shi_prediction_plots.py` yet. The Schmitz
reference has different labels, species choices, and expression-source caveats,
so it needs its own output namespace and notebook/module path.

## 2026-06-11 Bershteyn 2025 Author-vs-Our Shi Prediction Audit

Reason:

```text
The downloaded/provided Bershteyn 2025 Seurat object already contains
author/paper-style prediction metadata:
  predicted.GEtype
  predicted.GEcluster
  predicted.GEgws
  predicted.macaclass
  predicted.musclass
  type

Those author metadata columns are distinct from our from-scratch cross-study
Shi Seurat TransferData outputs:
  shi_seurat_full_*
  shi_seurat_ge_only_*

The audit asks whether our Shi-reference predictions look like the author
metadata and records the answer without using author labels as transfer inputs.
```

Execution rule:

```text
Do not load the Bershteyn 2025 RDS interactively on the login shell for this
audit. The Seurat object is ~1.1 GB and should be opened through Slurm with a
real memory allocation.
```

New audit code:

```text
scripts/15_bershteyn_author_vs_shi_prediction_audit.R
slurm_templates/32_bershteyn_author_vs_shi_prediction_audit.sbatch.template
```

Submitted Slurm job:

```text
job_id: 51679759
state: COMPLETED
exit: 0:0
elapsed: 00:00:39
MaxRSS: 6900280K
node: gl3255
memory requested: 160G

job file:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/32_bershteyn_author_vs_shi_prediction_audit.sbatch

logs:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/bersh-author-shi-audit-51679759.out
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/bersh-author-shi-audit-51679759.err
```

Follow-up score-slot inventory job:

```text
job_id: 51682246
state: COMPLETED
exit: 0:0
elapsed: 00:00:38
MaxRSS: 5931192K
node: gl3010
memory requested: 16G

logs:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/bersh-score-inventory-51682246.out
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/bersh-score-inventory-51682246.err
```

Follow-up full metadata schema job:

```text
job_id: 51685282
state: COMPLETED
exit: 0:0
elapsed: 00:00:33
MaxRSS: 5742920K
node: gl3010
memory requested: 16G

logs:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/bersh-meta-schema-51685282.out
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/bersh-meta-schema-51685282.err
```

Inputs:

```text
Author metadata source:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/bershteyn_2025/bershteyn_2025_seurat.rds

Our Shi-transfer source:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_shi_seurat_label_transfer/cross_study_shi_seurat_label_transfer_v2_ge_only_age/tables/per_study/bershteyn_2025_shi_seurat_label_transfer_obs.tsv.gz
```

Audit outputs:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_shi_seurat_label_transfer/cross_study_shi_seurat_label_transfer_v2_ge_only_age/tables/author_comparison/

bershteyn_2025_author_vs_our_shi_prediction_summary.tsv
bershteyn_2025_author_score_like_object_inventory.tsv
bershteyn_2025_author_metadata_schema.tsv
bershteyn_2025_author_object_top_level_slots.tsv
bershteyn_2025_author_vs_our_label_vocabulary.tsv
bershteyn_2025_author_vs_our_label_set_comparison.tsv
bershteyn_2025_author_vs_our_shi_prediction_distributions.tsv
bershteyn_2025_author_vs_our_shi_prediction_coarse_comparison.tsv
bershteyn_2025_author_vs_our_shi_prediction_cell_level_confusion.tsv
bershteyn_2025_our_shi_major_confidence_by_label.tsv
bershteyn_2025_our_shi_full_week_confidence_by_label.tsv
bershteyn_2025_our_shi_ge_only_week_confidence_by_label.tsv
bershteyn_2025_our_shi_full_week_by_sample.tsv
bershteyn_2025_author_vs_our_shi_prediction_joined_cells.tsv.gz
bershteyn_2025_author_vs_our_shi_prediction_notes.tsv
```

Join result:

```text
Joined cells: 124,583
Missing author metadata after join: 0
```

Author score-slot inventory:

```text
The follow-up audit searched score/prediction-like names in the downloaded
Bershteyn Seurat object's main metadata and object-level slots:
  meta.data column names
  misc names
  tools names
  commands names
  reductions names
  assays names

Only these author prediction-like entries were found:
  meta.data  predicted.GEtype     metadata_column  124583
  meta.data  predicted.GEcluster  metadata_column  124583
  meta.data  predicted.GEgws      metadata_column  124583
  meta.data  predicted.macaclass  metadata_column  124583
  meta.data  predicted.musclass   metadata_column  124583

No author score matrices or per-class author score columns were found in this
object-level inventory. The downloaded object appears to contain final author
labels, not the underlying author prediction-score vectors.

The full metadata-schema dump independently confirmed that the only numeric
metadata columns in the author object are QC fields:
  nCount_RNA
  nFeature_RNA
  percent.mt
  percent.ribo

The author prediction fields are character/factor labels, not numeric scores.
Top-level object slots are:
  assays: RNA
  reductions: pca, umap
  graphs: integrated_nn, integrated_snn
  tools: Integration
  commands: FindIntegrationAnchors, withCallingHandlers, ScaleData.integrated,
            RunPCA.integrated, RunTSNE, FindNeighbors.integrated.pca,
            FindClusters, RunUMAP.integrated.pca
```

Main result:

```text
Major GE identity broadly agrees:
  author predicted.GEtype MGE fraction:        0.987727
  our shi_seurat_full MGE fraction:            0.976538
  both author and our call MGE:              120,311 cells
  author-MGE recovered as our MGE:             97.77%
  our MGE calls that are author-MGE:           98.89%
  binary MGE Jaccard overlap:                  96.71%
  binary MGE agreement:                        96.72%
  median major-label prediction score:         0.913156
  fraction major-label score >= 0.75:          0.831654

Age/GW labels do not agree:
  author predicted.GEgws GW18 fraction:        0.902916
  our whole-Shi GW18 fraction:                 0.002841
  our GE-only GW18 fraction:                   0.000867
  our corrected whole-Shi expected GW mean:   10.428001
  our corrected GE-only expected GW mean:     10.148277
```

Exact MGE overlap:

```text
Binary definition:
  author MGE = predicted.GEtype == "MGE"
  our MGE = shi_seurat_full_predicted_shi_label == "MGE"

2x2 overlap:
  both MGE:                         120,311
  author MGE, our not-MGE:            2,743
  author not-MGE, our MGE:            1,349
  both not-MGE:                         180
  total cells:                      124,583

Metrics:
  recall of author-MGE by our MGE call:       120,311 / 123,054 = 97.77%
  precision of our MGE call vs author-MGE:    120,311 / 121,660 = 98.89%
  Jaccard overlap of MGE sets:                120,311 / 124,403 = 96.71%
  binary MGE/non-MGE agreement:              (120,311 + 180) / 124,583 = 96.72%
```

Exact GE-type vocabulary comparison:

```text
Author `predicted.GEtype` winner labels:
  MGE, CGE, POA, LGE, EN, NPC

Our `shi_seurat_full_predicted_shi_label` winner labels:
  MGE, Excitatory neuron, LGE, CGE, Thalamic neurons, progenitor

Shared exact label names:
  MGE, LGE, CGE

Author-only winner labels in this comparison:
  EN, NPC, POA

Our-only winner labels in this comparison:
  Excitatory neuron, Thalamic neurons, progenitor

Exact same-string agreement across the shared labels MGE/LGE/CGE is dominated
by MGE:
  MGE same-string matches: 120,311
  LGE same-string matches:       5
  CGE same-string matches:      12
  total same-string matches: 120,328 / 124,583 = 96.58%

This means the strong agreement is specifically a coarse MGE identity result,
not proof that the two pipelines have identical full label vocabularies.
```

Complete label vocabulary:

```text
Full machine-readable tables:
  bershteyn_2025_author_vs_our_label_vocabulary.tsv
  bershteyn_2025_author_vs_our_label_set_comparison.tsv

Author columns and labels:

  predicted.GEtype:
    MGE 123,054 (98.77%)
    CGE   1,059 (0.85%)
    POA     378 (0.30%)
    LGE      76 (0.06%)
    EN       15 (0.01%)
    NPC       1 (<0.01%)

  predicted.GEcluster:
    MGE0    122,182 (98.07%)
    CGE1      1,070 (0.86%)
    MGElhx8     766 (0.61%)
    POA        425 (0.34%)
    LGE7        88 (0.07%)
    ENs         25 (0.02%)
    MGE9        23 (0.02%)
    LGE8         2 (<0.01%)
    NPC          1 (<0.01%)
    LGE5         1 (<0.01%)

  predicted.GEgws, normalized:
    GW18 112,488 (90.29%)
    GW16   4,449 (3.57%)
    GW12   3,919 (3.15%)
    GW13   2,374 (1.91%)
    GW09   1,353 (1.09%)

  predicted.macaclass:
    MGE_LHX6/MAF     90,660 (72.77%)
    MGE_LHX6/NPY     31,255 (25.09%)
    VMF_ZIC1/ZIC2     2,657 (2.13%)
    MGE_CRABP1/MAF        7 (<0.01%)
    CGE_NR2F2/PROX1       4 (<0.01%)

  predicted.musclass:
    MGE_LHX6/MAF    107,623 (86.39%)
    MGE_LHX6/NPY     10,826 (8.69%)
    VMF_PEG10/DLK1    3,732 (3.00%)
    LGE_FOXP1/ISL1    1,458 (1.17%)
    VMF_TMEM163/OTP     473 (0.38%)
    VMF_CRABP1/LHX8     466 (0.37%)
    VMF_NR2F2/LHX6        3 (<0.01%)
    VMF_LHX1/POU6F2       2 (<0.01%)

  type:
    LHX6/MAF/ZEB2   114,958 (92.27%)
    LHX6/SST/NPY      7,955 (6.39%)
    LHX6/ENC1/LHX8    1,670 (1.34%)

Our columns and labels:

  shi_seurat_full_predicted_shi_label:
    MGE                121,660 (97.65%)
    Excitatory neuron    1,537 (1.23%)
    LGE                  1,354 (1.09%)
    CGE                     25 (0.02%)
    Thalamic neurons         4 (<0.01%)
    progenitor               3 (<0.01%)

  our_full_week_normalized:
    GW09 82,301 (66.06%)
    GW12 40,305 (32.35%)
    GW16    962 (0.77%)
    GW13    661 (0.53%)
    GW18    354 (0.28%)

  our_ge_only_week_normalized:
    GW09 92,355 (74.13%)
    GW12 30,733 (24.67%)
    GW16    761 (0.61%)
    GW13    626 (0.50%)
    GW18    108 (0.09%)
```

Coarse comparison:

```text
GE identity:
  MGE author 123,054 cells (98.77%) vs our 121,660 cells (97.65%)
  LGE author      76 cells (0.06%)  vs our   1,354 cells (1.09%)
  CGE author   1,059 cells (0.85%)  vs our      25 cells (0.02%)

GW/stage:
  GW09 author   1,353 cells (1.09%)  vs our 82,301 cells (66.06%)
  GW12 author   3,919 cells (3.15%)  vs our 40,305 cells (32.35%)
  GW13 author   2,374 cells (1.91%)  vs our    661 cells (0.53%)
  GW16 author   4,449 cells (3.57%)  vs our    962 cells (0.77%)
  GW18 author 112,488 cells (90.29%) vs our    354 cells (0.28%)
```

Age/GW cell-level disagreement:

```text
Exact author-GW vs our whole-Shi GW winner match:
  2,496 / 124,583 = 2.00%

Author GW18 cells:
  author GW18 total:                  112,488
  author GW18 -> our GW09:             74,126 (65.90% of author GW18)
  author GW18 -> our GW12:             36,515 (32.46% of author GW18)
  author GW18 -> our GW18:                346 (0.31% of author GW18)
  author GW18 -> our GW09/GW12:       110,641 (98.36% of author GW18)

Winner-label GW shift:
  mean absolute winner difference:       7.49 GW units
  mean signed difference, ours-author:  -7.48 GW units

This is a systematic early-shift in the Shi-reference age projection, not a
small random disagreement around the same age distribution.
```

How our scores were computed:

```text
The MGE/type score and GW/age score are two separate Seurat TransferData
applications, even though they use the same whole-Shi anchor set for the
`shi_seurat_full_*` outputs:

1. Major GE/cell-state score:
   refdata = reference@meta.data[[label_col]]
   output winner = shi_seurat_full_predicted_shi_label
   output max score = shi_seurat_full_prediction_score
   per-label scores = shi_seurat_full_prediction_score_<label>

2. Gestational-week score:
   refdata = reference@meta.data[[week_col]]
   output winner = shi_seurat_full_predicted_shi_week_label
   output max score = shi_seurat_full_week_prediction_score
   per-label scores = shi_seurat_full_week_prediction_score_<GW>

The GE-only week analysis then reruns FindTransferAnchors/TransferData after
subsetting the Shi reference to MGE/LGE/CGE cells. It still maps most Bershteyn
cells to GW09/GW12:
  our GE-only GW09: 92,355 cells (74.13%)
  our GE-only GW12: 30,733 cells (24.67%)
  our GE-only GW18:    108 cells (0.09%)

Therefore the age disagreement is not caused by accidentally including
non-GE labels in the week model. It persists when the age model is restricted
to GE reference cells.
```

Interpretation:

```text
Our Shi-reference transfer appears to recapitulate the broad MGE identity of
Bershteyn 2025, but it does not recapitulate the author `predicted.GEgws`
stage assignment. The author object marks most cells as GW18-like, whereas the
Shi-reference classifier maps most cells to GW09/GW12 and has a corrected
expected-GW mean near 10.4.

For methods/publication wording, phrase this as follows:
  We recomputed reference projection scores using Seurat TransferData against
  the Shi human GE/organoid reference. These scores strongly supported broad
  MGE identity and agreed with the deposited Bershteyn author MGE labels at the
  cell-set level. However, the recomputed Shi-reference age scores did not
  reproduce the deposited author GEgws labels; most author-GW18 cells projected
  to GW09/GW12 in the Shi reference. Because the deposited Bershteyn object
  contains author labels but not author per-class score vectors, author-score
  equality cannot be directly tested from the available object.

Do not present our Shi GW output as a reproduction of the Bershteyn author
`predicted.GEgws` metadata. Treat it as a separate reference-system projection.
```

Important expected-GW caveat:

```text
The per-study obs table contains stale/buggy expected-GW numeric columns from
the R transfer step:
  per_study_obs_full_expected_gw_mean_stale_buggy:    8.286551
  per_study_obs_ge_only_expected_gw_mean_stale_buggy: 8.430542

Cause:
  the R helper `parse_gw_numeric()` averaged all numbers found in duplicate
  score-column labels such as `GW12_01`, so `GW12_01` was not treated as GW12.

The Python finalizer and this audit recompute corrected expected-GW values by
parsing the `GWxx` prefix only. Use corrected/finalizer summary values for
expected-GW interpretation, not the stale per-study obs numeric expected-GW
columns. Winner-take-all week labels are unaffected except for display
normalization of duplicate labels such as `GW12_01` and `GW12_02` to `GW12`.

Patch:
  `scripts/13_run_cross_study_shi_seurat_label_transfer.R` now parses the
  `GWxx` prefix first and falls back to the first numeric token only when no
  GW-style token is present. Existing v2 per-study obs files are not rewritten
  by this patch; rerun the transfer array if refreshed per-study expected-GW
  columns are needed.
```

Re-run command:

```bash
cp slurm_templates/32_bershteyn_author_vs_shi_prediction_audit.sbatch.template \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/32_bershteyn_author_vs_shi_prediction_audit.sbatch
sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/32_bershteyn_author_vs_shi_prediction_audit.sbatch
```

# Handoff: Jia Program DIV30 Scoring

This is a standalone Slurm-executed notebook workflow. It is not Notebook 02
and does not continue the Notebook 01 pipeline. Its narrow goal is to score the
Jia et al. 2026 Science RGC/IPC programs on the converted DIV30 AnnData object
and compare those scores against the built-in Seurat clusters plus the existing
DIV30 resolution-sweep cluster assignments.

## Current Run Target

Default run label:

```text
jia_program_div30_scoring_v1
```

Default run directory:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/jia_program_div30_scoring/jia_program_div30_scoring_v1
```

Executed notebook path:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/jia_program_div30_scoring/executed/jia_program_div30_scoring.jia_program_div30_scoring_v1.executed.ipynb
```

## Code Added

Notebook:

```text
python_notebooks/notebooks/jia_program_div30_scoring.ipynb
```

Reusable methods module:

```text
python_notebooks/src/mge_organoid_python/gene_program_scoring.py
```

Slurm template:

```text
slurm_templates/19_execute_jia_program_div30_scoring.sbatch.template
```

## Inputs

Marker/program CSV:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/reference/Jia_et_al_2026_Science_3_progs.csv
```

Expected CSV columns:

```text
gene, cluster, p_val, avg_log2FC, pct.1, pct.2, p_val_adj
```

Default program order:

```text
IPC,RGC1,RGC2
```

DIV30 converted AnnData input:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/varela_div30.h5ad
```

DIV30 resolution-sweep assignments:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/seurat_cluster_resolution_sweep/seurat_cluster_resolution_sweep_v1/varela_div30/tables/resolution_sweep_cluster_assignments_wide.tsv
```

## Scoring Method

The first-pass method uses `scanpy.tl.score_genes` on the converted DIV30
AnnData `.X` matrix. This is the Scanpy analogue of Seurat-style
`AddModuleScore`: each program score is the average expression of the program
genes adjusted against expression-binned control genes.

The run preserves marker-gene audit tables:

```text
RUN_DIR/tables/jia_program_markers_full.tsv
RUN_DIR/tables/jia_program_markers_selected.tsv
RUN_DIR/tables/jia_program_gene_overlap_summary.tsv
RUN_DIR/tables/jia_program_gene_overlap_detail.tsv
```

By default the run uses all markers in the CSV. Optional environment controls:

```text
JIA_PROGRAM_TOP_N=
JIA_PROGRAM_MIN_AVG_LOG2FC=
JIA_PROGRAM_MAX_P_VAL_ADJ=
JIA_PROGRAM_CTRL_SIZE=50
JIA_PROGRAM_RANDOM_STATE=0
JIA_PROGRAM_HIGH_SCORE_QUANTILE=0.9
```

## Output Contract

Run outputs:

```text
RUN_DIR/tables/
RUN_DIR/plots/
RUN_DIR/logs/
```

Important tables:

```text
RUN_DIR/tables/jia_program_run_parameters.tsv
RUN_DIR/tables/jia_program_marker_selection_summary.tsv
RUN_DIR/tables/jia_program_gene_overlap_summary.tsv
RUN_DIR/tables/jia_program_gene_overlap_detail.tsv
RUN_DIR/tables/jia_program_score_thresholds.tsv
RUN_DIR/tables/div30_jia_program_scores_obs.tsv
RUN_DIR/tables/div30_jia_program_summary_by_seurat_clusters.tsv
RUN_DIR/tables/div30_jia_program_summary_by_resolution_sweep.tsv
RUN_DIR/tables/div30_jia_program_summary_all_groupings.tsv
RUN_DIR/tables/div30_jia_program_best_matches_by_grouping.tsv
RUN_DIR/tables/jia_program_output_manifest.tsv
RUN_DIR/tables/jia_program_div30_scoring_complete.tsv
```

Required plot:

```text
RUN_DIR/plots/div30_umap_seurat_clusters_jia_program_scores_panel.png
```

That panel is one row:

```text
[DIV30 built-in Seurat clusters] [IPC score overlay] [RGC1 score overlay] [RGC2 score overlay]
```

The first panel colors all cells by the built-in Seurat cluster column
(`seurat_clusters`, falling back to `RNA_snn_res.0.2`). The program panels use
the same UMAP coordinates with all cells drawn in light gray, then the program
score drawn over the top as a continuous overlay.

Additional summary plots:

```text
RUN_DIR/plots/div30_jia_program_mean_score_by_seurat_clusters_heatmap.png
RUN_DIR/plots/div30_jia_program_fraction_high_by_seurat_clusters_heatmap.png
```

## Slurm Command

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
cp slurm_templates/19_execute_jia_program_div30_scoring.sbatch.template \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/19_execute_jia_program_div30_scoring.sbatch
sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/19_execute_jia_program_div30_scoring.sbatch
```

Slurm logs:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/jia-program-div30-jia-prog-div30-<jobid>.out
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/jia-program-div30-jia-prog-div30-<jobid>.err
```

## Interpretation Notes

This run does not recompute UMAP. It uses the exported Seurat UMAP coordinates
stored in the DIV30 AnnData object.

This run does not change cluster resolution inside the AnnData object. It joins
the already-completed DIV30 Seurat resolution-sweep assignment table and
summarizes Jia scores by each sweep column.

The method answers whether the DIV30 Seurat clusters, or any resolution-sweep
clusters, are enriched for the Jia `IPC`, `RGC1`, or `RGC2` programs by looking
at mean program score and fraction of cells above the global high-score
threshold for each program.

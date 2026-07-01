# Handoff: Siebert Marker Coexpression UMAP Analysis

Date started: 2026-06-30

Purpose: non-final exploratory analysis for Siebert et al. 2026 cells. This is
not a final figure package.

## Goal

Make a compact two-row UMAP figure:

```text
Top row:
  clusters | DCX | PCDH19 | VIM | Ki67/MKI67

Bottom row:
  DCX+ Ki67+
  DCX+ Ki67+ SOX2+
  DCX+ Ki67+ SOX2+ PCDH19+
  sequential percent-positive bars
  final all-four coexpression percentage
```

Positive-cell cutoff:

```text
expression >= 0.5
```

Expression scale:

```text
log1p(CP10K), computed directly from RNA counts in the Siebert Seurat object.
```

`ki67m` in the user request is treated as `MKI67`, with aliases `KI67`,
`KI67M`, `Ki67`, and `Ki67m` accepted.

## Code

Renderer:

```text
scripts/16_siebert_marker_coexpression_umap_analysis.R
```

Slurm template:

```text
slurm_templates/61_siebert_marker_coexpression_umap_analysis.sbatch.template
```

Copied job script:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/61_siebert_marker_coexpression_umap_analysis.sbatch
```

## Inputs

Canonical Siebert Seurat object:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siebert_2026/siebert_2026_seurat.rds
```

Requested markers:

```text
DCX
PCDH19
VIM
MKI67
SOX2
```

UMAP reduction:

```text
umap
```

Cluster column:

```text
seurat_clusters
```

## Outputs

Default output root:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siebert_2026/analysis/siebert_marker_coexpression_umap_v1
```

Expected files:

```text
plots/siebert_marker_coexpression_umap_grid.png
plots/siebert_marker_coexpression_umap_grid.pdf
plots/siebert_marker_coexpression_umap_grid.svg
tables/siebert_marker_gene_matches.tsv
tables/siebert_marker_coexpression_summary.tsv
tables/siebert_marker_coexpression_by_cluster.tsv
tables/siebert_marker_expression_and_coexpression_per_cell.tsv.gz
siebert_marker_coexpression_manifest.tsv
```

## Run Command

Submitted on 2026-06-30:

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
mkdir -p /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs
cp slurm_templates/61_siebert_marker_coexpression_umap_analysis.sbatch.template \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/61_siebert_marker_coexpression_umap_analysis.sbatch
sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/61_siebert_marker_coexpression_umap_analysis.sbatch
```

Initial Slurm job:

```text
52603853
```

Initial job status:

```text
FAILED, 2026-06-30 11:57 EDT
```

Failure reason:

```text
The extraction and tables completed, but plot rendering failed because the
cluster-label `geom_text()` layer did not explicitly set x/y aesthetics.
```

Patch:

```text
scripts/16_siebert_marker_coexpression_umap_analysis.R
```

The cluster-label layer now uses:

```text
aes(x = umap_1, y = umap_2, label = cluster)
```

Resubmitted fixed job:

```text
52605895
```

Status at update:

```text
PENDING (Priority), 2026-06-30 12:02 EDT
```

## 2026-06-30 High-contrast Rerender and Object Audit

The first completed plot was visually too faint. The renderer was updated to:

```text
- Plot all cells as a grey background, then overlay nonzero marker-expression
  cells with a saturated high-contrast color scale.
- Increase UMAP point sizes for expression/coexpression overlays.
- Wrap the long DCX+ Ki67+ SOX2+ PCDH19+ title.
- Write explicit Seurat object and marker-expression QC tables.
```

Updated output root:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siebert_2026/analysis/siebert_marker_coexpression_umap_v2_high_contrast
```

Successful Slurm job:

```text
52609398
```

Runtime/resource:

```text
COMPLETED, 2026-06-30 13:07:48-13:09:09 EDT
Elapsed: 00:01:21
MaxRSS: 11469040K
Requested memory: 64G
```

Important audit output:

```text
tables/siebert_seurat_object_qc.tsv
tables/siebert_marker_gene_matches.tsv
tables/siebert_marker_object_expression_qc.tsv
```

The object audit confirms:

```text
Seurat object: results/siebert_2026/siebert_2026_seurat.rds
Assays: RNA,SCT
Assay used: RNA
Layer used: counts
Reduction used: umap
Cluster column: seurat_clusters
Cells: 64,676
Features: 27,823
```

All requested markers matched exact feature rows:

```text
DCX -> DCX
PCDH19 -> PCDH19
VIM -> VIM
MKI67 -> MKI67
SOX2 -> SOX2
```

Marker raw-count support from the Seurat object:

| Gene | Raw total counts | Raw nonzero cells | Raw nonzero % | Cells >= 0.5 log1p(CP10K) |
| --- | ---: | ---: | ---: | ---: |
| DCX | 148,519 | 40,338 | 62.37 | 62.07 |
| PCDH19 | 3,835 | 3,331 | 5.15 | 4.88 |
| VIM | 547,981 | 50,930 | 78.75 | 78.72 |
| MKI67 | 15,132 | 6,015 | 9.30 | 9.17 |
| SOX2 | 108,384 | 35,442 | 54.80 | 54.75 |

Sequential coexpression at cutoff expression >= 0.5:

| State | Cells | Percent |
| --- | ---: | ---: |
| DCX+ Ki67+ | 4,070 / 64,676 | 6.29 |
| DCX+ Ki67+ SOX2+ | 3,512 / 64,676 | 5.43 |
| DCX+ Ki67+ SOX2+ PCDH19+ | 464 / 64,676 | 0.72 |

Primary high-contrast figure:

```text
plots/siebert_marker_coexpression_umap_grid.png
plots/siebert_marker_coexpression_umap_grid.pdf
plots/siebert_marker_coexpression_umap_grid.svg
```

## 2026-06-30 v3 Sample Breakdown Rerender

The user requested the per-sample breakdown at the end of the figure. The
renderer was updated to append:

```text
- Top-right sample UMAP panel.
- Bottom-right sample x sequential-state percent-positive heatmap.
- Exact per-sample table:
  tables/siebert_marker_coexpression_by_sample.tsv
```

Updated output root:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siebert_2026/analysis/siebert_marker_coexpression_umap_v3_sample_breakdown
```

Successful Slurm job:

```text
52620046
```

Runtime/resource:

```text
COMPLETED, 2026-06-30 14:44:38-14:46:11 EDT
Elapsed: 00:01:33
MaxRSS: 15851680K
Requested memory: 64G
```

Per-sample sequential coexpression at cutoff expression >= 0.5:

| Sample | DCX+ Ki67+ | DCX+ Ki67+ SOX2+ | DCX+ Ki67+ SOX2+ PCDH19+ |
| --- | ---: | ---: | ---: |
| Young_1 | 1,134 / 15,620 = 7.26% | 923 / 15,620 = 5.91% | 95 / 15,620 = 0.61% |
| Young_2 | 1,368 / 16,377 = 8.35% | 1,164 / 16,377 = 7.11% | 121 / 16,377 = 0.74% |
| Old_1 | 819 / 16,606 = 4.93% | 736 / 16,606 = 4.43% | 116 / 16,606 = 0.70% |
| Old_2 | 749 / 16,073 = 4.66% | 689 / 16,073 = 4.29% | 132 / 16,073 = 0.82% |

## 2026-06-30 DIV30 Organoid Add-on

The user requested adding the paper's DIV30 organoids to the same folder using
the same marker/coexpression logic.

The renderer was generalized to support study-specific labels, sample columns,
and output prefixes. DIV30 writes `div30_marker_*` files into the same output
root as the Siebert v3 sample-breakdown render:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siebert_2026/analysis/siebert_marker_coexpression_umap_v3_sample_breakdown
```

DIV30-specific Slurm template:

```text
slurm_templates/62_div30_marker_coexpression_umap_analysis.sbatch.template
```

Successful Slurm job:

```text
52622543
```

Runtime/resource:

```text
COMPLETED, 2026-06-30 15:16:39-15:19:42 EDT
Elapsed: 00:03:03
MaxRSS: 50547596K
Requested memory: 64G
```

DIV30 object/run contract:

```text
Seurat object: results/varela_this_paper/varela_this_paper_seurat.rds
Study label: This Study, DIV 30
Assay used: RNA
Count layers used: counts.9583-MW-1 through counts.9583-MW-6
Reduction used: umap
Sample column: orig.ident
Cluster column: seurat_clusters
Cells: 90,631
Features: 18,082
```

All requested markers matched exact feature rows:

```text
DCX -> DCX
PCDH19 -> PCDH19
VIM -> VIM
MKI67 -> MKI67
SOX2 -> SOX2
```

DIV30 outputs:

```text
plots/div30_marker_coexpression_umap_grid.png
plots/div30_marker_coexpression_umap_grid.pdf
plots/div30_marker_coexpression_umap_grid.svg
tables/div30_marker_coexpression_summary.tsv
tables/div30_marker_coexpression_by_sample.tsv
tables/div30_marker_coexpression_by_cluster.tsv
tables/div30_marker_object_expression_qc.tsv
tables/div30_marker_expression_and_coexpression_per_cell.tsv.gz
```

DIV30 sequential coexpression at cutoff expression >= 0.5:

| State | Cells | Percent |
| --- | ---: | ---: |
| DCX+ Ki67+ | 17,468 / 90,631 | 19.27 |
| DCX+ Ki67+ SOX2+ | 15,279 / 90,631 | 16.86 |
| DCX+ Ki67+ SOX2+ PCDH19+ | 7,420 / 90,631 | 8.19 |

DIV30 per-organoid all-four coexpression:

| Sample | Cells | Percent |
| --- | ---: | ---: |
| 9583-MW-1 | 3,003 / 14,564 | 20.62 |
| 9583-MW-2 | 199 / 4,532 | 4.39 |
| 9583-MW-3 | 1,752 / 15,773 | 11.11 |
| 9583-MW-4 | 466 / 11,004 | 4.23 |
| 9583-MW-5 | 1,157 / 22,230 | 5.20 |
| 9583-MW-6 | 843 / 22,528 | 3.74 |

## 2026-06-30 Combined Cutoff Summary

The combined per-study/per-sample cutoff summary was generated from the
exported per-cell expression tables, so it can show both expression cutoffs
without rerunning the Seurat jobs.

Script:

```text
python_notebooks/scripts/summarize_marker_coexpression_cutoffs.py
```

Command:

```bash
/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python python_notebooks/scripts/summarize_marker_coexpression_cutoffs.py \
  --analysis-dir /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siebert_2026/analysis/siebert_marker_coexpression_umap_v3_sample_breakdown \
  --cutoffs 0.5,1.0
```

Combined outputs, all in the shared Siebert/DIV30 output root:

```text
tables/combined_marker_coexpression_by_study_sample_cutoff.tsv
tables/combined_marker_coexpression_overall_by_study_cutoff.tsv
plots/combined_marker_coexpression_by_sample_cutoff_heatmap.png
plots/combined_marker_coexpression_by_sample_cutoff_heatmap.pdf
plots/combined_marker_coexpression_by_sample_cutoff_heatmap.svg
```

Overall all-four coexpression:

| Study | Cutoff | Cells | Percent |
| --- | ---: | ---: | ---: |
| This Study, DIV 30 | 0.5 | 7,420 / 90,631 | 8.19 |
| This Study, DIV 30 | 1.0 | 3,283 / 90,631 | 3.62 |
| Siebert 2026 | 0.5 | 464 / 64,676 | 0.72 |
| Siebert 2026 | 1.0 | 122 / 64,676 | 0.19 |

Per-sample all-four coexpression:

| Study | Sample | Cutoff 0.5 | Cutoff 1.0 |
| --- | --- | ---: | ---: |
| This Study, DIV 30 | 9583-MW-1 | 3,003 / 14,564 = 20.62% | 1,598 / 14,564 = 10.97% |
| This Study, DIV 30 | 9583-MW-2 | 199 / 4,532 = 4.39% | 106 / 4,532 = 2.34% |
| This Study, DIV 30 | 9583-MW-3 | 1,752 / 15,773 = 11.11% | 763 / 15,773 = 4.84% |
| This Study, DIV 30 | 9583-MW-4 | 466 / 11,004 = 4.23% | 162 / 11,004 = 1.47% |
| This Study, DIV 30 | 9583-MW-5 | 1,157 / 22,230 = 5.20% | 361 / 22,230 = 1.62% |
| This Study, DIV 30 | 9583-MW-6 | 843 / 22,528 = 3.74% | 293 / 22,528 = 1.30% |
| Siebert 2026 | Young_1 | 95 / 15,620 = 0.61% | 35 / 15,620 = 0.22% |
| Siebert 2026 | Young_2 | 121 / 16,377 = 0.74% | 39 / 16,377 = 0.24% |
| Siebert 2026 | Old_1 | 116 / 16,606 = 0.70% | 24 / 16,606 = 0.14% |
| Siebert 2026 | Old_2 | 132 / 16,073 = 0.82% | 24 / 16,073 = 0.15% |

Check status:

```bash
squeue -j 52605895 -o '%.18i %.9P %.32j %.8u %.2t %.10M %.6D %R'
sacct -j 52605895 --format=JobID,JobName%30,State,ExitCode,Elapsed,MaxRSS,ReqMem,Start,End
```

Expected logs:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/siebert-marker-coexpr-siebert-marker-coexpr-52603853.out
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/siebert-marker-coexpr-siebert-marker-coexpr-52603853.err
```

## Notes

- Existing cross-study compact marker-expression tables include `DCX` and
  `MKI67` for Siebert, but not all requested markers (`PCDH19`, `VIM`, `SOX2`),
  so this analysis extracts directly from the Seurat object.
- The Siebert object previously had empty RNA `data`/`logcounts` layers in
  label-transfer work, so this renderer intentionally uses RNA counts and
  computes log1p(CP10K) for the five requested genes.
- This analysis writes exploratory outputs under `results/siebert_2026/analysis`
  and should not be treated as a final-figure package unless explicitly
  promoted later.

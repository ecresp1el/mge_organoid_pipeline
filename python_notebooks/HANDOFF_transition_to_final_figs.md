# Handoff: Transition to Final Figs

Living handoff for moving candidate analysis plots into publication-quality,
reproducible final figure folders.

This workflow will be reused many times and may eventually track 20-30+ figures.
Keep it organized as a registry: a compact tracker table at the top, then one
detail section per figure. Amend this file whenever a figure is found,
confirmed, modified, rerendered, or finalized.

## How to Keep This Handoff Scalable

Use two levels of documentation:

```text
1. Current Figure Tracker
   One compact row per candidate/final figure.

2. Figure Detail Sections
   One section per figure with paths, logs, decisions, modifications, and final
   output locations.
```

Do not put long narrative notes in the tracker table. Put only the current
status and the most important pointer there. Put details in the figure-specific
section.

Use stable figure IDs that can become folder names:

```text
fig_<short_topic>_<version_or_panel>
```

Examples:

```text
fig_cross_study_marker_expression_v12
fig_shi_prediction_gw_grid_v2
fig_div30_urd_branchpoint_v1
```

Recommended status vocabulary:

```text
Found
Log-audited
Confirmed
Modifying
Modified
Validated
Finalized
Parked
Rejected
```

Each figure should have exactly one detailed section named:

```text
## Figure: <figure_id>
```

Inside that section, keep the same subsections in the same order:

```text
Status
Candidate Paths
Source Code and Handoffs
Prepared Assets
Log Audit
Rerun Decision
Modification Log
Validation
Final Figure Package
Open Questions
```

If a figure is rejected or parked, still leave the section in place with the
reason. That prevents rediscovering the same dead end later.

## Core Rule

Prefer plot regeneration from existing prepared assets over rerunning upstream
analysis.

Before rerunning anything, check:

```text
1. What generated the candidate figure?
2. Are prepared tables, cached objects, or plotting-ready assets already present?
3. Do logs show the original job was short enough to rerun safely?
4. Is a plot-only rerender path available?
5. Can formatting be changed without recomputing upstream analysis?
```

If the full job was long, memory-heavy, or upstream-analysis heavy, do not rerun
the full workflow unless explicitly needed. Instead, isolate the plotting code
and point it at the existing assets.

## Standard Final Figure Location

Use a structured final figure area under the Great Lakes project root:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/final_figures/
```

For each figure, create a dedicated folder:

```text
final_figures/<figure_id>/
  README.md
  code/
  figures/
    png/
    pdf/
    svg/
  tables/
  logs/
  provenance/
```

Expected contents:

- `figures/png/`: high-resolution PNG output.
- `figures/pdf/`: vector or publication-quality PDF output.
- `figures/svg/`: SVG output whenever technically feasible.
- `code/`: exact scripts, wrappers, notebook extracts, or Slurm scripts used for
  the final render.
- `tables/`: copied or symlinked plotting inputs when small enough; otherwise
  include a manifest pointing to the original prepared assets.
- `logs/`: Slurm logs or local command logs from the final render.
- `provenance/`: source paths, run labels, git commit, environment notes, and
  checksums/manifests where useful.

Do not rely on memory or notebook state for final figures. The final folder
must explain how the figure was made.

## Figure Transition Checklist

For every candidate figure:

```text
[ ] Found: record candidate path(s), run label, and source handoff/code.
[ ] Confirmed: inspect the plot and decide it is the right biological figure.
[ ] Log audit: check Slurm logs/accounting and record runtime/resources.
[ ] Rerun decision: decide plot-only vs full rerun vs no rerun.
[ ] Modify: make formatting changes in an isolated test/final location.
[ ] Validate: inspect PNG/PDF/SVG outputs and confirm labels/layout are correct.
[ ] Finalize: copy exact code, outputs, logs, and provenance into final_figures.
```

Recommended log checks:

```bash
find /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs -maxdepth 1 -type f -iname '*<keyword>*' -printf '%TY-%Tm-%Td %TH:%TM %s %p\n' | sort -r | head
sacct -j <job_id> --format=JobID,JobName%28,State,ExitCode,Elapsed,MaxRSS,AllocCPUS,ReqMem,Start,End
tail -n 80 <log.out>
tail -n 80 <log.err>
```

Record whether the job was:

- `plot-only`: safe default for formatting changes.
- `full-but-short`: acceptable to rerun if it does not redo expensive upstream
  analysis or if the user explicitly approves.
- `long/upstream-heavy`: avoid rerun; extract plotting code and use existing
  assets.

## Reproducibility Requirements

Each finalized figure must include:

- Source candidate path.
- Final output path.
- Exact code used for final render.
- Input data/assets used for plotting.
- Run command or Slurm submit command.
- Log files from final render.
- Runtime and resource notes.
- Git commit or working-tree status at render time.
- Notes on any manual decisions: excluded studies, gene order, label edits,
  color scales, dimensions, DPI, file formats.

For publication-quality output, prefer:

```text
PNG: 300-600 dpi, depending on figure size and journal needs.
PDF: vector when possible; otherwise high-quality embedded raster.
SVG: generated whenever Matplotlib/ggplot output supports it cleanly.
```

## Current Figure Tracker

| Figure ID | Status | Candidate Path | Notes |
| --- | --- | --- | --- |
| `fig_cross_study_marker_expression_v12` | Found, Log-audited | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_marker_expression/cross_study_marker_expression_v12` | Multi-study ON/OFF-target and PV precursor marker-expression UMAP grids. Formatting changes pending. |

## Figure: fig_cross_study_marker_expression_v12

### Status

```text
Found: yes
Log-audited: yes
Confirmed: pending user visual review
Modified: no
Validated: no
Finalized: no
```

### Candidate Paths

Source handoff:

```text
python_notebooks/HANDOFF_shi_reference_div30_label_transfer.md
Section: Cross-Study Marker Expression Plot Workflow
```

Candidate output root:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_marker_expression/cross_study_marker_expression_v12
```

Plot files:

```text
plots/cross_study_marker_expression_on_target.png
plots/cross_study_marker_expression_on_target.pdf
plots/cross_study_marker_expression_off_target.png
plots/cross_study_marker_expression_off_target.pdf
plots/cross_study_marker_expression_on_off_target.png
plots/cross_study_marker_expression_on_off_target.pdf
plots/cross_study_marker_expression_pv_precursors_on_target.png
plots/cross_study_marker_expression_pv_precursors_on_target.pdf
plots/cross_study_marker_expression_pv_precursors_off_target.png
plots/cross_study_marker_expression_pv_precursors_off_target.pdf
plots/cross_study_marker_expression_pv_precursors_on_off_target.png
plots/cross_study_marker_expression_pv_precursors_on_off_target.pdf
```

### Source Code and Handoffs

Code entrypoints:

```text
python_notebooks/src/mge_organoid_python/cross_study_marker_expression.py
python_notebooks/scripts/run_cross_study_marker_expression.py
slurm_templates/25_cross_study_marker_expression_plot_only.sbatch.template
slurm_templates/25_cross_study_marker_expression.sbatch.template
```

### Prepared Assets

Prepared plotting assets:

```text
tables/per_study/*_marker_expression.tsv.gz
tables/cross_study_marker_expression_gene_panels.tsv
tables/cross_study_marker_expression_plot_manifest.tsv
tables/cross_study_marker_expression_distribution_audit.tsv
tables/cross_study_marker_expression_internal_plot_filter_summary.tsv
```

Cluster QC assets created during final-figure transition:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_marker_expression/cross_study_marker_expression_v12/cluster_qc/
```

Cluster QC source manifest:

```text
cluster_qc/cluster_umap_qc_source_manifest.tsv
```

This manifest records, per study:

```text
study_id
study_label
marker_table_path
source_seurat_path
source_h5ad_path
reduction
cluster_col_source
sample_col_source
n_cells_loaded
```

Cluster QC tables:

```text
cluster_qc/cluster_umap_qc_figure_filter_summary.tsv
cluster_qc/cluster_umap_qc_study_summary.tsv
cluster_qc/cluster_umap_qc_cluster_counts.tsv
cluster_qc/cluster_umap_qc_sample_cluster_counts.tsv
```

Cluster QC plots:

```text
cluster_qc/cross_study_marker_expression_v12_cluster_umap_qc_figure_default.png
cluster_qc/cross_study_marker_expression_v12_cluster_umap_qc_figure_default.pdf
cluster_qc/cross_study_marker_expression_v12_cluster_umap_qc_figure_default.svg
cluster_qc/cross_study_marker_expression_v12_cluster_umap_qc_all_prepared_cells.png
cluster_qc/cross_study_marker_expression_v12_cluster_umap_qc_all_prepared_cells.pdf
cluster_qc/cross_study_marker_expression_v12_cluster_umap_qc_all_prepared_cells.svg
```

### Log Audit

Log audit:

```text
Full v12 job:
  Job ID: 51527782
  Log: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/cross-marker-expression-cross-marker-expr-51527782.out
  State: COMPLETED
  Elapsed: 00:04:58
  MaxRSS: 14564532K
  AllocCPUS: 8
  ReqMem: 160G

Plot-only v12 rerender:
  Job ID: 51530179
  Log: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/cross-marker-plot-cross-marker-plot-51530179.out
  State: COMPLETED
  Elapsed: 00:01:42
  MaxRSS: 920412K
  AllocCPUS: 4
  ReqMem: 32G

Earlier plot-only rerender:
  Job ID: 51529553
  State: COMPLETED
  Elapsed: 00:01:47
```

### Rerun Decision

Decision:

```text
Formatting changes should use plot-only rerender from existing prepared tables.
Full rerun is short, but plot-only is the correct default unless genes/studies or
prepared marker tables need to change.
```

Cluster QC decision:

```text
Use existing v12 per-study marker-expression tables. These already contain
cell_id, study_id, study_label, sample, cluster, umap_1, and umap_2. No upstream
analysis or Seurat/H5AD extraction was rerun.
```

### Modification Log

```text
2026-06-25:
  - Candidate identified.
  - Logs checked.
  - No formatting modifications made yet.
  - Created quick cluster UMAP QC script:
    python_notebooks/scripts/plot_cross_study_cluster_umap_qc.py
  - Generated cluster QC outputs under:
    /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_marker_expression/cross_study_marker_expression_v12/cluster_qc
  - Figure-default QC applies the same study/filter logic as v12 plotting:
    exclude bershteyn_2025; Samarasinghe controls only.
  - Updated cluster QC labels for this-study panels:
    DIV30 uses collapsed paper/manual cluster classes.
    DIV90 uses the audited cluster_number_name mapping table.
```

Mapped label sources:

```text
DIV30:
  Source script: python_notebooks/scripts/map_div30_paper_cluster_annotations.py
  Source handoff: python_notebooks/HANDOFF_div30_paper_clusters_urd.md
  Mapping:
    raw clusters 0,3,7 -> 1 - Radial glia
    raw cluster 6     -> 2 - Inhibitory progenitors
    raw cluster 1     -> 3 - SST+ cIN
    raw cluster 4     -> 4 - PV neuron precursor
    raw cluster 2     -> 5 - MGE subpallial neurons

DIV90:
  Source audit:
    /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_umap_cluster_label_audit/div90_umap_cluster_label_audit_v1
  Mapping table:
    tables/div90_cluster_number_name_to_biology_mapping.tsv
  Metadata columns confirmed in source obs:
    cluster_id
    cluster_number_name
```

### Validation

```text
Pending user visual review of candidate v12 plots.
Cluster QC grid was visually spot-checked after render. The figure-default grid
uses study table order and shows cluster IDs over each study's UMAP.
Second QC render confirmed DIV30 collapsed to 5 mapped paper/manual classes and
DIV90 joined to the 13 audited cluster_number_name labels.
```

### Final Figure Package

```text
Not created yet.
Target root:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/final_figures/fig_cross_study_marker_expression_v12
```

### Open Questions

```text
Which exact panels should be finalized?
What formatting changes are needed before final render?
Should SVG export be added to the plotting function for this figure family?
```

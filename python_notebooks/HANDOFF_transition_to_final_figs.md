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

Publication export standards for this final-figure series:

```text
Text/editability:
  - SVG text must remain editable in Illustrator.
  - Matplotlib SVG export must use:
      svg.fonttype = none
      font.family = Arial
      font.sans-serif = Arial first, then available fallbacks
  - PDF export should use TrueType text when possible:
      pdf.fonttype = 42
      ps.fonttype = 42
  - Great Lakes may not have Microsoft Arial installed locally. This can emit
    findfont warnings, but the SVG must still be checked directly for editable
    <text> elements with Arial font-family declarations and no DejaVu glyph-path
    text definitions.

DPI/resolution:
  - Record DPI separately for PNG, PDF, and SVG because rasterized UMAP point
    layers embedded inside PDF/SVG depend on savefig dpi.
  - Keep DPI consistent within a figure package whenever the render can
    complete safely.
  - For cluster/QC UMAP outputs, 600 dpi export has completed successfully and
    should remain the current target.
  - For very large marker-expression grids, do not attempt high-DPI final
    renders interactively. Use Slurm and copy outputs into final_figures only
    after logs confirm the job completed and all PNG/PDF/SVG files are nonzero.

Slurm execution standard:
  - Follow the project's existing Slurm-template logic instead of inventing a
    new activation pattern.
  - Prefer copying/editing the closest existing template under
    slurm_templates/ into PROJECT_ROOT/jobs/ and submitting that file.
  - Use the established path variables:
      REPO_ROOT=/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
      PROJECT_ROOT=/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder
      CONDA_ENV_BIN=/home/elcrespo/miniconda3/envs/mge-organoid-python/bin
  - In Python Slurm jobs, use:
      export PATH="${CONDA_ENV_BIN}:${PATH}"
      export PYTHONPATH="${REPO_ROOT}/python_notebooks/src:${PYTHONPATH:-}"
      cd "${REPO_ROOT}"
      "${CONDA_ENV_BIN}/python" ...
    This matches the existing project templates such as
    25_cross_study_marker_expression_plot_only.sbatch.template,
    17b_execute_seurat_anndata_umap_inventory_merge.sbatch.template, and
    50_siletti_div90_harmony_integration_sensitivity.sbatch.template.
  - Log to PROJECT_ROOT/logs with Slurm %x/%j or equivalent job IDs.
  - Record the exact sbatch file, sbatch command, job ID, stdout/stderr log
    paths, runtime, memory request, and whether the job exited cleanly.
  - Do not sync final outputs from a Slurm render until all required PNG/PDF/SVG
    files exist, are nonzero, and the log shows successful completion.
```

## Current Figure Tracker

| Figure ID | Status | Candidate Path | Notes |
| --- | --- | --- | --- |
| `fig_cross_study_marker_expression_v12` | Found, Log-audited, Modifying, Validated, Final package started | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_marker_expression/cross_study_marker_expression_v12` | First final folder created for cluster UMAP QC / DIV90 published recode outputs. Rerendered with editable Arial SVG text and 600 dpi export for rasterized UMAP layers. Other marker-expression multi-grids still pending formatting/finalization. |
| `fig_cross_study_marker_expression_pv_precursors_on_off_target_v12` | Modified, Slurm-rendered, Packaged candidate | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_marker_expression/cross_study_marker_expression_v12_pv_precursors_final_candidate` | Plot-only rerender from v12 prepared marker tables. Added Bershteyn 2025, added PV ON-target LHX6/LHX8/NKX2.1, retained ERBB4, applied DIV90 vertical plotting-only UMAP orientation, removed DIV90 stressed clusters 6/7 from visualization only, drew expression values 0-1 as background gray with blue scale starting above 1, colorbars labeled from 0, exported PNG/PDF/SVG at 600 dpi through Slurm job 52370542 with editable Arial SVG text. |
| `fig_cross_study_shi_label_transfer_v1_umap_score_grids` | Modified, Slurm-rendered, Packaged candidate | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_shi_seurat_label_transfer/cross_study_shi_seurat_label_transfer_v1` | Plot-only rerender of Shi Seurat label-transfer UMAP score grids and matched sample-composition panels from saved v1 combined table. Updated score grids to fixed 0-1 grey-to-blue scaling, applied DIV90 visualization-only stressed-cluster 6/7 removal, applied DIV90 vertical plotting-only UMAP orientation, wrapped long all-label headers, exported requested PNG/PDF/SVG at 600 dpi through Slurm job 52371207 with editable Arial SVG text. Added sample-composition panels through Slurm job 52371553 using the same final denominator and Bershteyn 2023 shorthand sample labels. |

## Figure: fig_cross_study_marker_expression_v12

### Status

```text
Found: yes
Log-audited: yes
Confirmed: pending user visual review
Modified: no
Validated: partial; cluster UMAP QC and DIV90 published Fig. D recode outputs validated
Finalized: partial; first final folder created for cluster UMAP QC / DIV90 published recode outputs
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

DIV90 published Fig. D recode assets:

```text
cluster_qc/div90_published_fig_d_10_class_recode.tsv
cluster_qc/div90_published_fig_d_cluster_counts.tsv
cluster_qc/div90_published_fig_d_sample_composition.tsv
cluster_qc/div90_published_fig_d_10_class_umap.png
cluster_qc/div90_published_fig_d_10_class_umap.pdf
cluster_qc/div90_published_fig_d_10_class_umap.svg
cluster_qc/div90_published_fig_d_sample_composition.png
cluster_qc/div90_published_fig_d_sample_composition.pdf
cluster_qc/div90_published_fig_d_sample_composition.svg
```

Cross-study source/GEO/sample provenance table:

```text
cluster_qc/cross_study_marker_expression_v12_geo_sample_provenance.tsv
```

This table records each configured study's external accession/source, whether it
was used in the v12 figure-default UMAP QC plot, loaded/plotted cell counts,
sample IDs loaded/plotted, and any v12 subsetting/filtering decisions. DIV30 and
DIV90 intentionally leave GEO accession fields blank because they are this-study
internal datasets.

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
    DIV90 first joins the audited current/raw cluster_number_name mapping table.
  - Added DIV90 published Fig. D recode layer on top of current/raw cluster IDs.
    Raw/current DIV90 cluster IDs are provenance only and must not be treated as
    the published 10-class annotation.
  - DIV90 current clusters 6 and 7 are excluded from published-style outputs as
    Stressed Cells. The Stressed Cells label is allowed only in the recode audit
    table, not in published-style UMAP legends or sample composition plots.
  - DIV90 published-style visualizations now remove excluded cells only at
    plotting time. Audit/provenance tables retain the excluded current cluster
    rows and reasons.
  - Added plotting-only published orientation coordinates for the DIV90 Fig. D
    UMAP. Original umap_1/umap_2 values are unchanged. The transform is:
      UMAP1_published = umap_1
      UMAP2_published = -1 * umap_2
    The same vertical plotting-only DIV90 transform is now applied to the DIV90
    panel inside the cross-study UMAP QC grids as well as the standalone DIV90
    published Fig. D UMAP.
  - Applied UMAP formatting pass:
    all UMAP scatter plots use the same larger point size, UMAP axes are fully
    hidden, multi-study QC panels are arranged in a single row, and panel titles
    report cluster count and cell count.
  - Added source/GEO/sample provenance table for all configured v12 cross-study
    UMAP QC studies:
    cluster_qc/cross_study_marker_expression_v12_geo_sample_provenance.tsv
  - Created first structured final figure package:
    /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/final_figures/fig_cross_study_marker_expression_v12
    This package currently freezes the cluster UMAP QC grids, DIV90 published
    recode UMAP/composition outputs, tables, logs, code, and provenance. The
    original marker-expression multi-grid panels remain open for formatting and
    later final packaging.
  - Rerendered the first final package outputs with publication-editable SVG
    text:
      font.family = Arial
      svg.fonttype = none
      pdf.fonttype = 42
    All four packaged SVGs now contain editable <text> elements with Arial
    declarations and no DejaVu glyph-path text definitions.
  - Rerendered the first final package outputs with consistent 600 dpi export
    for PNG/PDF/SVG save calls. This matters because UMAP points are rasterized
    inside PDF/SVG; savefig dpi controls the embedded raster resolution.
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
  Current/raw cluster labels:
    0  - MGE Striatal/GP Fated
    1  - SST+, NPY +, Cortical Fated
    2  - CRABP1+/PV Precursors
    3  - PV precursors/Migrating cells/Cortical-fated
    4  - Pre-Astrocytes/Astrocytes 1
    5  - LHX8+ vMGE GABergic Striatal/GP fated 1
    6  - Stressed Cells
    7  - Stressed Cells
    8  - LHX8+ vMGE GABergic Striatal/GP fated 2
    9  - Pre-OPCs/OPCs
    10 - Pre-Astrocytes/Astrocytes 2
    11 - PV Precursors
    12 - Dividing cells
```

DIV90 published Fig. D recode:

```text
Current raw cluster 0  -> 3. MGE Striatal/GP fated
Current raw cluster 1  -> 1. SST+, NPY+ Cortical fated
Current raw cluster 2  -> 2. CRABP1+/PV Precursors
Current raw cluster 3  -> 7. PV Precursors/Migrating cells/Cortical fated
Current raw cluster 4  -> 8. Pre-Astrocytes/Astrocytes
Current raw cluster 5  -> 4. LHX8+ vMGE GABAergic Striatal/GP fated 1
Current raw cluster 6  -> excluded; Stressed Cells
Current raw cluster 7  -> excluded; Stressed Cells
Current raw cluster 8  -> 5. LHX8+ vMGE GABAergic Striatal/GP fated 2
Current raw cluster 9  -> 10. Pre-OPCs/OPCs
Current raw cluster 10 -> 8. Pre-Astrocytes/Astrocytes
Current raw cluster 11 -> 6. PV Precursors
Current raw cluster 12 -> 9. Dividing cells

Published order:
  1. SST+, NPY+ Cortical fated
  2. CRABP1+/PV Precursors
  3. MGE Striatal/GP fated
  4. LHX8+ vMGE GABAergic Striatal/GP fated 1
  5. LHX8+ vMGE GABAergic Striatal/GP fated 2
  6. PV Precursors
  7. PV Precursors/Migrating cells/Cortical fated
  8. Pre-Astrocytes/Astrocytes
  9. Dividing cells
  10. Pre-OPCs/OPCs
```

### Validation

```text
Pending user visual review of candidate v12 plots.
Cluster QC grid was visually spot-checked after render. The figure-default grid
uses study table order and shows cluster IDs over each study's UMAP.
Second QC render confirmed DIV30 collapsed to 5 mapped paper/manual classes and
DIV90 joined to the 13 audited cluster_number_name labels.
DIV90 published Fig. D outputs validated after rerender:
  - 22,338 DIV90 cells before published recode keep filter.
  - 20,049 DIV90 cells after excluding current clusters 6 and 7.
  - 2,289 cells removed as Stressed Cells.
  - Published output TSV/SVG files have no "Stressed" or "EXCLUDED" matches.
  - UMAP and sample-composition PNGs were visually spot-checked.
  - Final published-style legend/bar plot contains exactly the 10 published
    Fig. D classes.
DIV90 published Fig. D orientation validation:
  - Plotting-only coordinates are UMAP1_published = umap_1 and
    UMAP2_published = -1 * umap_2.
  - Vertical-only flip requested after visual review; this is a plotting-only
    transform and does not modify original UMAP embeddings.
  - Visualization SVGs and published-style composition/count TSVs have no
    "Stressed" or "EXCLUDED" matches.
UMAP formatting validation:
  - Multi-study figure-default QC plot regenerated as one row.
  - Multi-study all-prepared QC plot regenerated as one row.
  - DIV90 panel in both multi-study QC grids uses the same vertical
    plotting-only flip as the standalone DIV90 published Fig. D UMAP.
  - DIV90 published Fig. D UMAP title reports "10 clusters, n=20,049".
  - UMAP scatter point size is shared across the cross-study grids and DIV90
    published UMAP.
Editable SVG validation:
  - cross_study_marker_expression_v12_cluster_umap_qc_figure_default.svg:
    editable text present, Arial declarations present, no DejaVu glyph paths.
  - cross_study_marker_expression_v12_cluster_umap_qc_all_prepared_cells.svg:
    editable text present, Arial declarations present, no DejaVu glyph paths.
  - div90_published_fig_d_10_class_umap.svg:
    editable text present, Arial declarations present, no DejaVu glyph paths.
  - div90_published_fig_d_sample_composition.svg:
    editable text present, Arial declarations present, no DejaVu glyph paths.
Resolution validation:
  - First final package was rerendered and recopied at 600 dpi for all four
    figure families:
      cross_study_marker_expression_v12_cluster_umap_qc_figure_default
      cross_study_marker_expression_v12_cluster_umap_qc_all_prepared_cells
      div90_published_fig_d_10_class_umap
      div90_published_fig_d_sample_composition
  - Packaged SVG files are larger after 600 dpi rerender, consistent with
    higher-resolution embedded raster UMAP layers.
```

### Final Figure Package

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/final_figures/fig_cross_study_marker_expression_v12

Created: 2026-06-25
Status: initial final package created for cluster UMAP QC / DIV90 published Fig.
D recode outputs. This is the first final_figures folder and should be used as
the template for subsequent multi-grid plot finalization.

Included structure:
  README.md
  code/
  figures/png/
  figures/pdf/
  figures/svg/
  tables/
  logs/
  provenance/

Included figure families:
  cross_study_marker_expression_v12_cluster_umap_qc_figure_default
  cross_study_marker_expression_v12_cluster_umap_qc_all_prepared_cells
  div90_published_fig_d_10_class_umap
  div90_published_fig_d_sample_composition

Included provenance:
  code copies for QC and original v12 marker-expression workflow
  v12 source study config
  GEO/sample provenance table
  copied v12 Slurm logs
  final render command note
  git commit/status snapshot
  file manifest
  sha256 manifest
```

### Open Questions

```text
Which exact panels should be finalized?
What formatting changes are needed for the original marker-expression multi-grid
panels before final render?
Should SVG export be added to the original marker-expression plotting function
for this figure family?
```

## Figure: fig_cross_study_marker_expression_pv_precursors_on_off_target_v12

### Status

```text
Found: yes
Log-audited: yes; source v12 workflow logs were previously checked and this
  update used the already prepared v12 marker-expression tables.
Confirmed: pending user visual review
Modified: yes
Validated: partial; output manifest and PNG visual spot-check completed
Finalized: candidate package created
```

### Candidate Paths

Original source plot:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_marker_expression/cross_study_marker_expression_v12/plots/cross_study_marker_expression_pv_precursors_on_off_target.pdf
```

Modified candidate run:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_marker_expression/cross_study_marker_expression_v12_pv_precursors_final_candidate
```

Final-figure candidate package:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/final_figures/fig_cross_study_marker_expression_pv_precursors_on_off_target_v12_candidate
```

Target outputs:

```text
figures/png/cross_study_marker_expression_pv_precursors_on_off_target.png
figures/pdf/cross_study_marker_expression_pv_precursors_on_off_target.pdf
figures/svg/cross_study_marker_expression_pv_precursors_on_off_target.svg
```

### Modification Log

```text
2026-06-25:
  - Copied existing v12 prepared marker-expression tables into an isolated
    candidate run folder. No expression extraction, UMAP recomputation,
    reclustering, or annotation updates were run.
  - Added Bershteyn et al. 2025 to the plotted study set by running the plot
    command without the previous Bershteyn 2025 study exclusion.
  - Updated the PV precursor ON-target gene list:
      MAFB, MEF2C, ERBB4, LHX6, LHX8, NKX2-1, ETV1, CRABP1, TAC1, ST18, PVALB
    ERBB4 was already present and was retained. NKX2-1 is displayed as NKX2.1.
  - Added plotting-only coordinates for marker-expression grids:
      UMAP1_plot = umap_1
      UMAP2_plot = umap_2 for all studies except DIV90
      DIV90 UMAP2_plot = -1 * umap_2
    Original umap_1/umap_2 columns are unchanged.
  - Added DIV90 visualization-only stressed-cell filter for marker-expression
    grids, matching the earlier published-style DIV90 logic:
      current clusters 6 and 7 removed from plotted DIV90 row
      22,338 cells before filter
      20,049 cells after filter
      2,289 cells removed
    Source marker-expression tables and original cluster annotations are not
    modified.
  - Changed the marker-expression colormap and thresholding so expression
    values from 0 through 1 remain background gray (#d0d0d0). The blue overlay
    and color scale start above expression value 1.
  - Corrected marker-expression colorbar labeling:
      colorbar begins at 0
      values from 0 through 1 are drawn as the background gray floor
      blue overlay starts above expression value 1
      right tick is each gene's q99 positive-expression upper limit
  - Updated Matplotlib export settings for publication editing:
      font.family = Arial
      svg.fonttype = none
      pdf.fonttype = 42
    Great Lakes does not have Microsoft Arial installed locally, so render logs
    may show findfont warnings. SVG output is still written with editable
    <text> elements styled as Arial for Illustrator.
  - Added SVG export to the marker-expression plotting function.
  - Initial plot-only command used for earlier candidate versions:
      PYTHONPATH=python_notebooks/src /home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python python_notebooks/scripts/run_cross_study_marker_expression.py --project-root /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder --run-label cross_study_marker_expression_v12_pv_precursors_final_candidate plot
  - Refreshed lightweight setup metadata tables after plotting so
    cross_study_marker_expression_gene_panels.tsv and
    cross_study_marker_expression_genes.tsv match the expanded PV precursor
    gene panel.
  - Runtime was a few minutes; SVG/PDF writing was the slowest step.
  - 2026-06-25 high-DPI update attempt:
      Interactive 600 dpi marker-expression rerender was killed.
      Interactive targeted 600 dpi and 450 dpi rerenders were also killed.
      Do not continue high-DPI marker-expression final renders interactively.
      Submit through Slurm with sufficient memory/time, following the existing
      project Slurm environment logic: CONDA_ENV_BIN points to
      /home/elcrespo/miniconda3/envs/mge-organoid-python/bin, PYTHONPATH includes
      REPO_ROOT/python_notebooks/src, and the command runs from REPO_ROOT.
      Sync outputs only after PNG/PDF/SVG are all nonzero and logs show a clean
      exit.
      A failed interactive attempt left partial files in the candidate run
      plots directory; check and replace source outputs before packaging.
  - 2026-06-25 Slurm final render completed successfully:
      Job ID: 52370542
      Submitted sbatch:
        /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/render_pv_precursor_final_panel_600dpi.sbatch
      Logs:
        /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/pv-final-panel-pv-final-panel-52370542.out
        /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/pv-final-panel-pv-final-panel-52370542.err
      Requested resources:
        4 CPUs, 180G memory, 3 hours, standard partition, account parent0
      DPI:
        PNG = 600
        PDF = 600
        SVG = 600
      Slurm script validated PNG/PDF/SVG were nonzero before copying into
      final_figures.
```

### Validation

```text
Manifest check:
  - Target panel genes:
    MAFB, MEF2C, ERBB4, LHX6, LHX8, NKX2-1, ETV1, CRABP1, TAC1, ST18, PVALB,
    SP8, EBF1, NKX2-2, RAX, HMX3, DBH
  - Target panel studies:
    This Study DIV30, This Study DIV90, Siebert 2026, Walsh 2025,
    Bershteyn 2025, Bershteyn 2023, Samarasinghe 2021
  - Excluded study IDs: none
  - Target outputs exist as PNG, PDF, and SVG.
  - DIV90 filter summary confirms 22,338 cells before filter, 20,049 after
    filter, and 2,289 current-cluster 6/7 stressed cells removed.

Visual spot-check:
  - PNG opened successfully.
  - Bershteyn et al. 2025 appears as a new row.
  - DIV90 row uses the vertical-only plotting orientation.
  - DIV90 row label reports "Stressed cells removed" and n = 20,049.
  - Values from 0 through 1 are gray; blue expression overlay starts above 1.
  - Colorbars now start at 0 while retaining the 0-1 gray floor.
  - Packaged SVG has editable <text> elements, font-family Arial declarations,
    and no DejaVu glyph-path text definitions.
  - Gene header displays NKX2.1.
  - High-DPI final render completed via Slurm at 600 dpi for PNG/PDF/SVG.
  - Slurm output log records nonzero output validation:
      PNG 7.8M
      PDF 7.3M
      SVG 11M
```

### Final Figure Package

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/final_figures/fig_cross_study_marker_expression_pv_precursors_on_off_target_v12_candidate

Created: 2026-06-25
Status: candidate package created; awaiting user visual approval before calling
it final.

Included structure:
  README.md
  code/
  figures/png/
  figures/pdf/
  figures/svg/
  tables/
  logs/
  provenance/

Included provenance:
  copied plotting code and CLI wrapper
  targeted final-panel render script
  submitted 600 dpi Slurm script
  plot-only Slurm template
  plot manifest and audit tables
  Slurm stdout/stderr logs
  final render command note
  git commit/status snapshot
  file manifest
  sha256 manifest
```

## Figure: fig_cross_study_shi_label_transfer_v1_umap_score_grids

### Status

```text
Found: yes
Log-audited: yes
Confirmed: pending user visual review
Modified: yes
Validated: yes
Finalized: candidate package created; awaiting user visual approval
```

### Candidate Paths

Source run:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_shi_seurat_label_transfer/cross_study_shi_seurat_label_transfer_v1
```

Requested candidate plots:

```text
plots/umap_grids/cross_study_umap_shi_seurat_full_individual_gw_scores_grid.pdf
plots/umap_grids/cross_study_umap_shi_seurat_full_all_label_scores_grid.pdf
```

Final candidate package:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/final_figures/fig_cross_study_shi_label_transfer_v1_umap_score_grids_candidate
```

### Source Code and Handoffs

Primary plotting module:

```text
python_notebooks/src/mge_organoid_python/cross_study_shi_prediction_plots.py
```

CLI wrapper:

```text
python_notebooks/scripts/run_cross_study_shi_prediction_plots.py
```

Original Seurat label-transfer workflow that created the saved tables:

```text
scripts/13_run_cross_study_shi_seurat_label_transfer.R
```

Slurm render template added for final packaging:

```text
slurm_templates/49_render_cross_study_shi_label_transfer_final_umap_grids.sbatch.template
```

Relevant broader handoff:

```text
python_notebooks/HANDOFF_shi_reference_div30_label_transfer.md
```

### Prepared Assets

Primary plot input used by the final render:

```text
tables/cross_study_shi_seurat_label_transfer_obs.tsv.gz
```

Per-study source tables remain in:

```text
tables/per_study/*_shi_seurat_label_transfer_obs.tsv.gz
seurat/per_study/*_shi_seurat_full_prediction_scores.tsv.gz
seurat/per_study/*_shi_seurat_full_week_prediction_scores.tsv.gz
```

Data-processing path for the final plot-only render:

```text
plot_umap_only()
  load_combined_table()
    normalize_obs_table()
    compute_expected_gw()
    validate_combined()
  make_umap_grids()
    apply_internal_umap_plot_filters()
    add_plot_coordinates()
    downsample_by_study(max_cells_per_study=None)
    plot_continuous_umap_grid()
```

The two requested grids are generated by `plot_continuous_umap_grid()` using:

```text
all_label_scores_grid:
  features = label_score_columns(plot_data)
  vmin = 0
  vmax = 1

individual_gw_scores_grid:
  features = week_score_columns(plot_data)
  vmin = 0
  vmax = 1
```

### Log Audit

Earlier v1 plot-only timing:

```text
Job: cross-shi-plot 51537151
Start: 2026-06-08 20:02:27 EDT
End: 2026-06-08 20:05:19 EDT
Runtime: about 3 minutes
Result: completed plot-only UMAP grids from existing prediction tables
```

Final render:

```text
Job ID: 52371207
Submitted sbatch:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/49_render_cross_study_shi_label_transfer_final_umap_grids.sbatch
Logs:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/shi-final-umaps-52371207.out
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/shi-final-umaps-52371207.err
Runtime: 00:06:10
Max RSS: 12,191,916K
Exit code: 0:0
Resources requested: 4 CPUs, 64G memory, 2 hours, standard partition, account parent0
```

### Rerun Decision

```text
Rerun type: plot-only.

Reason:
  - Saved combined and per-study label-transfer tables already existed.
  - Earlier v1 plot-only Slurm logs showed the UMAP grids could be regenerated in
    a few minutes.
  - Requested changes were formatting/visualization-only:
      score colormap
      DIV90 visualization filter
      DIV90 plotting orientation
      output formats/DPI/SVG text

Not rerun:
  - Seurat label transfer
  - UMAP computation
  - clustering
  - annotation recomputation
```

### Modification Log

```text
2026-06-26:
  - Changed continuous Shi score UMAP grids from white-to-blue to grey-to-blue.
    The score range remains fixed at 0..1 for every panel.
  - Set background/zero-score grey to #d0d0d0 to match the final marker-expression
    figure style.
  - Added DIV90 visualization-only stressed-cell exclusion:
      current clusters 6 and 7 removed from plotted DIV90 row
      22,338 cells before filter
      20,049 cells after filter
      2,289 cells removed
    The saved source tables and original annotations are not modified.
  - Added plotting-only coordinates:
      UMAP1_plot = umap_1
      UMAP2_plot = umap_2 for all studies except DIV90
      DIV90 UMAP2_plot = -1 * umap_2
    Original umap_1/umap_2 columns are unchanged.
  - Added 600 dpi savefig export for PNG/PDF/SVG.
  - Added SVG export for continuous UMAP score grids.
  - Wrapped long all-label grid headers:
      Excitatory IPC
      Excitatory neuron
      Thalamic neurons
    and widened continuous score grids slightly to prevent column-title overlap.
  - Added Matplotlib text export settings:
      font.family = Arial
      svg.fonttype = none
      pdf.fonttype = 42
      ps.fonttype = 42
  - Added Slurm final render/package template:
      slurm_templates/49_render_cross_study_shi_label_transfer_final_umap_grids.sbatch.template
  - 2026-06-26 sample-composition addition:
      Added two matched all-study sample-composition panels to the same final
      package:
        all_studies_all_shi_major_labels_shi_predicted_age_sample_composition_stacked_bar
        all_studies_all_shi_major_labels_shi_major_class_sample_composition_stacked_bar
      These were recalculated from the v1 combined table rather than copied from
      the old source-run PDFs.
      Both use the same final visualization denominator as the UMAP score grids:
        DIV90 current clusters 6/7 excluded as stressed cells
        Samarasinghe controls only
        finite plotting coordinates required
        no downsampling
        no score cutoff
      Predicted-stage plot uses:
        shi_seurat_full_predicted_shi_week_label
      Major-class plot uses:
        shi_seurat_full_predicted_shi_label
      Slurm job:
        52371553
      Runtime:
        00:00:26
      Max RSS:
        1,264,748K
      2026-06-26 label update:
        Bershteyn et al. 2023 sample labels were updated to the requested
        shorthand style in both sample-composition figures and exported TSVs:
        D0 -> DIV0 hESC
        D14 -> DIV14 NPC
        MB460 -> DIV42 EOP L1 U
        MB461 -> DIV42 EOP L1 S
        MB279 -> DIV42 EOP L2 U
        MB280 -> DIV42 EOP L2 S
        MB527 -> DIV42 EOP L3 U
        MB528 -> DIV42 EOP L3 S
      Added code:
        python_notebooks/scripts/render_cross_study_shi_sample_composition_final_panels.py
      Added Slurm template:
        slurm_templates/49b_render_cross_study_shi_sample_composition_final_panels.sbatch.template
```

### Validation

```text
Final package outputs exist and are nonzero:
  figures/png/cross_study_umap_shi_seurat_full_all_label_scores_grid.png
  figures/pdf/cross_study_umap_shi_seurat_full_all_label_scores_grid.pdf
  figures/svg/cross_study_umap_shi_seurat_full_all_label_scores_grid.svg
  figures/png/cross_study_umap_shi_seurat_full_individual_gw_scores_grid.png
  figures/pdf/cross_study_umap_shi_seurat_full_individual_gw_scores_grid.pdf
  figures/svg/cross_study_umap_shi_seurat_full_individual_gw_scores_grid.svg

PNG dimensions:
  all_label_scores_grid.png = 6521 x 4975
  individual_gw_scores_grid.png = 4181 x 4989

DIV90 filter summary:
  before = 22,338
  after = 20,049
  removed = 2,289
  plot_filter = exclude_div90_stressed_clusters_6_7
  plot_filter_label = Stressed cells removed

Manifest:
  - DIV90 rows in both requested grids report n_cells_plotted = 20,049.
  - All requested score panels report vmin = 0 and vmax = 1.

SVG text:
  - all_label_scores_grid.svg has 50 editable <text> elements.
  - individual_gw_scores_grid.svg has 32 editable <text> elements.
  - Both SVGs contain font-family: 'Arial' text declarations.
  - No DejaVu/Nimbus/Arial glyph-path text definitions were found.

Sample-composition validation:
  - Predicted-stage PNG/PDF/SVG exist and are nonzero.
  - Major-class PNG/PDF/SVG exist and are nonzero.
  - Predicted-stage table:
      tables/cross_study_shi_final_all_shi_major_labels_predicted_age_sample_composition.tsv
  - Major-class table:
      tables/cross_study_shi_final_all_shi_major_labels_major_class_sample_composition.tsv
  - Filter summary:
      tables/cross_study_shi_sample_composition_final_filter_summary.tsv
  - DIV90 final denominator = 20,049.
  - Samarasinghe final denominator = 26,935 controls only.
  - Both sample-composition PNGs opened successfully.
  - Predicted-stage SVG has 71 editable Arial <text> elements and no glyph-path text.
  - Major-class SVG has 76 editable Arial <text> elements and no glyph-path text.
```

### Final Figure Package

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/final_figures/fig_cross_study_shi_label_transfer_v1_umap_score_grids_candidate

Included:
  README.md
  code/cross_study_shi_prediction_plots.py
  code/run_cross_study_shi_prediction_plots.py
  figures/png/
  figures/pdf/
  figures/svg/
  tables/cross_study_shi_umap_plot_manifest.tsv
  tables/cross_study_shi_umap_internal_plot_filter_summary.tsv
  tables/cross_study_shi_seurat_label_transfer_readiness.tsv
  tables/cross_study_shi_transfer_diagnostics_summary.tsv
  tables/cross_study_shi_final_all_shi_major_labels_predicted_age_sample_composition.tsv
  tables/cross_study_shi_final_all_shi_major_labels_major_class_sample_composition.tsv
  tables/cross_study_shi_sample_composition_final_filter_summary.tsv
  tables/cross_study_shi_final_sample_composition_manifest.tsv
  logs/shi-final-umaps-52371207.out
  logs/shi-final-umaps-52371207.err
  logs/shi-final-samples-52371553.out
  logs/shi-final-samples-52371553.err
  provenance/render_manifest.tsv
  provenance/code_diff.patch
  provenance/sample_composition_code_diff.patch
  provenance/git_status_short.txt
  provenance/sha256sums.txt
```

### Open Questions

```text
- User visual approval is still needed before marking this final.
- If the user wants only these two grids preserved in the source run folder, no
  further action is needed. The plot-only rerender also refreshed the other v1
  UMAP-grid outputs in the source run because the existing command renders all
  UMAP grids together.
```

## Exploratory Validation Diagnostic: Shi Query Routing Static

```text
Status:
  Exploratory diagnostic rendered; not a finalized publication figure.

Purpose:
  Visualize how DIV30 and DIV90 UMAP regions are being sorted by saved Shi
  Seurat TransferData predictions without rerunning Seurat transfer.

Output directory:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_shi_seurat_label_transfer/cross_study_shi_seurat_label_transfer_v1/plots/validation/shi_query_projection_routing_static

Script:
  python_notebooks/scripts/render_shi_query_routing_static_diagnostic.py

Slurm template:
  slurm_templates/49c_render_shi_query_routing_static_diagnostic.sbatch.template

Final Slurm job:
  52371684

Runtime:
  00:00:47

Max RSS:
  2,824,996K

Key outputs:
  shi_query_projection_routing_predicted_stage_static.png/pdf/svg
  shi_query_projection_routing_major_class_static.png/pdf/svg
  shi_query_projection_routing_predicted_stage_routes.tsv
  shi_query_projection_routing_major_class_routes.tsv

Interpretation:
  Query cells are plotted on their DIV30/DIV90 UMAP plane. Source nodes are
  recoded DIV30/DIV90 class centroids. Shi target nodes are plotted on a raised
  target plane. Route width represents fraction of the source class assigned to
  the target, and route opacity represents mean prediction score.

Important limitation:
  This is not a true Seurat anchor/reference-cell projection. The v1 workflow
  saved per-cell TransferData predictions and scores, not the Seurat anchor
  object or top reference-cell matches.
```

## Pause Point: True Shi Anchor Projection Assets

```text
Date:
  2026-06-26

User request:
  User wants the true Seurat anchor/reference-cell projection assets, not only
  the score-routing diagnostic. Specifically, rerun only what is needed to save
  the missing anchors/projection links for DIV30 and DIV90 against Shi et al.

Important clarification:
  The finalized figures DO have saved plot-ready assets:
    - Shi final UMAP grids/sample composition use saved v1 prediction tables.
    - Marker-expression figures use saved per-study marker-expression tables.
  The missing asset is narrower:
    - Seurat TransferAnchorSet / true query-to-Shi-reference-cell anchor links
      were not saved by the original v1 Shi workflow.

Why rerun is needed:
  scripts/13_run_cross_study_shi_seurat_label_transfer.R creates:
    anchors <- Seurat::FindTransferAnchors(...)
  then uses those anchors for TransferData and writes prediction-score TSVs.
  It does not save anchors. At the end of run_transfer_one(), it removes:
    rm(query, anchors, ge_anchors, ...)

Current implementation started:
  Added targeted helper script:
    scripts/14_save_shi_query_anchor_projection_assets.R

Current status:
  - Script has been written but not yet Slurm-submitted.
  - No anchor assets have been generated yet.
  - No Slurm template for this new anchor-asset job has been added yet.
  - Local/login-node R inspection failed/OOMed with exit 137 when trying to
    read large RDS objects, so this must be run on Slurm.

Intended output run label:
  cross_study_shi_seurat_anchor_projection_v1

Intended output directory:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_shi_seurat_label_transfer/cross_study_shi_seurat_anchor_projection_v1

Planned scope:
  - Only studies:
      varela_div30
      varela_div90
  - Whole-Shi FindTransferAnchors matching script 13 logic.
  - Save true anchor objects and tables.
  - Do not overwrite the finalized v1 prediction tables or final_figures.

Planned saved assets:
  seurat/anchors/varela_div30_shi_full_transfer_anchors.rds
  seurat/anchors/varela_div90_shi_full_transfer_anchors.rds
  tables/varela_div30_shi_full_anchor_pairs.tsv.gz
  tables/varela_div90_shi_full_anchor_pairs.tsv.gz
  tables/varela_div30_shi_full_top_anchor_per_query.tsv.gz
  tables/varela_div90_shi_full_top_anchor_per_query.tsv.gz
  tables/*_shi_reference_coordinates_for_anchor_plot.tsv.gz
  tables/*_query_coordinates_for_anchor_plot.tsv.gz
  plots/*_shi_full_anchor_projection_static.png/pdf/svg
  diagnostics/shi_anchor_projection_diagnostics.tsv
  README.md

Important behavior:
  - Full anchor objects/tables should include all query cells.
  - DIV90 static plots should exclude current query clusters 6/7 for visual
    consistency with finalized figures, but this visualization filter should
    not remove cells from saved anchor objects/full anchor-pair tables.
  - This output is distinct from the score-routing diagnostic:
      score-routing = saved TransferData predictions/scores
      true projection = saved Seurat anchor links/reference-cell matches

Next steps tomorrow:
  1. Read/review scripts/14_save_shi_query_anchor_projection_assets.R.
  2. Create a Slurm template, likely:
       slurm_templates/49d_save_shi_query_anchor_projection_assets.sbatch.template
     using the Seurat-capable env/module pattern from existing Seurat jobs.
     The conda Rscript exists at:
       /home/elcrespo/miniconda3/envs/mge-organoid-python/bin/Rscript
     Existing Seurat array templates also use:
       module load Bioinformatics
       module load r-seurat/5.1.0-R-4.4.1-c3m7yfq
  3. Run an R parse check in the chosen Slurm environment if possible.
  4. Submit via Slurm, not locally.
  5. Watch logs, timing, memory.
  6. Validate outputs:
       anchor RDS files exist and are nonzero
       anchor-pair TSVs exist and have cell1/cell2/query/reference IDs
       diagnostics report n_anchors and unique query/reference anchor cells
       static PNG/PDF/SVG render
       README/provenance/checksums are present
  7. Then decide whether to make a polished publication-style anchor plot from
     those true assets.

Git/worktree note at pause:
  Untracked:
    scripts/14_save_shi_query_anchor_projection_assets.R
```

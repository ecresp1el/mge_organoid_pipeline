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
| `fig_cross_study_marker_expression_pv_precursors_on_off_target_v12` | Updated, Slurm-rendered, Validated package | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_marker_expression/cross_study_marker_expression_v12_pv_precursors_final_candidate` | 2026-06-27 revision groups PV ON-target MGE genes first as NKX2.1/LHX6/LHX8/ERBB4 with a top `MGE` span, replaces ST18 with ZEB2, removes HMX3, and refreshed marker extraction because ZEB2 was absent from the prepared v12 tables. Exported PNG/PDF/SVG at 600 dpi through Slurm job 52439584. |
| `fig_cross_study_shi_label_transfer_v1_umap_score_grids` | Modified, Slurm-rendered, Packaged candidate | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_shi_seurat_label_transfer/cross_study_shi_seurat_label_transfer_v1` | Plot-only rerender of Shi Seurat label-transfer UMAP score grids and matched sample-composition panels from saved v1 combined table. Updated score grids to fixed 0-1 grey-to-blue scaling, applied DIV90 visualization-only stressed-cluster 6/7 removal, applied DIV90 vertical plotting-only UMAP orientation, wrapped long all-label headers, exported requested PNG/PDF/SVG at 600 dpi through Slurm job 52371207 with editable Arial SVG text. Added sample-composition panels through Slurm job 52371553 using the same final denominator and Bershteyn 2023 shorthand sample labels. |
| `maturation_scores` | Slurm-rendered, Validated package | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/maturation_scores/maturation_scores_v1` | Unified DIV30/DIV90 Jia RGC1/RGC2/IPC Scanpy `score_genes` UMAP overlays plus derived `IPC - mean(RGC1,RGC2)` maturation index. Final package is `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/final_figures/maturation_scores`; Slurm job 52441418 completed 2026-06-27 with PNG/PDF/SVG plus score and audit tables. |
| `fig_div90_jia_urd_marker_pseudotime_tree_v1_candidate` | Found, Confirmed, Final package started | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_allcells_jia_root10_neuron_s9_7tips_urd_resumable_v1` | New candidate package created for the DIV90 Jia all-cell URD marker-validation, UMAP pseudotime, and cluster-number-name tree pseudotime panels. Current assets were copied into `final_figures`; formatting should proceed by plot-only re-render from existing assets, not by recomputing URD pseudotime, the lineage tree, or marker validation. |
| `fig_siletti_div90_restricted_mge_llc_msn_emsn_chat_rivers_v6_no_cge_min10_ordered` | Generated, Validated, Packaged candidate | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_FINAL_FIGURE_CANDIDATE_restricted_mge_llc_msn_emsn_chat_river_plots_v6_no_cge_min10_ordered` | Current preferred Siletti river candidate. v6 removes CGE before bridge export and before unified fast-kNN transfer, keeps MGE, LAMP5-LHX6/chandelier, MSN, eccentric MSN, and Splatter-CHAT, filters plotted river edges with fewer than 10 DIV90 cells, recomputes rectangles/ribbons from the filtered edge tables so tiny populations do not leave boxes, and fixes the left DIV90 cluster order by descending Pallial/cortical fraction. Copied into `final_figures` as a clearly marked candidate package; not final-approved yet. |

## Figure: maturation_scores

### Status

```text
Found: yes
Renderer added: yes
Slurm render: completed as job 52441418 on 2026-06-27
Final validation: passed
```

### Candidate Paths

```text
Old DIV30 Jia candidate:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/jia_program_div30_scoring/jia_program_div30_scoring_v1/plots/div30_umap_seurat_clusters_jia_program_scores_panel.png

Old DIV30 Jia tables:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/jia_program_div30_scoring/jia_program_div30_scoring_v1/tables/

Older combined DIV30/DIV90 RGC/IPC progression path:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/mgeo_rgc_ipc_monocle3/

Later DIV90 Jia-root metadata path:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_allcells_jia_root10_neuron_s9_7tips_urd_resumable_v1/inputs/div90_jia_rootscore_scored_metadata.tsv
```

The old combined RGC/IPC progression workflow used Shi Table S5 RGC/IPC genes
and plain mean log-normalized expression, not Jia RGC1/RGC2/IPC Scanpy module
scores. The later DIV90 URD RootScore metadata also used mean log-normalized
Jia program expression for DIV90. The DIV30 Jia notebook used Scanpy
`score_genes`, so the new final render recalculates DIV30 and DIV90 together
with one scoring contract.

### Source Code and Handoffs

```text
python_notebooks/HANDOFF_jia_program_div30_scoring.md
python_notebooks/notebooks/jia_program_div30_scoring.ipynb
python_notebooks/src/mge_organoid_python/gene_program_scoring.py
python_notebooks/scripts/render_maturation_scores_final.py
slurm_templates/60_render_maturation_scores_final.sbatch.template
```

### Prepared Assets

```text
Marker CSV:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/reference/Jia_et_al_2026_Science_3_progs.csv

DIV30 AnnData:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/varela_div30.h5ad

DIV90 AnnData:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/varela_div90.h5ad
```

### Log Audit

```text
Submitted Slurm job:
sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/60_render_maturation_scores_final.sbatch

Job ID:
52441418

Logs:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/maturation_scores-52441418.out
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/maturation_scores-52441418.err

sacct:
State COMPLETED, ExitCode 0:0, Elapsed 00:03:00, batch MaxRSS 24357580K,
AllocCPUS 8, ReqMem 180G.
```

### Rerun Decision

```text
Rerun type: scoring rerun, not upstream reconstruction.
Reason: DIV30 and DIV90 prior score paths were not identical.
Method: Scanpy score_genes on each timepoint with the same Jia marker CSV,
ctrl_size=50, random_state=0, use_raw=False.
```

### Modification Log

```text
2026-06-27:
- Added a unified final renderer for DIV30 and DIV90 Jia RGC1/RGC2/IPC scores.
- Added derived display-only maturation index:
  jia_score_IPC - mean(jia_score_RGC1, jia_score_RGC2).
- Exports PNG/PDF/SVG at 600 dpi, with rasterized point layers and editable
  SVG text settings.
```

### Validation

```text
Validated 2026-06-27:
- Slurm job 52441418 completed with exit code 0.
- PNG, PDF, and SVG outputs are nonzero.
- PNG is 6908 x 3016 px.
- SVG contains editable `<text>` elements with Arial-family declarations.
- Visual inspection of PNG confirmed two rows (DIV30, DIV90) and four score
  columns (RGC1, RGC2, IPC, IPC-minus-mean-RGC index).
```

### Final Figure Package

```text
Target:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/final_figures/maturation_scores/

Expected outputs:
figures/png/maturation_scores_umap_grid.png
figures/pdf/maturation_scores_umap_grid.pdf
figures/svg/maturation_scores_umap_grid.svg
tables/div30_div90_jia_maturation_scores_obs.tsv.gz
tables/*gene_overlap*
tables/*scanpy_control_gene*
tables/*scanpy_program_gene_bins*
provenance/maturation_scores_provenance.json
```

### Open Questions

```text
- Whether to use only the three primary Jia program scores in the manuscript, or
  also include the derived IPC-minus-RGC maturation index as a compact summary.
- Whether DIV90 UMAP orientation should stay native or use a plotting-only flip
  for visual comparison with other final figures.
```

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
Log-audited: yes; source v12 workflow logs were previously checked, and the
  2026-06-27 marker revision refreshed marker extraction before rendering.
Confirmed: pending user visual review
Modified: yes
Validated: yes; manifest/table checks and PNG visual spot-check completed
Finalized: candidate package refreshed
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
  - 2026-06-27 requested marker revision:
      PV ON-target marker order now starts with the MGE block:
        NKX2-1, LHX6, LHX8, ERBB4
      NKX2-1 is displayed as NKX2.1.
      Added a horizontal top span labeled "MGE" over those four columns.
      Replaced ST18 with ZEB2.
      Removed HMX3 from the paired OFF-target marker list.
      Existing prepared v12 marker-expression tables contain ST18 and HMX3 but
      not ZEB2, so this revision requires fresh H5AD/Seurat marker extraction
      before the final 600 dpi render can be packaged.
  - 2026-06-27 Slurm refresh completed successfully:
      Job ID: 52439584
      Submitted sbatch:
        /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/render_pv_precursor_final_panel_mge_zeb2_600dpi.sbatch
      Logs:
        /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/pv-mge-zeb2-pv-mge-zeb2-52439584.out
        /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/pv-mge-zeb2-pv-mge-zeb2-52439584.err
      Runtime/resources:
        00:05:16 elapsed, 34,003,032K MaxRSS, 8 CPUs, 180G memory requested
      DPI:
        PNG = 600
        PDF = 600
        SVG = 600
      Outputs were copied to final_figures after PNG/PDF/SVG were confirmed
      nonzero.
```

### Validation

```text
Manifest check:
  - 2026-06-27 target panel genes:
    NKX2-1, LHX6, LHX8, ERBB4, MAFB, MEF2C, ETV1, CRABP1, TAC1, ZEB2, PVALB,
    SP8, EBF1, NKX2-2, RAX, DBH
  - Manifest confirms ZEB2 is present and ST18/HMX3 are absent.
  - Refreshed per-study marker-expression tables contain ZEB2 and do not
    contain ST18 or HMX3.
  - Target outputs exist as PNG, PDF, and SVG:
      PNG 6.8M, 8047 x 5026
      PDF 6.5M
      SVG 9.3M
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
  - 2026-06-27 PNG spot-check confirms the first four columns are NKX2.1,
    LHX6, LHX8, and ERBB4 with a horizontal MGE span above them; ZEB2 appears
    in place of ST18; HMX3 is absent; and the ON/OFF divider remains after
    PVALB.
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
  - Script has been written.
  - Slurm template has been written.
  - Combined Slurm job 52371983 was submitted, started, and was cancelled by
    user request because the user prefers separate DIV30 and DIV90 jobs.
  - Replacement separate Slurm jobs completed successfully:
      DIV30: job 52372021, COMPLETED, 00:14:29, MaxRSS 99904956K
      DIV90: job 52372022, COMPLETED, 00:07:07, MaxRSS 37740076K
  - True Seurat anchor assets were generated for both DIV30 and DIV90:
      anchor RDS files saved
      anchor-pair TSVs saved
      top-anchor-per-query TSVs saved
      query/reference coordinate TSVs saved
      diagnostics, README, provenance, logs, and checksums saved
  - Static anchor line plots did not render for either study. Diagnostics report
    static_plot_rendered = FALSE. Logs indicate no sampled line links survived
    the coordinate merge. This is a plotting/ID-merge issue to inspect from
    the saved TSV/RDS assets; it does not mean the anchor computation failed.
  - Follow-up plot-only renderer was added and run successfully from saved
    assets only:
      python_notebooks/scripts/render_shi_true_anchor_projection_from_assets.py
      slurm_templates/49e_render_shi_true_anchor_projection_from_assets.sbatch.template
    Final plot-only Slurm job:
      52391075, COMPLETED, 00:01:23, MaxRSS 1122144K, ExitCode 0:0
    Output directory:
      /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_shi_seurat_label_transfer/cross_study_shi_seurat_anchor_projection_v1/true_anchor_projection_plots
    This renderer fixes the plotting-only reference ID mismatch by stripping
    the anchor-table suffix "_reference" before joining to the Shi coordinate
    table. It does not rerun FindTransferAnchors.
    Current formatted figure adds lower reference subpanels under each
    projection panel:
      query UMAP colored by saved winner-take-all TransferData label, with
      query cluster numbers overlaid at centroids
      Shi reference UMAP colored by matching Shi major class or canonical
      GW/stage label, with label names printed on the lower reference UMAP
      top query UMAPs, top Shi reference UMAPs, and anchor links share the same
      palette for the selected label mode
    Current renderer uses full saved anchor-pair tables for lines, not the
    top-anchor-only table, and plots all cells from saved coordinate tables
    without point downsampling.
    Winner-take-all label source:
      /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_shi_seurat_label_transfer/cross_study_shi_seurat_label_transfer_v1/tables/cross_study_shi_seurat_label_transfer_obs.tsv.gz
    Plot-only GW canonicalization:
      raw Shi reference labels such as GW12_02 and GW16/12_01 are plotted as
      GW12 and GW16 so the reference UMAP, query WTA labels, and lines share
      the same canonical GW palette.
  - Final-figures candidate package created:
      /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/final_figures/fig_cross_study_shi_true_anchor_projection_v1_candidate
    Package contains:
      figures/ PNG, PDF, SVG outputs
      code/ render script, anchor asset script, original v1 transfer script,
            and Slurm templates 49d/49e
      logs/ final plot Slurm log plus DIV30/DIV90 anchor asset Slurm logs
      tables/ render diagnostics plus DIV30/DIV90 anchor diagnostics and
              selected-study manifest
      provenance/ handoff snapshot, git status/commit, manifests, checksums
      README.md with method-ready documentation and exact code-line references
  - Local/login-node R inspection failed/OOMed with exit 137 when trying to
    read large RDS objects, so any heavy R/Seurat inspection should use Slurm.
  - A first submission with 200G memory was rejected by Slurm because the node
    configuration was unavailable; the accepted separate submissions used 160G,
    matching the existing Shi transfer array memory request.

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
  For separate per-study output directories:
    <run_label>/varela_div30/
    <run_label>/varela_div90/
  Each directory should contain:
    seurat/anchors/<study>_shi_full_transfer_anchors.rds
    tables/<study>_shi_full_anchor_pairs.tsv.gz
    tables/<study>_shi_full_top_anchor_per_query.tsv.gz
    tables/<study>_shi_reference_coordinates_for_anchor_plot.tsv.gz
    tables/<study>_query_coordinates_for_anchor_plot.tsv.gz
    plots/<study>_shi_full_anchor_projection_static.png/pdf/svg
    diagnostics/shi_anchor_projection_diagnostics.tsv
    README.md

Important behavior:
  - Full anchor objects/tables should include all query cells.
  - One whole-Shi Seurat TransferAnchorSet is generated per query study.
  - That single anchor object is shared for both:
      Shi major cell-type interpretation
      Shi GW/stage interpretation
  - The exported anchor-pair table must clearly include both:
      reference_shi_label
      reference_shi_week_label
    This is the explicit bridge allowing the same anchors to be viewed by cell
    type or by GW/stage.
  - DIV90 static plots should exclude current query clusters 6/7 for visual
    consistency with finalized figures, but this visualization filter should
    not remove cells from saved anchor objects/full anchor-pair tables.
  - This output is distinct from the score-routing diagnostic:
      score-routing = saved TransferData predictions/scores
      true projection = saved Seurat anchor links/reference-cell matches

Completed Slurm details:
  Cancelled combined job:
    52371983
  Cancelled state:
    CANCELLED at 00:05:19
  Replacement separate jobs:
    DIV30: 52372021, COMPLETED, 00:14:29, MaxRSS 99904956K, ExitCode 0:0
    DIV90: 52372022, COMPLETED, 00:07:07, MaxRSS 37740076K, ExitCode 0:0
  Output directories:
    /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_shi_seurat_label_transfer/cross_study_shi_seurat_anchor_projection_v1/varela_div30
    /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_shi_seurat_label_transfer/cross_study_shi_seurat_anchor_projection_v1/varela_div90

Completed diagnostics:
  DIV30:
    n_query_cells = 90631
    n_reference_cells = 38831
    n_shared_features = 14354
    dims = 1:50
    n_anchors = 665
    n_unique_query_anchor_cells = 492
    n_unique_reference_anchor_cells = 219
    anchor_rds_saved = TRUE
    static_plot_rendered = FALSE
  DIV90:
    n_query_cells = 22338
    n_reference_cells = 38831
    n_shared_features = 14354
    dims = 1:50
    n_anchors = 1096
    n_unique_query_anchor_cells = 694
    n_unique_reference_anchor_cells = 424
    div90_plot_filter = static_plot_excludes_query_clusters_6_7
    anchor_rds_saved = TRUE
    static_plot_rendered = FALSE

Completed true-anchor plot-only diagnostics:
  Source:
    saved full anchor-pair tables and saved query/reference coordinate tables from
    cross_study_shi_seurat_anchor_projection_v1
  Output:
    true_anchor_projection_plots/figures
  DIV30:
    n_query_cells_plotted_no_downsampling = 90631
    n_reference_cells_plotted_no_downsampling = 38831
    n_anchor_links_after_link_flag_filter = 665
    n_query_cells_with_wta_major_label = 90631
    n_query_cells_with_wta_gw_label = 90631
    n_reference_ids_exact_match_before_suffix_fix = 0
    n_reference_ids_match_after_suffix_fix = 665
    n_query_ids_match_coordinates = 665
    n_links_after_coordinate_merge = 665
    n_anchor_links_removed_by_coordinate_visual_filter = 0
    n_unique_query_anchor_cells_after_merge = 492
    n_unique_reference_anchor_cells_after_merge = 219
  DIV90:
    n_query_cells_plotted_no_downsampling = 20049
    n_reference_cells_plotted_no_downsampling = 38831
    n_anchor_links_after_link_flag_filter = 1096
    n_query_cells_with_wta_major_label = 20049
    n_query_cells_with_wta_gw_label = 20049
    n_reference_ids_exact_match_before_suffix_fix = 0
    n_reference_ids_match_after_suffix_fix = 1096
    n_query_ids_match_coordinates = 939
    n_links_after_coordinate_merge = 939
    n_anchor_links_removed_by_coordinate_visual_filter = 157
    n_unique_query_anchor_cells_after_merge = 576
    n_unique_reference_anchor_cells_after_merge = 353
  Figures written for both PNG/PDF/SVG:
    div30_div90_shi_true_anchor_projection_major_class_side_by_side_grid
    div30_div90_shi_true_anchor_projection_gw_stage_side_by_side_grid
    varela_div30_shi_true_anchor_projection_major_class_side_by_side
    varela_div30_shi_true_anchor_projection_gw_stage_side_by_side
    varela_div90_shi_true_anchor_projection_major_class_side_by_side
    varela_div90_shi_true_anchor_projection_gw_stage_side_by_side
  SVG validation:
    Text remains editable as <text>.
    SVG font styles contain Arial and no DejaVu fallback text styles.
  Visual logic:
    Each panel places query UMAP next to corresponding Shi et al. reference UMAP.
    Lines connect all saved true Seurat anchor links that survive plotting
    coordinate/visual filters.
    A lower subpanel under each projection panel shows:
      query winner-take-all labels with query cluster numbers overlaid
      Shi major class or canonical Shi GW/stage labels with printed label names
    "Filter" in the true-anchor plot diagnostics means visualization-only
    coordinate/cluster exclusion, not winner-take-all assignment. Winner-take-all
    labels come from the saved Seurat TransferData prediction columns.
    DIV90 query uses visualization-only stressed-cell exclusion and vertical
    flip, matching finalized Shi UMAP logic.

Current next step:
  Review the rendered PNG/PDF/SVG in true_anchor_projection_plots/figures. If
  aesthetics need more tuning, rerun only the plot-only script/template; do not
  rerun FindTransferAnchors unless the anchor assets themselves are questioned.
  For methods writing, use the package README:
    /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/final_figures/fig_cross_study_shi_true_anchor_projection_v1_candidate/README.md
  It documents:
    original TransferData WTA label derivation
    true Seurat anchor asset generation
    plot-only rendering inputs and exact code line references
    visualization rules and filter/winner-take-all distinction

Git/worktree note at pause:
  Modified/tracked:
    python_notebooks/HANDOFF_transition_to_final_figs.md
    scripts/14_save_shi_query_anchor_projection_assets.R
  Untracked/new:
    slurm_templates/49d_save_shi_query_anchor_projection_assets.sbatch.template
```

## Figure: fig_div90_jia_urd_marker_pseudotime_tree_v1_candidate

### Status

```text
Found: yes
Confirmed: yes; user selected these three current DIV90 URD assets
Log-audited: pending
Modified: no
Validated: package-start validation only; source files exist and are nonzero
Finalized: no; candidate folder created for formatting/replotting
```

### Candidate Paths

Source run root:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_allcells_jia_root10_neuron_s9_7tips_urd_resumable_v1
```

User-facing mounted-path alias:

```text
/Volumes/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_allcells_jia_root10_neuron_s9_7tips_urd_resumable_v1
```

Candidate plots:

```text
jia_fig_s11_style_marker_validation_v1/plots/jia_fig_s11_style_urd_marker_validation.png
lineage_decision_report/plots/umap_pseudotime.png
lineage_tree_cluster_number_name_v1/plots/urd_tree_pseudotime.png
```

### Source Code and Handoffs

Relevant source scripts/templates copied into the candidate package:

```text
scripts/15_div30_urd_lineage_decision_report.R
scripts/17_div30_urd_finalize_lineage_tree_report.R
scripts/25_div30_urd_jia_fig_s11_marker_validation.R
slurm_templates/34_div90_jia_lineage_urd_smoke.sbatch.template
slurm_templates/36_refresh_urd_posthoc_plots.sbatch.template
```

Relevant source handoff:

```text
python_notebooks/HANDOFF_div90_jia_lineage_urd_plan.md
```

### Prepared Assets

Initial final-figure candidate folder:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/final_figures/fig_div90_jia_urd_marker_pseudotime_tree_v1_candidate
```

Frozen current assets copied into the package:

```text
figures/png/current_jia_fig_s11_style_urd_marker_validation.png
figures/png/current_umap_pseudotime.png
figures/png/current_urd_tree_pseudotime_cluster_number_name.png
figures/pdf/current_jia_fig_s11_style_urd_marker_validation.pdf
```

Copied plotting/provenance tables:

```text
tables/jia_fig_s11_marker_order.tsv
tables/jia_fig_s11_marker_expression_summary.tsv
tables/pseudotime_ordering_by_annotation.tsv
tables/root_annotation_composition.tsv
tables/cluster_number_name_tree_status.tsv
tables/cluster_number_name_tree_tip_composition.tsv
```

Package provenance:

```text
README.md
provenance/source_paths.tsv
provenance/file_manifest.txt
provenance/sha256_manifest.txt
provenance/git_commit.txt
provenance/git_status_short.txt
provenance/HANDOFF_div90_jia_lineage_urd_plan.md
provenance/HANDOFF_transition_to_final_figs_at_package_start.md
logs/render_status.txt
```

### Log Audit

```text
Pending. The package currently records copied source assets only.
No new Slurm or local render was launched for this package on 2026-06-26.
```

### Rerun Decision

```text
Do not recompute upstream URD analysis for formatting.

For this figure set, do not recompute:
  - URD pseudotime
  - lineage tree construction
  - marker-validation expression summaries

Formatting should proceed by plot-only re-render from the existing production
run assets and copied plotting tables. If script edits are needed, make them in
a dedicated final/reformat script and record the command/logs in this package.
```

### Modification Log

```text
2026-06-26:
  - User selected the current marker-validation, UMAP pseudotime, and
    cluster-number-name tree pseudotime assets for final-figure reformatting.
  - Confirmed the /nfs/turbo source files exist and are nonzero.
  - Created the candidate final-figure folder.
  - Copied current PNG/PDF assets, key tables, source scripts/templates, and
    provenance into the package.
  - Recorded explicitly that formatting should begin from current assets and
    should not rerun upstream URD computation.
```

### Validation

```text
Package-start validation:
  current_jia_fig_s11_style_urd_marker_validation.png: copied, nonzero
  current_umap_pseudotime.png: copied, nonzero
  current_urd_tree_pseudotime_cluster_number_name.png: copied, nonzero
  current_jia_fig_s11_style_urd_marker_validation.pdf: copied, nonzero

Visual formatting review pending user direction.
```

### Final Figure Package

```text
Started, not finalized:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/final_figures/fig_div90_jia_urd_marker_pseudotime_tree_v1_candidate
```

### Open Questions

```text
Awaiting user direction on formatting changes: dimensions, labels, panel
composition, color scales, fonts, export formats, and whether to keep the three
plots separate or assemble them into a multi-panel figure.
```

## Pending Reference Figure: Siletti All-Supercluster Transfer

Status: generated but rejected/diagnostic; do not package as final.

Primary handoff:

```text
python_notebooks/HANDOFF_siletti_2023_whb_reference_metadata.md
```

Important correction:

```text
The completed/available `mge_cge_llc_splatter` Siletti runs are not all Siletti
WHB superclusters. They include only four staged H5AD superclusters:
MGE interneuron, CGE interneuron, LAMP5-LHX6 and Chandelier, and Splatter.

The true all-supercluster run must use all 31 Siletti/CELLxGENE WHB
supercluster H5ADs and transfer a label that exists for every adult reference
cell. The first all-supercluster pass uses `source_supercluster`, not
`candidate_jia_group`.
```

Current true all-supercluster job chain:

```text
52396197  siletti-fetch-all    COMPLETED  00:20:36  max RSS 37247344K  gl3343
52396203  siletti-all-bridge   COMPLETED  00:04:32  max RSS 17670324K  gl3253
52396227  siletti-all-knn      COMPLETED  00:00:55  max RSS 7514800K   gl3470
52396231  siletti-all-plot     COMPLETED  00:03:09  max RSS 4406628K   gl3470
```

Download status:

```text
31 / 31 complete Siletti WHB supercluster H5ADs are staged.
```

Completed all-supercluster plot output:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_all_supercluster_source_supercluster_plots_v1/plots
```

Completed plot files:

```text
figure_B_all_supercluster_reference_div90_overlay.png/pdf
figure_B_div90_predicted_siletti_superclusters.png/pdf
figure_C_div90_class_to_siletti_supercluster_river.png/pdf
figure_D_sample_predicted_siletti_supercluster_proportions.png/pdf
```

Completed plotting tables/assets:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_all_supercluster_source_supercluster_plots_v1/tables/div90_query_cells_with_all_siletti_supercluster_assignments.tsv.gz
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_all_supercluster_source_supercluster_plots_v1/tables/siletti_reference_cells_with_all_supercluster_umap.tsv.gz
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_all_supercluster_source_supercluster_plots_v1/tables/figure_C_div90_class_to_siletti_supercluster_edges.tsv
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_all_supercluster_source_supercluster_plots_v1/tables/figure_D_sample_predicted_siletti_supercluster_proportions.tsv
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_all_supercluster_source_supercluster_plots_v1/tables/plot_config.json
```

Validation note:

```text
The four PNGs and four PDFs are nonzero and were written at 2026-06-26
12:53 EDT. PNG dimensions are:
  figure_B_all_supercluster_reference_div90_overlay.png: 3185 x 2096
  figure_B_div90_predicted_siletti_superclusters.png: 2894 x 2096
  figure_C_div90_class_to_siletti_supercluster_river.png: 3571 x 1817
  figure_D_sample_predicted_siletti_supercluster_proportions.png: 2635 x 1753
```

Post-generation issue:

```text
Do not use these candidate outputs for final figures.

The first plot output showed `Unassigned` because of a plotter merge bug:
predictions were merged using the barcode-like `cell_id` column from the
query-obs prediction table instead of `seurat_cell_id`. The merge bug was fixed
in:
  python_notebooks/scripts/plot_siletti_div90_all_supercluster_figure.py

A plot-only rerun was submitted as:
  52398318  siletti-all-plot-fix

However, the underlying kNN assignment is itself degenerate, not merely a plot
bug. The completed kNN table assigns all 16,206 DIV90 cells to:
  Upper-layer intratelencephalic

By DIV90 broad class:
  MGE Striatal/GP Fated: 3,601 -> Upper-layer intratelencephalic
  SST+, NPY +, Cortical Fated: 3,548 -> Upper-layer intratelencephalic
  CRABP1+/PV Precursors: 3,503 -> Upper-layer intratelencephalic
  PV precursors/Migrating cells/Cortical-fated: 2,283 -> Upper-layer intratelencephalic
  LHX8+ vMGE GABergic Striatal/GP fated 1: 1,924 -> Upper-layer intratelencephalic
  LHX8+ vMGE GABergic Striatal/GP fated 2: 922 -> Upper-layer intratelencephalic
  PV Precursors: 425 -> Upper-layer intratelencephalic

This should be treated as a failed all-supercluster diagnostic, not a final
biological result.
```

Root cause and corrected rerun:

```text
Root cause:
  The bridge exporter reused the first H5AD's gene-column index for every H5AD
  in a multi-H5AD scope. The 31 Siletti/CELLxGENE supercluster H5ADs have the
  same gene set but different `var["Gene"]` order. In the all-supercluster run,
  the first H5AD was Upper-layer intratelencephalic, so other supercluster
  reference blocks were gene-column scrambled and the kNN collapsed to the
  first supercluster.

Fix:
  python_notebooks/scripts/export_siletti_div90_seurat_bridge.py now builds
  per-H5AD gene indices and intersects unique genes across all selected
  reference H5ADs and DIV90 query.

Corrected v2 chain:
  52398418_0  siletti-all-bridge-v2  COMPLETED  00:05:35  max RSS 22772156K  gl3027
  52398421    siletti-all-knn-v2     COMPLETED  00:01:47  max RSS 14302360K  gl3478
  52398423    siletti-all-plot-v2    COMPLETED  00:03:41  max RSS 12706328K  gl3027

Corrected v2 plot output:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_all_supercluster_source_supercluster_plots_v2/plots

Canceled obsolete v1 plot-only fix:
  52398318
```

Corrected v2 assignment summary:

```text
The corrected v2 run no longer collapses to Upper-layer intratelencephalic.
Most cells assign to Siletti `Splatter`, with small MGE/CGE/Miscellaneous calls.
This is a corrected broad all-supercluster diagnostic, but it should still be
reviewed biologically before final packaging.

Predicted supercluster counts:
  Splatter: 16,122
  Miscellaneous: 31
  MGE interneuron: 26
  CGE interneuron: 24
  Committed oligodendrocyte precursor: 2
  Fibroblast: 1

Corrected v2 files:
  figure_B_all_supercluster_reference_div90_overlay.png/pdf
  figure_B_div90_predicted_siletti_superclusters.png/pdf
  figure_C_div90_class_to_siletti_supercluster_river.png/pdf
  figure_D_sample_predicted_siletti_supercluster_proportions.png/pdf
```

Scale-up audit:

```text
Detailed audit file:
  python_notebooks/HANDOFF_siletti_scaleup_transfer_audit_2026_06_26.md

Do not finalize the all-supercluster source-supercluster figure until the audit
question is resolved. The key issue is that the old MGE/Jia-style transfer used
a curated, restricted label space, while the scaled v2 run uses raw broad
`source_supercluster` labels and is dominated by Splatter.

Audit job:
  52413882  siletti-all-jia-audit-v2  COMPLETED  00:01:01

Audit conclusion:
  The corrected all-supercluster v2 bridge can produce MGE/Jia-style labels
  when the original curated `candidate_jia_group` label and exclusions are used.
  Therefore, the Splatter-dominant all-supercluster figure is not a remaining
  matrix-corruption problem; it is a label-space/reference-composition problem.
```

## Figure: fig_siletti_div90_restricted_mge_llc_msn_emsn_chat_rivers_v6_no_cge_min10_ordered

### Status

```text
Generated v6 final-candidate river plots from one unified leaf-label transfer.
CGE is excluded completely from the adult reference before bridge export and
before fast-kNN transfer; it is not merely hidden from the plots.
Pallial/subpallial and major labels are derived from that same prediction, so
the three river levels use one nested assignment logic. Vertical Pallial/cortical
vs Subpallial grouping bars, a combined 1x3 panel, and visual accounting audits
are included. v6 filters plotted river edges with fewer than 10 DIV90 cells and
uses one fixed left-side DIV90 class order across all river plots. A
final_figures candidate package has been created, but this is not
final-approved yet; visual review and final wording remain pending.
```

### Candidate Paths

```text
Final-candidate output:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_FINAL_FIGURE_CANDIDATE_restricted_mge_llc_msn_emsn_chat_river_plots_v6_no_cge_min10_ordered

Final-figures candidate package:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/final_figures/fig_siletti_div90_restricted_mge_llc_msn_emsn_chat_rivers_v6_no_cge_min10_ordered_candidate

Expected plots:
plots/river_div90_class_to_adult_pallial_subpallial_bin.png/pdf
plots/river_div90_class_to_adult_major_interneuron_subtypes.png/pdf
plots/river_div90_class_to_adult_final_fine_subtypes.png/pdf
plots/river_div90_class_to_adult_combined_1x3_pallial_major_fine.png/pdf

Audit plots in same folder:
plots/audit_reference_supercluster_by_pallial_subpallial_bin.png/pdf
plots/audit_reference_major_subtype_by_pallial_subpallial_bin.png/pdf
plots/audit_reference_major_subtype_totals_including_excluded.png/pdf
plots/audit_reference_final_fine_subtype_by_pallial_subpallial_bin.png/pdf
plots/audit_reference_final_fine_subtype_totals_including_excluded.png/pdf
plots/audit_div90_class_to_pallial_subpallial_counts.png/pdf
plots/audit_div90_class_to_major_subtype_counts.png/pdf
plots/audit_div90_class_to_final_fine_subtype_counts.png/pdf
plots/audit_transfer_score_distributions.png/pdf
tables/audit_reference_labels_with_zero_query_assignments.tsv
```

### Source Code and Handoffs

```text
python_notebooks/scripts/export_siletti_div90_seurat_bridge.py
python_notebooks/scripts/siletti_div90_fast_knn_label_transfer.py
python_notebooks/scripts/plot_siletti_div90_final_restricted_rivers.py
python_notebooks/scripts/plot_siletti_div90_jia_figure.py
python_notebooks/HANDOFF_siletti_2023_whb_reference_metadata.md
```

### Rerun Decision

```text
Rerun submitted with corrected per-H5AD gene indexing and restricted final
candidate reference.

Reference scope:
  MGE interneuron
  LAMP5-LHX6 and Chandelier
  Medium spiny neuron
  Eccentric medium spiny neuron
  Splatter restricted to subpallial cholinergic/CHAT-like rows only

Pallial/subpallial rule:
  CerebralCortex + Hippocampus -> Pallial/cortical
  all other primary ROI groups -> Subpallial
```

### Modification Log

```text
Bridge:
  52439510  siletti-final-restricted-bridge  COMPLETED  00:02:29  max RSS 17571972K  gl3018
  run_label: siletti_div90_final_candidate_restricted_ge_msn_chat_bridge_v1
  scope: mge_cge_llc_msn_emsn_cholinergic

Pallial/subpallial transfer:
  52439512  siletti-final-pallsub-transfer  COMPLETED  00:01:37  max RSS 13764068K  gl3018
  label_column: roi_pallial_subpallial_bin

Major subtype transfer:
  52439514  siletti-final-subtype-transfer  COMPLETED  00:01:34  max RSS 12580828K  gl3047
  label_column: major_interneuron_subtype
  excluded labels: Other selected reference

River plots:
  52439515  siletti-final-rivers  COMPLETED  00:00:18  max RSS 650588K  gl3300

Audit-layer rerender:
  plot_siletti_div90_final_restricted_rivers.py rerun after completion to add
  explicit visual/accounting audit outputs to the same final-candidate folder.
```

### Validation

```text
Generated plot files:
  plots/river_div90_class_to_adult_pallial_subpallial_bin.png/pdf
  plots/river_div90_class_to_adult_major_interneuron_subtypes.png/pdf
  plots/audit_reference_supercluster_by_pallial_subpallial_bin.png/pdf
  plots/audit_reference_major_subtype_by_pallial_subpallial_bin.png/pdf
  plots/audit_reference_major_subtype_totals_including_excluded.png/pdf
  plots/audit_div90_class_to_pallial_subpallial_counts.png/pdf
  plots/audit_div90_class_to_major_subtype_counts.png/pdf
  plots/audit_transfer_score_distributions.png/pdf

PNG dimensions:
  river_div90_class_to_adult_pallial_subpallial_bin.png: 3607 x 1817
  river_div90_class_to_adult_major_interneuron_subtypes.png: 3905 x 1817
  audit_div90_class_to_major_subtype_counts.png: 4085 x 1403
  audit_div90_class_to_pallial_subpallial_counts.png: 3787 x 1403
  audit_reference_major_subtype_by_pallial_subpallial_bin.png: 3593 x 1403
  audit_reference_major_subtype_totals_including_excluded.png: 2942 x 1172
  audit_reference_supercluster_by_pallial_subpallial_bin.png: 3593 x 1403
  audit_transfer_score_distributions.png: 3041 x 1172

Saved tables:
  tables/river_div90_class_to_adult_pallial_subpallial_bin_edges.tsv
  tables/river_div90_class_to_adult_major_interneuron_subtypes_edges.tsv
  tables/div90_query_with_pallial_subpallial_assignments.tsv.gz
  tables/div90_query_with_major_interneuron_subtype_assignments.tsv.gz
  tables/audit_accounting_totals.tsv
  tables/audit_reference_source_supercluster_by_pallial_subpallial_bin.tsv
  tables/audit_reference_major_subtype_by_pallial_subpallial_bin.tsv
  tables/audit_reference_major_subtype_totals.tsv
  tables/audit_div90_class_to_pallial_subpallial_counts.tsv
  tables/audit_div90_class_to_major_subtype_counts.tsv

Audit accounting totals:
  reference_cells_total: 60,000
  reference_cells_used_for_subtype_transfer: 51,195
  reference_cells_excluded_from_subtype_transfer_other_selected_reference: 8,805
  div90_cells_pallial_transfer: 16,206
  div90_cells_major_subtype_transfer: 16,206
  river_pallial_edges_total: 16,206
  river_major_subtype_edges_total: 16,206
```

Assignment counts:

```text
Pallial/subpallial:
  Subpallial: 11,442
  Pallial/cortical: 4,764

Major subtype:
  Medium spiny neuron: 4,549
  Eccentric medium spiny neuron: 4,111
  Vip: 3,527
  Pvalb: 1,819
  Lamp5 Lhx6: 1,150
  Sst: 742
  Subpallial Cholinergic neurons: 306
  Pax6: 2
```

### V2 Aligned River Update

```text
Status:
  Generated and validated. This is the current preferred candidate for visual
  review because it fixes river rectangle/ribbon alignment and adds the
  fine-subtype river.

Output:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_FINAL_FIGURE_CANDIDATE_restricted_ge_msn_chat_river_plots_v2

Main plots:
  plots/river_div90_class_to_adult_pallial_subpallial_bin.png/pdf
  plots/river_div90_class_to_adult_major_interneuron_subtypes.png/pdf
  plots/river_div90_class_to_adult_final_fine_subtypes.png/pdf

Geometry/rendering changes:
  Rectangles and ribbons now use the same scaled vertical segment positions.
  Right-side labels are ordered with pallial/cortical groups above subpallial
  groups where the label space carries that distinction.
  Small labels are spread with leader lines so tiny bins remain readable without
  moving the rectangles themselves.

Fine subtype label column:
  final_fine_subtype

Fine subtype options:
  Cortical PV+ basket neurons
  Cortical PV+ Chandelier neurons
  Cortical SST+ LRP neurons
  Cortical SST+ Mt neurons
  Cortical SST+ nMt neurons
  Pallial/cortical Medium spiny neuron
  Subpallial PV+ neurons
  Subpallial SST+ LRP neurons
  Subpallial SST+ neurons
  Subpallial Cholinergic neurons
  Subpallial Medium spiny neuron
  Subpallial Eccentric medium spiny neuron

Other selected reference is retained for accounting and excluded from the
fine-subtype transfer.

Jobs:
  52439770  siletti-final-restricted-bridge-v2   COMPLETED  00:02:27  max RSS 17589848K  gl3435
  52439771  siletti-final-pallsub-transfer-v2    COMPLETED  00:01:50  max RSS 13849648K  gl3009
  52439772  siletti-final-major-transfer-v2      COMPLETED  00:01:29  max RSS 12468616K  gl3018
  52439773  siletti-final-fine-transfer-v2       COMPLETED  00:01:13  max RSS 11204292K  gl3408
  52439774  siletti-final-rivers-v2              COMPLETED  00:00:26  max RSS 953068K    gl3378
```

V2 accounting totals:

```text
reference_cells_total: 60,000
reference_cells_used_for_subtype_transfer: 51,195
reference_cells_excluded_from_subtype_transfer_other_selected_reference: 8,805
reference_cells_used_for_final_fine_subtype_transfer: 37,417
reference_cells_excluded_from_final_fine_subtype_transfer_other_selected_reference: 22,583
div90_cells_pallial_transfer: 16,206
div90_cells_major_subtype_transfer: 16,206
div90_cells_final_fine_subtype_transfer: 16,206
river_pallial_edges_total: 16,206
river_major_subtype_edges_total: 16,206
river_final_fine_subtype_edges_total: 16,206
```

V2 fine subtype assignments:

```text
Subpallial Eccentric medium spiny neuron: 6,745
Subpallial Medium spiny neuron: 4,301
Cortical PV+ basket neurons: 3,550
Subpallial SST+ neurons: 763
Subpallial Cholinergic neurons: 394
Cortical SST+ Mt neurons: 286
Cortical SST+ LRP neurons: 145
Cortical SST+ nMt neurons: 22

Cortical PV+ Chandelier neurons are present in the reference
(613 reference cells) but received 0 DIV90 assignments in this transfer.
```

### V3 Grouped River Update

```text
Status:
  Generated and validated. This is the current preferred candidate for visual
  review.

Output:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_FINAL_FIGURE_CANDIDATE_restricted_ge_msn_chat_river_plots_v3

Main plots:
  plots/river_div90_class_to_adult_pallial_subpallial_bin.png/pdf
  plots/river_div90_class_to_adult_major_interneuron_subtypes.png/pdf
  plots/river_div90_class_to_adult_final_fine_subtypes.png/pdf
  plots/river_div90_class_to_adult_combined_1x3_pallial_major_fine.png/pdf

V3 changes:
  Added vertical Pallial/cortical and Subpallial grouping bars on all river
  plots.
  Added split major label column: major_interneuron_subtype_roi.
  Major labels are now ordered as Pallial/cortical labels first, then
  Subpallial labels.
  Zero-count ordered target labels remain visible so reference options such as
  PV/chandelier are not mistaken for absent reference classes.

Jobs:
  52440165  siletti-final-restricted-bridge-v3    COMPLETED  00:02:40  max RSS 17645408K  gl3028
  52440166  siletti-final-pallsub-transfer-v3     COMPLETED  00:01:52  max RSS 13868608K  gl3009
  52440167  siletti-final-majorroi-transfer-v3    COMPLETED  00:01:26  max RSS 12586348K  gl3408
  52440168  siletti-final-fine-transfer-v3        COMPLETED  00:01:39  max RSS 11204324K  gl3018
  52440169  siletti-final-rivers-v3               COMPLETED  00:00:29  max RSS 1096588K   gl3329
```

V3 accounting totals:

```text
reference_cells_total: 60,000
reference_cells_used_for_subtype_transfer: 51,195
reference_cells_excluded_from_subtype_transfer_other_selected_reference: 8,805
reference_cells_used_for_final_fine_subtype_transfer: 37,417
reference_cells_excluded_from_final_fine_subtype_transfer_other_selected_reference: 22,583
div90_cells_pallial_transfer: 16,206
div90_cells_major_subtype_transfer: 16,206
div90_cells_final_fine_subtype_transfer: 16,206
river_pallial_edges_total: 16,206
river_major_subtype_edges_total: 16,206
river_final_fine_subtype_edges_total: 16,206
```

V3 split-major assignments:

```text
Subpallial Eccentric medium spiny neuron: 4,614
Subpallial Medium spiny neuron: 4,587
Pallial/cortical Vip: 2,737
Pallial/cortical Pvalb: 1,840
Pallial/cortical Lamp5 Lhx6: 1,135
Pallial/cortical Sst: 887
Subpallial Cholinergic neurons: 337
Subpallial Vip: 65
Pallial/cortical Pax6: 4
```

PV chandelier status:

```text
Adult reference contains:
  Cortical PV+ Chandelier neurons: 613 reference cells
  Pallial/cortical Chandelier: 613 reference cells
  Subpallial Chandelier: 87 reference cells

DIV90 assignments to those chandelier labels in v3 are 0. They are shown as
zero-count target labels in the river plots, not as nonzero flows.
```

### V4 Unified Transfer Update

```text
Status:
  Generated and validated. Superseded for final-figure plotting by v5, which
  uses the same unified transfer but adds the <10-cell edge filter and fixed
  left-side ordering.

Output:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_FINAL_FIGURE_CANDIDATE_restricted_ge_msn_chat_river_plots_v4_unified

Final-figures candidate package:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/final_figures/fig_siletti_div90_restricted_ge_msn_chat_rivers_v4_unified_candidate

Main plots:
  plots/river_div90_class_to_adult_pallial_subpallial_bin.png/pdf
  plots/river_div90_class_to_adult_major_interneuron_subtypes.png/pdf
  plots/river_div90_class_to_adult_final_fine_subtypes.png/pdf
  plots/river_div90_class_to_adult_combined_1x3_pallial_major_fine.png/pdf

Logic:
  Transfer is performed once using `unified_leaf_subtype`.
  `unified_major_subtype_roi` and `unified_pallial_subpallial_bin` are derived
  from that same predicted leaf label for plotting.
  Zero-count reference labels are not drawn in the assignment rivers; they are
  tracked in `tables/audit_reference_labels_with_zero_query_assignments.tsv`.

Jobs:
  52440422  siletti-final-bridge-v4-unified    COMPLETED  00:02:23  max RSS 17672732K  gl3014
  52440423  siletti-final-unified-transfer-v4  COMPLETED  00:01:47  max RSS 12666716K  gl3009
  52440424  siletti-final-rivers-v4-unified    COMPLETED  00:00:29  max RSS 1073464K   gl3434
```

V4 accounting totals:

```text
reference_cells_total: 60,000
reference_cells_used_for_subtype_transfer: 51,580
reference_cells_excluded_from_subtype_transfer_other_selected_reference: 8,420
reference_cells_used_for_final_fine_subtype_transfer: 51,580
reference_cells_excluded_from_final_fine_subtype_transfer_other_selected_reference: 8,420
div90_cells_pallial_transfer: 16,206
div90_cells_major_subtype_transfer: 16,206
div90_cells_final_fine_subtype_transfer: 16,206
river_pallial_edges_total: 16,206
river_major_subtype_edges_total: 16,206
river_final_fine_subtype_edges_total: 16,206
```

V4 derived counts:

```text
Pallial/subpallial:
  Subpallial: 10,295
  Pallial/cortical: 5,911

Major:
  Subpallial Medium spiny neuron: 5,547
  Subpallial Eccentric medium spiny neuron: 3,906
  Pallial/cortical Vip: 2,460
  Pallial/cortical Pvalb: 1,954
  Pallial/cortical Lamp5 Lhx6: 1,226
  Subpallial Sst: 394
  Subpallial Cholinergic neurons: 367
  Pallial/cortical Sst: 266
  Subpallial Vip: 81
  Pallial/cortical Pax6: 5

Fine/leaf:
  Subpallial Medium spiny neuron: 5,547
  Subpallial Eccentric medium spiny neuron: 3,906
  Pallial/cortical Vip: 2,460
  Cortical PV+ basket neurons: 1,954
  Pallial/cortical Lamp5 Lhx6: 1,226
  Subpallial SST+ neurons: 394
  Subpallial Cholinergic neurons: 367
  Cortical SST+ Mt neurons: 160
  Cortical SST+ LRP neurons: 86
  Subpallial Vip: 81
  Cortical SST+ nMt neurons: 20
  Pallial/cortical Pax6: 5
```

PV chandelier isolated audit:

```text
Adult reference contains Cortical PV+ Chandelier neurons, but the unified v4
transfer assigns 0 DIV90 cells to that label.

Audit job 52440381 restricted the fine reference choices to:
  Cortical PV+ basket neurons
  Cortical PV+ Chandelier neurons
  Subpallial PV+ neurons

Result:
  Cortical PV+ basket neurons: 16,206
  Cortical PV+ Chandelier neurons: 0
  Subpallial PV+ neurons: 0

Soft-score/non-winner audit:
  tables/audit_chandelier_soft_score_summary.tsv
  tables/audit_chandelier_soft_score_by_div90_class.tsv
  tables/audit_chandelier_soft_score_candidate_cells.tsv

  In the unified v4 transfer, Cortical PV+ Chandelier is never rank 1.
  Only 3 DIV90 cells have any nonzero Cortical PV+ Chandelier score, all in
  CRABP1+/PV Precursors. This is a tiny chandelier-neighbor signal, not a
  robust alternate chandelier assignment.
```

### V5 Min-10 Ordered Plot Update

```text
Status:
  Generated and visually spot-checked. Superseded by v6, which reruns the
  bridge and unified transfer after removing CGE from the reference.

Output:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_FINAL_FIGURE_CANDIDATE_restricted_ge_msn_chat_river_plots_v5_min10_ordered

Final-figures candidate package:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/final_figures/fig_siletti_div90_restricted_ge_msn_chat_rivers_v5_min10_ordered_candidate

Computation:
  Plot-only rerender from the saved v4 bridge and unified transfer. No bridge
  export and no fast-kNN transfer were rerun.

Plotting changes:
  Edges with fewer than 10 DIV90 cells are not drawn.
  River rectangles and ribbons are recomputed from the filtered plotted-edge
  tables, so sub-10-cell populations do not leave tiny visible target boxes or
  spacing artifacts.
  Full unfiltered edge tables remain in `tables/` for accounting.
  Filtered edges are recorded in explicit audit tables.

Left-side DIV90 class order:
  PV Precursors
  CRABP1+/PV Precursors
  SST+, NPY +, Cortical Fated
  PV precursors/Migrating cells/Cortical-fated
  LHX8+ vMGE GABergic Striatal/GP fated 1
  MGE Striatal/GP Fated
  LHX8+ vMGE GABergic Striatal/GP fated 2

Filter audit:
  pallial_subpallial: 0 edges filtered, 16,206 / 16,206 cells plotted
  major_subtype: 14 edges filtered, 16,146 / 16,206 cells plotted
  final_fine_subtype: 20 edges filtered, 16,128 / 16,206 cells plotted

New tables:
  tables/river_plot_edge_filter_audit.tsv
  tables/river_major_subtype_filtered_edges_lt10.tsv
  tables/river_final_fine_subtype_filtered_edges_lt10.tsv
  tables/river_div90_class_to_adult_pallial_subpallial_bin_plotted_edges_min10.tsv
  tables/river_div90_class_to_adult_major_interneuron_subtypes_plotted_edges_min10.tsv
  tables/river_div90_class_to_adult_final_fine_subtypes_plotted_edges_min10.tsv

Validation:
  Main PNG/PDF files are present and nonzero.
  Plotted-edge minimums are >=10:
    pallial/subpallial min plotted edge: 18
    major subtype min plotted edge: 12
    final fine subtype min plotted edge: 11
  Combined 1x3 PNG dimensions: 11794 x 3559.
```

### V6 No-CGE Min-10 Ordered Update

```text
Status:
  Generated, visually spot-checked, and copied into final_figures as the current
  preferred candidate package. Not final-approved yet; visual review and final
  wording remain pending.

Output:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_FINAL_FIGURE_CANDIDATE_restricted_mge_llc_msn_emsn_chat_river_plots_v6_no_cge_min10_ordered

Final-figures candidate package:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/final_figures/fig_siletti_div90_restricted_mge_llc_msn_emsn_chat_rivers_v6_no_cge_min10_ordered_candidate

Computation:
  CGE was removed from the bridge scope before transfer. This is not a
  plot-only hide-CGE update.

Reference scope:
  MGE interneuron
  LAMP5-LHX6 and Chandelier
  Medium spiny neuron
  Eccentric medium spiny neuron
  Splatter restricted to subpallial cholinergic/CHAT-like rows only

Reference after export:
  reference_cells_total: 50,190
  MGE interneuron: 20,887
  Medium spiny neuron: 14,030
  Eccentric medium spiny neuron: 9,236
  LAMP5-LHX6 and Chandelier: 5,403
  Splatter CHAT/cholinergic subset: 634

Jobs:
  52440897  siletti-final-bridge-v6-nocge            COMPLETED  00:02:06  max RSS 15431728K
  52440898  siletti-final-unified-transfer-v6-nocge  COMPLETED  00:01:36  max RSS 11414312K
  52440899  siletti-final-rivers-v6-nocge            COMPLETED  00:00:28  max RSS 1068412K

Filter audit:
  pallial_subpallial: 1 edge filtered, 16,205 / 16,206 cells plotted
  major_subtype: 11 edges filtered, 16,163 / 16,206 cells plotted
  final_fine_subtype: 17 edges filtered, 16,144 / 16,206 cells plotted

Validation:
  Main PNG/PDF files are present and nonzero.
  Combined 1x3 PNG dimensions: 11794 x 3559.
  Visual spot-check confirms CGE-derived target labels are absent.
```

### Final Figure Package

```text
Status:
  Packaged candidate only. This folder is for final-figure development and
  provenance capture; it does not mean the Siletti river figure is final
  approved.

Package:
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/final_figures/fig_siletti_div90_restricted_mge_llc_msn_emsn_chat_rivers_v6_no_cge_min10_ordered_candidate

Contents:
  README.md
  figures/png/river_div90_class_to_adult_pallial_subpallial_bin.png
  figures/png/river_div90_class_to_adult_major_interneuron_subtypes.png
  figures/png/river_div90_class_to_adult_final_fine_subtypes.png
  figures/png/river_div90_class_to_adult_combined_1x3_pallial_major_fine.png
  figures/pdf/ matching PDF files
  figures/png/audit/ and figures/pdf/audit/ visual accounting checks
  tables/ all river edge, assignment, accounting, zero-label, and chandelier
    audit tables from the v4 output
  code/ exact scripts copied for the candidate package
  logs/ Slurm logs for bridge, unified transfer, river rendering, and PV-only
    chandelier audit
  provenance/ git commit, git status, relevant uncommitted diff, source README,
    file manifest, and sha256 manifest

Note:
  No SVG files were generated for this candidate. The package keeps an empty
  figures/svg/ directory to match the final_figures standard layout.
```

### Open Questions

```text
Visual review pending before final approval.

Interpretation caveat to review visually:
  Pallial/subpallial is an adult primary-ROI bin, not a developmental origin
  label. Under this rule, many adult MGE/CGE/LLC interneurons are counted as
  Pallial/cortical because their sampled ROI is CerebralCortex or Hippocampus.
  MSN, eccentric MSN, and Splatter-CHAT/cholinergic rows are mostly Subpallial.
  Use the audit plots/tables above to decide whether this rule is acceptable for
  the final figure wording.
```

Assignment/methods summary:

```text
This run transfers `source_supercluster`, not `candidate_jia_group`.
The reference is the all-31-supercluster Siletti/CELLxGENE WHB set, downsampled
at bridge export to max 100 cells per Siletti subcluster and max 60,000 total
adult reference cells. DIV90 query cells are uncapped.

The assignment is not Seurat anchors. It is the repo fast-kNN path:
shared unique genes -> library-size normalize to 10,000 counts/cell -> log1p ->
top 3,000 variable genes -> TruncatedSVD with 50 components -> cosine kNN with
k = 50 -> normalized cosine-similarity weighted vote -> highest-weight
`source_supercluster` becomes `predicted.id`.
```

Detailed methods/input/output record:

```text
python_notebooks/HANDOFF_siletti_2023_whb_reference_metadata.md
section: All-Supercluster Assignment Workflow Details
```

## 2026-06-27 DIV90 URD Marker/Pseudotime Tree Final-Figure Package Update

Final-figure package:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/final_figures/fig_div90_jia_urd_marker_pseudotime_tree_v1_candidate

User-facing alias:
/Volumes/umms-parent/mgeo_neuron_scrnaseq_projectfolder/final_figures/fig_div90_jia_urd_marker_pseudotime_tree_v1_candidate
```

Source run:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_allcells_jia_root10_neuron_s9_7tips_urd_resumable_v1
```

Plot-only update completed and submitted as Slurm job `52440501`
(`div90_urd_6marker`), which completed successfully with exit code `0:0`
on `gl3411` in 19 seconds. This did not recompute URD pseudotime, random
walks, or the lineage tree; it rerendered marker overlays from the saved
tree object:

```text
lineage_tree_jia_endpoint_tips_v1/div30_urd_lineage_tree_object.rds
```

Current marker panel order:

```text
Hes1 | Nkx2.1 | Lhx6 | Lhx8 | Crabp1 | Kcnc1
```

Matrix features used:

```text
HES1 | NKX2-1 | LHX6 | LHX8 | CRABP1 | KCNC1
```

`Nkx2.1` is the display label for feature `NKX2-1`. The earlier requested
`Stmn5` panel was removed from the current packaged figure because `STMN5`
was not present in the saved expression matrix.

Current packaged figures:

```text
figures/png/current_jia_fig_s11_style_urd_marker_validation.png
figures/pdf/current_jia_fig_s11_style_urd_marker_validation.pdf
figures/png/current_div90_urd_marker_overlays_hes1_nkx21_lhx6_lhx8_crabp1_kcnc1.png
figures/pdf/current_div90_urd_marker_overlays_hes1_nkx21_lhx6_lhx8_crabp1_kcnc1.pdf
figures/png/current_urd_tree_cluster_number_name.png
figures/png/current_urd_tree_cluster_number_name_grid.png
figures/png/current_urd_tree_pseudotime.png
figures/png/current_urd_tree_pseudotime_cluster_number_name.png
figures/png/current_umap_pseudotime.png
figures/png/current_pseudotime_overlay_umap.png
```

Package metadata/provenance refreshed:

```text
README.md
logs/render_status.txt
tables/jia_fig_s11_marker_order.tsv
tables/jia_fig_s11_marker_expression_summary.tsv
provenance/source_paths.tsv
provenance/file_manifest.txt
provenance/sha256_manifest.txt
code/25_div30_urd_jia_fig_s11_marker_validation.R
```

The marker renderer now supports optional `--genes`, `--gene-labels`, and
`--panel-title` arguments while preserving the original Jia Fig. S11 default
panel when those arguments are omitted.

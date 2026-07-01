# Handoff: Cross-study MGE/PV-marker synthesis concept

Date started: 2026-06-29

Purpose: capture the methods and conceptual frame for a new synthesis figure
that summarizes the repeated single-gene cross-study marker panels before
building the final integrated plot.

## Concept

The single-gene panels were created to answer a focused biological question:

```text
Across studies and samples, how consistently are cells with an MGE/PV-precursor
marker profile being made?
```

The planned synthesis figure should roll up the per-gene panels into a compact
study-by-marker and sample-by-marker view. It should preserve the key message
from the UMAP/violin/sample-positive panels without requiring one full panel
per gene in the final figure.

Primary interpretation axis:

```text
Early/identity markers:
  NKX2.1/NKX2-1
  LHX6
  LHX8

Maturation / interneuron-lineage-associated markers:
  MAFB
  MEF2C
  CRABP1
  TAC1
  VIPR2

Late PV marker:
  PV/PVALB
```

Working biological readout:

```text
High NKX2.1 plus LHX6/LHX8 supports an MGE-like precursor/interneuron identity.
Low PVALB across studies indicates that canonical PV maturation is scarce under
this threshold, even where upstream MGE markers are present.
```

## Source Data

Use the prepared compact marker-expression tables from the PV precursor final
candidate run. `VIPR2` was not present in the original candidate compact
tables, so a copied/enriched run was created by appending a freshly exported
`VIPR2` column:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/cross_study_marker_expression/cross_study_marker_expression_v12_pv_precursors_final_candidate_plus_vipr2/tables/per_study/*_marker_expression.tsv.gz
```

Do not rerun upstream extraction, UMAP, clustering, label transfer, or Seurat
work for the first synthesis pass. The needed per-cell marker tables already
exist.

## Generated Single-gene Package

The synthesis figure should treat this consolidated package as the immediate
plotting input:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/final_figures/fig_cross_study_single_gene_expression_v12_candidate
```

The previous separate per-gene final-figure folders were deleted on 2026-06-29
and replaced by this fresh consolidated package.

The package contains:

```text
figures/png|pdf|svg/cross_study_<gene>_expression_umap_violin_sample_positive.*
tables/cross_study_<gene>_expression_study_summary.tsv
tables/cross_study_<gene>_expression_sample_positive_fraction.tsv
tables/cross_study_<gene>_expression_plot_filter_summary.tsv
tables/cross_study_<gene>_expression_output_manifest.tsv
tables/cross_study_single_gene_expression_combined_study_summary.tsv
tables/cross_study_single_gene_expression_combined_output_manifest.tsv
provenance/batch_gene_manifest.tsv
provenance/batch_render_manifest.tsv
```

## Shared Methods

The repeated panels use one method contract:

```text
Expression values:
  log1p(CP10K), carried from the prepared marker-expression tables.

Positive-cell threshold:
  expression >= 0.5 log expression.

Study denominator:
  visualization-filtered cells with finite marker expression.

Sample denominator:
  visualization-filtered cells with finite marker expression within each sample.

Visualization filters:
  Samarasinghe et al. 2021: controls only.
  This Study, DIV 90: remove current clusters 6 and 7 as stressed cells.
  Other studies: no additional single-gene synthesis filter.

DIV90 UMAP orientation:
  display-only UMAP2 flip is used for UMAP panels. Original coordinates are
  not overwritten.

Aliases:
  NKX2.1 maps to marker table column NKX2-1.
  PV maps to marker table column PVALB.

Renderer:
  python_notebooks/scripts/render_cross_study_lhx6_expression_final_panel.py
  python_notebooks/scripts/render_cross_study_single_gene_expression_batch.py
```

## Study-level Percent-positive Matrix

All values are percent of cells with marker expression >= 0.5 after the shared
visualization filters.

| Study | NKX2.1 | LHX6 | LHX8 | MAFB | MEF2C | CRABP1 | TAC1 | VIPR2 | PV/PVALB |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| This Study, DIV 30 | 94.96 | 42.20 | 58.23 | 20.57 | 12.35 | 40.84 | 13.57 | 1.94 | 0.21 |
| This Study, DIV 90 | 61.58 | 35.73 | 17.72 | 12.46 | 14.96 | 27.24 | 7.81 | 3.03 | 0.43 |
| Siebert et al. 2026 | 19.98 | 7.07 | 0.53 | 5.96 | 27.22 | 10.77 | 1.62 | 0.67 | 0.13 |
| Walsh et al. 2025 | 43.64 | 36.22 | 4.05 | 29.96 | 54.30 | 12.95 | 8.23 | 0.73 | 0.44 |
| Bershteyn et al. 2025 | 17.34 | 39.49 | 5.65 | 52.78 | 61.22 | 0.24 | 5.11 | 0.22 | 0.001 |
| Bershteyn et al. 2023 | 10.99 | 36.42 | 14.71 | 28.58 | 55.32 | 12.68 | 8.56 | 0.86 | 0.10 |
| Samarasinghe et al. 2021 controls | 0.31 | 0.05 | 0.04 | 4.70 | 22.22 | 14.26 | 1.34 | 0.82 | 0.73 |

## Recommended Synthesis Figure

Start with a plot-only figure package:

```text
final_figures/fig_cross_study_mge_pv_marker_synthesis_v1_candidate
```

Suggested layout:

```text
Panel A:
  Study x marker heatmap of percent-positive cells.
  Rows = studies.
  Columns = NKX2.1, LHX6, LHX8, MAFB, MEF2C, CRABP1, TAC1, VIPR2, PV/PVALB.
  Cell labels can show rounded percentages.

Panel B:
  Sample x marker heatmap or dot plot.
  Rows grouped by study and ordered with the same sample ordering used in the
  single-gene panels.
  Dot color = percent-positive cells.
  Dot size = sample cell count or percent-positive, whichever reads better.

Panel C:
  Marker-state summary for each study, such as:
    MGE precursor signal = mean or geometric-style average of NKX2.1/LHX6/LHX8
    PV maturation signal = PV/PVALB percent-positive
  Use this as an interpretive aid, not as a replacement for the direct marker
  matrix.
```

Conservative first version:

```text
Build Panel A and Panel B only.
Avoid inventing a composite score until the direct marker patterns are reviewed.
```

## Rendered Publication-style Follow-up

The user next requested a more specific integrated publication-style figure,
which superseded rendering the generic synthesis concept above as a standalone
package.

Rendered package:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/final_figures/fig_cross_study_integrated_mge_marker_expression_v1_candidate
```

Renderer:

```text
python_notebooks/scripts/render_cross_study_integrated_marker_publication_figures.py
```

Figure set:

```text
1. Cross-study canonical MGE marker comparison:
   NKX2.1, LHX6, LHX8, CRABP1.
   Paired gene-block layout: each gene has a compact sample-level expression
   violin row immediately followed by the matched sample-level percent-positive
   bar row.
   Sample names are shown once at the bottom only.
   Gene labels are neutral black; color encodes study/group only.

2. This Study DIV30 vs DIV90 replicate reproducibility:
   NKX2.1, LHX6, LHX8.
   Unfilled DIV30/DIV90 boxplots with replicate/sample dots.
   Three variants are rendered:
     no per-cell cutoff
     per-cell inclusion cutoff expression >= 0.5 log1p(CP10K)
     per-cell inclusion cutoff expression >= 1.0 log1p(CP10K)
   For each replicate/sample and gene, the cutoff is applied at the cell level
   first; the dot is then the mean expression across included cells.
   Gene titles are neutral black; DIV30/DIV90 colors encode group only.
```

Important scale note:

```text
All values remain the existing log1p(CP10K) values from the prepared marker
tables. No log10 conversion or rescaling is applied. The dashed reference line
and positive-cell threshold are expression >= 0.5 log1p(CP10K).
```

## Open Design Decisions

Before final synthesis rendering, decide:

```text
1. Should the marker threshold remain 0.5 for every gene, or should the final
   figure also show sensitivity at 1.0?
2. Should the synthesis figure include violin/UMAP thumbnails, or only the
   summary matrix?
3. Should PVALB be labeled as PV in the figure while preserving PVALB in table
   column names?
4. Should the study-level heatmap use one color scale for all markers, or split
   early identity markers and late PV marker into separate visual bands?
5. Should the sample-level panel include all samples, or collapse large
   Bershteyn sample sets into study-level summaries for readability?
```

## Guardrails

```text
- Keep synthesis as a plot-only operation from existing tables.
- Preserve the shared denominator and filters unless explicitly changing them.
- Write any new synthesis tables under the new final_figures package.
- Include PNG, PDF, and SVG exports.
- Keep SVG text editable.
- Copy the exact synthesis renderer and this handoff into package provenance.
- Record whether any composite score is used; direct percent-positive marker
  values must remain available in the tables.
```

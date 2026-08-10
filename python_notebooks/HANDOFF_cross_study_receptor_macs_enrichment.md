# Handoff: cross-study receptor-guided enrichment framework

## Objective

Build a new final-figure candidate that asks what cell populations would be
theoretically enriched by positive selection for ERBB4, CXCR4, or PLXNA2 across
the existing organoid studies/protocol cohorts.

The package combines:

1. Cross-study UMAPs colored by sample.
2. Receptor-specific expected target recovery, target purity, and enrichment.
3. Expected post-MACS cell-state composition.

## Canonical output

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/final_figures/fig_cross_study_receptor_macs_enrichment_v1_candidate
```

Primary figure:

```text
figures/png/cross_study_receptor_macs_enrichment.png
figures/pdf/cross_study_receptor_macs_enrichment.pdf
figures/svg/cross_study_receptor_macs_enrichment.svg
```

Expanded context-plus-MACS figure:

```text
figures/png/cross_study_marker_context_and_receptor_macs_enrichment.png
figures/pdf/cross_study_marker_context_and_receptor_macs_enrichment.pdf
```

The expanded figure adds the full PV-precursor ON/OFF-target UMAP atlas as
panel a and the canonical sample-level MGE expression matrix as panel b. The
receptor/MACS framework follows as panels c-e. The two upstream context panels
are raster source figures, so the expanded PDF is a high-resolution raster
assembly; use the standalone source packages for editable vector versions.

## Source data

Existing marker/state genes and UMAP/sample metadata:

```text
results/cross_study_marker_expression/cross_study_marker_expression_v12_pv_precursors_final_candidate_plus_vipr2/tables/per_study/
```

New CXCR4 and PLXNA2 extraction:

```text
results/cross_study_marker_expression/cross_study_marker_expression_v13_receptor_macs/tables/per_study/
```

ERBB4 is loaded from the established v12 tables. CXCR4 and PLXNA2 are exported
from the same canonical Seurat objects used by the cross-study framework, then
joined by study ID and cell ID.

## Model contract

Expression positivity:

```text
expression >= 0.5 log1p(CP10K)
```

Independent target definition:

```text
GAD2 positive
AND (LHX6 positive OR LHX8 positive)
AND MKI67 negative
```

This definition is deliberately independent of ERBB4, CXCR4, and PLXNA2 so
that no candidate capture receptor is advantaged by being part of the target
definition.

Working expected-MACS model:

```text
marker-positive retention probability = 0.80
marker-negative nonspecific retention probability = 0.05
```

For each cell, the expected retained weight is 0.80 or 0.05. Post-sort purity
and composition are the normalized retained weights. Target recovery is the
retained target weight divided by the number of target cells in the input.
Fold enrichment is post-sort target purity divided by input target fraction.

Across-study median ranking excludes cohorts with an input target fraction
below 1%. The cohort-specific table retains every study, but the figure labels
rows with fewer than 25 target cells as `low n` rather than emphasizing an
unstable fold estimate.

Sensitivity tables cover marker-positive retention values of 0.60, 0.80, and
0.95 crossed with background retention values of 0.01, 0.05, and 0.10.

## Cell-state composition categories

Categories are assigned in the following priority order, with cycling applied
last and therefore overriding the other states:

```text
LHX6+ interneuron-like: GAD2+ and LHX6+
LHX8+ ventral-like: GAD2+ and LHX8+ without the LHX6 assignment
NKX2.1+ progenitor-like: NKX2.1+ without the interneuron-like assignments
Dorsal/neurogenic-like: SP8+, PAX6+, or NEUROD2+
Cycling: MKI67+
Other / undetermined
```

## Visualization filters

The figure preserves the established cross-study rules:

```text
Samarasinghe et al. 2021: controls only.
This Study, DIV90: exclude current clusters 6 and 7 as stressed cells.
Other studies: no additional visualization filter.
DIV90 UMAP2: vertically flipped for display only.
```

## Interpretation guardrail

This is an RNA-expression proxy and a computational ranking framework. It does
not demonstrate surface protein abundance, antibody specificity, epitope
accessibility, bead recovery, or live-cell compatibility. A receptor that
ranks well computationally still requires flow cytometry/protein validation
and a MACS reagent titration before experimental use.

The study objects are treated as empirical protocol-cohort proxies. The figure
does not claim that study effects are caused only by protocol; age, cell line,
sampling, QC, and study-specific preprocessing also differ.

## Rerun

Submit the complete extraction-and-render job:

```bash
cp slurm_templates/63_cross_study_receptor_macs_enrichment.sbatch.template \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/63_cross_study_receptor_macs_enrichment.sbatch

sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/63_cross_study_receptor_macs_enrichment.sbatch
```

The renderer can be rerun plot-only after receptor extraction exists:

```bash
export PYTHONPATH="$PWD/python_notebooks/src:${PYTHONPATH:-}"
/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python \
  python_notebooks/scripts/render_cross_study_receptor_macs_enrichment.py \
  --project-root /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder \
  --fresh
```

## Tables

```text
cross_study_receptor_capture_summary.tsv
cross_study_receptor_post_macs_composition.tsv
cross_study_receptor_sample_capture_summary.tsv
cross_study_receptor_macs_sensitivity_summary.tsv
cross_study_receptor_rank_summary.tsv
cross_study_sample_key.tsv
cross_study_receptor_expression_availability.tsv
cross_study_receptor_plot_filter_summary.tsv
```

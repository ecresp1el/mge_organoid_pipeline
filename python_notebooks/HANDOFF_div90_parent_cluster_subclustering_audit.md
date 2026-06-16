# Handoff: DIV90 Parent Cluster Subclustering Audit

Date: 2026-06-14

Status: Slurm array v2 completed successfully. This is a screening audit for
whether existing DIV90 parent clusters should be further subclustered. It is not
a final cell-type annotation update.

## Question

We already have major DIV90 clusters. The goal was not to recluster the whole
object from scratch, but to ask:

```text
Within each existing parent cluster, is there marker-supported substructure that
is not obviously explained by nCount_RNA, nFeature_RNA, percent.mt, sample,
cell line, or cell cycle?
```

Dividing cells were intentionally excluded.

## Inputs

Cached DIV90 AnnData:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/varela_div90.h5ad
```

Parent cluster metadata:

```text
cluster_id
cluster_number_name
```

Technical / experimental covariates tested:

```text
orig.ident
cell_line
Phase
nCount_RNA
nFeature_RNA
percent.mt
S.Score
G2M.Score
CC.Difference
```

`cell_line` was derived from:

```text
metadata/div30_div90_sample_id_to_biolabel_map.tsv
```

## Code Added

Main audit script:

```text
python_notebooks/scripts/div90_parent_cluster_subclustering_audit.py
```

Merge helper:

```text
python_notebooks/scripts/merge_div90_parent_cluster_subclustering_audit.py
```

Slurm array template:

```text
slurm_templates/37_div90_parent_cluster_subclustering_array.sbatch.template
```

Copied job file:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/37_div90_parent_cluster_subclustering_array.sbatch
```

## Slurm Runs

### v1 Historical / Superseded

Array job:

```text
51778536
```

This ran successfully on compute nodes, but the marker-support table was not
valid because marker ranking was performed after scaling `X`. Scanpy's
log-fold-change calculations became unstable on scaled values, and marker
support was incorrectly zero.

Do not interpret v1 recommendations.

### v2 Current Run

Array job:

```text
51778695
```

All 12 array tasks completed:

```text
51778695_0 through 51778695_11: COMPLETED
```

Each task processed one parent cluster:

```text
0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
```

Parent cluster `12 - Dividing cells` was excluded.

Run directory:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_parent_cluster_subclustering_audit/div90_parent_cluster_subclustering_audit_v2
```

Per-parent outputs:

```text
RUN_DIR/per_parent/cluster_<cluster_id>/
```

Merged report:

```text
RUN_DIR/div90_parent_cluster_subclustering_audit_merged_report.md
```

Merged tables:

```text
RUN_DIR/tables/div90_parent_cluster_subclustering_recommendations.tsv
RUN_DIR/tables/div90_parent_cluster_subclustering_resolution_summary.tsv
RUN_DIR/tables/div90_parent_cluster_subclustering_marker_support.tsv
RUN_DIR/tables/div90_parent_cluster_subclustering_top_markers.tsv.gz
RUN_DIR/tables/div90_parent_cluster_subclustering_decision_audit_with_fail_reasons.tsv
RUN_DIR/tables/div90_failed_interesting_clusters_marker_support_top12.tsv
```

Per-parent local UMAP audits:

```text
RUN_DIR/per_parent/cluster_<cluster_id>/plots/parent_cluster_<cluster_id>_local_umap_audit.{png,pdf}
```

## Statistical / Decision Procedure

For each parent cluster, the script:

1. Subsets cells belonging to that parent `cluster_id`.
2. Stores the original log-expression matrix in `layers["log_expr"]`.
3. Selects highly variable genes within the parent cluster.
4. Recomputes PCA, nearest-neighbor graph, and local UMAP within the parent.
5. Runs Leiden clustering over a resolution sweep.
6. Ranks genes for each candidate subcluster using Wilcoxon tests on
   `layers["log_expr"]`, not scaled expression.
7. Counts marker-supported subclusters.
8. Quantifies technical / experimental association of candidate subclusters.
9. Picks the lowest resolution that passes all filters.

Resolution sweep used by v2:

```text
0.1,0.2,0.3,0.4,0.6,0.8,1.0,1.2,1.5
```

Clustering parameters:

```text
n_top_genes = 2000
n_pcs = 30
n_neighbors = 30
min_parent_cells = 100
min_subcluster_cells = 40
```

Marker-support rule:

```text
For each candidate subcluster:
  enriched marker if adjusted P <= 0.05 and log fold-change >= 0.25

A subcluster is marker-supported if:
  it has at least 2 enriched markers

A resolution is marker-supported if:
  at least 75% of its subclusters are marker-supported
```

Technical / experimental confounding metrics:

```text
Categorical covariates:
  orig.ident, cell_line, Phase
  metric = Cramer's V between candidate subcluster and covariate

Numeric covariates:
  nCount_RNA, nFeature_RNA, percent.mt, S.Score, G2M.Score, CC.Difference
  metric = eta-squared from one-way grouping by candidate subcluster

Sample/cell-line domination:
  max fraction of any subcluster made up by one sample or one cell line
```

Technical-confounding flag:

```text
flag if any of:
  max Cramer's V >= 0.65
  max eta-squared >= 0.25
  max dominant sample fraction >= 0.80
  max dominant cell-line fraction >= 0.80
```

Recommendation rule:

```text
candidate_subcluster if a tested resolution has:
  n_subclusters >= 2
  min_subcluster_cells >= 40
  marker_supported_fraction >= 0.75
  technical_confounding_flag == False

Choose the lowest-resolution passing split.

Otherwise:
  keep_parent_only
```

Important interpretation: this is a conservative screening rule. A
`keep_parent_only` call means "do not accept a subcluster split under these
filters yet." It does not prove there is no biological continuum or finer state.

## Current Findings

v2 recommendations with parent cluster names:

| Parent cluster | Parent label | Recommendation | Resolution | Candidate subclusters | Notes |
|---:|---|---|---:|---:|---|
| 0 | MGE Striatal/GP Fated | keep parent only | NA | 1 | marker-supported split was sample/cell-line dominated |
| 1 | SST+, NPY+, Cortical Fated | keep parent only | NA | 1 | marker-supported split failed covariate/confounding filters |
| 2 | CRABP1+/PV Precursors | keep parent only | NA | 1 | marker-supported split strongly sample/cell-line dominated |
| 3 | PV precursors/Migrating cells/Cortical-fated | keep parent only | NA | 1 | marker-supported split failed covariate/confounding filters |
| 4 | Pre-Astrocytes/Astrocytes 1 | keep parent only | NA | 1 | marker-supported split failed covariate/confounding filters |
| 5 | LHX8+ vMGE GABAergic Striatal/GP fated 1 | candidate subcluster | 0.2 | 3 | passes current marker and confounding filters |
| 6 | Stressed Cells | candidate subcluster | 0.3 | 2 | passes filters, but biological value is limited because parent is stressed |
| 7 | Stressed Cells | keep parent only | NA | 1 | marker-supported split failed covariate/confounding filters |
| 8 | LHX8+ vMGE GABAergic Striatal/GP fated 2 | candidate subcluster | 0.2 | 2 | passes current marker and confounding filters |
| 9 | Pre-OPCs/OPCs | candidate subcluster | 0.1 | 2 | passes current marker and confounding filters |
| 10 | Pre-Astrocytes/Astrocytes 2 | keep parent only | NA | 1 | marker-supported split failed covariate/confounding filters |
| 11 | PV Precursors | keep parent only | NA | 1 | low-level split exists but fails current covariate thresholds |

## Decision Values And Failure Reasons

The most useful table for auditing the calls is:

```text
RUN_DIR/tables/div90_parent_cluster_subclustering_decision_audit_with_fail_reasons.tsv
```

This table records, for the audited/best resolution per parent:

```text
parent_cluster_id
parent_cluster_name
recommendation
audited_resolution
n_parent_cells
n_subclusters
min_subcluster_cells
marker_supported_subclusters
marker_supported_fraction
technical_confounding_flag
max_cramers_v
max_eta_squared
max_dominant_sample_fraction
max_dominant_cell_line_fraction
max_technical_effect
top_categorical_metric
top_categorical_value
top_numeric_metric
top_numeric_value
top_dominance_metric
top_dominance_value
explicit_fail_reasons
```

Summary of the exact v2 decision values:

| Parent | Parent label | Audited resolution | n subclusters | Marker-supported fraction | Max Cramer's V | Max eta-squared | Max dominant sample | Max dominant cell line | Decision | Explicit fail reason |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---|---|
| 0 | MGE Striatal/GP Fated | 0.1 | 2 | 1.00 | 0.446 | 0.022 | 1.000 | 1.000 | keep parent only | max dominant sample and cell-line fractions >= 0.80 |
| 1 | SST+, NPY+, Cortical Fated | 0.2 | 2 | 1.00 | 0.821 | 0.019 | 0.434 | 0.594 | keep parent only | max Cramer's V >= 0.65 |
| 2 | CRABP1+/PV Precursors | 0.2 | 3 | 1.00 | 0.466 | 0.020 | 0.973 | 0.973 | keep parent only | max dominant sample and cell-line fractions >= 0.80 |
| 3 | PV precursors/Migrating cells/Cortical-fated | 0.2 | 2 | 1.00 | 0.910 | 0.133 | 0.812 | 0.838 | keep parent only | max Cramer's V >= 0.65 plus dominant sample/cell-line fractions >= 0.80 |
| 4 | Pre-Astrocytes/Astrocytes 1 | 0.1 | 2 | 1.00 | 0.844 | 0.091 | 0.530 | 0.566 | keep parent only | max Cramer's V >= 0.65 |
| 5 | LHX8+ vMGE GABAergic Striatal/GP fated 1 | 0.2 | 3 | 1.00 | 0.336 | 0.012 | 0.630 | 0.675 | candidate subcluster | passed all audit filters |
| 6 | Stressed Cells | 0.3 | 2 | 1.00 | 0.357 | 0.100 | 0.293 | 0.469 | candidate subcluster | passed all audit filters |
| 7 | Stressed Cells | 0.2 | 2 | 1.00 | 0.868 | 0.126 | 0.963 | 0.969 | keep parent only | max Cramer's V >= 0.65 plus dominant sample/cell-line fractions >= 0.80 |
| 8 | LHX8+ vMGE GABAergic Striatal/GP fated 2 | 0.2 | 2 | 1.00 | 0.484 | 0.005 | 0.534 | 0.540 | candidate subcluster | passed all audit filters |
| 9 | Pre-OPCs/OPCs | 0.1 | 2 | 1.00 | 0.536 | 0.185 | 0.321 | 0.459 | candidate subcluster | passed all audit filters |
| 10 | Pre-Astrocytes/Astrocytes 2 | 0.2 | 3 | 1.00 | 0.710 | 0.110 | 0.958 | 0.989 | keep parent only | min subcluster cells < 40 plus Cramer's V/sample/cell-line failures |
| 11 | PV Precursors | 0.4 | 2 | 1.00 | 0.582 | 0.306 | 0.544 | 0.544 | keep parent only | max eta-squared >= 0.25 |

How to read these values:

- `marker_supported_fraction = 1.00` means every candidate subcluster at that
  resolution had at least two enriched genes by the marker rule.
- A parent can still fail if the candidate split is too aligned with sample,
  cell line, phase, QC, or cell-cycle covariates.
- `max_cramers_v` is the strongest categorical covariate association observed
  among `orig.ident`, `cell_line`, and `Phase`.
- `max_eta_squared` is the strongest numeric covariate association observed
  among `nCount_RNA`, `nFeature_RNA`, `percent.mt`, `S.Score`, `G2M.Score`,
  and `CC.Difference`.
- `max_dominant_sample_fraction` asks whether any candidate subcluster is mostly
  one sample.
- `max_dominant_cell_line_fraction` asks whether any candidate subcluster is
  mostly one cell line.

## How Marker Support Is Captured

Markers are captured in two linked tables:

```text
RUN_DIR/tables/div90_parent_cluster_subclustering_marker_support.tsv
RUN_DIR/tables/div90_parent_cluster_subclustering_top_markers.tsv.gz
```

For each parent cluster and each tested Leiden resolution:

1. The script performs `scanpy.tl.rank_genes_groups(..., method="wilcoxon")`.
2. The test is run on `layers["log_expr"]`, which stores the original log-like
   DIV90 H5AD expression before scaling.
3. A gene is counted as an enriched marker for a candidate subcluster if:

```text
adjusted P <= 0.05
log fold-change >= 0.25
```

4. A candidate subcluster is considered marker-supported if it has at least two
   enriched genes.
5. A resolution is marker-supported if at least 75% of its candidate subclusters
   are marker-supported.

Important caveat: marker support by itself is not sufficient for the conservative
audit call. A split can have strong markers and still be flagged if those
subclusters align strongly with sample, cell line, QC depth, or cell-cycle
covariates.

Because cell line may be part of the biology in this experiment, cell-line
association should be interpreted carefully. The current audit treats it as a
conservative warning, not as proof that the split is artifactual.

Candidate splits and top markers at recommended resolution:

```text
Parent 5, resolution 0.2, 3 subclusters:
  subcluster 0: n=1028, top marker GRIK2, 174 enriched markers
  subcluster 1: n=631,  top marker SOX4,  146 enriched markers
  subcluster 2: n=265,  top marker PEG10, 159 enriched markers

Parent 6, resolution 0.3, 2 subclusters:
  subcluster 0: n=744, top marker GRIK2,  68 enriched markers
  subcluster 1: n=482, top marker POLR2A, 101 enriched markers

Parent 8, resolution 0.2, 2 subclusters:
  subcluster 0: n=494, top marker ARHGAP36, 146 enriched markers
  subcluster 1: n=428, top marker DPYD,    151 enriched markers

Parent 9, resolution 0.1, 2 subclusters:
  subcluster 0: n=560, top marker SMOC1, 374 enriched markers
  subcluster 1: n=355, top marker DCX,   127 enriched markers
```

## Important Biological Caveats

Clusters `1`, `3`, and `11` are biologically interesting from marker-panel UMAPs,
but the subclustering audit did not accept their Leiden splits under the current
filters.

This matters:

- Cluster `3` showed marker-supported substructure, but the best split had
  strong covariate association:
  `max_technical_effect = 0.9097`.
- Cluster `11` showed a two-way split at resolution `0.4`, but it failed because
  `eta_squared` for numeric technical/cell-cycle-like covariates reached
  `0.3064`, above the `0.25` threshold.
- Cluster `1` showed a two-way split at resolution `0.2`, but it failed because
  Cramer's V reached `0.8209`, above the `0.65` threshold.

So the conservative interpretation is:

```text
Do not promote clusters 1, 3, or 11 to official subclusters yet.
Review their local UMAP audit plots and marker tables if we want a more
biology-forward, less conservative call.
```

### Marker-Supported But Failed Splits: Clusters 1, 3, 11

The table below records the actual marker evidence for the three biologically
interesting failed splits:

```text
RUN_DIR/tables/div90_failed_interesting_clusters_marker_support_top12.tsv
```

Cluster `1 - SST+, NPY +, Cortical Fated`

Audited split:

```text
resolution = 0.2
n_subclusters = 2
marker_supported_fraction = 1.00
failure = max_cramers_v >= 0.65
```

Exact failure values:

```text
orig.ident_cramers_v = 0.8209
cell_line_cramers_v = 0.6491
Phase_cramers_v = 0.0462
max_eta_squared = 0.0194
max_dominant_sample_fraction = 0.4335
max_dominant_cell_line_fraction = 0.5943
```

Marker support:

| Subcluster | Cells | n enriched markers | Top marker | Top enriched markers |
|---:|---:|---:|---|---|
| 0 | 2,249 | 173 | `EDIL3` | `EDIL3`, `VWC2`, `LRFN5`, `ITGB8`, `CRABP1`, `NRCAM`, `RAB3IP`, `ADRA1A`, `ADGRB1`, `PCDH10`, `CDH4`, `RYR3` |
| 1 | 1,299 | 116 | `ACKR3` | `ACKR3`, `ZEB2`, `LCN9`, `BTG1`, `GPR173`, `SLAIN1`, `PLS3`, `CEP85L`, `PDZRN3`, `CSAD`, `CCSER1`, `WLS` |

Interpretation: this split has plausible biological markers, including an
`ACKR3/ZEB2`-rich side versus an `EDIL3/VWC2/LRFN5/ITGB8/CRABP1`-rich side.
It failed because the split was strongly associated with sample identity. Since
sample and cell line are coupled in this dataset, this should be treated as a
warning rather than an automatic rejection.

Cluster `3 - PV precursors/Migrating cells/Cortical-fated`

Audited split:

```text
resolution = 0.2
n_subclusters = 2
marker_supported_fraction = 1.00
failure = max_cramers_v >= 0.65; max dominant sample/cell-line fractions >= 0.80
```

Exact failure values:

```text
orig.ident_cramers_v = 0.9097
cell_line_cramers_v = 0.7395
Phase_cramers_v = 0.0209
nFeature_RNA_eta_squared = 0.1326
nCount_RNA_eta_squared = 0.1273
max_dominant_sample_fraction = 0.8117
max_dominant_cell_line_fraction = 0.8377
```

Marker support:

| Subcluster | Cells | n enriched markers | Top marker | Top enriched markers |
|---:|---:|---:|---|---|
| 0 | 1,975 | 76 | `L3MBTL1` | `L3MBTL1`, `PEG3`, `SYT4`, `LRFN5`, `RUNX1T1`, `MYO16`, `MAGEH1`, `EPHA7`, `EPHA6`, `TMEFF2`, `PLXNA4`, `TENT5A` |
| 1 | 308 | 168 | `BEX3` | `BEX3`, `USP11`, `TMSB4X`, `NAP1L3`, `CRABP1`, `MEST`, `DCX`, `ACTB`, `H1FX`, `TXNIP`, `TUBB2B`, `NNAT` |

Interpretation: this is the strongest "marker-supported but confounded" case.
The smaller subcluster has strong developmental/neuronal markers including
`BEX3`, `CRABP1`, `DCX`, `TUBB2B`, and `NNAT`, while the larger subcluster has
a different gene set including `L3MBTL1`, `PEG3`, `SYT4`, `LRFN5`, and
`RUNX1T1`. It failed because the split was highly aligned with sample and cell
line. This may still be biological if sample/cell-line captures real genetic or
differentiation-state differences, but it is not cleanly separable from those
experimental variables in the current audit.

Cluster `11 - PV Precursors`

Audited split:

```text
resolution = 0.4
n_subclusters = 2
marker_supported_fraction = 1.00
failure = max_eta_squared >= 0.25
```

Exact failure values:

```text
orig.ident_cramers_v = 0.5819
cell_line_cramers_v = 0.4863
Phase_cramers_v = 0.0442
nCount_RNA_eta_squared = 0.3064
nFeature_RNA_eta_squared = 0.2986
percent.mt_eta_squared = 0.1325
max_dominant_sample_fraction = 0.5444
max_dominant_cell_line_fraction = 0.5444
```

Marker support:

| Subcluster | Cells | n enriched markers | Top marker | Top enriched markers |
|---:|---:|---:|---|---|
| 0 | 245 | 16 | `CCSER1` | `CCSER1`, `ZEB2`, `ACKR3`, `SLITRK2`, `MAF`, `RASA1`, `PRKG1`, `HCN1`, `NPAS3`, `WSB1`, `PTCHD4`, `BCLAF3` |
| 1 | 180 | 159 | `CPE` | `CPE`, `TUBB2A`, `PEG10`, `CNTN1`, `PCDH10`, `MT-CO3`, `DKK3`, `MT-ND4L`, `PCDH17`, `MEST`, `MT-ATP6`, `ATP1B1` |

Interpretation: unlike clusters `1` and `3`, cluster `11` did not fail mainly
because of sample or cell line. It failed because the split was associated with
RNA depth / detected genes (`nCount_RNA` and `nFeature_RNA`). The marker split
is still biologically interesting: one side is `ACKR3/ZEB2/MAF`-leaning, while
the other is `CPE/TUBB2A/PEG10/CNTN1`-leaning, but the depth association needs
manual review before promoting it.

## What "Reasonable Number Of Clusters" Means Here

Under the current conservative rules, the reasonable number is:

```text
Parent clusters retained as-is:
  0,1,2,3,4,7,10,11 = 8 parent states

Accepted candidate subcluster parents:
  5 -> 3 subclusters
  6 -> 2 subclusters
  8 -> 2 subclusters
  9 -> 2 subclusters
```

This gives a provisional screened structure of:

```text
8 retained parent clusters
+ 3 + 2 + 2 + 2 accepted subclusters replacing parents 5,6,8,9
= 17 screened states among non-dividing DIV90 cells
```

However, cluster `6` is a stressed-cell cluster, so the biologically useful
number may be:

```text
15 non-dividing, non-stress screened states
```

depending on whether stressed-cell substructure is retained or tracked
separately.

## Recommended Next Review

## Follow-up Validation Of Clusters 1, 3, And 11

After the parent-cluster audit, a targeted validation run was added for the
three biologically interesting but audit-failed neuronal splits:

```text
python_notebooks/scripts/div90_validate_failed_interesting_subclusters.py
slurm_templates/38_div90_validate_failed_interesting_subclusters.sbatch.template
```

Successful Slurm run:

```text
job = 51779685
run label = div90_failed_interesting_subcluster_validation_v2
status = COMPLETED
elapsed = 00:01:05
max RSS = 3.9 GB
```

Output directory:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_failed_interesting_subcluster_validation/div90_failed_interesting_subcluster_validation_v2
```

Main interpretation file:

```text
OUTDIR/div90_failed_interesting_subcluster_validation_interpretation.md
```

Cluster `1 - SST+, NPY+, Cortical Fated`

Deliverables:

```text
plots/cluster1_local_umap_marker_axis.png
plots/cluster1_module_scores_by_sample.png
plots/cluster1_module_scores_by_cell_line.png
tables/cluster1_sample_composition_table.tsv
tables/cluster1_cell_line_composition_table.tsv
tables/cluster1_within_sample_module_correlation_table.tsv
```

Module axes:

```text
cluster1_EDIL3_CRABP1_axis:
  EDIL3, VWC2, LRFN5, ITGB8, CRABP1, NRCAM, PCDH10, CDH4, RYR3

cluster1_ACKR3_ZEB2_axis:
  ACKR3, ZEB2, BTG1, GPR173, CSAD, WLS
```

Key values:

```text
Total cells = 3,548
subcluster 0 = 2,249 cells, EDIL3/CRABP1 axis
subcluster 1 = 1,299 cells, ACKR3/ZEB2 axis

By cell line:
  2E:  90.1% subcluster 0,  9.9% subcluster 1
  79B: 74.0% subcluster 0, 26.0% subcluster 1
  H9:   9.0% subcluster 0, 91.0% subcluster 1

Strongest within-sample anti-correlation:
  10496-MW-3 Pearson = -0.327, Spearman = -0.335
```

Interpretation:

```text
The marker axis is biologically plausible, and both sides are seen in more than
one sample. However, the split is strongly sample/cell-line structured. H9 is
mostly ACKR3/ZEB2-high, while 2E and 79B are mostly EDIL3/CRABP1-high. Keep
cluster 1 as parent-only for official labels, but preserve the axis as a
candidate cortical interneuron maturation/migration hypothesis.
```

Cluster `3 - PV precursors/Migrating cells/Cortical-fated`

Deliverables:

```text
plots/cluster3_local_umap_marker_axis.png
plots/cluster3_neighbor_cluster_marker_dotplot.png
plots/cluster3_qc_marker_review.png
tables/cluster3_sample_composition_table.tsv
tables/cluster3_cell_line_composition_table.tsv
tables/cluster3_neighbor_cluster_marker_dotplot_table.tsv
tables/cluster3_qc_marker_review_summary.tsv
```

Module axes:

```text
cluster3_L3MBTL1_SYT4_axis:
  L3MBTL1, PEG3, SYT4, LRFN5, RUNX1T1, MYO16, MAGEH1, TMEFF2, PLXNA4, TENT5A

cluster3_BEX3_CRABP1_DCX_axis:
  BEX3, USP11, TMSB4X, NAP1L3, CRABP1, MEST, DCX, TUBB2B, NNAT
```

Key values:

```text
Total cells = 2,283
subcluster 0 = 1,975 cells, L3MBTL1/SYT4 axis
subcluster 1 = 308 cells, BEX3/CRABP1/DCX axis

Subcluster 1 by sample:
  10496-MW-2: 250 cells, 91.6% of that sample's cluster-3 cells
  10496-MW-5:  35 cells, 92.1% of that sample's cluster-3 cells
  other samples: 4-8 cells each

Subcluster 1 by cell line:
  H9: 71.7% of H9 cluster-3 cells
  2E:  5.7% of 2E cluster-3 cells
  79B: 0.8% of 79B cluster-3 cells
```

QC/stress values:

```text
subcluster 0:
  percent.mt = 0.298
  stress module = 0.314
  nCount_RNA = 2819
  nFeature_RNA = 1994

subcluster 1:
  percent.mt = 0.343
  stress module = 0.506
  nCount_RNA = 4057
  nFeature_RNA = 2631
```

The H5AD var names contained no `RPL*` or `RPS*` symbols, so the ribosomal
module could not be scored from this export.

Interpretation:

```text
The BEX3/CRABP1/DCX-high state has real neuroblast-like marker support
including BEX3, CRABP1, DCX, TUBB2B, NNAT, and MEST. It is not purely absent in
other samples, but it is overwhelmingly concentrated in H9, especially
10496-MW-2 and the small 10496-MW-5 sample. It also has higher nCount/nFeature
and higher stress score. Keep cluster 3 parent-only for official labels; record
this as a sample/cell-line-enriched candidate neuroblast-like state requiring
replication.
```

Cluster `11 - PV Precursors`

Deliverables:

```text
plots/cluster11_local_umap_marker_axis.png
plots/cluster11_nfeature_tertile_marker_scores.png
plots/cluster11_ncount_tertile_marker_scores.png
plots/cluster11_cluster3_comparison_marker_dotplot.png
plots/cluster11_regressed_ncount_nfeature_local_umap.png
tables/cluster11_nfeature_tertile_marker_scores.tsv
tables/cluster11_ncount_tertile_marker_scores.tsv
tables/cluster11_covariate_adjusted_model_results.tsv
tables/cluster11_cluster3_comparison_marker_dotplot_table.tsv
tables/cluster11_regressed_reclustering_comparison.tsv
```

Module axes:

```text
cluster11_ACKR3_ZEB2_MAF_axis:
  CCSER1, ZEB2, ACKR3, SLITRK2, MAF, RASA1, PRKG1, HCN1, WSB1, PTCHD4

cluster11_CPE_TUBB2A_CNTN1_axis:
  CPE, TUBB2A, PEG10, CNTN1, PCDH10, DKK3, MEST, ATP1B1
```

nFeature/nCount tertile values:

```text
nFeature low:    85.9% subcluster 0, 14.1% subcluster 1
nFeature middle: 63.8% subcluster 0, 36.2% subcluster 1
nFeature high:   23.2% subcluster 0, 76.8% subcluster 1

nCount low:      85.9% subcluster 0, 14.1% subcluster 1
nCount middle:   63.8% subcluster 0, 36.2% subcluster 1
nCount high:     23.2% subcluster 0, 76.8% subcluster 1
```

Covariate-adjusted model:

```text
module_score ~ audit_subcluster + nFeature_RNA + nCount_RNA + percent.mt + orig.ident + cell_line

ACKR3/ZEB2/MAF module:
  audit_subcluster 1 coefficient = -0.447, p = 3.81e-16, R2 = 0.394

CPE/TUBB2A/CNTN1 module:
  audit_subcluster 1 coefficient =  0.462, p = 3.26e-17, R2 = 0.639
```

Regression sensitivity:

```text
ARI original audit split vs nCount/nFeature-regressed reclustering = 0.033
```

Interpretation:

```text
The module-score difference remains statistically associated with audit
subcluster after covariate adjustment. However, the split is still strongly
organized by transcriptome complexity: low nFeature/nCount cells are mostly
ACKR3/ZEB2/MAF-high subcluster 0, while high nFeature/nCount cells are mostly
CPE/TUBB2A/CNTN1-high subcluster 1. The nCount/nFeature-regressed reclustering
does not reproduce the original split. Keep cluster 11 parent-only for official
labels, but track it as a candidate PV-precursor maturation/neuroblast-to-
precursor gradient requiring matched-depth validation.
```

Overall validation call:

```text
Do not promote clusters 1, 3, or 11 to official subclusters yet.

Preserve these as biological hypotheses:
  cluster 1: EDIL3/CRABP1 vs ACKR3/ZEB2 cortical interneuron axis
  cluster 3: H9-enriched BEX3/CRABP1/DCX neuroblast-like candidate state
  cluster 11: ACKR3/ZEB2/MAF to CPE/TUBB2A/CNTN1 PV-precursor maturation gradient
```

## DIV90 Neuron-Only Reclustering

After the per-parent audit and targeted validation, a separate neuron-only
reclustering was run. This is important because it asks a different question:
instead of testing substructure inside each original cluster, it reclusters all
neuronal DIV90 cells together as one neuronal compartment.

Code:

```text
python_notebooks/scripts/div90_neuron_only_reclustering.py
slurm_templates/43_div90_neuron_only_reclustering.sbatch.template
```

Successful Slurm run:

```text
job = 51780411
run label = div90_neuron_only_reclustering_v1
status = COMPLETED
elapsed = 00:04:33
max RSS = 5.3 GB
```

Output directory:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_neuron_only_reclustering/div90_neuron_only_reclustering_v1
```

Report:

```text
OUTDIR/neuron_only_reclustering_report.md
```

Important terminology:

```text
"Parent cluster" in this section means the original DIV90 cluster labels already
present in the input H5AD: cluster_id and cluster_number_name.

It does not mean a URD parent, lineage-tree parent, or developmental ancestor.
The neuron-only clusters are new Leiden clusters computed after subsetting to
the neuronal parent clusters only.
```

Neuron-only subset:

```text
Included original DIV90 clusters:
  0 = MGE Striatal/GP Fated
  1 = SST+, NPY+, Cortical Fated
  2 = CRABP1+/PV Precursors
  3 = PV Precursors/Migrating Cells/Cortical-Fated
  5 = LHX8+ vMGE GABAergic Striatal/GP fated 1
  8 = LHX8+ vMGE GABAergic Striatal/GP fated 2
  11 = PV Precursors

Excluded original DIV90 clusters:
  4 = Pre-Astrocytes/Astrocytes 1
  6 = Stressed Cells
  7 = Stressed Cells
  9 = Pre-OPCs/OPCs
  10 = Pre-Astrocytes/Astrocytes 2
  12 = Dividing cells

Neuron-only cells = 16,206
```

Resolution sweep:

```text
0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0
```

Recommended neuron-only resolution:

```text
resolution = 0.6
subcluster key = neuron_leiden_r0_6
number of neuron-only states = 9
minimum subcluster size = 896
marker-supported fraction = 1.0
```

How the resolution was selected:

```text
score = 4 * marker_supported_fraction
      + cluster_count_bonus
      + min(min_subcluster_size / 150, 1)
      + 0.8 * mean_parent_purity
      - 2 * max(0, max_cell_line_dominance - 0.85)
      - 2 * max(0, max_sample_dominance - 0.75)

cluster_count_bonus = +1 if 7-12 clusters, otherwise -0.25 * distance from 9 clusters
```

A neuron-only subcluster was considered marker-supported if it had at least two
positive markers among the ranked genes:

```text
adjusted P <= 0.05
log fold-change >= 0.25
```

So "good enough" here means:

```text
marker-supported;
not too many tiny clusters;
enough clusters to separate major neuronal programs;
not so many clusters that the map fragments;
no single sample/cell line completely dominates the overall resolution.
```

Resolution summary:

| Resolution | n neuron clusters | Min size | Marker-supported fraction | Mean parent purity | Max cell-line dominance | Recommendation |
|---:|---:|---:|---:|---:|---:|---|
| 0.1 | 3 | 1861 | 1.00 | 0.679 | 0.591 | no, too coarse |
| 0.2 | 5 | 920 | 1.00 | 0.726 | 0.600 | no, still coarse |
| 0.3 | 6 | 922 | 1.00 | 0.750 | 0.598 | no, still modestly coarse |
| 0.4 | 8 | 102 | 1.00 | 0.747 | 0.971 | no, tiny/line-dominated cluster |
| 0.6 | 9 | 896 | 1.00 | 0.731 | 0.949 | recommended |
| 0.8 | 14 | 330 | 1.00 | 0.752 | 0.949 | no, fragmented |
| 1.0 | 16 | 282 | 1.00 | 0.771 | 0.961 | no, fragmented |

### Working Labels For The 9 Neuron-Only Clusters

These are working labels, not final official annotations. They should be used as
interpretation handles for plots, URD/tip selection, and follow-up validation.

| Neuron-only cluster | Cells | Working interpretation | Main original source | Marker/support clues |
|---:|---:|---|---|---|
| 0 | 3,461 | Immature MGE striatal/GP-fated neuroblast-like state | mostly original cluster 0 | `SOX5`, `SLIT2`, `PEG10`, `STMN2`, `TUBB2B`, `DCX`, `CPE`, `NKX2-1` |
| 1 | 2,581 | SST/NPY cortical-fated interneuron state | mostly original cluster 3, plus cluster 1 | `SST`, `CHODL`, `LRFN5`, `EDIL3`, `NPY`, `ERBB4`, `GAD1` |
| 2 | 2,175 | ACKR3/ZEB2/MAF PV/cortical precursor-like state | mostly original cluster 2, plus clusters 11/1 | `ACKR3`, `ZEB2`, `MAF`, `ERBB4`, `GAD1`, `DCX` |
| 3 | 1,927 | CRABP1+/CXCR4+ PV precursor / migratory state | mostly original cluster 2, plus cluster 1 | `CRABP1`, `NFIB`, `ERBB4`, `CXCR4`, `LHX6`, `GAD1` |
| 4 | 1,866 | LHX8+ ventral MGE / striatal-GP GABAergic state 1 | very strongly original cluster 5 | `LHX8`, `GABRA1`, `GRIK2`, `TMEFF2`, `CNTN1`, `CPE` |
| 5 | 1,206 | 79B-biased cortical/PV-like EDIL3/VWC2/ERBB4 state | mostly original clusters 1 and 2 | `VWC2`, `VSTM2A`, `ERBB4`, `EDIL3`, `ZEB2`, `SOX6`; 94% 79B |
| 6 | 1,153 | 2E-biased SST/ERBB4/GAD1 cortical interneuron-like state | mostly original cluster 1 | `VSTM2A`, `VWC2`, `ERBB4`, `GAD1`, `SST`; 95% 2E |
| 7 | 941 | LHX8+ ventral MGE / striatal-GP GABAergic state 2 | mostly original cluster 8 | `LHX8`, `CDH23`, `ZFHX3`, `TENM1`, `ACHE`, `NKX2-1` |
| 8 | 896 | H9-biased ACKR3/ZEB2-related immature/cortical-PV mixed state | mixed original clusters 1, 2, and 0 | `LCN9`, `CSAD`, `PDZRN3`, `NR2F1`, `MEIS2`, `DLX5`, `ZEB2`; 84.5% H9 |

Broad grouping:

```text
Striatal/GP / ventral MGE side:
  neuron-only clusters 0, 4, 7

Cortical/SST/PV precursor side:
  neuron-only clusters 1, 2, 3, 5, 6, 8

Line-biased but biologically plausible:
  neuron-only cluster 5 = 79B-biased
  neuron-only cluster 6 = 2E-biased
  neuron-only cluster 8 = H9-biased
```

### Candidate Tips / Most Differentiated States

A small marker-score table was added:

```text
OUTDIR/neuron_only_candidate_tip_maturity_scores.tsv
```

This is not a URD pseudotime result. It is a marker-based prioritization of
which neuron-only clusters look more differentiated versus more immature.

Marker sets used:

```text
immature/migration:
  DCX, STMN2, DCLK1, BEX3, TUBB2A, CPE, CNTN1

pan-GABA maturation:
  GAD1, GAD2, ERBB4

SST/cortical maturation:
  SST, NPY, SATB1, MAF, MAFB, ERBB4, GAD1, GAD2

PV maturation:
  SOX6, ERBB4, KCNC1, KCNC2, GAD1, GAD2

striatal/GP specification:
  LHX8, ISL1, GBX2, TAC1, NKX2-1, GAD1, GAD2
```

Tip-like score:

```text
maturity_minus_immature =
  average(pan-GABA maturation, strongest lineage maturation score)
  - immature/migration score
```

Current ranking from the canonical marker panel:

| Rank | Neuron-only cluster | Tip-like score | Interpretation |
|---:|---:|---:|---|
| 1 | 2 | 0.743 | strongest candidate mature PV/cortical precursor-like tip; `ACKR3/ZEB2/MAF/ERBB4/GAD1` |
| 2 | 3 | 0.697 | strong candidate mature CRABP1+/CXCR4+ PV/migratory tip; `CRABP1/ERBB4/LHX6/GAD1` |
| 3 | 5 | 0.523 | possible differentiated cortical/PV-like tip, but strongly 79B-biased |
| 4 | 6 | 0.484 | possible differentiated SST/ERBB4/GAD1 tip, but strongly 2E-biased |
| 5 | 1 | 0.034 | SST/NPY cortical-fated candidate; has terminal markers but still high immature markers |
| 6 | 8 | -0.090 | H9-biased mixed immature/cortical-PV state; not a clean tip |
| 7 | 7 | -0.765 | LHX8+ striatal/GP state 2; biologically coherent but immature/intermediate by marker panel |
| 8 | 4 | -1.029 | LHX8+ striatal/GP state 1; biologically coherent but immature/intermediate by marker panel |
| 9 | 0 | -1.225 | immature MGE striatal/GP neuroblast-like state |

Interpretation for URD/tip choice:

```text
Best candidate differentiated tips from the neuron-only data:
  cluster 2 = ACKR3/ZEB2/MAF PV/cortical precursor-like tip
  cluster 3 = CRABP1+/CXCR4+ PV precursor / migratory tip

Secondary tip candidates, but line-biased:
  cluster 5 = 79B-biased EDIL3/VWC2/ERBB4 cortical/PV-like state
  cluster 6 = 2E-biased SST/ERBB4/GAD1 cortical interneuron-like state

Candidate cortical/SST tip with caution:
  cluster 1 = SST/NPY cortical-fated state, but still high DCX/STMN2/CPE,
              so it may be a maturing/intermediate cortical branch rather than
              the most differentiated terminal state.

Not good terminal tips yet:
  clusters 0, 4, 7 = strong immature/migration scores, likely intermediates
                    despite coherent MGE/striatal-GP identities.
  cluster 8 = H9-biased mixed immature state, not a clean terminal tip.
```

### Science Data S9 / Five-Class Reference Module Scoring

Module-scoring index terms for future search:

```text
Jia, S9 module, Science Data S9, Module 1, Module 5, TF-only, TFs,
science.adw1803_data_s9.xlsx, div90_neuron_only_reference_module_scoring
```

The user provided the Science Data S9 module workbook:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/reference/science.adw1803_data_s9.xlsx
```

This workbook contains five full module sheets and five TF-only module sheets:

```text
Module 1, Module 2, Module 3, Module 4, Module 5
Module 1（TFs）, Module 2（TFs）, Module 3（TFs）, Module 4（TFs）, Module 5（TFs）
```

The user also identified five major class anchors:

```text
EPHA5/MEF2C
LHX6/NFIA
CRABP1/ANGPT2
NR2F1/NR2F2
LHX8/ISL1
```

New scoring script:

```text
python_notebooks/scripts/div90_neuron_only_reference_module_scoring.py
```

Slurm template:

```text
slurm_templates/44_div90_neuron_only_reference_module_scoring.sbatch.template
```

Submitted Slurm job:

```text
job = 51793516
run label = div90_neuron_only_reference_modules_v1
status = COMPLETED
elapsed = 00:01:43
max RSS = 4.5 GB
```

Planned output directory:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_neuron_only_reference_modules/div90_neuron_only_reference_modules_v1
```

Outputs:

```text
div90_neuron_only_s9_reference_module_report.md
div90_neuron_only_s9_gene_set_coverage.tsv
div90_neuron_only_s9_module_scores_by_neuron_cluster.tsv
div90_neuron_only_s9_module_scores_by_parent_cluster.tsv
div90_neuron_only_s9_module_scores_by_cell_line.tsv
div90_neuron_only_s9_module_scores_by_sample.tsv
div90_neuron_only_s9_best_module_calls.tsv
div90_neuron_only_anchor_class_heatmap.png/pdf
div90_neuron_only_s9_full_module_heatmap.png/pdf
div90_neuron_only_s9_tf_module_heatmap.png/pdf
div90_neuron_only_anchor_class_umaps.png/pdf
div90_neuron_only_s9_full_module_umaps.png/pdf
div90_neuron_only_s9_tf_module_umaps.png/pdf
```

Important design choice:

```text
The Excel sheet numbers are not being forced to equal the five biological class
anchors. The script scores both:

1. The workbook modules exactly as given.
2. The five user-provided anchor classes separately.

This avoids prematurely forcing a S9 module number onto a biological class if
the anchor genes and full-module sheets do not align perfectly.
```

Gene-set coverage:

```text
EPHA5/MEF2C:      2/2 genes present
LHX6/NFIA:        2/2 genes present
CRABP1/ANGPT2:    1/2 genes present; ANGPT2 missing, so this is CRABP1-driven
NR2F1/NR2F2:      2/2 genes present
LHX8/ISL1:        2/2 genes present

S9 full module coverage:
  module 1: 236/298 genes present
  module 2: 103/118 genes present
  module 3: 216/228 genes present
  module 4: 336/362 genes present
  module 5: 264/393 genes present

S9 TF-module coverage:
  module 1 TFs:  6/6 genes present
  module 2 TFs:  7/7 genes present
  module 3 TFs: 12/13 genes present
  module 4 TFs: 14/14 genes present
  module 5 TFs: 13/13 genes present
```

Best five-class anchor calls for the recommended neuron-only clusters:

| Neuron cluster | Best anchor class | Best z | Margin vs second | Interpretation |
|---:|---|---:|---:|---|
| 0 | LHX8/ISL1 | -0.29 | 0.07 | weak anchor match; better treated as immature MGE/striatal-GP neuroblast-like |
| 1 | LHX6/NFIA | 0.80 | 0.37 | moderate LHX6/NFIA-biased cortical/PV/SST state |
| 2 | EPHA5/MEF2C | 1.78 | 1.18 | strong EPHA5/MEF2C-like state |
| 3 | CRABP1/ANGPT2 | 2.41 | 0.84 | strong CRABP1-driven state; ANGPT2 missing |
| 4 | LHX8/ISL1 | 1.43 | 1.69 | strong LHX8/ISL1 striatal/GP-like state |
| 5 | NR2F1/NR2F2 | 1.08 | 0.12 | mixed NR2F1/2 vs LHX6/NFIA; 79B-biased |
| 6 | NR2F1/NR2F2 | 0.72 | 0.10 | mixed NR2F1/2 vs EPHA5/MEF2C; 2E-biased |
| 7 | LHX8/ISL1 | 2.23 | 1.52 | strongest LHX8/ISL1 striatal/GP-like state |
| 8 | NR2F1/NR2F2 | 1.56 | 1.55 | strong NR2F1/NR2F2-like state; H9-biased |

What appears present in DIV90 neurons:

```text
Clear EPHA5/MEF2C-like class:
  neuron-only cluster 2

Clear CRABP1/ANGPT2-like class:
  neuron-only cluster 3, but call is CRABP1-driven because ANGPT2 is missing
  from the H5AD.

Clear LHX8/ISL1-like class:
  neuron-only clusters 4 and 7

Clear NR2F1/NR2F2-like class:
  neuron-only cluster 8

Moderate / mixed LHX6/NFIA-like class:
  neuron-only cluster 1 is the cleanest LHX6/NFIA-biased state, but it is not
  as sharply separated as the other anchor classes.

Line-biased mixed states:
  cluster 5 = NR2F1/2 vs LHX6/NFIA mixed, 79B-biased
  cluster 6 = NR2F1/2 vs EPHA5/MEF2C mixed, 2E-biased

Weakly captured by the five anchor classes:
  cluster 0 = immature MGE/striatal-GP neuroblast-like state
```

Revised candidate-tip interpretation after S9 scoring:

```text
Best endpoint/tip candidates:
  cluster 2 = EPHA5/MEF2C-like
  cluster 3 = CRABP1-driven CRABP1/ANGPT2-like
  cluster 4 = LHX8/ISL1-like striatal/GP state
  cluster 7 = LHX8/ISL1-like striatal/GP state
  cluster 8 = NR2F1/NR2F2-like, H9-biased

Secondary / mixed tip candidates:
  cluster 1 = LHX6/NFIA-biased but moderate
  cluster 5 = mixed NR2F1/2 and LHX6/NFIA, 79B-biased
  cluster 6 = mixed NR2F1/2 and EPHA5/MEF2C, 2E-biased

Not a good terminal tip:
  cluster 0 = immature state, weakly matched to all five anchor classes
```

### Final Working Branch Model For Next DIV90 All-Cell URD

Recap of what was done in this thread:

```text
1. Per-parent DIV90 subclustering audit:
   - tested substructure within original parent clusters
   - kept clusters 1, 3, and 11 parent-only under conservative confounding gates
   - documented marker-supported but sample/depth-associated splits

2. Targeted validation of clusters 1, 3, and 11:
   - cluster 1: EDIL3/CRABP1 vs ACKR3/ZEB2 axis, sample/cell-line structured
   - cluster 3: BEX3/CRABP1/DCX neuroblast-like state, H9-enriched
   - cluster 11: ACKR3/ZEB2/MAF vs CPE/TUBB2A/CNTN1, nCount/nFeature-associated

3. Neuron-only reclustering:
   - subset original neuronal clusters 0,1,2,3,5,8,11
   - reclustered as one neuron-only object
   - selected resolution 0.6, giving 9 neuron-only states

4. Science Data S9 / five-class module scoring:
   - scored the Excel full modules and TF-only modules
   - separately scored user-provided anchor classes:
     EPHA5/MEF2C, LHX6/NFIA, CRABP1/ANGPT2, NR2F1/NR2F2, LHX8/ISL1
   - used these to refine candidate branches and tips
```

The next all-cell DIV90 URD should use the same all-cell inclusion logic as the
previous DIV90 glia-tip run:

```text
Retain all non-stressed cells:
  original clusters 0,1,2,3,4,5,8,9,10,11,12

Exclude:
  original clusters 6,7 = stressed cells

Root:
  original cluster 12 = Dividing cells
  top 10% RootScore within cluster 12, matching the later DIV90 root convention
```

For neuronal cells, tips should now be assigned from the neuron-only
`neuron_leiden_r0_6` labels rather than from original parent clusters. Glial
tips should keep the same definitions as before.

The resulting 7-tip model is:

| Tip | Cells used as tip | Meaning |
|---|---|---|
| `tip_lhx8_isl1_state1` | neuron-only cluster 4 | LHX8/ISL1-like ventral MGE branch state 1 |
| `tip_lhx8_isl1_state2` | neuron-only cluster 7 | LHX8/ISL1-like ventral MGE branch state 2 |
| `tip_crabp1_angpt2_fetal_precursor` | neuron-only cluster 3 | CRABP1/ANGPT2-like fetal precursor branch |
| `tip_lhx6_nfia_epha5_mef2c_cortical` | neuron-only clusters 1 + 2 | LHX6/NFIA + EPHA5/MEF2C cortical interneuron branch |
| `tip_nr2f1_nr2f2` | neuron-only clusters 5 + 6 + 8 | NR2F1/NR2F2-like branch |
| `tip_astrocytes` | original clusters 4 + 10 | combined astrocyte endpoint, same as prior DIV90 v4 |
| `tip_opc` | original cluster 9 | OPC endpoint, same as prior DIV90 v4 |

Important: neuron-only cluster `0` should not be passed as a tip. It should be
retained as an upstream/intermediate state, likely upstream of the LHX8/ISL1-like
ventral MGE branch.

Machine-readable plan:

```text
metadata/div90_jia_lineage_urd_plan_neuron_s9_7tips.tsv
```

Code added for the all-cell 7-tip setup:

```text
python_notebooks/scripts/div90_assign_neuron_s9_7tip_groups.py
```

Exporter/template updates:

```text
python_notebooks/scripts/export_div90_jia_lineage_urd_inputs.py
slurm_templates/34_div90_jia_lineage_urd_smoke.sbatch.template
```

The template now supports:

```text
DIV90_TIP_MODE=neuron_s9_7tip
TREE_TIP_LABELS=tip_lhx8_isl1_state1,tip_lhx8_isl1_state2,tip_crabp1_angpt2_fetal_precursor,tip_lhx6_nfia_epha5_mef2c_cortical,tip_nr2f1_nr2f2,tip_astrocytes,tip_opc
```

Before changing official labels:

1. Open local UMAP audit plots for accepted parents `5`, `8`, and `9`.
2. Review `top_markers.tsv.gz` for whether markers are biologically coherent.
3. Treat parent `6` as stress/QC structure unless there is a specific reason to
   interpret stressed-cell subclusters biologically.
4. For clusters `1`, `3`, and `11`, inspect marker-panel UMAPs and local UMAP
   audit plots manually before deciding whether to loosen the covariate gates.

## Re-run Command

Submit the array:

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
cp slurm_templates/37_div90_parent_cluster_subclustering_array.sbatch.template \
  /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/37_div90_parent_cluster_subclustering_array.sbatch

RUN_LABEL=div90_parent_cluster_subclustering_audit_v2 \
sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/37_div90_parent_cluster_subclustering_array.sbatch
```

Merge outputs:

```bash
/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python \
  python_notebooks/scripts/merge_div90_parent_cluster_subclustering_audit.py \
  --run-dir /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_parent_cluster_subclustering_audit/div90_parent_cluster_subclustering_audit_v2 \
  --job-id 51778695
```

## Link To Resumable URD Rerun

The neuron-only S9 branch model from this handoff is now wired into the
production-style DIV90 all-cell URD rerun as seven tips:

```text
div90_allcells_jia_root10_neuron_s9_7tips_urd_resumable_v1
```

Full expected output directory:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_jia_lineage_urd/div90_allcells_jia_root10_neuron_s9_7tips_urd_resumable_v1
```

The run keeps the same seven-tip definition:

```text
tip_lhx8_isl1_state1 = neuron 4
tip_lhx8_isl1_state2 = neuron 7
tip_crabp1_angpt2_fetal_precursor = neuron 3
tip_lhx6_nfia_epha5_mef2c_cortical = neuron 1 + neuron 2
tip_nr2f1_nr2f2 = neuron 5 + neuron 6 + neuron 8
tip_astrocytes = original clusters 4 + 10
tip_opc = original cluster 9
```

Neuron `0` remains retained as an upstream/intermediate state and is not passed
as a tip. The shared URD runner is now resumable, but parameters stay aligned
with the previous all-cell runs: `MAX_CELLS=0`, `URD_KNN=100`,
`URD_N_FLOODS=20`, `URD_NUM_VARIABLE_GENES=3000`, `URD_PCA_MP_FACTOR=2`, and
`URD_SIGMA=local`.

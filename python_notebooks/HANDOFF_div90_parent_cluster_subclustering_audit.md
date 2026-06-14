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

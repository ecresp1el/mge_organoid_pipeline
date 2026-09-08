# PCDH19 expanded coarse-to-fine pilot: inputs, outputs and tuning

## Authorized scope and review boundary

The user requested a substantially larger per-sample subset, explicit coarse
clustering followed by fine clustering within those groups, and visible
DEGs/top-20 genes, heatmaps and dendrograms. This checkpoint selects **1,000 cells
per sample, 12,000 total**, ten times the previous pilot, across all 12 samples.
It retains all 19,071 genes. This is an expanded pilot, not the full 446,349-cell
analysis; more cells improve coverage but do not guarantee more or better clusters.

The two levels carry numerical IDs without biological annotation. Expression
from all samples is pooled. Genotype/sex/sample labels do not drive the fits.
No batch correction, integration, QC regression or additional cell filtering
is performed. Successful execution stops **IN_REVIEW**, without approving the
old pilot, changing preprocessing statuses, or claiming a locked full-data K.

## Inputs and representations

| Input | Exact role |
|---|---|
| Approved Step 02 run `02_qc_filtering_20260830_124611_97e1bb5` | Raw integer expression; SHA-256 `fadba4a25a7b6b7320219b21c189b6325687493519ba0fe1bd27efd79606b103`. |
| Balanced random subset, seed 20260908 | Exactly 1,000 cells/sample. Selected IDs/source rows are saved. Selection ignores clusters and genotype; this is a new subset, not guaranteed to contain the previous 1,200 cells. |
| Step 06 run `06_technical_sample_batch_diagnostics_20260903_114110_6087eeb` | Existing UMAP/PCA coordinates copied for display only after checking cell alignment. No Step 06 cluster or annotation enters fitting. |
| Frozen `config/hicat_pilot.json` | Snapshot of the editable `config/hicat_coarse_fine.json`; same frozen filename lets the existing input preparer be reused. |

**Normalization occurs here, once:** `ln(1 + raw_count / all-gene cell total × 1,000,000)`.
The identical normalized expression feeds both levels; normalization is not
refitted separately within groups. HVGs/PCA/graphs are refitted locally. The
main H5AD keeps raw sparse integer `.X` and exact sparse float64
`.layers['log1p_cpm']`. There is no corrected expression matrix. Dense working
expression is bounded to the configured subset and checked against a memory
budget before fitting; full-data memory scalability remains unqualified.

## Coarse → fine algorithm

1. Fit one coarse Allen feature/PCA/graph/DE-merge pass with stringent criteria.
   Stop its children at the explicit coarse boundary. Complete the coarse
   cross-branch merge and pairwise audit; save the resulting coarse assignments.
2. For each saved coarse group, subset those exact cells and run recursive fine
   clustering using relaxed separation criteria. Recompute local HVGs, PCA and
   graphs at each eligible branch. No other coarse parent's cells enter that fit.
3. Complete fine merging across recursive branches **within each coarse parent**.
   Fine IDs encode their parent, e.g. `C0002.F0003`. Save the actual cell-level
   mapping, group sizes, models, split/merge records and terminal reasons.
4. Audit **every final fine pair globally**, including pairs from different coarse
   parents. Report insufficient separation without silently merging across
   parents. This explicit nesting constraint differs from a freely merging global
   final partition. Any flagged cross-parent pair requires review before locking
   a hierarchy; the coarse boundary is not proof of final separation.
5. Repeat the entire coarse and fine sequence with seed 20260909 on the same
   selected cells. Coarse parents are refitted in the repeat. Compare partitions
   by cell membership (ARI/NMI and contingency tables), not by matching ID names.
   This is not bootstrap consensus or a resampling stability assessment.

The [user's tutorial](https://taxonomy.shinyapps.io/scrattch_tutorial/#section-coarse-level-clustering)
explicitly demonstrates stringent coarse clustering followed by iterative finer
splits with relaxed parameters. The
[Allen methods vignette](https://alleninstitute.github.io/scrattch.hicat/articles/scrattch.hicat_release.html)
also discusses marker heatmaps and expression dendrograms, as well as lower
foreground detection requirements for lower-depth assays. These guide the staged
design; the R numerical thresholds are not asserted equivalent to Python units.

The frozen Python source remains Allen commit
`99154957c74023235763025fda9dc4eb3ca952c6`. See the
[original adapter protocol](HICAT_CLUSTER_DISCOVERY_PROTOCOL.md) for source
edge cases and instrumentation. The adapter retains its guarded single-group
return, mutable merge memberships, explicit PCA rows, elbow filter, effective
neighbor count capped by retained PCs, and uncapped Python DE score. It does
not silently implement R behavior or remove dimensions associated with sample/QC.

## Tuning controls

Edit [config/hicat_coarse_fine.json](config/hicat_coarse_fine.json) for a future
versioned run. Do not edit a frozen run's configuration.

| Control | Coarse | Fine | Meaning |
|---|---:|---:|---|
| DE score minimum | 500 | 150 | Higher requires stronger aggregate evidence to keep groups separate. |
| Foreground detection q1 minimum | 0.5 | 0.3 | Fine permits evidence expressed in fewer foreground cells. |
| Normalized detection difference qdiff minimum | 0.7 | 0.5 | Fine permits less binary expression differences. |
| Mean natural-log expression difference minimum | 1 | 1 | Python's effect threshold; not a log2 fold change. |
| Pairwise adjusted p maximum | 0.05 | 0.05 | Allen eBayes/within-pair Holm filtering. |
| Minimum evidence genes | 5 | 5 | Independently checked in final pairwise audits. |
| Small-cluster merge threshold | 20 | 20 | Also used in upstream DE detection filtering. |
| Minimum cells to attempt recursion | Coarse boundary | 40 | Small fine branches have an explicit size-stop reason. |

Common controls: up to 3,000 HVGs; 50 fitted PCs, elbow selection up to 20;
Annoy 50 trees; requested graph k=15 (actual k recorded per node); Euclidean
neighbors, Jaccard weights, taynaud Louvain resolution=1.0; two merge neighbors;
depth guard=20 (hitting it fails rather than asserting convergence); one worker
for reproducibility. SLURM requests 4 CPUs, 64 GiB, four hours; the allocation is
a limit, not a runtime estimate. A single dense expression matrix is about
1.71 GiB; the preflight reserves ten such matrices within the 48 GiB working budget.

Fine relaxation is a declared hypothesis for this review checkpoint, not a
validated biological resolution. Three criteria change between coarse and fine;
this is a staged comparison, not a controlled one-parameter sensitivity experiment.
The aim is supported, stable subdivisions, not tuning to reach a desired K.

## DEGs and top-20 genes: two distinct evidence products

**HiCAT pairwise merge evidence** remains in each fitted node and endpoint audit:
pair scores, number of qualifying genes, up/down gene lists and memberships.
Those are the actual criteria used to retain or merge clusters.

**Inspection gene tables** are new and independent of merge decisions. For every
baseline coarse and fine cluster, save all 19,071 gene statistics, plus up to
20 qualifying positive genes. There are three explicit comparison contexts:

- `coarse_global`: each coarse group versus all other selected cells.
- `fine_global`: each fine group versus all other selected cells.
- `fine_within_Cxxxx`: each fine group versus its siblings inside that parent.

A parent with only one fine group has no sibling comparison; its absence is
recorded, not filled with a synthetic DEG result. All candidates retain HiCAT
pairwise evidence; detailed top-20 reporting and figures cover the baseline.

Inspection tests use two-sided Welch t-tests on natural log1p(CPM), then
Benjamini-Hochberg correction across **all genes in that one contrast**. Tables
retain gene ID/symbol, foreground/background cell counts, means, mean difference,
detection fractions (expression >0), t statistic, p-value and BH-adjusted p-value.
Positive genes require BH <=0.05, mean difference >=0.25, foreground detection
>=25%; rank by descending t, then mean difference, then gene ID. All these
reporting controls are in `marker_reporting`, separate from clustering controls.
Fewer than 20 qualifying genes are reported without padding. All-gene files are
compressed TSVs; top-20 tables and coverage files are readable TSVs.

These are exploratory cell-level markers of clusters inferred from the same
expression data. They are not independent confirmation, sample-replicated
WT/KO DE, or biological annotations. Reported natural-log mean differences must
not be renamed log2 fold changes. HiCAT uses eBayes; the inspection tables use
Welch/BH explicitly and must not be mistaken for HiCAT's merge tests.

## Output contract

Run root: `PAPER3_ROOT/results/hicat/01_coarse_fine_expanded_pilot/<RUN_ID>/`.

| Asset | Contents |
|---|---|
| `inputs/pilot_raw_counts.h5ad`, `pilot_cell_selection.tsv`, `input_identity.json` | Raw subset, exact source rows/barcodes/sample counts, source/input hashes. |
| `outputs/pcdh19_hicat_coarse_fine.h5ad` | Raw `.X`, normalized layer, coarse/fine assignments for both seeds, baseline coarse PCA/graph, existing display coordinates. |
| `outputs/anndata_uns_inventory.json` | Structured exact mirror of `.uns['hicat']`: identities, scope, normalization, no-correction flag, hierarchy/marker definitions, settings and counts. |
| `outputs/parent_child_mapping.tsv`, `cell_assignments.tsv.gz` | Baseline fine-to-coarse mapping with sizes; all cells with both seeds' assignments. |
| `outputs/cluster_count_summary.tsv`, `partition_agreement.tsv`, `*_seed_overlap.tsv` | Coarse/fine counts, unresolved pairs, seed ARI/NMI and membership overlap. |
| `outputs/baseline/`, `seed_repeat/` | Per-candidate hierarchy; `coarse/` and `fine_Cxxxx/` hold per-node features, PCA loadings/centers/projections, graph/index, cells, split/merge and eBayes evidence. Global fine audit saved at candidate root. |
| `outputs/top20_genes_all_contexts.tsv` | Combined baseline top genes with cluster, rank, comparison context and statistics. |
| `outputs/markers/<context>/*_all_genes.tsv.gz` | Complete all-gene inspection statistics per tested cluster. |
| `outputs/markers/<context>/top20_genes.tsv`, `marker_coverage.tsv`, `comparison.json` | Qualified ranked genes, tested/untestable status and exact methods/thresholds. |
| `outputs/markers/<context>/cluster_mean_log1p_cpm.tsv.gz`, `cluster_detection_fraction.tsv.gz` | Means and detection fractions for every gene and cluster. |
| `outputs/markers/<context>/expression_linkage.npy`, `linkage_input_labels.tsv`, `dendrogram_input_cluster_means.tsv` | Reusable expression dendrogram model/input; explicit absence JSON if one cluster or no marker genes. |
| `outputs/markers/<context>/*_heatmap_*.tsv`, `heatmap_cluster_order.tsv` | Exact plotted row z-scores, detection percentages and column order. |
| `outputs/figures/coarse_fine_review.pdf`, PNGs, `figure_index.tsv` | Paginated overview, parent-child hierarchy, coarse/fine/sibling expression dendrograms, top-20 mean/detection heatmaps for every baseline cluster, sample composition. |
| `outputs/validation_checks.tsv`, `output_manifest.tsv`, `software_versions.json` | Saved-data checks, per-file hashes/sizes and package versions. |
| `outputs/REVIEW_README.md`, `STEP_STATUS.json` | Completed package status, counts and review entry points. |
| `code/`, `config/`, `upstream/`, `provenance/`, `logs/` | Frozen code/reference/config, package freeze, authorization scope, hashes, scheduler logs and printed/file-backed progress (including top genes). |

H5AD `.uns` has only `hicat`; JSON strings inside it hold the full resolved
config and candidate summaries for H5AD compatibility. `.raw`, `.varm` and
`.varp` are absent/empty: raw counts are `.X`; fitted loadings are per-node files.
The H5AD alone does not contain all fitted models, marker statistics or figures.
Keep the full run package with it. No new UMAP, full-data assignments, bootstrap
consensus, reference mapping or single universal recursive model is produced.

## Reading the figures

Expression dendrograms use average linkage/Euclidean distance on cluster means
across the union of qualified top-20 genes for the stated comparison context.
They summarize similarity; branch heights are not confidence or merge thresholds.
The separate parent-child diagram shows the actual coarse/fine membership tree.
Heatmaps show per-gene z-scores of mean expression and adjacent detection
percentages. Every cluster gets its own ranked top-gene page; genes are labeled
with both symbols and IDs. The reused Step 06 UMAP is a common display and is
not used for clustering or for choosing K.

## Execution and validation

```bash
/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python \
  paper3_pcdh19/bin/submit_hicat_coarse_fine.py \
  --upstream /tmp/pcdh19_hicat_reference_dev
```

Submission snapshots the exact source and configuration before SLURM runs them.
Failures retain staging/logs; successful outputs publish by atomic rename.
Validation checks raw and normalized expression, identities, exact assignments,
parent nesting, sample counts, display coordinates and `.uns` inventory after
H5AD reopening. Synthetic tests exercise direct Welch equivalence, BH values,
invalid nesting and the full hierarchy/marker/plot export path. No previous
scientific package is overwritten.

## Function navigation

### hierarchy.py

| Definition | Responsibility |
|---|---|
| [`CoarseEngine`](scripts/hicat/hierarchy.py#L20) | Reuse the instrumented Allen pass but stop before fitting child branches. |
| [`CoarseEngine._visit`](scripts/hicat/hierarchy.py#L23) | Fit the root normally; record its children as coarse terminal groups. |
| [`HierarchyEngine`](scripts/hicat/hierarchy.py#L40) | Fit and save a paired coarse/fine hierarchy for one random seed. |
| [`HierarchyEngine.__init__`](scripts/hicat/hierarchy.py#L55) | Resolve candidate state without fitting or modifying input data. |
| [`HierarchyEngine._config`](scripts/hicat/hierarchy.py#L63) | Return an independent level config so relaxed settings cannot leak. |
| [`HierarchyEngine.run`](scripts/hicat/hierarchy.py#L69) | Return coarse/fine labels, parent table, and global final DE audit. |
| [`HierarchyEngine.validate_nesting`](scripts/hicat/hierarchy.py#L132) | Reject missing assignments or any fine cluster spanning two parents. |

### markers.py

| Definition | Responsibility |
|---|---|
| [`MarkerReporter`](scripts/hicat/markers.py#L17) | Rank positive separation genes while retaining all tested gene statistics. |
| [`MarkerReporter.__init__`](scripts/hicat/markers.py#L20) | Store explicit marker filters, separate from clustering thresholds. |
| [`MarkerReporter.bh`](scripts/hicat/markers.py#L28) | Benjamini-Hochberg adjusted p-values in original gene order. |
| [`MarkerReporter.moments`](scripts/hicat/markers.py#L37) | Return n, sums, squared sums and nonzero counts without dense copies. |
| [`MarkerReporter.stats`](scripts/hicat/markers.py#L45) | Convert sufficient statistics to means, unbiased variance, detection. |
| [`MarkerReporter.run`](scripts/hicat/markers.py#L52) | Save every contrast and return means, detection fractions and top genes. |

### hierarchy_report.py

| Definition | Responsibility |
|---|---|
| [`HierarchyReport`](scripts/hicat/hierarchy_report.py#L18) | Build a readable review book without silently truncating cluster coverage. |
| [`HierarchyReport.__init__`](scripts/hicat/hierarchy_report.py#L21) | Open the report and prepare a page manifest for every exported figure. |
| [`HierarchyReport.save`](scripts/hicat/hierarchy_report.py#L29) | Save one PDF page and a matching PNG, recording its interpretation. |
| [`HierarchyReport.overview`](scripts/hicat/hierarchy_report.py#L38) | Show both levels on fixed coordinates and their exact parent-child counts. |
| [`HierarchyReport.expression`](scripts/hicat/hierarchy_report.py#L81) | Save a mean-expression dendrogram and one top-20 page per cluster. |
| [`HierarchyReport.composition`](scripts/hicat/hierarchy_report.py#L139) | Save counts/fractions and show sample contributions at either level. |
| [`HierarchyReport.close`](scripts/hicat/hierarchy_report.py#L154) | Finalize the PDF and write a complete PNG/page description index. |

### hierarchy_workflow.py

| Definition | Responsibility |
|---|---|
| [`HierarchyWorkflow`](scripts/hicat/hierarchy_workflow.py#L29) | Own a versioned larger pilot from verified raw input through review assets. |
| [`HierarchyWorkflow.__init__`](scripts/hicat/hierarchy_workflow.py#L32) | Load the exact frozen settings and create a new staging directory. |
| [`HierarchyWorkflow.run`](scripts/hicat/hierarchy_workflow.py#L40) | Fit hierarchies, export complete marker/plot assets, and validate I/O. |
| [`main`](scripts/hicat/hierarchy_workflow.py#L198) | Run only the specified frozen expanded-pilot package. |


## Additional location review

The [completed handoff](HICAT_COARSE_FINE_HANDOFF.md) links an additional five-page
PDF highlighting each fine cluster on the fixed Step 06 UMAP. It reads the
completed object without recomputing clusters or coordinates. Rendering source,
coordinates, page/cluster mapping, input hash and file manifest are stored in
the run's separate `review/` directory. The original 89-page scientific report
and output manifest are unchanged. The helper is
[`location_report.py`](scripts/hicat/location_report.py); `render_locations`
accepts an input H5AD and a new output directory.

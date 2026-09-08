# PCDH19 expanded coarse-to-fine pilot — IN_REVIEW

Completed run `hicat_coarse_fine_20260908_203708_7e050f23`, SLURM job **60563267**. This is the user's requested
larger coarse/fine checkpoint: **12,000 cells, 1,000 per sample, all 19,071 genes**.
No full-data count is locked; no batch correction or biological annotation was
performed. Successful computation does not approve the scientific resolution.

## Start here: plots and top genes

- [Coarse/fine review PDF](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/01_coarse_fine_expanded_pilot/hicat_coarse_fine_20260908_203708_7e050f23/outputs/figures/coarse_fine_review.pdf) — **89 pages**,
  including coarse/fine UMAP displays, parent-child mapping diagrams, expression
  dendrograms, top-20 mean-expression/detection heatmaps and sample composition.
- [Figure index](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/01_coarse_fine_expanded_pilot/hicat_coarse_fine_20260908_203708_7e050f23/outputs/figures/figure_index.tsv) — page descriptions and PNG filenames.
- [All ranked top-20 gene tables](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/01_coarse_fine_expanded_pilot/hicat_coarse_fine_20260908_203708_7e050f23/outputs/top20_genes_all_contexts.tsv) — cluster ID,
  comparison context, rank, gene ID/symbol, effect, detection and test statistics.
- [Coarse top genes](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/01_coarse_fine_expanded_pilot/hicat_coarse_fine_20260908_203708_7e050f23/outputs/markers/coarse_global/top20_genes.tsv) and
  [fine top genes versus all other cells](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/01_coarse_fine_expanded_pilot/hicat_coarse_fine_20260908_203708_7e050f23/outputs/markers/fine_global/top20_genes.tsv).
- [All marker contexts and complete per-gene statistics](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/01_coarse_fine_expanded_pilot/hicat_coarse_fine_20260908_203708_7e050f23/outputs/markers) —
  `fine_within_Cxxxx/top20_genes.tsv` gives each parent's sibling comparisons.
- [Main AnnData](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/01_coarse_fine_expanded_pilot/hicat_coarse_fine_20260908_203708_7e050f23/outputs/pcdh19_hicat_coarse_fine.h5ad) and
  [exact structured `.uns` inventory](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/01_coarse_fine_expanded_pilot/hicat_coarse_fine_20260908_203708_7e050f23/outputs/anndata_uns_inventory.json).

The global fine UMAP has crowded labels. Use the separate
[individual fine-cluster location PDF](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/01_coarse_fine_expanded_pilot/hicat_coarse_fine_20260908_203708_7e050f23/review/fine_cluster_locations/fine_cluster_locations.pdf)
for five pages showing all 38 baseline groups one at a time, with exact membership
and coordinate tables. This review supplement reads the completed H5AD; it does
not recluster or change the original output package. Its rendering source,
input hash and separate file manifest are saved in `review/`.

## Observed counts and stability

| candidate | coarse_clusters | fine_clusters | coarse_unresolved_pairs | fine_unresolved_pairs |
|---|---|---|---|---|
| baseline | 4 | 38 | 1 | 2 |
| seed_repeat | 5 | 39 | 0 | 4 |

| level | adjusted_rand | normalized_mutual_information |
|---|---|---|
| coarse | 0.6759 | 0.555 |
| fine | 0.3358 | 0.5989 |

The baseline parent-child counts are:

| coarse_cluster | cells | fine_clusters |
|---|---|---|
| C0001 | 3854 | 13 |
| C0002 | 7295 | 21 |
| C0003 | 824 | 3 |
| C0004 | 27 | 1 |

**Fine-cluster membership is not stable enough to lock K:** adjusted Rand index
0.336, despite similar 38/39 counts. Baseline has two unresolved fine pairs;
the repeat has four. Baseline coarse groups also retain one pair with only
three qualifying genes against the declared minimum of five.

Counts changed with a tenfold larger subset **and** an explicit relaxation of fine
thresholds. This cannot isolate a cell-count effect or prove the new count is
better. Repeat-seed agreement is measured on these same cells, not bootstrap
resampling. Every fine cluster remains nested in its coarse parent. The global
fine audit includes cross-parent pairs but does not silently merge them and erase
that hierarchy. Unresolved pairs require review before locking any resolution.

## What the gene evidence means

The ranked inspection tables compare each cluster with the rest using two-sided
Welch tests, BH adjustment across all genes per contrast, and positive-expression
and detection filters. They are exploratory cell-level markers of groups defined
from the same data, not independent confirmation or genotype DE with biological
replicates. Fine-sibling comparisons restrict the background to the same parent.
All 19,071 genes are retained in each testable all-gene file, and at most 20
qualifying positive genes are displayed per cluster/context. No padding is used.

HiCAT's actual eBayes merge evidence is preserved separately under
[baseline](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/01_coarse_fine_expanded_pilot/hicat_coarse_fine_20260908_203708_7e050f23/outputs/baseline) and [seed repeat](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/01_coarse_fine_expanded_pilot/hicat_coarse_fine_20260908_203708_7e050f23/outputs/seed_repeat), including all
visited nodes and final pairwise audits. These inspection Welch/BH tables do not
replace or change the HiCAT merge rules.

Expression dendrograms summarize mean-expression similarity on the union of
reported marker genes. Their branch heights are distances, not statistical
support. The separate coarse-to-fine diagrams show the actual membership
hierarchy. Every heatmap matrix, dendrogram input, linkage and cluster ordering
is saved alongside the corresponding marker context.

## Explicit marker-coverage exceptions

- `fine_global` / `C0001.F0007`: 0 reported genes (tested).
- `fine_global` / `C0001.F0011`: 0 reported genes (tested).
- `fine_global` / `C0001.F0012`: 18 reported genes (tested).
- `fine_within_C0001` / `C0001.F0011`: 0 reported genes (tested).
- `fine_within_C0003` / `C0003.F0002`: 4 reported genes (tested).
- `fine_within_C0004` / `C0004.F0001`: 0 reported genes (not_testable_no_rest_or_insufficient_cells).

A missing positive-marker list is retained as an inspection result; all testable
contrasts still have complete all-gene statistics. No gene list was padded.

## Input/output boundary

Input is a new balanced subset from the exact approved Step 02 raw counts.
Step 06 coordinates are copied only for display. Normalization is performed
once as `ln(1 + raw_count / all-gene cell total × 1,000,000)` and shared by both
levels. Features/PCA/graphs are recomputed locally inside each group.

The H5AD saves raw sparse integer `.X`, exact normalized `.layers['log1p_cpm']`,
coarse/fine labels for both seeds, baseline coarse PCA/graph, inherited display
coordinates and `.uns['hicat']`. Models/graphs for each fine branch, per-gene
statistics, dendrogram linkages, figures and full logs live outside the H5AD.
Keep the [complete run directory](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/01_coarse_fine_expanded_pilot/hicat_coarse_fine_20260908_203708_7e050f23) with the object.

Not produced: full-data assignments/K, biological names, batch-corrected
expression or coordinates, new UMAP, bootstrap consensus, reference mapping,
new-cell classifier or one universal recursive PCA. `.raw` is absent because
raw counts are already in `.X`.

## Validation and reproducibility

All **18 saved-data checks passed**. Independent verification checked
**2746 output manifest entries**, exact sample sizes and parent nesting,
structured inventory equality, all-gene row counts, and top-gene agreement with
the complete tables. Eight pre-submission tests passed and also reproduced against frozen code using
the saved `provenance/test_bundle/` fixtures, including a real
synthetic end-to-end fit/export and direct Welch/BH checks. The source, upstream
commit, exact config, package freeze, input identities, logs and tests used are
preserved with this run. Figures were inspected separately for readability.

[Method, tuning, assets and function guide](HICAT_COARSE_FINE_PROTOCOL.md) explains
the coarse/fine thresholds and separate marker-reporting controls. Future changes
must produce a new versioned package. The
[HiCAT ledger](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/APPROVAL_LEDGER.tsv)
records this checkpoint IN_REVIEW without changing prior approvals.

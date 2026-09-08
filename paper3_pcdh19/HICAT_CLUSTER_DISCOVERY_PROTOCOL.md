# PCDH19 HiCAT cluster discovery: method, assets and review contract

## Objective and current authorized checkpoint

Determine a defensible number of expression-defined clusters in the pooled
12-sample PCDH19 dataset. Cluster IDs are arbitrary identifiers, **not biological
annotations**. There is no target K and no label-transfer/reference-annotation
step. A cluster count is conditional on separation criteria and data coverage;
a visually appealing UMAP or one successful run does not lock it in.

The user explicitly selected **“Prepare method/assets and run a small technical
pilot”** on 2026-09-08. This checkpoint therefore comprises the method contract,
source/environment qualification, a bounded real-data pilot, exact saved assets
and review plots. It stops **IN_REVIEW**. Full-data clustering, a consensus
partition, and locking the final K are later review checkpoints, not implied by
a successful pilot. The earlier major-step review model remains in effect.

Pilot expression input: the exact **approved Step 02** raw-count file, SHA-256
`fadba4a25a7b6b7320219b21c189b6325687493519ba0fe1bd27efd79606b103`.
All 446,349 source cells and 19,071 genes remain in that immutable checkpoint.
Step 06 supplies fixed display PCA/UMAP coordinates only. Its provisional
biological labels, Leiden partition and HVGs do not enter HiCAT fitting.
Rejected scDblFinder fields and skipped ambient/contaminant results are not used.
This pilot does not alter Step 06's existing review status or endorse its headline.

## Sources and what is actually being used

- User-supplied [Python repository, dev branch](https://github.com/AllenInstitute/transcriptomic_clustering/tree/dev).
  The inspected/frozen commit is **`99154957c74023235763025fda9dc4eb3ca952c6`**,
  package version `1.0.0`. The moving branch name is not the reproducibility pin.
- User-supplied [interactive scrattch tutorial](https://taxonomy.shinyapps.io/scrattch_tutorial/#section-overview).
  Its HTML was fetched successfully, including the overview. The tutorial explains
  iterative subdivision and DE-based stopping, with resampling/consensus as a
  separate robustness procedure.
- [Allen's HiCAT methods vignette](https://alleninstitute.github.io/scrattch.hicat/articles/scrattch.hicat_release.html)
  provides the R-method context. It is not a guarantee of numerical equivalence
  with the Python implementation. In particular, the R normalization uses log2,
  while the inspected Python normalizer uses natural log.

This is **Python HiCAT primitive orchestration with an explicitly documented
adapter**, not a claim that unmodified R HiCAT or the Python convenience script
was executed verbatim. The core mean/variance, feature selection, PCA, PC
filtering, projection, Jaccard graph/community clustering, empirical-Bayes DE,
and merge functions come from the pinned Allen source. Our classes add explicit
stage control, bug-path guards, persistence, progress, validation and plotting.
Upstream source is frozen unchanged in each pilot package.

## Algorithm, in order

1. Select exactly **100 cells per registered technical sample**, without
   replacement, using seed `20260908`: **1,200 cells**, all 19,071 genes. This
   balances pilot coverage and is not a representative estimate of rare-cluster
   frequency or the full dataset's K. Save every selected barcode and source row.
2. Compute **`ln(1 + 1,000,000 × raw_count / all-gene cell total)`**. Save this
   expression as `.layers['log1p_cpm']`, while preserving raw integer `.X`.
   Do not reuse Step 06's `ln(1 + counts per 10,000)` representation.
3. Within each eligible branch, use the Allen mean/variance and HVG procedures,
   then PCA and its elbow-based component selection. Recompute locally at every
   subdivision; no predefined marker panels determine groups.
4. Build the Allen approximate-neighbor/Jaccard graph and run its `taynaud`
   Louvain backend. This backend is chosen because its resolution argument is
   explicitly wired through in the inspected source. It is not Step 06 Leiden.
5. Merge small communities and communities with insufficient expression
   separation using the Allen empirical-Bayes DE merge. Save the tested pairs,
   scores, evidence genes and local memberships at each evaluation.
6. Recurse into surviving groups with at least 40 cells. A sole merged group
   terminates normally; a group below 40 stops by the stated size rule. Hitting
   the depth guard or encountering invalid/zero-dimensional input is a **failure**,
   never quietly labeled convergence or “one cluster.”
7. Across terminal branches, fit a final PCA using the union of data-derived
   separation genes and **all pilot cells**, then perform a cross-branch DE
   merge. Save counts before and after this step; recursive leaves alone are
   not the final reported candidate K.
8. Audit every final cluster pair for the declared score and minimum-gene
   criteria. Weak pairs are flagged in a table rather than silently changing
   the partition. A candidate with unresolved pairs is not ready for locking K.
9. Repeat the same selected cells with a second seed and a stricter merge score.
   Save all three partitions and compare them. This is a small sensitivity
   check, **not bootstrapped consensus clustering**.

No biological names are generated. Data-derived genes are saved/displayed only
to explain why a split is supported; gene evidence is essential to HiCAT even
when cell-type annotation is outside scope. DE here is algorithm-internal cell
separation evidence, not sample-replicated genotype differential expression.
Sample, genotype and sex metadata are available for review; only expression
enters clustering. No batch integration, QC-mode regression or cell deletion
is introduced.

## Pilot controls and where to tune them

The complete machine-readable configuration is
[`config/hicat_pilot.json`](config/hicat_pilot.json). It is frozen per run.
Identity fields and the pilot cap are guards, not knobs to bypass before review.

| Setting | Pilot value | Effect / interpretation |
| --- | --- | --- |
| Normalization | CPM target 1,000,000; natural log1p | Must match the Python input representation and threshold units. |
| Gene eligibility | expression >1 in at least 4 cells | Evaluated within each branch before HVG selection. |
| HVGs | Up to 3,000 | Local features for each branch; source may return fewer. |
| PCA | Up to 50 fitted PCs; all pilot cells | Centered PCA on log-CPM, without Step 06's additional SD scaling/clipping. |
| PC filter | Allen elbow; at most 20 retained | Recomputed locally; zero retained components fail. |
| Graph | k=15, Euclidean, Annoy 50 trees, Jaccard weights | Effective k also follows the pinned one-step convention of capping by retained PC dimensions; each node records requested/effective k. |
| Community method | `taynaud` Louvain; resolution 1.0 | Creates candidate groups before DE merge, not the final K by itself. |
| q1 | >0.5 | A qualifying gene must be detected above low_thresh in more than half the foreground cells. |
| q2 | None | No additional background detection ceiling. |
| qdiff | >0.7 | Normalized detection difference: abs(q1-q2)/max(q1,q2). |
| p-adjusted | <0.05 | Python eBayes code uses Holm correction within each pair. |
| lfc | absolute value >1.0 | Difference of mean **natural-log normalized expression** in this implementation; do not label it log2 fold change. |
| DE score | 150 baseline; 300 sensitivity | Upstream sum of -log10(adjusted p), without the R tutorial's per-gene cap. Higher threshold generally makes splits harder to retain; monotonic K is not guaranteed by the entire recursive procedure. |
| Minimum genes | 5 | Declared final separation audit; upstream merge control flow has a caveat below. |
| Cluster-size threshold | 20 | Both the small-cluster merge threshold and an expression-detection cell-count condition inside DE filtering; not a guarantee every final group must have 20 cells. |
| Recursive eligibility | at least 40 cells | Avoid attempting meaningful subdivision of a branch too small to support two 20-cell groups. |
| Merge neighbors | 2 | Nearby cluster comparisons during iterative merging; final all-pair audit detects any untested weak pairs. |
| Seed | 20260908, plus 20260909 repeat | Same cell subset for comparisons; randomness changes clustering, not sample selection. |
| Workers | 1 | Keeps small pilot resource use/reproducibility controlled. Full-data performance is not established by this pilot. |
| Depth guard | 20 | Safety guard raises; it is not a legitimate scientific stopping rule. |

These are an explicitly specified **pilot baseline**, not validated PCDH19
thresholds. Native Python threshold values are retained except for the stated
20-cell minimum and 40-cell recursion eligibility. Lower-depth/Flex-specific
parameter choices require full-data sensitivity review; R numerical thresholds
must not be copied without accounting for log base and implementation differences.

## Source qualification and explicit adapter behavior

Source inspection and runtime probes found issues relevant to a reproducible run:

| Finding at pinned commit | Our bounded pilot handling |
| --- | --- |
| `merge_clusters()` returns a dict rather than the documented tuple when only one community exists. | Detect one community and return that exact unchanged membership before invoking the broken convenience return path. |
| `merge_two_clusters()` requires mutable list memberships, although orchestration may carry NumPy index arrays. | Convert memberships to native integer lists before invoking upstream merging; an identical-distribution merge test exercises this path. |
| `filter_components(method=None)` actually returns zero PCs. | Use the explicit elbow filter; test and document this defect. Never interpret an empty projection as one biological cluster. |
| `pca(cell_select=None)` can attempt `len(slice)`. | Supply explicit selected row indices at every PCA call. |
| `onestep_clust()` caps graph k by the number of retained dimensions. | Preserve that convention and record effective k per node; do not silently correct a potentially influential rule. |
| The `vtraag` branch does not forward the caller's resolution argument. | Use the explicit `taynaud` branch for this pilot. |
| R tutorial score caps do not appear in Python `calc_de_score`. | Preserve Python scores and clearly label their definition; possible infinite scores are retained in tables and only color scales are capped. |
| The merge loop can break on score before enforcing minimum-gene criteria. | Preserve upstream merging, then independently audit **all final pairs** for score AND minimum-gene evidence. Flag unresolved pairs; do not claim full separation automatically. |
| Convenience final-merge code samples cells and uses implicit masks. | Explicit all-pilot-cell marker-union PCA plus the same Allen merge primitives; no sampling mask ambiguity. This is a declared orchestration difference. |
| Old AnnData/Leiden pins and private Scanpy APIs differ from the preprocessing environment. | Separate Python 3.8 environment; exact installed package freeze; formula/import/I/O checks. Preprocessing environment remains unchanged. |

The adapter does not fabricate partitions if upstream numerical code fails.
The pilot's **dense** normalized working matrix is bounded to 1,200 cells; a
successful pilot does not qualify memory/time or sparse/backed execution for
446,349 cells. Full-data preparation must validate those execution paths.

## What will be saved at this checkpoint

All paths are relative to a new run under
`PAPER3_ROOT/results/hicat/00_method_and_technical_pilot/<RUN_ID>/`.
Run preparation freezes source/config first; SLURM executes those copies.
`outputs/` appears only after validation and one atomic directory rename.
If execution fails, `staging/`, logs and failure evidence remain for diagnosis.

| Asset | Contents and purpose |
| --- | --- |
| `inputs/pilot_raw_counts.h5ad` | Minimal raw pilot input, sample/QC metadata, full gene set, fixed Step 06 display coordinates. |
| `inputs/pilot_cell_selection.tsv` | Cell ID, exact source row, registered sample and selection seed. |
| `inputs/input_identity.json` | Verified Step 02/Step 06 file hashes and pilot input hash/dimensions. |
| `outputs/pcdh19_hicat_pilot.h5ad` | Main reusable pilot object: raw `.X`; saved `.layers['log1p_cpm']`; baseline/repeat/stricter final and recursive cluster columns; root PCA/graph; fixed display coordinates; explicit `.uns['hicat']`. |
| `outputs/anndata_uns_inventory.json` | Exact human-readable metadata dictionary saved with the H5AD; no stale latest_step fields. |
| `outputs/cell_assignments.tsv.gz` | Every pilot cell's metadata and all candidate cluster assignments, easy to join elsewhere. |
| `outputs/cluster_count_summary.tsv` | For each candidate: root communities, post-merge groups, recursive leaves, final K, node/depth counts and unresolved final pairs. |
| `outputs/partition_agreement.tsv`, `cluster_stability.tsv`, `baseline_seed_overlap.tsv` | ARI/NMI, per-baseline-cluster best Jaccard matches, and exact membership overlap. Numeric cluster IDs themselves are not compared across runs. |
| `outputs/sample_cluster_*.tsv` | Absolute sample contributions, within-sample fractions and within-cluster fractions. |
| `outputs/<candidate>/recursion_tree.tsv` | Every visited node, its parent/depth/cell count, graph/merged K, and explicit terminal reason. |
| `outputs/<candidate>/<node>/cell_membership.tsv` | Local-to-pilot cell identity mapping; essential to decode node-local labels. |
| `outputs/<candidate>/<node>/feature_selection.tsv`, `eligible_gene_statistics.tsv` | Actual feature mask and eligible linear-CPM means/variances. Final merge uses the marker union rather than fresh HVG eligibility statistics. |
| `outputs/<candidate>/<node>/pca_components_all.tsv`, `pca_components_retained.tsv`, `pca_center.tsv`, `pca_variance.tsv`, `projection.npz` | Reusable per-node PCA loadings, centering, variances and cell projections. There is no single PCA model for the entire recursive hierarchy. |
| `outputs/<candidate>/<node>/graph.npz`, `annoy.index`, `graph_membership.tsv`, `merged_membership.json` | Exact pre/post-merge memberships and the graph/index for each recursively fitted node. Final cross-branch merge does not build a new cell graph. |
| `outputs/<candidate>/<node>/de_evaluation_*.tsv` and matching membership JSON | Each tested pair's score/gene evidence plus the cell memberships in effect at that evaluation. |
| `outputs/<candidate>/final_pairwise_evidence.tsv` | All final cluster pairs, evidence genes/scores, and whether declared score/gene criteria are met. |
| `outputs/<candidate>/cluster_mean_log1p_cpm.tsv.gz`, `cluster_detection_fraction.tsv.gz`, `marker_union.tsv` | Cluster expression/detection summaries and discovered evidence-gene union. |
| `outputs/baseline_cluster_mean_linkage.npy`, `baseline_linkage_labels.tsv` | Numerical dendrogram assets when at least two clusters exist; absent with a documented one-cluster case. |
| `outputs/qualification_checks.tsv`, `validation_checks.tsv`, `software_versions.json`, `output_manifest.tsv` | Runtime/source checks, exact partition/raw-count/normalization round-trip checks, versions and content hashes. |
| `outputs/figures/hicat_pilot_review.pdf` plus three PNG previews | The coherent A–L report described below. |
| `outputs/PILOT_REPORT.md`, `outputs/STEP_STATUS.json` | Findings, limitations, assets and explicit `IN_REVIEW` boundary. |
| `code/`, `upstream/`, `config/`, `provenance/` | Frozen adapter and upstream source, exact settings/dependency freeze, source-file hashes, original repository identity/diff, authorization record and live events. |

### H5AD schema and what is NOT saved

- `.X`: raw sparse integer counts for the **pilot cells only**.
- `.layers['log1p_cpm']`: exact normalized expression used for fitting, saved as
  sparse float64. No separate `.raw` copy is needed.
- `.obs['hicat_baseline']`, `hicat_seed_repeat`, `hicat_stricter_merge`: final
  candidate IDs; matching `*_recursive` columns preserve pre-global-merge leaves.
- `.obsm['X_pca_hicat_root']`: baseline root PCA projection. Per-node models
  and projections are separate files.
- `.obsm['X_umap_step06_display']`, `X_pca_step06_display`: selected rows of
  the already-existing unintegrated Step 06 coordinates, explicitly display-only.
- `.obsp['hicat_root_connectivities']`: baseline root graph; recursive graphs
  are per-node assets. This is not a single final HiCAT graph for all depths.
- `.uns['hicat']`: schema/stage/status/output-run ID, scope, exact upstream/input
  identity, dimensions, count/annotation flags, normalization/slot meanings,
  complete resolved-config JSON, candidate summaries JSON, adapter/model/stability
  descriptions. `anndata_uns_inventory.json` mirrors it exactly.

Not produced: full-data cluster assignments or full-data K; cell-type names;
reference mapping; bootstrap consensus; an N×N full-data co-clustering matrix;
a new HiCAT UMAP; normalized/scaled duplicates per recursive node; regression
or integration models; a single universal fitted PCA; a classifier for assigning
new cells. Per-gene p-values/effect sizes for every tested pair are not retained
by the upstream summary API—the pilot saves its returned pair scores and gene
lists, not unreturned statistics. Figures are not evidence of assets absent
from this inventory. Keep the full run directory with the main H5AD.

## What the review plots will look like

One three-page A–L report, all from the observed pilot data:

| Page / panel | Display | What it answers |
| --- | --- | --- |
| 1A | Lines tracing K from root communities → root merge → recursive leaves → final merge, for each candidate | Where did the cluster count change? |
| 1B | Final K bars for baseline, repeat seed and stricter merge score | Is the pilot count sensitive to those choices? |
| 1C | Fixed Step 06 UMAP with numeric HiCAT IDs | Where do the discovered groups lie on a common display? No UMAP islands are counted as clusters. |
| 1D | Identical coordinates/cells colored by sample | Are groups plausibly aligned with sample identity? This is descriptive, not automatic correction. |
| 2E | Cluster-size bars | Are there tiny groups or very uneven partitions? |
| 2F | Sample × cluster cell-count heatmap | Which samples contribute to each group? Fractions are saved alongside counts. |
| 2G | Data-derived separation-gene × cluster expression heatmap | Is there readable gene evidence for the groups, without naming cell types? |
| 2H | All-final-pair DE-score heatmap | Are weakly separated pairs left after merging? Plot color is capped; exact scores and pass/fail live in tables. |
| 3I | Actual visited split tree, nodes labeled with n and post-merge K | Which branches split, and at what depth? Not a taxonomy. |
| 3J | Baseline vs repeat-seed contingency heatmap | Are the same cells grouped together when IDs change? |
| 3K | Best-match Jaccard overlap per baseline cluster | Which groups are stable or fragmented across seeds? |
| 3L | Dendrogram of cluster-mean expression | Which groups are close in expression? This dendrogram is distinct from the recorded split tree. |

There is no invented mock result. One-cluster cases explicitly show that pair
comparisons/dendrograms are unavailable. PNG previews and the PDF are inspected
before reporting the pilot as complete.

## Prints, logs, saves and review gates

Each operation prints and fsyncs structured START/COMPLETE/FAILED events with
candidate/node identity, exact settings, dimensions, K before/after, output
references, duration and lifetime peak RSS. The same events are stored in
`provenance/hicat_progress_events.jsonl`; the latest event is atomically replaced
in `hicat_progress_latest.json`. Allen's detailed messages are retained in
`logs/upstream.log` and scheduler stdout/stderr. Failure is recorded and raised;
there is no silent fallback to a successful partition.

Validation must establish: exact selected-cell membership, no duplicated/lost
cells in any candidate, raw count preservation, saved normalized-layer equality,
H5AD read-back, metadata/coordinate alignment, all 12 samples, and correctly
labeled pilot scope. Scientific sensitivity and unresolved pairs are reported
separately; a technically valid file does not validate a final K.

Later, after pilot review, prepare/freeze the full input and qualify the
sparse/backed runtime. A full baseline should be followed by reviewed seed,
separation-threshold and sample/resampling checks. Use partition overlap,
per-cluster stability, all-pair separation evidence, sample composition and
explicit stopping records to choose a supported partition. A stable range of K
may be more honest than forcing one number. No full-data count is locked by
this pilot or by an automatic default rule.

## Source-code navigation

The object-oriented pilot code is under [`scripts/hicat/`](scripts/hicat/):

- `prepare.py`: approved-source identity, label-independent sampling, input asset.
- `engine.py`: Allen primitive adapter, per-node feature/PCA/graph/merge persistence,
  recursion, final cross-branch merge and pairwise endpoint audit.
- `report.py`: panels A–L using saved observed values.
- `provenance.py`: immutable input hashes and durable printed/file events.
- `workflow.py`: qualification, normalization, candidate comparisons, H5AD assembly,
  validation and atomic publication.
- `cli.py`: thin explicit-run entry point.

The pilot submitter freezes exact files even when the repository has pending
documentation edits: it records the Git identity, working diff, and content
hashes. It does not pretend a commit alone identifies a dirty working tree.
All scientific execution uses the frozen copies. No historical run is replaced.

## Function index

Each link opens the relevant definition; module docstrings explain the overall role.

### cli.py

| Definition | Responsibility |
|---|---|
| [`main`](scripts/hicat/cli.py#L7) | Parse the explicit run directory, run the pilot, and print its asset path. |

### engine.py

| Definition | Responsibility |
|---|---|
| [`PilotEngine`](scripts/hicat/engine.py#L20) | Own one candidate run's recursive clustering and saved node artifacts. |
| [`PilotEngine.__init__`](scripts/hicat/engine.py#L23) | Store frozen controls and initialize trace collections; no fit yet. |
| [`PilotEngine._table`](scripts/hicat/engine.py#L41) | Save a numerical TSV without silently dropping its identifying index. |
| [`PilotEngine._capture_de`](scripts/hicat/engine.py#L45) | Call upstream DE unchanged and persist each evaluated pair's evidence. |
| [`PilotEngine._representation`](scripts/hicat/engine.py#L70) | Fit/persist local HVGs, PCA model and projected cells using Allen code. |
| [`PilotEngine._merge`](scripts/hicat/engine.py#L128) | Merge using unchanged Allen criteria, handling its single-group API bug. |
| [`PilotEngine._visit`](scripts/hicat/engine.py#L182) | Visit one recursive branch and record why it splits or terminates. |
| [`PilotEngine.labels`](scripts/hicat/engine.py#L239) | Return deterministic C0001-style labels after checking complete partition. |
| [`PilotEngine.run`](scripts/hicat/engine.py#L254) | Run recursion, cross-branch merge and all-pair endpoint evidence. |

### prepare.py

| Definition | Responsibility |
|---|---|
| [`PilotInputPreparer`](scripts/hicat/prepare.py#L17) | Create a small raw-count input while binding it to approved source bytes. |
| [`PilotInputPreparer.__init__`](scripts/hicat/prepare.py#L20) | Load frozen config and resolve exact Step 02 and Step 06 paths. |
| [`PilotInputPreparer.run`](scripts/hicat/prepare.py#L26) | Validate source identity, sample rows, and write a minimal pilot H5AD. |
| [`main`](scripts/hicat/prepare.py#L100) | Parse run/project paths and prepare only the explicitly bounded pilot. |

### provenance.py

| Definition | Responsibility |
|---|---|
| [`json_value`](scripts/hicat/provenance.py#L19) | Convert Paths/NumPy scalars/arrays for provenance JSON, not expression data. |
| [`write_json`](scripts/hicat/provenance.py#L34) | Atomically write a JSON metadata file at path; replace only that file. |
| [`sha256`](scripts/hicat/provenance.py#L42) | Return the streaming SHA-256 of a file without changing its contents. |
| [`manifest`](scripts/hicat/provenance.py#L51) | Return relative paths, byte sizes and hashes of files below root. |
| [`Progress`](scripts/hicat/provenance.py#L65) | Record one execution's operations and their resolved inputs/results. |
| [`Progress.__init__`](scripts/hicat/provenance.py#L68) | Create a new single-writer stream; refuse to mix with an older stream. |
| [`Progress.note`](scripts/hicat/provenance.py#L78) | Print and fsync one metadata event; RSS is lifetime peak GiB on Linux. |
| [`Progress.track`](scripts/hicat/provenance.py#L93) | Yield an output-metadata dict, then log completion or re-raise failure. |

### report.py

| Definition | Responsibility |
|---|---|
| [`PilotReport`](scripts/hicat/report.py#L16) | Render a reviewable report with companion numerical tables. |
| [`PilotReport.__init__`](scripts/hicat/report.py#L19) | Create the figure directory and store display-only configuration. |
| [`PilotReport.scatter`](scripts/hicat/report.py#L26) | Draw categorical labels on fixed coordinates; no embedding is fitted. |
| [`PilotReport.heatmap`](scripts/hicat/report.py#L37) | Display a labeled numerical matrix; full values remain in saved TSVs. |
| [`PilotReport.publish`](scripts/hicat/report.py#L45) | Write three PDF pages and matching PNGs from all pilot cells. |

### workflow.py

| Definition | Responsibility |
|---|---|
| [`PilotWorkflow`](scripts/hicat/workflow.py#L27) | Keep normalization, candidate fitting, validation and publication explicit. |
| [`PilotWorkflow.__init__`](scripts/hicat/workflow.py#L30) | Resolve frozen config and create new staging/progress state. |
| [`PilotWorkflow.qualify`](scripts/hicat/workflow.py#L44) | Record source edge-case probes before any PCDH19 fitting. |
| [`PilotWorkflow.run`](scripts/hicat/workflow.py#L74) | Execute three bounded candidate runs and stop at IN_REVIEW. |

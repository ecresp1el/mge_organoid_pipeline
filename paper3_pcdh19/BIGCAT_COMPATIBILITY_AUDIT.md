# BigCAT compatibility audit — contingency only

**2026-09-09 UTC. Production remains frozen HiCAT, array `60617037`, `0-99%10`,
q1=0.4, qdiff=0.7, DEscore=150. No BigCAT job was submitted or package installed.**
Population: 446,349 cells from 12 dissected E14.5 mouse MGE samples.

Inspected exact Allen commits:
[scrattch.bigcat `4f2818a`](https://github.com/AllenInstitute/scrattch.bigcat/tree/4f2818a88fdc232099d031a54bc9fb6f95951f83),
[pybigcat `d9e2a94`](https://github.com/AllenInstitute/pybigcat/tree/d9e2a946e08a6c4b75991af9e2fd87d2d426fd78).
Current fitting reference is `transcriptomic_clustering` `99154957`; our frozen
adapter is `305ce00711da7829c168ca7d53c669d0220f9587`. Source snapshots, SHA256s,
and function comparisons are saved in
[the audit references](references/bigcat_compatibility_audit_20260909/SOURCE_MANIFEST.json).

BigCAT shares the HiCAT scheme, but its aligned graph and DE/merge behavior
are substantive method differences. A package substitution would change the
scientific run, even with the same q1/qdiff/DEscore numbers.

## What matches, and what differs

| Component | Frozen production HiCAT | Inspected BigCAT path and compatibility |
| --- | --- | --- |
| Feature selection / PCA | Up to 3,000 HVGs; 50 randomized PCs fitted on every node's cells; elbow retains at most 20; projection in blocks. | Python `highly_variable_genes`, mean/variance feature functions, `project`, and elbow/component filters have identical function ASTs to our upstream reference. Our sparse allocation adapters add memory controls. R `onestep_clust_big` uses count-based variable-gene selection, typically 50,000 PCA-training cells, and projection of remaining cells. Those defaults differ. |
| Approximate KNN | Annoy, Euclidean, 50 trees; Annoy seed follows iteration fit seed; k=min(15, retained PCs, cells−1). | Python supports the same Annoy metric/tree machinery and PC cap, but its graph entrypoint fixes `annoy_seed=1` independently of community-detection seed. R uses `get_knn_batch`/BiocNeighbors Annoy; the big-data entrypoint can restrict graph-reference cells to 300,000 and predict the rest by KNN. These are additional sampling/seed differences. |
| Jaccard / SNN graph | Jaccard weights only on direct KNN edges; no BigCAT pruning rule. | R `jaccard_big` and Python `jaccard_snn` connect all pairs sharing neighbors, using shared/(2k−shared). R computes block products and writes Parquet edges; Python forms the full sparse `B @ B.T` before pruning. For >50,000 cells, the aligned path removes weights ≤0.05. A sparse product can still expand substantially in memory. |
| Community detection | Python-Louvain (`taynaud`), resolution 1.0, at recursive nodes. | Python offers Louvain and Leiden; R offers both. Selecting Leiden, changing graph topology, or adopting R's automatic Ward clustering below 2,000 cells changes the method. |
| DE score | Empirical-Bayes tests; sum of −log10(adjusted p) over qualifying genes, without a per-gene cap. | R `fast_limma` and Python aligned eBayes cap each gene's contribution at **20**. Therefore score=150 is not the same criterion across implementations. R DE sample-size/statistic caps also need explicit matching; Python's aligned cap is disabled by default. |
| DE merging | Two neighboring clusters; correlation of cluster profiles when dimensions >2; later rounds focus on merged destinations; multiple eligible disjoint merges. The upstream condition also merges a pair with exactly five DE genes even at high score. | Python aligned mode uses Euclidean candidate distances, rechecks all clusters' neighbors, merges the weakest pair plus eligible extras below half the score threshold, and uses `num < min_genes`. R `merge_cl_big` defaults to four candidates. These changes can alter surviving boundaries. |
| Recursion / disk storage | Fresh representation per branch, minimum 40 cells, depth guard 20, saved per-node models/graphs/DE evidence; independent fit/mapping/final checkpoints. Expression is sparse H5AD loaded into memory. | Python `iter_clust` shares recursive splitting and has backed-H5AD child writers. R `iter_clust_big` has FBM/Parquet expression/statistics access and different stopping/representation defaults. Neither consumes our completion manifests directly; an adapter is needed. |
| Held-out assignment / consensus | Maximum Pearson correlation to sampled-cluster mean marker profiles; every cell receives one label per iteration; sparse membership blocks and denominator=completed iterations. | R `collect_subsample_cl_matrix` uses the same mean/Pearson mapping and accepts saved `result$test.cl`, avoiding remapping. `compile_cl_mat` and `get_cell.cl.co.ratio` implement the same membership-factorized affinity. R's separate `map_cells_knn_big` offers approximate cosine/correlation search, which is not automatically identical. No consensus/held-out driver was found in the inspected Python package. |

Source entrypoints:
[Python graph](https://github.com/AllenInstitute/pybigcat/blob/d9e2a946e08a6c4b75991af9e2fd87d2d426fd78/pybigcat/clustering.py),
[Python merge](https://github.com/AllenInstitute/pybigcat/blob/d9e2a946e08a6c4b75991af9e2fd87d2d426fd78/pybigcat/merging.py),
[R big-data clustering](https://github.com/AllenInstitute/scrattch.bigcat/blob/4f2818a88fdc232099d031a54bc9fb6f95951f83/R/cluster_big.R),
[R graph](https://github.com/AllenInstitute/scrattch.bigcat/blob/4f2818a88fdc232099d031a54bc9fb6f95951f83/R/harmonize.R#L1013),
[R DE definition](https://github.com/AllenInstitute/scrattch.bigcat/blob/4f2818a88fdc232099d031a54bc9fb6f95951f83/R/de.genes.R#L581),
[R consensus](https://github.com/AllenInstitute/scrattch.bigcat/blob/4f2818a88fdc232099d031a54bc9fb6f95951f83/R/consensusCluster.R#L245).

## Reusing inputs and completed work

**Step 02 preprocessing does not need to be repeated.** The existing
`consensus_benchmark_20260908_235116/inputs/full_log1p_cpm.h5ad` contains the
approved cell/gene universe and sparse float64 **ln(1+CPM)** expression. Python
can read this AnnData format. R needs a streamed transpose/export with explicit
gene/cell IDs into its `big.dat` storage. If using R's log2 conventions, convert
`x/ln(2)` and translate expression-scale thresholds; do not normalize twice.
Count-based HVG routines can read the already-approved Step 02 raw-count asset.
This is input conversion, not another QC/filtering run.

**Completed iterations can be reused for storage or consensus optimization.**
Their all-cell labels, sampled masks, held-out labels/scores, marker prototypes,
and `membership_B.npz` fully specify each iteration's consensus contribution.
An importer can reconstruct R named `cl`/`test.cl` objects or `cl.mat`, preserving
the frozen universe, cluster namespaces, denominator, and sampled assignments.
Replacing only storage or implementing the same affinity/refinement mathematics
does not require fitting those iterations again. Exact-output qualification
would still be required before adopting a new implementation.

**Completed HiCAT fits cannot be relabeled as BigCAT-aligned fits.** Adopting
SNN edges, fixed Annoy seeds, capped DE scoring, or aligned merging requires
new fits for a BigCAT ensemble. Existing labels could support a separately
named aggregation/final-DE sensitivity analysis, but that would be a hybrid
analysis. Never mix new BigCAT fits into this frozen 100-iteration ensemble.
Saved PCA/Annoy assets remain useful for controlled comparisons when cell IDs,
features, metric and seed match; the current weighted graph is not an SNN graph.

## Contingency implementation and qualification

Prefer a separate adapter that reads the saved membership blocks and existing
normalized asset, retaining the qualified consensus algebra and checkpoints.
Use R's Parquet/FBM patterns for bounded expression/statistic access if needed.
The inspected Python package implements backed H5AD paths, but **no Parquet
loader or consensus implementation was found in its package source**, despite
broader README wording; these interfaces must be implemented and tested.

A static compatibility blocker also exists in
[Python `pca`](https://github.com/AllenInstitute/pybigcat/blob/d9e2a946e08a6c4b75991af9e2fd87d2d426fd78/pybigcat/dimension_reduction.py#L119):
it sets the selected-gene count to `sum(vidx)`, unlike our corrected boolean
count. Integer gene indices can therefore produce an incorrect allocation width.
Retain the qualified PCA adapter or fix and qualify this in a separate future
environment; do not install or patch it into production.

Before any future switch, compare saved-cell/seed fixtures for representation,
graph edges, DE scores, merge decisions, membership affinities and refinement;
measure real-data RSS, disk use and runtime; then seek explicit authorization
for the documented scope. This audit performs source inspection only and makes
no unmeasured BigCAT speed or full-scale equivalence claim.

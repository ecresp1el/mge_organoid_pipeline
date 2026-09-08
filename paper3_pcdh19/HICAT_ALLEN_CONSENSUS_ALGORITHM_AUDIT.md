# Allen consensus algorithm audit

**The scalable Allen collector classifies held-out cells before consensus.**
The earlier conditional co-sampling denominator is superseded. The source
also has a large-graph branch that **bypasses `iter_consensus_clust`**.

Scope: **446,349 cells × 19,071 genes; 12 dissected E14.5 mouse MGE samples**,
from approved Step 02. The pilot's 38/17 clusters do not determine the full-data
hierarchy or K. This audit initially created documentation and pinned
references only. See [the execution handoff](HICAT_VALIDATION_HANDOFF.md)
for subsequent benchmark submission and current status.

The trace below uses Allen R commit
[`9af2f04`](https://github.com/AllenInstitute/scrattch.hicat/blob/9af2f04cb837a7b53b005dee88a92eae337afead/R/consensusCluster.R).
“Unchanged” means algorithmically usable at this scale; R functions still need
a Python port or an explicit R interface.

| Function, in execution order | Source behavior | Full-scale disposition |
| --- | --- | --- |
| `run_consensus_clust` → `sample` | Independently selects 80% without replacement and saves selected IDs. | **Unchanged:** 357,079 sampled cells per iteration; eventually 100 iterations. `init.result=NULL`, so no frozen pilot parents. |
| `iter_clust` → `onestep_clust` | Recursively fits expression-based clusters and DE merges; returns labels and fit-derived markers. | **Qualified scalable execution required:** the existing Python fitter covers this stage but has dense expression allocations and documented differences from R. |
| `collect_subsample_cl_matrix` → `map_by_cor` | Assigns the 89,270 held-out cells to the maximum-Pearson-correlation sampled-cluster mean, using that iteration's marker genes. Retains sampled labels. | **Faithful blockwise port:** save prototypes, marker IDs, scores, and fitted/inferred flags. Every cell then has a label in every iteration. |
| `compile_cl_mat` | Builds sparse, iteration-specific one-hot membership matrix B. | **Unchanged mathematics:** B has sum(K per iteration) rows and 446,349 columns. Never match numeric cluster IDs across iterations. |
| Wrapper dispatch → `iter_consensus_clust` **or** full-data `iter_clust` → `merge_cl_by_co` | Computes G = sum(first iteration's full-population cluster sizes²). If G < 1e9, recursively resolves a representative consensus graph. Otherwise obtains a fresh full-population expression fit and merges it using consensus affinities. | **Preserve the measured branch.** Port representative selection/cutting for the smaller branch; use exact factorized affinities for the larger branch. Sparse graph fill-in needs an explicit memory guard. |
| `refine_cl` → `get_cl_co_stats` | Reassigns by mean consensus affinity, then redistributes cells from confused/small clusters. Uses median cluster confusion and specific stopping/update order. | **Faithful N×K port:** wrapper tolerance .01, confusion .6, 50 inner iterations, configured minimum size20. No cells discarded. |
| `merge_cl` | Final DE merging, with neighboring clusters prioritized in marker-expression space. | **Scalable, but not interchangeable with any DE routine:** retain requested q1=.4/qdiff=.7/score150 and explicitly name the Python DE backend. Audit residual separations before accepting K. |

**The historical sparse collector also used held-out assignment.** At the
[2018 tutorial-era source](https://github.com/AllenInstitute/scrattch.hicat/blob/91d0934e7eec31794b5945a5570ad92d7a3a9521/R/consensusCluster.R#L48),
`collect_co_matrix_sparseM()` calls that classifier, selects representatives,
and divides co-clustering counts by the number of iterations. The conditional
co-sampled denominator belongs to the separate legacy `collect_co_matrix()`.
The vignette's intuitive prose should not be used to override its actual
sparse implementation.

That legacy function allocates dense numerator and denominator matrices.
One N×N float64 matrix would require **1.59 TB**; two matrices plus the ratio
require about **4.78 TB**. This path is infeasible for the intended allocation.
Sparse membership does not imply that its cell–cell cross-product stays sparse.

**Recommended implementation:** port the pinned R consensus control flow in
Python, retain the validated Python expression-fitting settings, and compute
cell-to-cluster affinities exactly as:

```text
A = B.T @ (B @ H) / (R × target_cluster_sizes)
```

Here H is current cell-to-cluster membership and R is the number of successful
iterations—100 in production. A is N×K; the full N×N matrix is never needed
for Allen's large-graph branch. The estimator measures reproducibility of
fitting **and held-out classification**, not conditional co-sampling alone.

The existing Python package has **no consensus implementation**. Its natural-log
normalization, all-node-cell PCA, Louvain backend, eBayes calculation and
two-neighbor DE merging also differ from R defaults. Preserve the user's
documented settings, disclose these differences, and do not call this an
unchanged R reproduction. The pilot's marker-PCA final merge cannot silently
replace the R wrapper's marker-expression final merge.

Before benchmarking, qualify the port against pinned R on small numerical
fixtures and qualify sparse statistics against the existing Python backend.
Then benchmark one complete real 80% iteration, held-out mapping and the
measured full-population branch. B for 100 iterations needs approximately
**537 MB** as float64/int32 CSC, excluding labels and working arrays;
**runtime, peak RAM and total disk remain unmeasured**. Report these before
launching the 100-fit production analysis. Benchmark outputs stay IN_REVIEW;
neither annotations nor final K are locked.

Exact line references, branching diagram, parameter differences, resource
arithmetic, planned outputs/plots and AnnData omissions are in
[the technical notes](HICAT_ALLEN_CONSENSUS_TECHNICAL_NOTES.md).
[Pinned sources and proposed method contract](references/allen_consensus_audit_20260908/README.md)
are saved locally with checksums and licenses.

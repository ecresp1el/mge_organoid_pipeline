# Current checkpoint: Step 07 complete; full-data consensus source audit

**Experiment:** dissected E14.5 mouse medial ganglionic eminence (MGE).

Run `hicat_validation_20260908_221625_d5d06a5d` **completed**, Great Lakes job
**60591567**, elapsed 43:01. All **1,529 independent checks passed**, covering
1,025 manifest files. Its outputs are **IN_REVIEW**, not approved. No full-data
clustering was performed in Step 07. The main PDF has five pages and the
detailed PDF 76 pages.

**Latest authorized work:** prepare an Allen-faithful full-data consensus
benchmark using the approved Step 02 population (446,349 cells × 19,071 genes),
not the 12,000-cell pilot. First complete the
[Allen consensus algorithm audit](HICAT_ALLEN_CONSENSUS_ALGORITHM_AUDIT.md).
That source audit is now complete; the implementation still needs qualification
before timing. **No benchmark or production consensus job has been launched.**
The audit preserves the requested Python fitting settings while identifying
their differences from R. Held-out assignment supersedes the earlier
conditional co-sampling denominator. The full-data hierarchy and K are learned;
the four pilot parents and K=38/17 are not fixed. Benchmarking is authorized;
the definitive 100-fit production launch is not. Stop after a measured resource
report with outputs IN_REVIEW. See the audit for exact proposed assets and
what has not yet been saved.

The requested focused **38→17 boundary review also completed**, without another
HiCAT fit or parameter change. Start with
[the focused review](HICAT_MERGE_COLLAPSE_REVIEW.md). Its final table contains
84 substantially lost/eroded sibling pairs: five `keep distinction`, three
`merge`, and 76 `unresolved`. The comparison has 76 nonzero flows and 19 original
clusters split across destinations. These are pairwise review recommendations,
not a new cluster count. Overall disposition: **STILL UNRESOLVED; do not lock 17**.

Run directory:
`/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/02_hierarchy_validation/hicat_validation_20260908_221625_d5d06a5d/`

The first attempt, `hicat_validation_20260908_220825_c8ef1930` (job 60589514),
completed the controlled comparison: 38 baseline fine clusters → 17 candidate
clusters, with parent counts 13→8, 21→7, 3→1, 1→1. It then failed because the
old AnnData reader could not parse an unused modern `uns` encoding in Step 06.
That failed run is preserved. The continuation reads required source metadata
selectively with h5py and reuses the 665 completed candidate assets only after
checking their manifest, identical inputs, fitting code and controls. It does
not rerun or retune the completed comparison. The corrected 17-test suite and
actual full-source metadata/UMAP checks pass.

Completed inputs: Step 02 approved raw counts; Step 06 unintegrated display;
expanded 12,000-cell pilot `hicat_coarse_fine_20260908_203708_7e050f23`, with
baseline 4 coarse/38 fine and repeat 5 coarse/39 fine clusters. Existing seed
fits are reused. The one new fine comparison holds baseline coarse parents
fixed and changes q1 0.3→0.4 and qdiff 0.5→0.7, retaining score 150.

Step 07 evaluates canonical biological programs, cell cycle/maturation,
sample representation, seed stability, Allen-parameter sensitivity, hierarchy
relationships, and full-data marker/display context. User annotations remain
hypotheses; missing positive top markers are not missing expression data.

The `02_hierarchy_validation` directory is the HiCAT subworkflow asset stage;
**Step 07** is its primary-workflow checkpoint name. This is one analysis.

Live state: `provenance/hicat_progress_latest.json`; detailed function trace:
`provenance/hicat_progress_events.jsonl`; scheduler logs: `logs/scheduler.out`
and `logs/scheduler.err`. Results publish to `outputs/` only after workflow
checks. `COMPUTATION_SUCCESS.txt` additionally requires independent saved-asset
verification. Failures retain a failure marker and logs.

Start reviewing `REVIEW_README.md`, the main/detailed PDFs, decision tables,
`pcdh19_hicat_hierarchy_validation.h5ad`, the selected-gene full-data sidecar,
and complete `uns`/slot inventories. See
[the protocol](HICAT_VALIDATION_PROTOCOL.md) for exact content and omissions.

# Current checkpoint: Step 07 complete; full-data HiCAT production submitted

**Production is authorized and running. Current restartable array: 60617037,
`0-99%10`, submitted 2026-09-09 00:43:55 UTC (September 8, 8:43 PM Detroit).**
Each task requests **8 CPUs, 150 GiB and 48 hours**. All frozen scientific
settings, seeds and sampled-cell lists remain unchanged.

The initial array **60614879** exposed an Annoy Python seed-conversion error in
iterations **2, 3, 4, 5, 6, 8, 9** after PCA. These failed attempts are preserved.
The qualified [exact-seed compatibility bridge](ANNOY_SEED_BINDING_RECOVERY.md)
forwards the original positive seed to the same installed native Annoy engine.
It changes neither seed values nor the scientific method or installed package.
All 100 seed fixtures passed; safe-seed graphs/labels matched the old binding;
checkpoint provenance and reuse passed.

The ten already-active iterations continue unchanged: benchmark iterations 0/1
and original-array iterations 7, 12, 13, 14, 15, 20, 26, 30. The replacement array
waits for these writers and then verifies/reuses completed checkpoints. Only
incomplete stages are computed. Pending superseded tasks were cancelled; no
running fit was cancelled or resized. This is one 100-iteration ensemble.

Exact production directory:
`/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/03_full_data_consensus_benchmark/consensus_restart_20260909_000256`.
`PRODUCTION_SUBMISSION.json` records current and previous arrays;
`production/STATUS.json` and `production/iteration_status.tsv` distinguish live
scientific iterations from replacement tasks waiting to reuse them. Refresh with
[the status command](bin/consensus_production_status.py).

Separately rerunnable production aggregation is **60617041**, and final DE is
**60617044**. Aggregation also waits for the real two-iteration benchmark/report
and requires all 100 sealed membership outputs. The benchmark is still running;
complete-iteration runtime, held-out mapping and real aggregation remain unmeasured.

The [BigCAT compatibility audit](BIGCAT_COMPATIBILITY_AUDIT.md) is saved. No
BigCAT package was installed and no BigCAT production jobs were submitted.
The user's later authorization supersedes the earlier benchmark approval gate.

**Experiment:** dissected E14.5 mouse medial ganglionic eminence (MGE).

Run `hicat_validation_20260908_221625_d5d06a5d` **completed**, Great Lakes job
**60591567**, elapsed 43:01. All **1,529 independent checks passed**, covering
1,025 manifest files. Its outputs are **IN_REVIEW**, not approved. No full-data
clustering was performed in Step 07. The main PDF has five pages and the
detailed PDF 76 pages.

**Earlier benchmark preparation:** prepare an Allen-faithful full-data consensus
benchmark using the approved Step 02 population (446,349 cells × 19,071 genes),
not the 12,000-cell pilot. First complete the
[Allen consensus algorithm audit](HICAT_ALLEN_CONSENSUS_ALGORITHM_AUDIT.md).
That source audit is complete. Seven numerical qualification tests passed,
including pinned-R comparisons, sparse/dense DE and fitting, and saved AnnData.
**The real full-data benchmark continues alongside the production submission above.**
Run `consensus_benchmark_20260908_235116` lives under
`results/hicat/03_full_data_consensus_benchmark/`. Preparation is job **60609367**,
the first 80% fit **60609369**, and full-data initialization **60609372**.
The user subsequently requested restart-safe independent iterations and a test
aggregator using 2–5 real completed iterations. The pending one-iteration
consensus/report jobs **60609374 / 60609377** were cancelled to replace those
stages with the expanded two-real-iteration validation. Completed/running
preparation and fit outputs are preserved for reuse. See
[the execution and asset contract](FULL_DATA_CONSENSUS_BENCHMARK.md).
The detailed independent-checkpoint design and real failure/retry test are in
[the restart protocol](HICAT_CONSENSUS_RESTART_PROTOCOL.md).
The restart package `consensus_restart_20260909_000256` has now frozen all
100 seed/sample-ID sets and adapter commit
`305ce00711da7829c168ca7d53c669d0220f9587`. The initial two-real-iteration benchmark
chain was submitted: adopt first fit **60612803**, second real fit **60612805**,
first mapping **60612807**, intentional second-mapping failure **60612809**,
mapping-only retry **60612811**, aggregation **60612813**, final DE **60612815**,
and validation/report **60612817**. The all-cell initializer remains **60609372**.
These job IDs establish submission, not completion; inspect scheduler state
and sealed checkpoints for actual progress. `AUTHORIZATION.json` now sets `production_allowed=true` under the later user
authorization. The frozen scientific config retains historical approval metadata;
the separate authorization record controls execution.
Actual states are in scheduler accounting and `progress/STAGE/` logs.
The user subsequently authorized production before the completed resource report;
the benchmark still continues to provide measured resource and restart evidence.
The audit preserves the requested Python fitting settings while identifying
their differences from R. Held-out assignment supersedes the earlier
conditional co-sampling denominator. The full-data hierarchy and K are learned;
the four pilot parents and K=38/17 are not fixed. Benchmarking is authorized;
the definitive 100-fit production launch is now explicitly authorized. Scientific
outputs remain IN_REVIEW. See the audit for exact proposed assets and
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

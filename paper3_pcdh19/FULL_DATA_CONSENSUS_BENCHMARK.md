# Full-data consensus benchmark — execution and saved assets

**Execution update, 2026-09-09 01:07 UTC:** the first full-data 80% fit and
held-out assignment completed in **01:03:22**. Its fit was adopted and its
all-cell membership was independently published/validated. One of 100
iterations is complete, nine are running, and the separate full-data consensus
initializer is running. The recovery array and two-real-iteration aggregation
test remain pending. See [current status, measurements and saved assets](HICAT_PRODUCTION_STATUS.md).

**Extension requested after submission:** the benchmark now requires two real
independent iterations and explicit restart/reuse validation. See the
[restart protocol](HICAT_CONSENSUS_RESTART_PROTOCOL.md). The original preparation,
first fit and full-data initialization are preserved; pending single-iteration
consensus/report jobs were cancelled. The original one-iteration plan below
records that initial submission, not the final report's expanded scope.

**Authorization update:** the user subsequently authorized production before
benchmark completion. Array **60617037**, `0-99%10`, is submitted with unchanged
settings and independently reusable checkpoints. See the
[current handoff](HICAT_VALIDATION_HANDOFF.md). The original benchmark plan below
is historical; its earlier stop-for-approval condition is superseded.
The production target is 100 independent 80% fits of approved Step02's
446,349 dissected E14.5 mouse MGE cells, with q1=.4, qdiff=.7, DE score150.
The source audit and historical reference files remain preserved.

## Historical initial benchmark plan

The original [launcher](bin/submit_consensus_benchmark.py) submitted a sequential
benchmark dependency chain. The table preserves that original R=1 plan; the
active benchmark now tests two real iterations, and production uses the
separate restartable 100-task array described above.

| Stage | Actual input and operation | Primary saved outputs |
| --- | --- | --- |
| `prepare` | Verify the approved Step02 SHA256, cell/gene dimensions and 12 sample IDs; normalize all 446,349 cells in sparse blocks | Shared `inputs/full_log1p_cpm.h5ad`, cell index, sample counts, source/hash provenance and complete `uns`/slot inventories |
| `fit80` | Uniformly sample exactly 357,079 source rows without replacement; fit the existing recursive Python HiCAT settings; classify all 89,270 held-out cells by fit-marker mean/Pearson correlation | Exact sampled rows/IDs, sampled assignments, all-cell assignments with fitted/inferred flags and correlations, prototype means, mapping marker IDs, every node's PCA/features/graph/DE evidence and recursion tree |
| `fullfit` | Measure Allen's G=sum(first iteration full-cell cluster sizes²). When G≥1e9, fit all 446,349 cells to initialize Allen's large-graph consensus branch | Actual G and chosen branch, full-population fit models/labels/DE evidence; no pilot parents or fixed K |
| `consensus` | Execute Allen's factorized membership, co-merging, refinement and final marker-expression DE merge using the one real iteration | Sparse R=1 B and index map, merge/refinement traces, all-cell diagnostic labels, full pairwise endpoint DE audit, expression means/detection, resource tables and review PDF |
| `report` | Collect completed stage wall time, process CPU time, peak RSS, exact on-disk bytes and cautious projections | `outputs/RESOURCE_REPORT.md`, JSON, measured stage/operation TSVs, manifest, `BENCHMARK_COMPUTATION_SUCCESS.txt`; state IN_REVIEW |

The measured branch is never chosen from pilot K=38/17. The current qualified
Python consensus implementation covers Allen's **large-graph branch**. If the
real first iteration selects the smaller representative-graph branch, the
benchmark stops explicitly for implementation of that exact branch; it does
not substitute a different inference method or call the result complete.

One real iteration cannot estimate 100-iteration stability or establish final
K. A separately labeled **synthetic capacity exercise** repeats the first
label vector into 100 disjoint iteration blocks to measure B storage and one
exact affinity pass at 44,634,900 memberships. It is not 100 biological fits;
its timings cannot establish production convergence or compression behavior.

## Scientific settings and fidelity

Use the [completed algorithm audit](HICAT_ALLEN_CONSENSUS_ALGORITHM_AUDIT.md)
and [detailed source trace](HICAT_ALLEN_CONSENSUS_TECHNICAL_NOTES.md).
Consensus membership, mean/Pearson mapping, co-merging and refinement are
qualified against Allen R commit `9af2f04`. Expression fitting and DE retain
the validated Allen Python `99154957` settings, including ln(1+CPM), all-node-cell
randomized PCA, Louvain, and two-neighbor eBayes merging. These are disclosed
differences from an unchanged R/limma run. The pilot adapter's pre-consensus
cross-branch marker-PCA merge remains part of the existing Python fitter;
the **post-consensus** merge uses marker expression directly, as the R wrapper
does, with the configured Python DE backend.

Memory changes preserve calculations: sparse sample variance/detection,
blocked projection, all-cell randomized PCA explicitly selected, and release
of parent expression before recursive children. No full N×N consensus object
is created. B is iteration-clusters × cells and the required affinity is N×K.
Numerical fixtures check factorization against explicit pair frequencies,
arbitrary-ID invariance, pinned R operations, sparse versus dense statistics/DE,
dense versus sparse fitting and the block-written normalized AnnData.

Mapping markers are the explicitly saved first20 up/down genes from each
persisted fitter DE call, matching the existing adapter's marker policy.
Canonical biological marker panels, provisional annotations, genotype, sample
identity, cell-cycle scores and previous clusters do not enter fitting.
No batch correction, integration or cell-cycle regression is performed.

## AnnData: exactly what is and is not saved

`inputs/full_log1p_cpm.h5ad` saves all source rows/genes, X as **float64 CSR
ln(1+CPM)**, approved-source metadata in `obs`, and gene symbols in `var`.
`uns['consensus_benchmark']` records the representation, linked raw-count asset
and hash, tissue context and explicit analysis omissions. Full inventories
are beside the object. `layers`, `raw`, `obsm`, `obsp`, `varm`, and `varp` are
empty. The raw counts remain in the unchanged approved Step02 H5AD.

Benchmark clustering assignments and consensus arrays are saved separately as
TSV/NPY/NPZ assets with source-row/cell-ID mappings. **No final consensus AnnData,
full pair matrix, annotation column or accepted final K is produced.** Expression
is saved once, not duplicated for each fitting stage.

## Plots and progress

`outputs/benchmark_review.pdf` contains measured stage-time/peak-RSS bars,
R=1 diagnostic cluster sizes, all-cell within-cluster sample-composition
heatmap, and an expression-mean display dendrogram. The dendrogram is explicitly
labeled a display summary, not the iterative fitting tree or definitive
consensus. The saved recursion tree and node memberships provide the actual
fit hierarchy. Complete DE evidence is numerical; biological annotation review
and definitive top-marker interpretation follow the eventual consensus.

Each stage writes `progress/STAGE/hicat_progress_events.jsonl` and
`hicat_progress_latest.json`, printing the same START/COMPLETE/FAILED events.
Stage resources distinguish process CPU time, wall time, and lifetime peak RSS.
Scheduler stdout/stderr are `logs/STAGE.out` and `.err`; stderr includes
`/usr/bin/time -v`. `job_ids.json` records actual IDs. Code/config/references
are copied before submission and checksum-verified by every stage.

## Historical resource plan and current measurements

The first fit/mapping job measured **01:03:22**, process peak **54.71 GiB**,
SLURM batch peak **61.24 GiB**. Subsequent adoption and membership publication
jobs took **01:05** and **00:49**. They reused the completed fit; their small
resource records must not be interpreted as fitting costs. Total ensemble and
final consensus costs remain provisional. The later user authorization fixes
production at **8 CPUs / 150 GiB / 48 hours / concurrency 10**; the older
resource recommendation below does not change that authorization.

The initial benchmark allocation is **8 CPUs, 150 GiB, 48-hour time limit**
for each expression-fit stage; these are limits, not measured consumption or
production recommendations. Preparation gets 1 CPU/32 GiB; consensus gets
8 CPUs/150 GiB; reporting gets 1 CPU/8 GiB. Stages execute sequentially. No
GPU is used. The report replaces estimates with observed costs and recommends
memory/time margins and initially at most two concurrent production fits.
Parallel completion-time estimates exclude queue time and filesystem contention.

Production total CPU cost is projected from the measured full 80% fit/mapping
cost, plus the full-data initialization and consensus costs. The report clearly
separates measured R=1 consensus work, synthetic R100 capacity timing, and the
unknown number of R100 refinement passes. It does not invent a measured total
for an analysis that has not run. Failed stages retain their artifacts and
prevent downstream success. Surviving weak DE pairs are reported, never turned
into an accepted final K.

The original benchmark launcher has no production-array interface. The separately
authorized [production launcher](bin/submit_consensus_production.py) submits the
array using unchanged frozen scientific code; it does not rerun completed fits.

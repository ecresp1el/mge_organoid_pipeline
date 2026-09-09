# Full-data HiCAT consensus: current execution status

For the **exact executed scRNA-seq parameters and how consensus is formed**,
read [the parameter and consensus guide](HICAT_PARAMETERS_AND_CONSENSUS.md).

**Latest snapshot: 2026-09-09 14:24:32 UTC / September 9, 10:24:32 AM Detroit (EDT).**
**99 of 100 production iterations are complete and sealed; iteration 99 is running.**
There are no pending iteration tasks or current failed tasks. Production
aggregation **60617041** and final DE **60617044** remain pending on dependencies.

**The requested 98-iteration consensus is complete and IN_REVIEW:** 34 clusters,
446,349 cells. Its corrected stage-03 publication is
`hicat_consensus_review_20260909_142348_efd9e13d`, job **60712671**, with all
35 validation checks passing. Use the [review handoff](HICAT_CONSENSUS_REVIEW_HANDOFF.md)
for canonical Turbo assets and numbered figures.

## Earlier 9:52 AM snapshot and measured timing

All clock times below are Detroit (EDT). Durations come from SLURM accounting;
completion counts also check fit/mapping completion manifest seals, not a fresh
hash of every scientific artifact.

| Work | Current state | Timing / evidence |
| --- | --- | --- |
| Shared normalized input; 100 frozen seeds and sample selections | Complete | Existing normalized H5AD, `FROZEN.json`, frozen per-iteration inputs |
| Iterations 0–97 | Complete; fit and mapping seals verified | 98 completed tasks in `60617037`; includes checkpoint reuse |
| Iteration 98 | Running | `60617037_98`; started Sep 9 **9:02:11 AM**, about **50 min** elapsed |
| Iteration 99 | Running | `60617037_99`; started Sep 9 **9:08:44 AM**, about **44 min** elapsed |
| Full-data expression fit used to initialize consensus | Complete | `60609372`; Sep 8 **8:58:56–10:36:20 PM**, **1 h 37 min 24 s** |
| Deliberate mapping failure and independent retry | Validation complete | `60612809` deliberately failed; `60612811` completed Sep 8 **9:43:54 PM** in **1 min 54 s**, reusing the fit |
| Two-iteration benchmark aggregation | Complete and sealed | `60612813`; Sep 8 **10:36:27–10:38:24 PM**, **1 min 57 s** |
| Two-iteration benchmark final DE | Complete and sealed | `60612815`; Sep 8 **10:38:32–10:42:42 PM**, **4 min 10 s** |
| Benchmark report and restart checks | Complete | `60612817`; finished Sep 8 **10:46:54 PM**, **3 min 42 s** |
| Production consensus aggregation/refinement | Pending: Dependency | `60617041`; awaits the full iteration array |
| Production final DE | Pending: Dependency | `60617044`; awaits production aggregation |

The recovery array started Sep 8 at **9:59:29 PM**, approximately **11 h 53 min**
before this snapshot. Of its completed tasks, **88 performed fresh iterations**
and ten reused earlier checkpoints. Fresh-task durations have a **74 min 39 s
median**, **89 min 18 s 90th percentile**, and **58 min 43 s–2 h 52 min 50 s
range**. Reuse tasks lasting about a minute are excluded from these statistics.

The seven earlier Annoy seed failures (**2, 3, 4, 5, 6, 8, 9**) have now all
completed successfully through the recovery array. Historical failed attempts
remain preserved; no retry is currently needed.

## Timing estimate recorded at 9:52 AM

- **Last two iterations: likely September 9, 10:15–10:40 AM Detroit.** This is
  an estimate from their start times and the observed median/upper-tail fresh
  iteration durations, not a scheduler guarantee. Applying the longest observed
  runtime to the later start gives approximately **12:02 PM**; that is an
  observed slow-run scenario, not a hard upper bound.
- **Production aggregation and final DE: already queued to follow.** The
  two-iteration benchmark took roughly six minutes of combined job runtime,
  but that does not establish the runtime for 100 memberships. Queue delay,
  filesystem load and consensus refinement can add time. A reliable completion
  clock time for the full run is not yet available; reassess once aggregation
  starts. The iteration ETA above is not the final-consensus ETA.
- **Biological review follows computation:** inspect markers, stability and
  sample composition before accepting cluster identities or a taxonomy.

## Current requested deliverable: 98-iteration consensus

The user requested the consensus of frozen iterations **0–97**, denominator
**98**, without waiting for the remaining two iterations. Computation job
**60707938** completed at **10:06:56 AM Detroit**, elapsed **6 min 12 s**, exit 0.
Aggregation/refinement produced 35 clusters; final DE produced **34 clusters**.
The 561-pair DE audit leaves **three pairs below the separation criteria**.

The validated review publication is within existing HiCAT stage **03**,
not a new numbered major step. It follows frozen object-oriented code, thin
entry points, numbered figures with source tables, manifested atomic outputs,
and `IN_REVIEW` status. See the [review handoff](HICAT_CONSENSUS_REVIEW_HANDOFF.md)
and [architecture/numbering contract](HICAT_CONSENSUS_REVIEW_PROTOCOL.md).

The source computation and cancelled early job **60707748** remain preserved.
The ad hoc script and generated result copies previously placed in the GitHub
checkout were relocated with checksum verification into the review run's
`provenance/replaced_repository_artifacts/`. Canonical scientific results and
review figures are on Turbo. The full 100-iteration production jobs continue.

## Inputs and the first completed iteration's saved assets

Population: **446,349 cells × 19,071 genes from 12 dissected E14.5 mouse MGE
samples**. Each iteration uses its frozen **357,079 sampled cells** and maps
**89,270 held-out cells**. The shared expression input is sparse float64
ln(1+CPM), saved in the original benchmark's `inputs/full_log1p_cpm.h5ad`.

Production run directory:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/03_full_data_consensus_benchmark/consensus_restart_20260909_000256
```

Paths below are relative to that directory. Resolve each stage's `CURRENT.json`
to its immutable successful attempt before reading assets.

| Saved asset | Meaning / validation |
| --- | --- |
| `iterations/000/fit/CURRENT.json` | Points to a validated fit attempt containing **2,127 manifest-listed artifacts**: models, node memberships, PCA/features, graphs, DE and provenance. Saved sampled partition and cell-ID alignment passed. |
| `iterations/000/mapping/CURRENT.json` | Points to a validated mapping attempt containing **12 manifest-listed artifacts**. |
| Mapping attempt: `all_cell_assignments.tsv.gz` | **446,349 rows** with `cell_id`, `cluster`, `sampled_for_fit`, and `heldout_correlation`. |
| Mapping attempt: `all_cell_labels.npy` | Integer cluster assignments in the frozen full-cell order. |
| Mapping attempt: `membership_B.npz` | Sparse **67 × 446,349** matrix, **446,349 nonzeros**; every cell assigned once and sampled labels preserved. No full cell–cell matrix was created. |
| Both attempts: `VALIDATION.json`, `ARTIFACTS.json`, `COMPLETE.json` | Reopen/alignment evidence, SHA256 manifests and completion seals; verified inputs for later reuse. |

The 67 membership rows describe **one iteration's output**, not an accepted
consensus K. Its fitter reports 48 unresolved final DE pairs, so completion is
an execution/checkpoint result, not evidence that every boundary is settled.
The two-iteration review PDF is available at `report/benchmark_review.pdf`,
alongside `RESOURCE_REPORT.md` and `RESTART_VALIDATION.json`. Final production
consensus labels and an accepted taxonomy are not yet available. Assignments are separate TSV/NPY/NPZ assets; no final
annotated consensus AnnData has been saved.

## Resources, recovery and remaining validation

The first fit/mapping job measured **1 hour 3 minutes 22 seconds**. Adoption and
membership publication were separate jobs; their short `resources.json` times
must not be mistaken for a fresh full-data fit. First-fit peak RSS was **54.71
GiB** from the process record and **61.24 GiB** from SLURM batch accounting.
These first-fit measurements are historical. The current production runtime
distribution is reported above; final consensus cost remains unmeasured.

The authorized request stays **8 CPUs, 150 GiB, 48 hours per task, maximum
concurrency 10**, q1=0.4, qdiff=0.7, DEscore=150. Normalization, seeds, sampled
cell lists and all other scientific settings remain frozen. Completed fit and
mapping stages are independently reusable; final aggregation and DE remain
separately rerunnable from saved outputs.

The [Annoy compatibility fix](ANNOY_SEED_BINDING_RECOVERY.md) preserves exact
positive seed values. All 100 seed fixtures, checkpoint reuse and Python/NumPy
RNG-state checks passed; 364 frozen scientific/input files matched. The seven previously failing full-data iterations have now succeeded in the
recovery array, providing full-data execution evidence beyond the seed fixtures.
Startup hash checks and bridge provenance are implemented. The real benchmark
mapping retry and successful-stage replay checks have also completed.

The [BigCAT audit](BIGCAT_COMPATIBILITY_AUDIT.md) is complete. No BigCAT production
job or package substitution was made. No batch correction, cycle regression,
or locked biological annotations were introduced.

Refresh the scheduler/checkpoint snapshot with the command in the
[restart protocol](HICAT_CONSENSUS_RESTART_PROTOCOL.md). `production/STATUS.json`
and `production/iteration_status.tsv` are timestamped operational snapshots;
they are not continuously updated displays. This Markdown records the stated
snapshot. Frozen run copies and sealed scientific artifacts are unchanged.

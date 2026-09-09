# Full-data HiCAT consensus: current execution status

**Snapshot: 2026-09-09 01:07 UTC / September 8, 9:07 PM Detroit.**
**1 of 100 iterations is complete, 9 are running, and 90 await execution/retry.**
The separate full-data consensus initializer is also running. Recovery array
**60617037** has **100 pending scheduler tasks**: it waits for the original
writers, then verifies/reuses their completed checkpoints. Scheduler tasks and
scientific iterations are different counts; this remains one 100-iteration run.

## What has actually run

| Work | Observed state | Job / evidence |
| --- | --- | --- |
| Shared normalized input; 100 frozen seeds and sample selections | Complete | Existing normalized H5AD; `FROZEN.json`; `config/seeds.tsv`; per-iteration input files |
| Iteration 0: HVG/PCA, graph, recursive HiCAT, DE merging and held-out assignment | Complete | `60609369`, elapsed **01:03:22**, exit 0 |
| Iteration 0: adopt/reopen fit and publish consensus membership | Complete and sealed | `60612803` (01:05), `60612807` (00:49); both exit 0 |
| Iteration 1 | Running, about 53 minutes elapsed | `60612805` |
| Iterations 7, 12, 13, 14, 15, 20, 26, 30 | Running, about 28–39 minutes elapsed | Corresponding `60614879_INDEX` tasks; no new failures recorded |
| Full-data expression fit used to initialize consensus | Running, about 8 minutes elapsed | `60609372`; separate from the 100 subsample iterations |
| Recovery array and seed-fix retries | Pending on existing writers | `60617037`, `0-99%10`; no retried full-data task has started |
| Iteration 1 mapping failure/retry test; two-real-iteration aggregation and final DE test | Not started | `60612809/60612811`, `60612813`, `60612815`, report `60612817` |
| Production consensus aggregation/refinement and final DE | Not started | Separate jobs `60617041` and `60617044` |

Seven earlier attempts failed at Annoy seed initialization: **2, 3, 4, 5, 6,
8, 9**. Those failures are preserved and queued for retry; they are not new
failures in the nine currently running attempts. Pending superseded tasks were
cancelled during recovery, without cancelling or resizing running fits.

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
Final consensus labels, the two-iteration review PDF, and an accepted taxonomy
are not available. Assignments are separate TSV/NPY/NPZ assets; no final
annotated consensus AnnData has been saved.

## Resources, recovery and remaining validation

The first fit/mapping job measured **1 hour 3 minutes 22 seconds**. Adoption and
membership publication were separate jobs; their short `resources.json` times
must not be mistaken for a fresh full-data fit. First-fit peak RSS was **54.71
GiB** from the process record and **61.24 GiB** from SLURM batch accounting.
These measurements cover one iteration; total ensemble wall time, variability,
and final consensus cost remain provisional.

The authorized request stays **8 CPUs, 150 GiB, 48 hours per task, maximum
concurrency 10**, q1=0.4, qdiff=0.7, DEscore=150. Normalization, seeds, sampled
cell lists and all other scientific settings remain frozen. Completed fit and
mapping stages are independently reusable; final aggregation and DE remain
separately rerunnable from saved outputs.

The [Annoy compatibility fix](ANNOY_SEED_BINDING_RECOVERY.md) preserves exact
positive seed values. All 100 seed fixtures, checkpoint reuse and Python/NumPy
RNG-state checks passed; 364 frozen scientific/input files matched. Iteration 0
used the original successful path, so its completion does **not** qualify a
full-data run through the new bridge. That check remains pending until a retry
starts. Startup hash checks and bridge provenance are already implemented.

The [BigCAT audit](BIGCAT_COMPATIBILITY_AUDIT.md) is complete. No BigCAT production
job or package substitution was made. No batch correction, cycle regression,
or locked biological annotations were introduced.

Refresh the scheduler/checkpoint snapshot with the command in the
[restart protocol](HICAT_CONSENSUS_RESTART_PROTOCOL.md). `production/STATUS.json`
and `production/iteration_status.tsv` are timestamped operational snapshots;
they are not continuously updated displays. This Markdown records the stated
snapshot. Frozen run copies and sealed scientific artifacts are unchanged.

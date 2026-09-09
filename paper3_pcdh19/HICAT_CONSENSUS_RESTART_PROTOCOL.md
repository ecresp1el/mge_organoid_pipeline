# Independent HiCAT iteration and consensus restart contract

**Current scope: production is explicitly authorized and array 60617037
(`0-99%10`) is submitted.** At **2026-09-09 01:07 UTC**, iteration 0 has sealed
fit/mapping outputs, nine iterations are running and 90 await execution/retry.
The recovery array itself is pending; full-data initialization is running.
The second benchmark iteration and remaining restart/aggregation tests are
incomplete. See [current status and saved assets](HICAT_PRODUCTION_STATUS.md).
All 100 frozen seed/sample-ID sets are used unchanged; benchmark
iterations 0/1 are reused by their corresponding array tasks. Each task requests
8 CPUs, 150 GiB and 48 hours. See [the live handoff](HICAT_VALIDATION_HANDOFF.md)
and `PRODUCTION_SUBMISSION.json` for the exact run and separately queued final
stages. The user superseded the earlier wait-for-benchmark approval condition.

The first array's seven Annoy signed-integer seed failures are preserved. The
replacement uses a separately hashed, qualified native seed bridge; all 100
original positive seeds, sample lists and scientific settings are unchanged.
See [the binding recovery](ANNOY_SEED_BINDING_RECOVERY.md). Existing active
iterations finish before the replacement verifies/reuses their checkpoints.

Population: approved Step02, **446,349 cells × 19,071 genes from 12 dissected
E14.5 mouse MGE samples**. Each iteration uniformly samples **357,079 cells**
without replacement and maps **89,270 held-out cells**. No pilot parents or K,
biological annotations, sample/genotype, integration, batch correction or cycle
regression enter fitting. Use q1=.4, qdiff=.7, DEscore150.

## Numerical method

The expression backend remains Allen Python `99154957`; consensus operations
follow the audited Allen R `9af2f04` large-graph branch. Use **PCA**, not WGCNA:
up to 3,000 HVGs; randomized PCA on all current-node cells; 50 components
computed and an elbow filter retaining at most 20. Each recursive branch fits
its own representation. The first real root fit retained six PCs; that is an
observed result, not a new fixed parameter. Graph settings retain the existing
PC-dimension cap on requested k=15.

See [the source audit](HICAT_ALLEN_CONSENSUS_ALGORITHM_AUDIT.md) for the exact
functions and the disclosed differences between Python fitting/DE and R.
Held-out labels come from maximum mean-prototype Pearson correlation on saved
fit-derived marker genes. Every successful iteration supplies one assignment
for every cell. The denominator counts complete iterations after mapping.

## Independent checkpoints

| Checkpoint | Reads | Saves | Retry behavior |
| --- | --- | --- | --- |
| `iterations/NNN/fit` | Frozen seed, selected IDs/rows, shared normalized expression, fitting config/code | Sampled labels; per-node features, PCA models, graph/Annoy index, memberships, merge/DE traces; resources; verified artifact manifest | Retry only that iteration's fit if no valid completion exists. Other iterations remain untouched. |
| `iterations/NNN/mapping` | That iteration's sealed fit and shared normalized expression | Held-out scores/prototypes; all-cell labels and fitted/inferred flags; exact CSC `membership_B.npz`; cell-universe/index contract; resources and manifest | Reuse the completed fit. A mapping failure never requires refitting that iteration. |
| `final/benchmark/aggregation` | Two **real completed** membership blocks plus the full-data initializer | Combined sparse B, co-merge/refinement traces, all-cell diagnostic partition and numerical validation | Rerun aggregation from saved blocks; no iteration recomputation. |
| `final/benchmark/DE` | Sealed aggregation, initializer marker IDs and shared expression | Final DE merge/audit, cluster means/detection and all-cell diagnostic assignments | Rerun DE only. Iterations and successful aggregation remain reusable. |
| `final/production/aggregation` | All 100 sealed real membership blocks and the full-data initializer | Combined sparse B, co-merge/refinement traces and candidate consensus assignments | Rerun this stage from saved iteration outputs; no iteration recomputation. |
| `final/production/DE` | Sealed production aggregation and shared marker expression | Final DE merge/audit and candidate all-cell assignments | Rerun DE only; accepted K still requires review. |

Each stage has append-only `attempts/TIMESTAMP_ID/` directories. A successful
attempt contains `VALIDATION.json`, `resources.json`, `ARTIFACTS.json` and
`COMPLETE.json`. `CURRENT.json` points only to a sealed successful attempt.
Failure leaves `FAILURE.json` and its partial artifacts; it does not publish
completion. Runtime lock ownership is checked against the scheduler before a
stale lock can be archived. A second live process cannot take the same stage.

Every artifact is reopened and SHA256-checked before completion. Validation
checks sampled/full cell IDs, sample-row alignment, partition lengths, PCA
dimensions/gene IDs, projection IDs, sparse graph structure, and membership
columns. All TSV/JSON/NPY/NPZ and binary Annoy artifacts have explicit readers.
Consumers revalidate the input seals; missing/corrupt/mismatched outputs cannot
silently enter the aggregate. No successful output directory is overwritten.

## What is frozen

- `adapter/` is an independent local Git snapshot with its own commit. It
  includes the actual worker code and oracle/restart tests; nothing is pushed.
- `reference/allen_python/` and the pinned R files record exact Allen commits.
- `config/run.json` binds normalization, dimensions, graph/PCA/DE settings,
  held-out assignment, storage format, final merge rules and input hashes.
- `config/environment.txt` freezes package versions; workers compare the actual
  environment before running. `FROZEN.json` checks code/config/reference bytes.
- `config/cell_universe.tsv.gz` and `gene_ids.npy` freeze the common row/column
  order. `config/seeds.tsv` lists **all 100** sampling/fitting seeds.
- Each `iterations/NNN/config.json` stores its seed/config/input hashes.
  `sampled_source_rows.npy` and `sampled_cell_ids.tsv.gz` freeze its exact draw.

Iteration 0 preserves the first completed benchmark seed and outputs.
Iterations 1–99 use separate NumPy SeedSequence children with distinct
sampling/fitting seeds. Seed derivation is frozen by the saved schedule; retries
read that schedule rather than generating a new seed. Adoption of the completed
first fit records the original code/config hashes and validates every copied
artifact. The second iteration runs the new production-like worker end to end.

`AUTHORIZATION.json` is separate from frozen scientific inputs. It now permits
all 100 iterations under the user's explicit production authorization, with a
maximum concurrency of ten. The prior authorization record is archived in the
production submission directory. Approval changes neither seeds nor completed
scientific contracts. BigCAT production remains unauthorized.

## Real restart and aggregation tests

These real-data restart/aggregation tests have not yet run as of the status
snapshot above. The second iteration's first mapping attempt is configured to deliberately raise an error
**after** writing the membership output and **before** publishing completion.
A separately queued retry reads the successful fit and completes mapping.
This validates recovery from actual full-data output, beyond small fixtures.
The failed attempt remains available for inspection.

The aggregator then reads both saved real membership files and verifies that
the concatenated B equals the saved inputs. Each cell contributes twice, and
the denominator is 2. It executes Allen co-merging/refinement and the configured
final DE stage. This is an end-to-end I/O/algorithm test; **K is not interpreted**.
No synthetic repeated partition substitutes for either real iteration.

The final report replays every successful fit/mapping/aggregation/DE request.
Each must return the same completed attempt and unchanged artifact hashes.
Thus restart tests do not create another fit. Small qualification fixtures also
test an injected failure, corruption detection and combination of saved B files.

## How to inspect or resubmit

Use [the production status command](bin/consensus_production_status.py) to
refresh scheduler counts and all 100 fit/mapping checkpoint states:

```bash
python paper3_pcdh19/bin/consensus_production_status.py --run-dir RUN_DIRECTORY
```

For failed iterations, use the status JSON's
`recommended_single_retry_array_command`: it supplies a fresh `sbatch` array
containing only failed indices, with cap 10 and the current script/resources.
For example, failures 17, 43 and 88 use `--array=17,43,88%10` in that full
command. It waits for the current array to terminate. **Do not use same-job-ID
`scontrol requeue` after a killed stage:** a surviving execution lock would
appear owned by that active job ID. Fresh job IDs let the worker recognize and
archive stale locks while preserving previous attempts.

The array runs one `iteration` job for each ID 0–99. `iteration` executes fit then
mapping, with separate completion
checkpoints. Resubmitting jobs 17, 43 and 88 would inspect/reuse their successful
stages and recompute only missing/failed stages in those three iterations.
It would not launch or recompute the other 97 iterations. No monolithic
100-iteration fitting job exists. Benchmark and production aggregation use
`final/benchmark/` and `final/production/` respectively. Do not use the older
benchmark submitter's aggregate command for production: it omits `--production`.

Production `aggregate.sbatch` and `merge.sbatch` are saved in the operations
directory recorded in `PRODUCTION_SUBMISSION.json`. They explicitly select
`--production`, consume all 100 saved memberships, and can be submitted separately
again with `sbatch --account=parent0 --partition=standard --cpus-per-task=8
--mem=150G --time=2-00:00:00 SCRIPT_PATH`. Run merge after aggregation completes;
if only final DE fails, rerun only merge. If an array element fails, the queued
aggregation detects its missing checkpoint rather than averaging fewer iterations.
After repairing that element, resubmit any failed/invalid-dependent final jobs;
successful iterations remain untouched. Completed task replays only verify/reuse
their sealed fit and mapping stages.

## Assets, plots and resource report

The normalized H5AD remains a shared input in the initial benchmark package:
X is float64 sparse ln(1+CPM), raw counts remain in approved Step02, and `uns`
and slot inventories are explicit. It is not copied per iteration. Labels,
models and membership blocks are the reusable iteration assets. No full
446,349×446,349 matrix or final annotated AnnData is produced by this benchmark.

`report/RESOURCE_REPORT.md` and JSON report measured complete-iteration wall
time, CPU time, peak process RSS, disk per completed iteration, projected
100-iteration storage/core-hours/wall time, and recommended CPU/RAM/concurrency.
`job_metrics/`, scheduler logs and `/usr/bin/time -v` provide individual-job
metrics, including failed attempts. `RESTART_VALIDATION.json` records unchanged
fit/output seals and the actual failed/retried mapping job IDs.

`report/benchmark_review.pdf` shows resource bars, a diagnostic sample-fraction
heatmap and an expression-mean display dendrogram from the two-iteration
endpoint. It is explicitly not a final taxonomy or proof of stability. Resource
projections must distinguish measured two-iteration operations from unknown
100-iteration consensus convergence, queue time and filesystem contention.

**Production is already explicitly authorized.** The benchmark report continues
to supply measurements; its original approval wording is reconciled by the
separate operational report-status job, without modifying scientific checkpoints.
No running job is resized in response to later resource recommendations.

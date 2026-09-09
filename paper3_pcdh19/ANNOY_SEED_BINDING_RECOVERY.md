# Exact frozen Annoy seed compatibility recovery

The initial production array `60614879` failed in iterations 2, 3, 4, 5, 6, 8,
9 at root Annoy initialization, after PCA. This was an integer-conversion error,
not an out-of-memory failure. Failed attempts and their PCA/log assets remain.

Annoy 1.17.3's Python `set_seed` parses a signed 32-bit integer. Fifty-three of
our already-frozen NumPy uint32 fit seeds exceed that limit. Its native
`AnnoyIndexInterface::set_seed` accepts uint64, and the Euclidean index uses
Kiss64Random. Wrapping seeds to negative int32 or reducing them modulo a smaller
range would change their native values; neither is done.

Source: [Annoy v1.17.3 binding](https://github.com/spotify/annoy/blob/75429e5dc930754698f1d37c44ea189a7521c7a3/src/annoymodule.cc),
[native interface](https://github.com/spotify/annoy/blob/75429e5dc930754698f1d37c44ea189a7521c7a3/src/annoylib.h),
[RNG](https://github.com/spotify/annoy/blob/75429e5dc930754698f1d37c44ea189a7521c7a3/src/kissrandom.h).

## Narrow runtime change

[seed64.cpp](scripts/annoy_seed64/seed64.cpp) validates the exact installed Annoy
object type/layout, then calls its existing native virtual setter with the
unchanged uint64 value. It compiles no distance, tree-building or RNG algorithm.
[worker.py](scripts/annoy_seed64/worker.py) verifies runtime and installed-binary
hashes, wraps only Allen's Annoy constructor and otherwise invokes the unchanged
frozen HiCAT worker. Representable seeds still use the original Python setter.
No installed package, scientific adapter file, sample list or parameter is edited.

Runtime package: `production/annoy_seed64_v1/` inside the production run. Its
`MANIFEST.json` freezes source, headers, compiled bridge, native Annoy binary
hash and build command. Manifest SHA256:
`0c70208b8ed376d98d161c69d404770c901243da2aeca3499f3e90a4511d9d62`.
Each new stage attempt saves `annoy_seed_binding_provenance.json` inside its
sealed artifact manifest. Existing valid stage contracts remain reusable.

## Qualification and submission

[The numerical qualification](scripts/annoy_seed64/qualify.py) passed all 100
frozen fit seeds: bridge-built indexes and KNN results matched direct calls to
the native setter; representable seeds therefore also matched the original
Python binding. The actual Allen on-disk graph/Jaccard/Louvain path produced
identical safe-seed graphs and labels. A high-seed eight-process path passed.
The frozen production graph worker count remains one; eight requested CPUs
also cover BLAS work. A separate checkpoint fixture confirmed provenance is
sealed and the second call reuses completed output without computing again.
These are implementation fixtures, not new biological parameter combinations.

Sign-off verification at 2026-09-09 00:53:45 UTC additionally rehashed **364**
frozen adapter/config/iteration files, including the seed schedule and saved
sample selections. All matched. After normal scientific imports, installing
the bridge and building/querying an index with exact seed 3245997830 left both
NumPy's global RNG state and Python's `random` state unchanged. Evidence:
`production/annoy_seed64_rng_scope_verification.json`. The sampled/held-out cell
sets are unchanged; final cluster assignments remain outputs to be computed.
The first recovery-array task was still pending, so this is source/fixture and
input-integrity confirmation, not a claim of completed full-data execution.
Production startup verifies the frozen hashes and saves bridge provenance
automatically. No production job was changed during this sign-off check.

Saved evidence: `production/annoy_seed64_qualification.json`,
`production/annoy_seed64_checkpoint_qualification.json` and the recovery
operations directory. The bridge has not yet completed a full-data iteration.

Replacement array **60617037**, `0-99%10`, keeps **8 CPUs / 150 GiB / 48 hours**.
The eight running original-array fits and two benchmark fits continue unchanged.
Only pending superseded tasks and final-stage wrappers were cancelled. The
replacement waits for the original writers, then checks all 100 iterations and
reuses every valid completed stage; this does not create a second ensemble.
Seven failed fits retry individually, as do any other incomplete iterations.
Completed PCA in an unsuccessful fit is retained for inspection; the current
restart unit is the completed fit, so that failed fit restarts from its start.

Use `bin/consensus_production_status.py --run-dir RUN` to list missing/failed
iterations and fresh `sbatch --array=INDEX%10 ...` retry commands. Use the
reported combined retry array to resubmit all failed indices at concurrency 10.
Fresh job IDs allow stale locks from killed tasks to be identified correctly;
do not use same-job-ID requeue after a killed stage.
Final aggregation **60617041** and DE **60617044** are separate jobs using
saved iteration outputs. A final-stage retry never refits completed iterations.
No BigCAT implementation or scientific setting was substituted.

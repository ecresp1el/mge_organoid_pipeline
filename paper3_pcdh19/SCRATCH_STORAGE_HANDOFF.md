# Shared-scratch migration handoff — 4 November 2026

The user requested scratch-based execution and persistent project memory for
the changed paths. This supersedes prior job IDs/storage instructions for new
execution, without rewriting historical records.

- Destination: `/scratch/parent_root/parent0/elcrespo/mge_organoid_pipeline/paper3_pcdh19/results/rebuild_20261103`.
- Preserved source/archive: `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/rebuild_20261103`.
- Migration **65352824**: bounded depth12, no symlink traversal; copies the rebuild, compares SHA256 of every copied file, then writes `SCRATCH_MIGRATION_SUCCESS.json`. Originals retained.
- HiCAT **65352825**: depends on migration, indices0–99, **20 maximum concurrent**, existing sealed fits skipped.
- Consensus/DE **65352826**: depends on successful array.
- Prior temporary recovery 65352448/65352449 cancelled to stop writers before copying. This does not erase completed fits.

This is a submitted migration, not a claim of completed transfer or resumed
computation. Inspect the migration seal and jobs before reporting completion.
`SCRATCH_LOCATION.json` is the new location/job pointer. Earlier `JOBS.json`
inside the copied run is historical until reconciled with this pointer.

## Path behavior

Run-relative inputs, outputs, iteration checkpoints, progress, and final
consensus now resolve under the scratch `HICAT_CYCLE_RUN`. Shared scratch job
temporary directories are created beneath that run and removed on normal exit.
The reference import root is the copied scratch baseline reference. The entire
rebuild is copied so sibling comparison/baseline assets are available.

Frozen configuration and provenance retain their original absolute paths.
The renderer currently reads the configured comparison from Turbo; raw counts,
installed Python environments and reference-hash validation may also read
existing preserved paths. Those remain valid intentionally; no source deletion
or global search-and-replace occurs. Execution/results/temp writes are scratch
based; this is not a fully self-contained relocation of all software/raw data.
Any future removal of Turbo source requires explicit path/schema migration and
revalidation first. Original FROZEN.json remains byte-identical.

ARC documents scratch per-user/account directories, 10TB and one million files
per root account by default, no backups, and purge after 60 days without use.
Validated final outputs must be copied to durable Turbo storage; do not rely
on scratch as the only long-term copy.

Sources:
- https://documentation.its.umich.edu/arc-hpc/greatlakes/user-guide/defaults-limits
- https://documentation.its.umich.edu/arc-hpc/greatlakes/policies


## Restart correction — 4 November 2026, 10:05 AM EST

Scratch migration65352824 completed at 09:23:57 EST; all908 files checksum
verified. Array65352825 then failed on existing interrupted progress logs:
`Progress` refuses to overwrite `hicat_progress_events.jsonl`. No new fits
completed;32 sealed fits preserved. Failed/interrupted iteration directories
are moved intact to `hicat_min500_cycle/interrupted_attempts_before_retry_20261104/`
by preparation65353943 after old writers stop. Old array/consensus cancelled.
Qualification iteration032:65353944. Remaining033–099:65353945, max20,
dependent on successful qualification. Final consensus:65353946. This gate
prevents repeating the failure across the full array before one full fit passes.
No scientific settings changed. Scratch `SCRATCH_LOCATION.json` and HiCAT
`JOBS.json` record current IDs. Prior IDs above are historical. Completion
ETA withdrawn pending successful qualification and resumed throughput.


### Immediate parallel release — 4 November 2026

Removed single-fit dependency from array65353945 at user request.
Iteration032 (65353944) continues independently;033–099 may run up to20
at once (up to21 total while032 is active). Final65353946 explicitly awaits
both65353944 and65353945. No cancelled/restarted numerical work in this change.
Scratch outputs and sealed fits preserved. Current numerical progress already
passed original scratch path and progress-file errors. See runtime
`PARALLEL_RELEASE_20261104.json`.


### User concurrency correction — 4 November 2026

User reports prior testing supports only **10 jobs simultaneously**. This
supersedes requests/records for20. Array65353945 throttle reduced to**9**,
plus independent iteration032/65353944 gives at most10. No single-fit
dependency on array; jobs eligible immediately. Keep total active HiCAT fits
<=10. After032 completes, array may be raised to10 after verification.

Scheduler had already launched20 array tasks before lowering throttle. Stopped
excess tasks[42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52], keeping9 plus032. Displaced tasks queued as
array65358814 after the original array ends and032 succeeds, max10.
Final65353946 dependencies adjusted for intentional cancellations; it still
requires all100 SUCCESS seals and label hashes. See USER_CONCURRENCY_LIMIT_20261104.json.


At12:37 EST:36 completed; independent032 succeeded. Array65353945
throttle increased9→10 so total active fit limit remains10. No new failures.

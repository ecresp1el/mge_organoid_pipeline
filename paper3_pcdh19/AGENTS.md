# PCDH19 active-run memory and storage contract

Updated 2026-11-04. Read `SCRATCH_STORAGE_HANDOFF.md` before submitting jobs or
reporting status. User requested scratch-based execution after temporary files
filled home storage. Do not inherit Codex TMPDIR for scientific jobs.

Active rebuild destination:
`/scratch/parent_root/parent0/elcrespo/mge_organoid_pipeline/paper3_pcdh19/results/rebuild_20261103`

Migration job 65352824 must finish and `SCRATCH_MIGRATION_SUCCESS.json` must
exist before this destination is treated as a completed migration. HiCAT array
65352825 depends on migration; consensus 65352826 depends on that array.
Concurrency maximum is 20, not four. Check actual running counts separately.
The array covers 0–99 and skips existing sealed SUCCESS.json checkpoints.

Turbo source is retained read-only as source/archive. Original counts, software
environments, and frozen absolute reference/provenance paths may intentionally
point there; do not mass-replace strings or modify frozen source/hash manifests.
Do not move directories beneath running writers. Scratch is not an archive.
Preserve durable final outputs on Turbo after successful validation.

Use explicit HICAT_CYCLE_RUN and HICAT_CYCLE_REFERENCE, unique scratch job temp
directories, and exit cleanup limited to the current job's temporary directory.
Never delete the shared run root in a cleanup trap. Never recursively search
scratch allocations or their parents; bound work to this exact project subtree.
No background polling. Scheduler checks no more than once per minute.


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

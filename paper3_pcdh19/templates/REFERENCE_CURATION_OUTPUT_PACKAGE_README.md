# Developing-mouse MGE reference-curation checkpoint package

This is a run-scoped package for the independent audit of La Manno 2021,
Bandler 2022, and Mayer 2018 as candidate WT developing-mouse MGE references.
It stops at the required early processed-object checkpoint. It does not map
Paper 3 cells, integrate references, reprocess raw sequencing, recluster large
objects, invent annotations, or choose the final reference.

Directory roles:

- `code/`: exact Python, R, shell, and SLURM files frozen before submission.
  The SLURM jobs execute these copies, not the repository files.
- `config/`: submitted source registry, environment/resource settings, and
  resolved run paths.
- `LaManno2021/`, `Bandler2022/`, `Mayer2018/`: per-study metadata and audit
  outputs. Large downloaded P0 objects are not duplicated here.
- `tables/`: source audit, raw-access registry, and combined early checkpoint.
- `logs/`: scheduler output and errors.
- `provenance/`: submitted command, job IDs, and runtime environment records.
- `SUCCESS.txt` or `FAILED.txt`: unambiguous checkpoint status.

Original P0 objects are cached once under the resolved
`PAPER3_CURATION_SOURCE_ROOT` inside `PAPER3_ROOT/inputs/`. Existing cached
objects are reused without being overwritten. A new default submission creates
a new versioned run. An existing inactive run can be regenerated only through
the guarded `--replace-run RUN_ID` submission option.

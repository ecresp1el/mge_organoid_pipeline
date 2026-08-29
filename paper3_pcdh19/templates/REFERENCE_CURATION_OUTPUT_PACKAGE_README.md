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
- `LaManno2021/`, `Bandler2022/`, `Mayer2018/`: per-study published metadata,
  object-linked sample inventory, standardized sample summary, and object
  audit outputs. Large downloaded P0 objects are not duplicated here.
- `tables/`: source/raw-access audit, all-candidate per-sample inventory,
  study-level sample/library summary, sample-field data dictionary, curation
  requirements ledger, combined early object checkpoint, and explicit
  candidate E15/MGE evidence status.
- `figures/`: cross-study E15/MGE evidence visual. La Manno's study-specific
  directory also contains the observed author annotation hierarchy/composition
  plot and its underlying TSVs.
- `logs/`: scheduler output and errors.
- `provenance/`: submitted command, job IDs, and runtime environment records.
- `SUCCESS.txt` or `FAILED.txt`: unambiguous checkpoint status.
- `REFERENCE_CURATION_REPORT.md`: observed-object comparison generated only
  after all three inspections and the combined checkpoint pass; it embeds the
  hierarchy/evidence visuals when they are present.

Original P0 objects are cached once under the resolved
`PAPER3_CURATION_SOURCE_ROOT` inside `PAPER3_ROOT/inputs/`. Existing cached
objects are reused without being overwritten. A new default submission creates
a new versioned run. An existing inactive run can be regenerated only through
the guarded `--replace-run RUN_ID` submission option.

Published sample metadata and membership in a downloaded P0 object are
different evidence states in this package. Blank, unresolved, partial, and
`NOT_ASSESSED` values are intentional; they prevent paper-level facts from
being presented as cell-level joins before those joins are proven.

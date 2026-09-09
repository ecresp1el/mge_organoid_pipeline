# HiCAT 03: saved-consensus review and publication

The requested deliverable is the consensus of completed iterations 0–97,
with denominator 98. This is a review of the existing full-data consensus work
within **HiCAT stage `03_full_data_consensus_benchmark`**. It does not create
primary-processing Step 08, rename the rejected primary Step 03, or advance
an unapproved analysis to a new major checkpoint.

The numbering follows the HiCAT sequence: `00_method_and_technical_pilot`,
`01_coarse_fine_expanded_pilot`, `02_hierarchy_validation`, then
`03_full_data_consensus_benchmark`. Stage 02 corresponds to primary Step 07;
that documented relationship must not be confused with primary Step 03.

The source computation is preserved at
`consensus_restart_20260909_000256/final/consensus_98_20260909/` under stage 03.
Its saved aggregation/refinement and final-DE outputs supply this reporting
stage. No clustering, normalization, aggregation or DE is rerun. Existing
100-iteration jobs remain independent of the requested 98-run result.

## Architecture and publication contract

- Canonical analysis classes live in `scripts/hicat/consensus_review_*.py`:
  immutable settings and verified readers in `consensus_review_io.py`, plotting
  in `consensus_review_report.py`, and workflow/atomic publication in
  `consensus_review_workflow.py`.
- `consensus_review_submission.py` freezes a new versioned package;
  `bin/submit_hicat_consensus_review.py` and
  `slurm/hicat_consensus_review.sbatch` are thin entry points.
- Runs use `PAPER3_ROOT/results/hicat/03_full_data_consensus_benchmark/` with
  IDs `hicat_consensus_review_<UTC timestamp>_<config hash>`.
- `code/`, `config/`, `upstream/`, `inputs/` and `provenance/` preserve executable
  source, environment, input identity, method history and the authorization.
  SLURM runs frozen copies. Jobs live in `PAPER3_ROOT/jobs/`, scheduler logs in
  `PAPER3_ROOT/logs/hicat/<RUN_ID>/`, linked from the run's `logs/`.
- Output is built in `staging/outputs/`, validated and manifested, then renamed
  atomically to `outputs/`. Computation leaves `IN_REVIEW`, never `APPROVED`.
- No scientific results are published inside the GitHub checkout.

## Validation and scope

Verify frozen source/configuration/environment; pin source identities; fully
hash and reopen every saved aggregation and DE artifact through the existing
checkpoint reader. Check parent seals, R=98, exact iteration IDs 0–97, cell
order/uniqueness, 446,349 cells, 19,071 genes, 12 samples, gene/cluster alignment,
and exact published table copies. Hash the completed publication and retain
its validation ledger. Reporting does not need to reread original per-iteration
PCA/graph models. The earlier consensus job's focused validation scope remains
explicit in provenance; publication does not rewrite that history.

## Review assets and figure numbering

`outputs/figures/hicat_consensus_review.pdf` contains numbered pages/PNGs
using the existing `HierarchyReport.save` convention:

1. `001_consensus_composition.png`: panels A–B, sizes and sample fractions.
2. `002_consensus_marker_expression.png`: panel C, the declared display-only
   marker panel, gene-wise z scores across cluster means.
3. `003_consensus_expression_dendrogram.png`: panel D, average linkage with
   correlation distance across all-gene cluster means; a display summary,
   not the recursive fitting tree.

`figure_index.tsv` maps each page to exact saved tables. `outputs/tables/`
contains assignments, full pairwise DE evidence, all-gene means/detection,
consensus diagnostics, plotted matrices, marker IDs/symbols, linkage and labels.
`validation/validation_checks.tsv`, `output_manifest.tsv`, `STEP_STATUS.json`
and `REVIEW_README.md` complete the review package. No new H5AD, UMAP, ranked
marker analysis, cell-type assignment or accepted taxonomy is implied.

This report is a same-stage publication correction of the requested result.
The ad hoc source run and its frozen script remain as historical provenance;
the canonical workflow entry point is now the stage-03 review launcher.

# Step 07: HiCAT hierarchy validation and parameter sensitivity

This checkpoint studies **dissected E14.5 mouse MGE**. It tests the existing
12,000-cell pilot (1,000 cells from each of 12 samples). Provisional annotation
text supplied by the user is a hypothesis reference, never a fitting input or
ground truth. Success means a complete review package marked `IN_REVIEW`.

## Execution and input contract

The frozen pilot is `hicat_coarse_fine_20260908_203708_7e050f23`, with 19,071
genes, four baseline coarse parents and 38 baseline fine clusters. Its repeat
seed yielded five coarse and 39 fine clusters. The exact H5AD hash, source
paths, cell counts, seed and bounded comparison are in
[`config/hicat_validation.json`](config/hicat_validation.json).

1. Verify source hashes, all cell/gene IDs, sample balance, raw integer counts,
   normalized-expression formula, stored clustering configuration, and display
   coordinates. Reuse both saved seed fits; do not recompute them.
2. Fit canonical program score controls on pilot expression independently of
   cluster labels and DEGs. Save every signature, mapping, missing gene, control
   gene and weight. Score cell cycle without regression. Assess biological
   combinations and alternatives appropriate to E14.5, including non-neural
   vascular and erythroid possibilities.
3. Quantify seed overlap, sample representation and phase-stratified identity
   differences. These are descriptive diagnostics: a seed repeat on identical
   cells is not a bootstrap, sample imbalance is not necessarily batch, and
   phase stratification cannot prove a causal cell-cycle effect.
4. Run one additional fine fit using the **same baseline coarse parents** and
   baseline seed 20260908. Change only fine `q1_thresh` 0.3→0.4 and
   `qdiff_thresh` 0.5→0.7; retain DE score 150 and all other engine controls.
   Candidate IDs use `.A` to distinguish this partition from original `.F`
   labels. Candidate fitting/merging does not modify original labels.
5. Apply the frozen pilot scoring model to the full 446,349-cell raw-count
   object in chunks. Display genes/programs on the existing Step 06 UMAP.
   Compare pilot versus full marker distributions and sampling coverage.
6. Generate the reports, all 38 individual fine-cluster review pages, source
   tables, machine summaries and saved-object checks. Publish atomically only
   after checks pass. Leave biological labels unlocked.

The Python engine retains the previous pilot's pinned Allen commit
`99154957c74023235763025fda9dc4eb3ca952c6`. See the
[Allen Python source](https://github.com/AllenInstitute/transcriptomic_clustering/tree/dev)
and the [scrattch tutorial](https://taxonomy.shinyapps.io/scrattch_tutorial/#section-overview).
The controlled settings are the user's specified Allen-reference comparison;
they do not establish a universal optimum or a predetermined cluster count.
The earlier [adapter protocol](HICAT_COARSE_FINE_PROTOCOL.md) still explains
the Python orchestration and differences from upstream convenience functions.

## Normalization, saved objects and exclusions

The pilot already contains `layers['log1p_cpm'] = ln(1 + raw / cell_total × 1e6)`.
Step 07 verifies and reuses it. `X` remains raw sparse integer counts. Full-data
display computes the same formula using totals across all 19,071 genes before
selecting canonical genes. Neither scoring nor plotting changes raw counts.

The output root is `results/hicat/02_hierarchy_validation/<RUN_ID>/outputs/`.

| Asset | Saved content and intended use |
|---|---|
| `pcdh19_hicat_hierarchy_validation.h5ad` | All 12,000 pilot cells × 19,071 genes; original raw X, normalized layer, existing labels/coordinates/models, separate Allen candidate labels, program scores and cell-cycle diagnostic fields. Review categories/hypotheses remain external tables. |
| `anndata_uns_inventory.json` | Complete readable inventory of saved `uns`, including original HiCAT provenance and Step 07 scoring/configuration/asset references. No hidden locked annotations. |
| `marker_programs/` | Canonical gene sets, sources, coverage, exact controls/weights, per-cell absolute scores and expression; coarse/fine mean, median and positive fractions in report source tables. |
| `stability/` | Both-level seed raw and normalized overlaps, per-cluster matches and global ARI/NMI. |
| `sample_composition/` | Cluster-by-sample counts, both fraction directions, expected-contribution enrichment and entropy/effective-sample diagnostics. |
| `cell_cycle/` | Per-cell scores/phases, phase composition and every within-parent phase-matched identity comparison. |
| `parameter_sensitivity/` | Exact config diff, candidate per-cell assignments, parent counts, all candidate engine models/merge and DE evidence, comparison overlaps and boundary-survival evidence. |
| `full_data_projection/` | Selected-canonical-gene normalized H5AD with existing UMAP, marker scores and pilot membership; per-cell source tables and coverage summaries. This is a display/scoring sidecar, not a full-gene HiCAT input. |
| `annotation_review/` and `tables/` | Hypothesis comparisons, fine and coarse review categories, technical/state relationships, developmental-order evidence and numerical plot sources. |
| `figures/` | Main and detailed PDFs, individual PNG/PDF panels, explicit figure/source index. Full-data vector panels rasterize point clouds for manageable files. |
| `validation/` | Input, partition, saved-object and scope checks. |
| `REVIEW_README.md`, `summary.json` | Direct answers to the requested review questions, limitations, readiness recommendation and final coarse/fine review-category table. |

Not produced: a full-data HiCAT partition, inferred full-data cluster labels,
an integrated/corrected expression matrix, regressed counts, removed cells,
locked annotations, a consensus clustering, a target K, or a newly fitted UMAP.
The original full raw-count object remains at its verified input path. The
full-data sidecar does not duplicate all 19,071 raw genes. Existing full marker
statistics remain in the earlier immutable pilot; Step 07 links them and copies
the top-marker source table for review. Fewer than 20 qualifying genes is
reported faithfully; C0001.F0011 has 197 cells and complete statistics but no
qualifying positive marker genes in those prior reporting contexts.

## Parameters to inspect and tune in a later checkpoint

`hicat_validation.json` controls the exact source identity, one DE comparison,
sample flags, seed-match thresholds, phase minimum cell counts and plot DPI.
`hicat_validation_programs.json` controls canonical sets, expression-bin
controls, positivity rules and sources. `hicat_validation_hypotheses.*`
contains review-only user hypotheses. Every run freezes all three; changing a
live config does not change a submitted job. Additional parameter tuning would
require a new versioned checkpoint, rather than silently changing this run.

High seed stability requires Jaccard ≥0.75 and precision/recall ≥0.8; moderate
requires Jaccard ≥0.5. Allen boundary survival is strong under the high rule,
partial for Jaccard ≥0.35 and recall ≥0.5, otherwise no. These are overlap-based
review rules, not biological truth. A target receiving >10% of a source group
counts toward fragmentation. Split/merge evidence retains source and target
IDs and fractions. Sample flags use maximum contribution >0.25 or effective
sample count <6 for clusters with at least 30 cells; small groups retain
descriptive values without treating a noisy fraction as decisive.

Phase-matched identity comparisons require at least 20 cells from each group
within a phase and a standardized program effect of at least 0.5. Canonical
identity scores exclude explicit cycle genes. Report direction, effect size,
phase cell counts and the number of evaluable phases; insufficient phase
coverage never becomes evidence of no identity difference. Categorical phase
scores are expression diagnostics, not direct measurements of cell-cycle time.
The final report records additional review-category rules and limitations.

## Code and operational trace

`bin/submit_hicat_validation.py` freezes code/config/tests/reference sources,
environment and repository state, then submits a bounded SLURM job.
`scripts/hicat/validation_workflow.py` orchestrates verification, the independent
scoring/metrics modules, `FixedParentSensitivity`, full-data display scoring,
review and reporting. Functions document their representations and outputs.
`provenance/hicat_progress_events.jsonl` and scheduler logs contain durable
START/COMPLETE/FAILED messages, parameters, dimensions and elapsed time.
Failed runs retain evidence and never receive an `IN_REVIEW` success record.

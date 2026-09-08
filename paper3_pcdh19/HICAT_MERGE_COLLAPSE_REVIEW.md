# Step 07 extension: explain the existing 38 → 17 comparison

**Completed; IN_REVIEW.** Final decision-table job **60601782** passed, including
630 reconstructed phase effects matching saved Step 07 values. The final
recommendations are five `keep distinction`, three `merge`, and 76 `unresolved`
pairwise boundaries. No new HiCAT fits or clustering parameters were run.

Start with the [concise searchable decision table](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/02_hierarchy_validation/hicat_validation_20260908_221625_d5d06a5d/review/merge_collapse_20260908_225259/decision_review_v2/outputs/decision_table.html)
and [interpretation with exact decision blockers](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/02_hierarchy_validation/hicat_validation_20260908_221625_d5d06a5d/review/merge_collapse_20260908_225259/REVIEW_README.md).
The [merge-flow tree](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/02_hierarchy_validation/hicat_validation_20260908_221625_d5d06a5d/review/merge_collapse_20260908_225259/outputs/merge_tree_38_to_17.pdf)
and [all 17 contribution summaries](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/02_hierarchy_validation/hicat_validation_20260908_221625_d5d06a5d/review/merge_collapse_20260908_225259/outputs/comparison_cluster_summary.tsv)
are saved alongside it. Overall judgment: **still unresolved; do not lock 17**.

`scripts/hicat/merge_decision_table.py` adds destination-specific program/marker
summaries and matched phases from already saved cell scores. It measures repeat
boundary separation independently of whole-cluster seed stability. The final
11-column table supersedes the broader original-group association screen for
recommendations; earlier versions remain marked as superseded. All methods and
operational review thresholds are documented in the final table's expandable
methods section. They do not change the clustering.

This review uses dissected E14.5 mouse MGE and the existing 12,000-cell pilot.
It introduces **no new HiCAT fit or parameter combination**. The original 38
labels and comparison 17 labels remain unchanged. No annotations are adopted.

The comparison is a refit under q1=0.4/qdiff=0.7 with fixed coarse parents, not
a pure agglomeration of the original 38 labels. An original cluster may split.
The requested merge tree therefore displays exact cell flows, including splits,
rather than suggesting an internal chronological merge history.

## Inputs and provenance

`bin/submit_merge_collapse_review.py --source-run <Step07 run>` copies completed
numeric evidence into a new version under `<Step07 run>/review/merge_collapse_*`.
It freezes the Python code, review thresholds, environment, source paths and
input hashes before submitting a small read-only analysis job. All output
checks remain IN_REVIEW. The parent report job can finish independently.

Inputs are saved memberships, the 38×17 overlap table, canonical expression,
phase-matched program effects, phase-pair summaries, sample fractions **within
each cluster**, seed diagnostics, provisional hypothesis review, canonical gene
mapping and the original final pairwise DE audit. No source expression matrix
is modified or newly normalized. Baseline H5AD access extracts metadata only.

The job independently reconstructs the overlap table from exact cell IDs,
checks 12,000 unique cells and fixed parent nesting, and requires each cluster's
sample fractions to sum to one. Input and output byte hashes are recorded.

## Outputs

- `all_cell_contributions.tsv`: every nonzero source/destination flow; count,
  fraction of original cluster and fraction of comparison cluster.
- `comparison_cluster_summary.tsv`: all 17 destinations and exact contributors.
- `cell_memberships.tsv.gz`: all 12,000 original/comparison/sample assignments.
- `merge_tree_38_to_17.pdf` and parent PNGs: cell-count-proportional flow diagrams.
- `lost_boundaries.tsv`: one row per lost/eroded original sibling boundary;
  `all_original_sibling_boundaries.tsv` retains all 291 pairs, including minor
  leakage and retained boundaries. Multiple pairwise boundaries may belong to
  one merger; their number is not simply 38 minus 17.
- `canonical_boundary_gene_contrasts.tsv`, cluster means/detection tables, and
  `matched_phase_program_evidence.tsv`: numerical biological support, including
  contrasts that are too weak or insufficiently sampled to classify.
- `decision_blocking_boundaries.tsv` and `merges_resolving_weak_subdivisions.tsv`:
  operational review screens requiring interpretation, never automatic merging.
- `README.md`, `summary.json`, `checks.json`, `output_manifest.tsv`: explanation,
  disposition, data-integrity checks and exact saved-file inventory.

## Review parameters, not clustering parameters

All cutoffs are explicit in each package's `review_config.json` and explained
in its README. A common destination receiving ≥90% of both originals defines
near-complete coalescence; ≥50% each defines majority coalescence; >10% each
defines partial erosion. Every smaller contribution remains in the full tables.

Canonical regional/developmental contrasts combine fixed-pilot-SD program
contrasts, independently specified genes and matched phase evidence. Missing
phase coverage is not a negative biological result. Cell-cycle and sample
associations are descriptive, not causal. Multiple axes can coexist; the
predominant-association column is a conservative review synthesis and does not
replace the numerical evidence or establish an anatomical lineage.

Original DE scores establish original separation under the original criteria.
They do not prove that a specific changed threshold caused a specific merger
in the comparison refit. No new original-pair DE tests under altered thresholds
are performed. No final cluster count is chosen by this extension.

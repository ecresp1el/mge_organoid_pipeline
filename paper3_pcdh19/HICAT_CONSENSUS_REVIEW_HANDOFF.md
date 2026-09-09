# HiCAT 03: current 98-iteration consensus review

**IN_REVIEW.** Job **60712671** published the corrected review package on
September 9, 2026 at approximately **10:24:31 AM Detroit**. All **35
checks passed**. It reuses the completed R=98 consensus from job **60707938**,
which finished at **10:06:56 AM**, elapsed **6 min 12 s**. No clustering,
aggregation, refinement or final DE was recomputed for this publication.

**Result: 34 clusters, all 446,349 cells, 19,071 genes and 12 samples.**
The consensus uses exactly frozen iterations 0–97 and denominator 98.
Aggregation/refinement gave 35 groups; final DE gave 34. Three of 561 audited
cluster pairs remain below the configured separation criteria. Cluster IDs
are not biological annotations or an approved final taxonomy.

## Review entry points

- [Review README](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/03_full_data_consensus_benchmark/hicat_consensus_review_20260909_142348_efd9e13d/outputs/REVIEW_README.md)
- [Numbered review PDF](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/03_full_data_consensus_benchmark/hicat_consensus_review_20260909_142348_efd9e13d/outputs/figures/hicat_consensus_review.pdf)
- [001: cluster sizes and sample composition](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/03_full_data_consensus_benchmark/hicat_consensus_review_20260909_142348_efd9e13d/outputs/figures/001_consensus_composition.png)
- [002: marker expression](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/03_full_data_consensus_benchmark/hicat_consensus_review_20260909_142348_efd9e13d/outputs/figures/002_consensus_marker_expression.png)
- [003: expression display dendrogram](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/03_full_data_consensus_benchmark/hicat_consensus_review_20260909_142348_efd9e13d/outputs/figures/003_consensus_expression_dendrogram.png)
- [Cell-level assignments](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/03_full_data_consensus_benchmark/hicat_consensus_review_20260909_142348_efd9e13d/outputs/tables/all_cell_assignments.tsv.gz)
- [Figure/source index](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/03_full_data_consensus_benchmark/hicat_consensus_review_20260909_142348_efd9e13d/outputs/figures/figure_index.tsv)
- [Full run package](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/hicat/03_full_data_consensus_benchmark/hicat_consensus_review_20260909_142348_efd9e13d)

This is a versioned review within **HiCAT stage 03**, following stages 00
(method pilot), 01 (expanded pilot) and 02 (hierarchy validation). The separate
primary-processing numbering remains unchanged. See the
[review protocol and architecture](HICAT_CONSENSUS_REVIEW_PROTOCOL.md).

Source and new output checks include full saved aggregation/DE artifact hashes
and reopening, exact cell/gene alignment, source identities, parent seals,
R=98, sample count, unchanged copied tables and figure/source availability.
Validation details and the output manifest are under `outputs/`.
No new H5AD, UMAP, ranked marker analysis, correction or annotation is claimed.

## Correction of the earlier publication

The unstructured result remains preserved at
`consensus_restart_20260909_000256/final/consensus_98_20260909/` as the source
computation. The ad hoc launcher and the two result folders created in the
GitHub checkout were relocated, with SHA256 verification, to this review run's
`provenance/replaced_repository_artifacts/`. They are historical provenance,
not the canonical publication or workflow entry point. No source computation
was deleted or silently rewritten.

Scientific code now follows the `scripts/hicat/consensus_review_*.py` classes;
the thin canonical launcher is `bin/submit_hicat_consensus_review.py`.
Frozen executed copies, inputs, configuration and environment are in this run.
The scheduler script is in `PAPER3_ROOT/jobs/`; scheduler logs are in
`PAPER3_ROOT/logs/hicat/hicat_consensus_review_20260909_142348_efd9e13d/`, linked from the run.
The HiCAT approval ledger records this run as `IN_REVIEW`.

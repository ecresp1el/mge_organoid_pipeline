# CODEX HANDOFF: FIRST E15 MGE MAPPING DIAGNOSTIC

## Purpose

Test two fast annotation approaches on the same Paper 3 cells while preserving
the existing analysis. This is downstream of reference curation and independent
of PCDH19 genotype classification.

## Immutable query contract

- 12 Cell Ranger samples, 450,788 total cells, 19,071 Ensembl features.
- Stable query ID: `<technical_sample_id>_<10x_barcode>`.
- Original raw UMI counts only; no integration or transformed matrix is used as
  the MapMyCells input.
- Original Cell Ranger graph clusters are retained as
  `<sample_id>:<graph_cluster>` because cluster numbers were computed separately
  and are not globally comparable.
- Original Cell Ranger UMAP coordinates are retained. They were also computed
  separately, so every figure must facet by sample and must not imply a common
  coordinate space.
- No reclustering, integration, manual annotation, or UMAP recomputation.

## Branch 1: Bandler/MIND

- Reference: the 4,481 definitively barcode-resolved cells from CA301 /
  GSM5684876 / WT E15.5 MGE.
- Labels: later-MIND `MIND_class` (Mitotic or Inhibitory Neuron Precursor) and
  all 12 `MIND_cluster` states.
- These later-MIND labels must never be described as the unrecovered original
  Bandler 2022 21-state taxonomy.
- Seurat reference-PCA label transfer; preserve the raw winning label, every
  label score, and maximum score for both hierarchy levels.
- Fixed, predeclared `Unassigned` rule: maximum score < 0.50. Do not optimize
  this cutoff from the results.
- Save a derived Seurat RDS with original counts, existing metadata/UMAP, and
  added MIND fields. Never alter the source delivery.

## Branch 2: MapMyCells

- Pinned `cell_type_mapper` 1.7.4 and official Allen WMB-10X CCN20230722 assets
  dated 20240831.
- Hierarchical, not flat, raw-count mapping; 100 bootstrap iterations, default
  0.5 marker downsampling, five runner-ups, deterministic seed 20260829.
- `CCN20230722_SUPT` is not directly traversed in the official workflow;
  MapMyCells retains it and backfills it from the lower assignment.
- Preserve every returned hierarchy level and all IDs, names, bootstrapping
  probabilities, correlations, aggregate probabilities, direct-assignment
  flags, and runner-up fields. Preserve the official CSV and extended JSON.
- Treat all calls as an adult WMB comparator, not developmental E15 ground truth.

## Concordance/report contract

Publish both a per-cell table and an aggregated table spanning:

`sample-qualified existing cluster × MIND class × MIND 12-state × MIND scores × all MapMyCells hierarchy levels × MapMyCells confidence fields`.

Plot both label sets and confidence on the unchanged, sample-faceted UMAPs.
Explicitly summarize MapMyCells names consistent with microglial, vascular,
oligodendroglial, astroglial, or other obvious broad out-of-reference cells.

## Execution

Run directory:
`/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/01_e15_mge_mapping_diagnostic/01_e15_mge_mapping_diagnostic_20260829_172214_a3f7ca7`

The run is complete and has top-level `SUCCESS.txt`. Successful SLURM stages:

- 59181267: query preparation, `COMPLETED 0:0` in 00:05:03
- 59181383: MapMyCells, `COMPLETED 0:0` in 00:15:18
- 59182039: Bandler/MIND successful retry, `COMPLETED 0:0` in 00:20:39
- 59182287: exact MapMyCells-to-Seurat attachment, concordance, plots, and report,
  `COMPLETED 0:0` in 00:06:52
- 59185031: report-only interpretation republish, `COMPLETED 0:0` in 00:02:20

The initial branch submission requested 192G, exceeding the 180G standard-node
capacity; no branch job was created. The frozen continuation records this and
uses 176G without changing any scientific setting. MIND job 59181382 then
failed before loading data because the Seurat module lacked optional `hdf5r`.
Retry 59181926 failed before loading data because `hdf5r`'s `bit64` dependency
was absent. The final isolated add-on environment contains `hdf5r` 1.3.12,
`bit64` 4.8.4, and `bit` 4.6.0 compiled against Great Lakes HDF5 1.10.8.
Both failed job scripts were archived under their job IDs. Dependency-invalid
report jobs 59181384 and 59181927 were cancelled automatically; pending report
59182040 was intentionally replaced by the combined-object job 59182287.

## Completed diagnostic result

- Query contract: 450,788 cells x 19,071 Ensembl features; all stable cell IDs,
  Cell Ranger graph clusters, and UMAP rows match exactly.
- Bandler/MIND: 332,740 Mitotic and 118,048 Inhibitory Neuron Precursor calls.
  At the 12-state level 427,589 (94.9%) are assigned and 23,199 are
  `Unassigned`; 247,985 (55.0% of all cells) are `Hist1h1b_Top2a`. Eleven of
  the 12 reference states receive assigned query cells. This concentration is
  diagnostic evidence, not a finalized annotation.
- The two-label broad winner is mathematically at least 0.50; its observed
  minimum is 0.500007, so the fixed 0.50 broad cutoff abstains on no cells.
  This limitation is explicit in the report.
- MapMyCells: 446,365 (99.0%) map to adult-taxonomy `05 OB-IMN GABA`. This is
  consistent with a broad immature GABAergic program but is not MGE-specific
  developmental ground truth.
- MapMyCells flags 3,116 explicitly named non-neural cells, including 1,195
  endothelial, 731 pericyte, 208 microglia, 145 OPC, and smaller glial/immune
  groups. The full result spans 25 classes, 114 subclasses, and 397 clusters.
- MapMyCells warned that 593 Allen marker genes were absent from the query and
  ignored them. The exact warning and log are preserved.

Primary outputs are `E15_MGE_MAPPING_DIAGNOSTIC_REPORT.md`,
`tables/per_cell_mapping_concordance.tsv.gz`,
`tables/cluster_label_concordance.tsv.gz`, figures 01--07, and the single
validated 31.84-GB combined object
`Bandler_MIND/paper3_query_with_mind_and_mapmycells_labels.rds`. The combined
object contains the original counts/existing UMAP, all MIND scores/calls, and
all 53 `mmc_*` hierarchy/confidence fields; the temporary MIND-only giant RDS
was removed after validation to avoid duplication.

## Reproducibility and reruns

Each run contains exact copies of every Python, R, shell, and SLURM file that
was actually executed, plus the exact sample key, resolved paths, full pinned
Python environment, and the copied CA301 barcode/label join. Large immutable
source matrices and Allen assets are cached once under `PAPER3_ROOT/inputs` and
recorded rather than duplicated.

Default submission creates a new versioned Step 01 run. Intentional overwriting
uses `submit_e15_mge_mapping_diagnostic.sh --replace-run <RUN_ID>` and is limited
to that exact inactive Step 01 directory; it refuses active jobs and never
deletes shared inputs or other steps.

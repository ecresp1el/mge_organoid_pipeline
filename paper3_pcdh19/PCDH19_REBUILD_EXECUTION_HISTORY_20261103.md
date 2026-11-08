# Historical execution notes — superseded by the current workflow guide

These are chronological records, including failed attempts and earlier data universes. They are not the current Methods recipe.

# PCDH19 Python rebuild — 3 November 2026

## User request and deletion boundary

User requested a clean-slate Python workflow: before/after cell-cycle correction
UMAPs with sample overlays, reclustering, then Allen MapMyCells. They also
requested replacing/removing old outputs to reduce Turbo usage. User explicitly chose to preserve the final archive and essential inputs while
removing replaceable outputs. Cleanup job 65311673 implements that scope, with
pre/post checks of the final archive and a SHA-256 check of the approved QC
input. It preserves raw Step00 and approved Step02, inputs, software, probe
audits and compact baseline UMAP/PCA/HiCAT assignments. All other replacement
analysis results, obsolete full matrices and draft figure packages are removed
only through the recorded exact-directory cleanup manifest.
Code in Git is not the storage target; generated PCDH19 runtime outputs are.

A read-only bounded storage inventory is running under job 65311425. Script:
`scripts/audit_pcdh19_storage.py`. Runtime report destination:
`results/rebuild_20261103/storage_audit/` under the study runtime. The scan
is restricted to the PCDH19 runtime, max depth 12, same filesystem, no symlinks.
It measures file sizes, not remaining allocation quota. A deletion candidate
list must enumerate exact directories and estimated reclaimable storage.

## Confirmed current data lineage

The figure UMAP is Step 06's shared, unintegrated embedding. The analysis code
uses raw-count sample-aware seurat_v3 HVGs (3000), normalize_total(10000),
log1p, HVG scaling, centered PCA (50), cosine neighbors (30), UMAP (min_dist
0.3, spread 1, seed 20260903). There is no S/G2M regression. Sample-aware HVG
selection is not batch correction. See `scripts/primary_processing/step06_analysis.py`.

The figure cluster labels come from a separate Python HiCAT consensus. It
also has no cell-cycle regression. Existing Step 06 UMAP coordinates were
only the display for these labels; HiCAT does not cluster the 2D UMAP.
See `HICAT_PARAMETERS_AND_CONSENSUS.md`. The user explicitly selected Python HiCAT consensus, matching the previous
clustering.
Do not silently describe Leiden as rerunning HiCAT.

Approved reusable input: primary Step 02 run
`02_qc_filtering_20260830_124611_97e1bb5`, checkpoint
`objects/pcdh19_step02_qc_filtered.h5ad`: 446349 cells × 19071 genes,
12 samples; 6554152927 bytes; SHA256
`fadba4a25a7b6b7320219b21c189b6325687493519ba0fe1bd27efd79606b103`.
Root: `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19`.
Python environment is configured in `config/primary_processing_step06.env`.
Rejected scDblFinder output is not a restart input.

## Planned computation and compact storage

1. Preserve one authoritative raw-count input; use the same cell IDs, feature
   selection, seed, normalization, scaling and graph parameters in both arms.
   Score S/G2M using unregressed log-normalized expression and audit mouse
   symbol coverage/aliases. Fit a baseline and an arm regressing S and G2M.
   Retain cell-cycle scores and lineage-marker context to judge biology lost
   as well as cycle dependence reduced. Do not regress genotype or sex.
2. Refit PCA/neighbors/UMAP after regression; a UMAP cannot be corrected just
   by recoloring it. Show combined 12-sample overlays, individual sample
   overlays and cycle scores before/after. Coordinate orientation can change
   between fits, so compare neighborhood structure rather than axis positions.
3. Recluster the corresponding high-dimensional representation with the
   user's selected Python algorithm. HiCAT and Scanpy/Leiden are distinct
   methods. Preserve original counts for markers/DE; residualized expression
   is an analysis representation, not measured counts.
4. Map the same raw-count cells with local Python Allen cell_type_mapper,
   keeping stable IDs and joining results to both embeddings/new clusters.
   Cell-level mapping does not require UMAP or cluster IDs, and should not use
   negative regression residuals as raw input. Current comparator entrypoint:
   `scripts/reference_mapping/run_mapmycells_comparator.py`; existing local
   environment and assets are configured in
   `config/e15_mge_mapping_diagnostic.env`.
5. Avoid repeated full count H5ADs. Keep compact cell-indexed coordinates,
   scores, labels, graphs/model artifacts only as needed, parameters and
   source snapshots. No full cell-by-cell dense coassignment matrix. Reuse
   existing software/reference resources. Execute numerical work via Slurm.

## MapMyCells interpretation

Existing reference is Allen adult whole-mouse-brain WMB-10X, CCN20230722.
The earlier comparator assigned approximately 99% of the original 450788-cell
query to adult OB-IMN GABA; it used an older pre-Step02 cell universe, so results
cannot simply be relabeled as a new fit on 446349 cells. Adult fine labels
are not validated embryonic progenitor identities. Preserve hierarchy,
confidence, missing-feature warnings and disagreement with developmental
marker context. No external upload is needed for local Python mapping.

Official implementation documentation consulted:
- https://scanpy.readthedocs.io/en/stable/how-to/cell-cycle.html
- https://github.com/AllenInstitute/cell_type_mapper
- https://github.com/AllenInstitute/cell_type_mapper/blob/main/docs/mapping_cells.md

This handoff records preparation and findings, not completed corrected UMAP,
reclustering or new Allen mapping. No scientific output has been overwritten.

## Execution state

Storage audit job 65311425 completed with zero read errors; approximately
247.5 GiB of path sizes before cleanup (hard-link duplicates may be counted).
Cleanup job 65311673 preserves the final frozen archive and essential inputs.
Comparison job 65311773 depends on successful cleanup. Its executed code and
configuration are frozen under `results/rebuild_20261103/code/`; parameters,
function-stage progress and source identity are retained. Tests confirm blocked
OLS against least-squares residuals and audit alias/missing/ambiguous marker
resolution. Pytest is absent in the scientific environment, so the two test
functions were executed directly and passed. No R is invoked.

Comparison output is `results/rebuild_20261103/cell_cycle_comparison/`: paired
sample overlays, 12 individual sample overlays, cycle overlays, all-cell
coordinates, scores, compact regression/PCA model and corrected graphs. It
retains all approved cells and the exact original 3000 HVGs. No duplicate
full-count H5AD is written. PCA/neighbors/UMAP are refit after intercept+S+G2M
OLS regression on the log-normalized HVG representation.

HiCAT consensus and Allen mapping have not yet been submitted. The next
scientific step is inspecting this requested before/after comparison, then
rerunning Python HiCAT with the selected correction made explicit. HiCAT
recomputes local representations at recursive nodes, so simply supplying
a corrected 2D UMAP would NOT apply the correction to HiCAT. Preserve raw
expression for DE/detection and mapping; do not feed negative residuals to
MapMyCells. A storage-conscious HiCAT rerun should write large per-iteration
graphs/indexes/model intermediates to job-local scratch and retain compact
assignments, fitted-model/prototype essentials, configuration and consensus
artifacts on NFS, rather than rebuilding the former 273494-file archive.

## Completed cleanup and corrected comparison submission

Cleanup job 65311673 completed successfully at 18:34:41 UTC on 3 November
2026. Removed 47 directories / 275712 files: 220.52 GiB logical path sizes
(189.86 GiB allocated path sizes; hard-link sharing means this is not a quota
measurement). Final archive hashes and approved QC input hash verified after
deletion; Step00 raw, Step02 approved QC, probe inputs, software/reference
resources and compact original embedding/HiCAT assignments were preserved.

Initial comparison job 65311773 stopped because the published Step06 raw-count
H5AD omits its transient highly_variable flags. Its outputs are preserved at
`cell_cycle_comparison_failed_65311773/`. The fixed code rederives the original
3000 sample-aware seurat_v3 HVGs from the same approved raw counts before
normalization. This reproduces the original Step06 feature-selection call;
it does not choose new genes based on corrected results. Source/config hashes
were refreshed before resubmission as job 65313839. Runtime `job_ids.json`
and `status.json` are authoritative for current execution. No HiCAT or Allen
rerun is claimed yet.

## Single QC / cell-cycle supplementary figure

User requested one assembled supplementary figure with QC before/after by
sample, counts entering/retained, and UMAPs before/after S+G2M regression and
scaling. Renderer: `scripts/render_pcdh19_qc_cycle_supplement.py`; output
`results/rebuild_20261103/supplementary_QC_cell_cycle/` (PNG/PDF/SVG, exact
retention TSV, legend and provenance). This does not imply a before-QC UMAP: 
both UMAPs use the identical 446349 retained cells; QC distributions/counts
show all 450788 versus the retained subset separately.

Comparison job 65313839 completed corrected PCA but failed at graph creation: 
Scanpy defaulted to the empty placeholder X rather than X_pca. Fixed by
explicit `use_rep="X_pca"`; corrected PCA is now checkpointed before graph
construction. Failed attempt is retained at
`cell_cycle_comparison_failed_65313839/`. New comparison job 65314443;
supplement job 65314444 depends on its success. Prior supplement dependency
65314368 was cancelled when its parent failed. Runtime job_ids.json is current.

## Completed cell-cycle comparison

Job 65314443 completed successfully (all446349cells, finite coordinates; no
new cell exclusion). Before/after PCA, graphs and UMAPs are saved in
`cell_cycle_comparison/`, with sample/cycle/lineage overlay figures. Maximum
absolute PC-score correlations across 50 PCs decreased from 0.390643 to
0.025803 (S) and 0.806257 to 0.035687 (G2M). These diagnose the intended
regression and are not proof of complete removal of biological cycle effects.

The single supplement rendered successfully as job 65314444. Final cosmetic
revision/PC-correlation annotation is job 65316378, overwriting that same
supplement figure package rather than adding duplicate data directories.
The output folder is `supplementary_QC_cell_cycle/`; all three image formats,
per-sample retention TSV, exact legend and input/code provenance are saved.
No clustering or Allen mapping is part of this completed correction-only step.


Supplement completed (2026-11-03): supplementary_QC_cell_cycle/PCDH19_QC_cell_cycle_supplement.{png,pdf,svg}. Input 450,788; excluded 4,439; retained 446,349. Matched post-QC UMAP comparison uses only S_score/G2M_score regression and rescaling; original counts preserved. Final render job 65316601 succeeded; PNG visually reviewed. See FIGURE_LEGEND.txt, PROVENANCE.json and VALIDATION.json alongside exports. HiCAT and Allen mapping have not been rerun.


## Minimum detected genes update — 3 November 2026

User selected an additional minimum detected-gene filter and authorized repeating
cell-cycle correction followed by the prior Python HiCAT consensus procedure.
The chosen explicit threshold is **at least 500 detected genes per cell**.
This removes **1,095 additional cells** from approved Step02 (0.2453%), leaving
**445,254 cells**, all 19,071 genes and all 12 samples. Relative to the original
450,788 delivered filtered-barcode cells, total exclusions are **5,534** and
retention is **98.7724%**. An earlier conversational total of 1,195 was an
arithmetic error, corrected against the actual filtered AnnData execution.
A 200-gene threshold would remove no additional approved cells. The 500-gene
threshold is an explicit QC choice, not an optimized or validated biological cutoff.

The new comparison is `cell_cycle_comparison_min500/`: both unregressed and
S/G2M-regressed PCA/neighbor/UMAP arms are refitted on the identical 445,254
retained cells with the same newly selected 3,000 HVGs. Reuses CP10K normalization,
scaling/clip10, PCA50, cosine k30, UMAP min_dist0.3/spread1 and seed20260903.
The supplement at `supplementary_QC_cell_cycle/` will be replaced on successful
completion, including per-cell final QC dispositions and per-sample retained totals.

HiCAT runtime is `hicat_min500_cycle/`; `JOBS.json` is the submitted dependency
chain and `config.json`/`FROZEN.json`/`REFERENCE.json` freeze the computation.
Uses the previous Allen Python reference and recursive feature/graph/DE rules,
100 independently sampled 80% fits and the prior exact seed schedule, plus
full-population initializer, membership coassignment merging/refinement and
final expression-DE merge/audit. No target cluster count or biological labels.
This is a modified rerun: local HVGs are selected with the original Allen rule
on unregressed ln(1+CPM); only PCA/graph and held-out Pearson marker mapping use
global S/G2M OLS residuals, SD scaling and upper clipping at 10. The global model
is fitted on the retained population before sampling. Original ln(1+CPM) gene
expression/detection remains the input for Allen DE; residuals are not counts.
These changes prevent reclustering on the old uncorrected HiCAT representation.
No genotype, sex, sample, library-size or mitochondrial covariates are regressed.

The count checkpoint is unchanged. One derived sparse ln(1+CPM) asset is shared
by jobs. Large node graphs/projections use job-local temporary storage and are
removed at exit; compressed numerical feature/PCA/DE audit archives, labels,
models and membership matrices remain on NFS. The array runs at most four
iterations concurrently after iteration0 succeeds. No Allen mapping submitted.
Numerical OLS/SD and block mapping checks passed, the old environment's count
reader matches QC metrics on a real 100-cell slice, and all 100 exact Annoy seeds
passed native forwarding/graph qualification. Execution status must be read
from individual completion files and logs; submission does not imply completion.


MapMyCells method locked by user (2026-11-03): use local Allen cell_type_mapper with a custom developmentally appropriate reference, not a replacement label-transfer algorithm. Official mapping_cells.md explicitly supports arbitrary taxonomies. Metadata validation uses the actual installed MapMyCells taxonomy API. La Manno published Class/Subclass/CellType need path-qualified machine IDs because generic Undefined labels occur beneath multiple parents; original published labels are preserved verbatim. Bandler CA301 later-MIND labels validate for 4,481 verified WT MGE cells, but one leaf has only one reference cell. These reference-coverage/support limitations remain to resolve before marker fitting and final mapping. See mapmycells_developmental_reference/REFERENCE_VALIDATION.json. Existing QC/HiCAT chain continues; no mapping job has been submitted.


## 4 November 2026 — concurrency increased

HiCAT array 65323908 was increased from four to 20 concurrent iteration jobs
at user request. At 08:39 EST, the scheduler confirmed 20 actually running
(tasks032–051), 48 pending; 32/100 iterations had completed at the preceding
checkpoint check. Previous September concurrency was ten. Scientific settings
and completed work were preserved. Final consensus/DE 65323909 remains
dependent on completion. See the current workflow guide for the operational
record. Cancellation target is awaiting clarification; no cancellation claimed.


### 2026-11-04 09:12 EST — temporary-storage failure and recovery

At 09:11 EST, five tasks failed (037, 038, 040, 042, 046). Logs show
Errno28/no space in inherited `/home/elcrespo/.codex/tmp`, contrary to the
intended node-local scratch behavior. Increasing concurrency exposed this
storage problem. The prior afternoon completion estimate is withdrawn.
Stopped original unfinished array 65323908 and dependent consensus 65323909;
all 32 completed iteration seals (000–031) preserved. Submitted remaining
68 iterations as recovery array **65352448**, maximum **20 concurrent**,
using a unique node-local `/tmp` directory per job with exit cleanup.
Replacement consensus/DE **65352449** depends on successful recovery.
Original frozen scientific source/config/seeds and wrapper are unchanged;
a separate recovery wrapper sets temporary-storage environment explicitly.
Runtime evidence: `hicat_min500_cycle/NODE_LOCAL_RECOVERY_20261104.json`;
`JOBS.json` now records replacement jobs and preserves original job IDs.
Recovery submitted; completion and actual concurrency require subsequent
scheduler/checkpoint verification. No final consensus count available.


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

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

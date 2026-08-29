CODEX HANDOFF: CURATE THREE DEVELOPING MOUSE MGE scRNA-seq REFERENCES
=======================================================================

PROJECT CONTEXT
---------------

I have single-cell RNA-seq data from DISSECTED EMBRYONIC MOUSE MEDIAL GANGLIONIC EMINENCE (MGE), approximately E15.

I need to select ONE published WT developing-mouse reference for downstream reference mapping / label transfer.

The three candidate references to curate are:

1. La Manno et al., Nature 2021
   "Molecular architecture of the developing mouse brain"
   DOI: 10.1038/s41586-021-03775-x
   Main code/data repository:
   https://github.com/linnarsson-lab/development-mouse
   Raw sequencing:
   NCBI BioProject PRJNA637987

2. Bandler et al., Nature 2022
   "Single-cell delineation of lineage and genetic identity in the mouse brain"
   DOI: 10.1038/s41586-021-04237-0
   GEO:
   GSE188528
   Critical WT MGE sample:
   GSM5684876 / CA301 = WT MGE E15
   Analysis repository:
   https://github.com/mayer-lab/Bandler-et-al_lineage

3. Mayer et al., Nature 2018
   "Developmental diversification of cortical inhibitory interneurons"
   DOI: 10.1038/nature25999
   GEO SuperSeries:
   GSE104158
   10x series:
   GSE104156
   Important MGE sample:
   GSM2790898 = MGE_E13.5_Lhx6pos
   Smart-seq2 MGE series:
   GSE104157

CURRENT EXECUTION STATUS
------------------------

The early processed-object and sample-metadata checkpoint completed on
2026-08-29 in run
`00_developing_mouse_mge_reference_curation_20260829_141944_3b2ad52`.
The successful replacement chain used jobs 59168348, 59168349, and 59168350.
The visual/report stage was added with job 59169776; layout and wording
revisions completed as jobs 59170025 and 59170165.
The Bandler artifact recovery/inspection first succeeded as job 59171169; the
official-supplement, hierarchy-figure, and report publication completed as job
59171983 after a schema-whitespace correction; final exact-age wording was
republished successfully as job 59172069. Failed attempts 59171052 and
59171918 remain in `provenance/job_ids.tsv` rather than being hidden.
The public MIND-atlas package was finalized by job 59175410. Network-only
attempts 59175246 and 59175304 failed because the compute-node proxy rejects
the public WebSocket; they also remain recorded. The corrected rerunnable
design captures the intended public PDFs on the submission host using the
already-frozen run copy, then SLURM validates those exact files, performs all
vector parsing/count reconciliation, and rebuilds the reports.
The subsequent deposited-barcode workflow completed successfully as job
59179410 (7,420/7,420 E15 joins definitive); its environment-only failed
attempt 59177809 remains recorded rather than hidden.
La Manno contains embedded author annotation columns and `X_UMAP`/`X_tSNE`;
the inspected Bandler and Mayer P0 files themselves are counts-only. Bandler's
separate later-MIND barcode recovery now displays later-atlas labels for 4,481
CA301 cells; Mayer still requires a stable barcode-to-label artifact, and the
original Bandler 2022 per-cell labels remain unavailable.
The evidence-based comparison is in that run's
`REFERENCE_CURATION_REPORT.md`. No Paper 3 cell was loaded or mapped.

The observed E15/MGE answer is now explicit. Bandler is the only candidate
with a directly identified exact WT E15.5 MGE sample (GEO shorthand E15) in the inspected P0
(`CA301`/`GSM5684876`, 4,516 cells). Mayer has a directly identified MGE
sample, but it is `MGE_E13.5_Lhx6pos` (6,515 cells), not E15 and not an
unbiased whole-MGE census. La Manno contains exact E15.0/E15.5 cells and
ventral-forebrain dissections, but no explicit author MGE label was found in
the inspected annotation vocabularies; exact E15 MGE membership therefore
remains unproven and ventral forebrain must not be relabeled as MGE.

The auxiliary visual audit publishes the La Manno author-label hierarchy at
`LaManno2021/figures/01_author_annotation_hierarchy_composition.{png,pdf}` and
the three-study evidence matrix at
`figures/02_candidate_e15_mge_annotation_evidence.{png,pdf}`. Its supporting
tables are `LaManno2021/metadata/class_by_subclass_composition.tsv`,
`LaManno2021/metadata/forebrain_gabaergic_celltype_composition.tsv`, and
`tables/candidate_e15_mge_status.tsv`. Those early plots show only observed
P0 labels/evidence; the later Bandler barcode-recovery package is a distinct
subsequent result and Mayer annotations remain unrecovered.

The Bandler recovery stage subsequently found the historical author link to
`STICR.seuratobject.RDS` in repository commit
`a8fa139a5ed6d8832b07b61f384982d630893c93`. The live 1.02-GB artifact is a
21,051-feature x 65,700-cell Seurat object with author broad/refined labels and
PCA, Harmony, and UMAP reductions. It is the postnatal STICR reference: its 18
sample IDs exclude CA298--CA303, so it does not label CA301. Official
Supplementary Data 1--4 and author-posted TrackerSeq metadata are cached
separately. They establish the seven-dataset embryonic integration design and
the complete 21-label embryonic vocabulary, but provide no barcode-to-label or
UMAP table. A later Mayer-lab interactive atlas names a local integrated object
(`EXCIT_INHIBIT_cleaned_sub.rds`) containing Bandler embryonic cells, labels,
stage, study, cell IDs, and UMAP; that object and its generated TSV/HDF5 files
are not included in the public atlas repository or returned by an intended
public object/table endpoint. However, the running Shiny app's intended public
vector-PDF handlers expose more evidence than the repository alone: exact
circle reconciliation yields 18,424 retained Bandler cells, split into 11,004
E13 and 7,420 E15 cells, with 12 later-atlas clusters, 2,877 `Mitotic` cells,
and 15,547 `Inhibitory Neuron Precursor` cells. Although sample/region is not a
discrete public field, 24 intended public RNA-expression vectors preserve both
expression fingerprints and cell order. Joining those fingerprints to the
deposited matrices definitively resolves every E15 row as 4,481 CA301 MGE,
2,937 CA302 CGE, or 2 CA303 LGE cells. Keep these later MIND-atlas labels
distinct from the original 2022 21-cluster taxonomy.

IMPORTANT BIOLOGICAL GOAL
-------------------------

The reference must be useful for annotating cells obtained from a whole E15 MGE dissection.

At MINIMUM I need broad developmental/cellular classes such as:

- radial glia / apical progenitor / neural stem-progenitor
- cycling progenitor
- intermediate progenitor / basal progenitor, if distinguishable
- neuroblast / newborn neuron
- immature neuronal / GABAergic precursor
- more differentiated MGE-derived neuronal states, if available
- microglia / macrophage
- oligodendrocyte lineage / OPC, if present
- astroglial lineage, if present
- endothelial / vascular / pericyte, if present
- other non-neural populations recovered from the dissection

Even better would be hierarchical subannotations within these broad classes, especially MGE lineage states.

For example, if a reference has:

level_1 = "Progenitor"
level_2 = "Radial glia"
level_3 = "MGE radial glia / NKX2-1+ state"

or:

level_1 = "Neuron"
level_2 = "Immature GABAergic neuron"
level_3 = "i_Nkx2-1", "i_Lhx6/Npy", "i_Nxph1", etc.

I want ALL published annotation levels preserved and documented.

PRIMARY PRIORITY
----------------

DO NOT start by reprocessing raw FASTQs.

The highest priority is:

1. obtain the MINIMAL published processed object(s) needed to reproduce the paper's annotated UMAP(s);
2. determine exactly what cell-level annotations are already present;
3. make UMAPs using the AUTHORS' OWN embeddings and AUTHORS' OWN labels whenever those are available;
4. determine whether MGE cells and approximately E15 cells can be directly selected;
5. determine how useful the existing annotation hierarchy is for my E15 MGE data.

I want to see the paper's annotation structure BEFORE deciding which reference to use.

DO NOT manually invent a new annotation scheme during this audit.

DO NOT recluster unless required only as a diagnostic.

DO NOT perform de novo annotation.

DO NOT run integration against my data yet.

DO NOT download raw FASTQs.

The goal of this stage is REFERENCE CURATION AND VERIFICATION.


GENERAL RULE: PROVE, DO NOT ASSUME
----------------------------------

For every dataset, distinguish these possibilities explicitly:

A. The downloadable processed object ALREADY contains cell-level annotations and UMAP coordinates.

B. The downloadable expression object contains cells/counts but annotations are in a separate metadata table that can be joined by a stable cell ID/barcode.

C. Published annotations exist, but only inside analysis code or intermediate files that must be reconstructed.

D. Published annotations appear only in the paper/supplement and cannot presently be mapped back to individual deposited cells.

These cases are NOT equivalent.

Do not say "annotated dataset available" unless you have inspected the actual file or directly verified the metadata schema.

Likewise, a file ending in .RDS is NOT automatically a Seurat object.
Use R:
    class(readRDS(...))
and inspect the object.


DIRECTORY / PIPELINE ORGANIZATION
---------------------------------

Build this as a small, isolated reference-curation module that can later plug into the existing pipeline.

Do not disturb existing analysis outputs.

Implemented physical structure:

PAPER3_ROOT/
├── inputs/developing_mouse_mge/
│   ├── LaManno2021/{source,metadata_sources}/
│   ├── Bandler2022/{source,metadata_sources}/
│   └── Mayer2018/{source,metadata_sources}/
└── results/00_developing_mouse_mge_reference_curation/<run_id>/
    ├── README.md
    ├── code/
    ├── config/
    ├── tables/
    ├── logs/
    ├── provenance/
    ├── LaManno2021/{metadata,figures,audit}/
    ├── Bandler2022/{metadata,figures,audit}/
    ├── Mayer2018/{metadata,figures,audit}/
    └── SUCCESS.txt or FAILED.txt

Large P0 objects are cached once under `inputs/` and are not duplicated into
every versioned result package. Each result package contains the exact Python,
R, shell-submission, and SLURM files copied before submission. The SLURM jobs
must execute the copies under that run's `code/` directory, not mutable files
from the repository.

Default submission creates a new timestamped run. Intentional reruns within
this step use `--replace-run RUN_ID`. Replacement is permitted only for a safe
existing run basename beneath this exact step, is refused while any recorded SLURM job
is active, and clears only that named run package. It never clears or
redownloads the shared P0 source cache. Frozen PCDH19 Steps 03–07 remain outside
this replacement policy.

Only put minimally necessary processed files into source/.

Do not duplicate giant objects unnecessarily.

Record:
- URL
- accession
- original filename
- download date
- byte size
- checksum
- file type
- source repository/database
- paper
- whether it was actually downloaded

SAMPLE / LIBRARY METADATA CONTRACT
----------------------------------

The first checkpoint must also prepare the evidence needed for the next
reference-curation stage. For every published or author-registered sample,
retain the exact sample ID and, when supplied, age/stage, tissue, anatomical
region, genotype/condition, sex, strain, organism, pool, replicate, QC state,
project/cohort, single-cell modality, capture technology, library strategy,
source, selection, layout, instrument, platform, chemistry, sequencing
protocol, reference genome, isolation/dissociation method, experimental
selection or enrichment, raw accession, and BioSample accession.

Published sample metadata and membership in the downloaded P0 object are two
different claims. Record both. Never claim that a published sample is present
in the P0 cells unless an author sample field, stable barcode/cell identifier,
or other explicit join proves it. Preserve all original GEO/author fields in a
long-form table even when they do not map into the standardized columns.

Required first-checkpoint outputs are:

- `<paper>/metadata/published_sample_inventory.tsv`
- `<paper>/metadata/published_sample_metadata_long.tsv`
- `<paper>/metadata/object_sample_inventory.tsv`
- `<paper>/metadata/sample_inventory.tsv`
- `<paper>/metadata/study_sample_summary.tsv`
- `tables/all_candidate_reference_samples.tsv`
- `tables/study_sample_summary.tsv`
- `tables/sample_metadata_data_dictionary.tsv`
- `tables/reference_curation_requirements_ledger.tsv`
- `tables/reference_annotation_availability.tsv`
- `<paper>/metadata/annotation_dictionary.tsv` (observed P0 labels only;
  header-only when the P0 has no labels)
- `LaManno2021/metadata/class_by_subclass_composition.tsv`
- `LaManno2021/metadata/forebrain_gabaergic_celltype_composition.tsv`
- `LaManno2021/figures/01_author_annotation_hierarchy_composition.{png,pdf}`
- `tables/candidate_e15_mge_status.tsv`
- `figures/02_candidate_e15_mge_annotation_evidence.{png,pdf}`
- `REFERENCE_CURATION_REPORT.md`
- `Bandler2022/BANDLER_AUTHOR_OBJECT_RECOVERY_REPORT.md`
- `Bandler2022/metadata/published_embryonic_sample_design.tsv`
- `Bandler2022/metadata/published_{embryonic_cluster,postnatal_broad_class,postnatal_refined_cluster}_{inventory,markers}.tsv`
- `Bandler2022/author_object_audit/author_artifact_scope.tsv`
- `Bandler2022/author_object_audit/author_seurat_{class_by_cluster,sample_inventory}.tsv`
- `Bandler2022/figures/01_bandler_recovered_and_published_annotation_structure.{png,pdf}`
- `Bandler2022/interactive_atlas/MIND_PUBLIC_ATLAS_CAPTURE_REPORT.md`
- `Bandler2022/interactive_atlas/metadata/atlas_{cluster,class,stage}_by_study.tsv`
- `Bandler2022/interactive_atlas/metadata/bandler_cluster_by_stage.tsv`
- `Bandler2022/interactive_atlas/figures/public_umap_{cluster,class,stage}_by_study.{png,pdf}`
- `Bandler2022/interactive_atlas/audit/{public_endpoint_scope,vector_count_validation}.tsv`
- `Bandler2022/interactive_atlas/barcode_recovery/BANDLER_E15_BARCODE_RECOVERY_REPORT.md`
- `Bandler2022/interactive_atlas/barcode_recovery/metadata/{bandler_e15,CA301}_later_atlas_barcode_join.tsv`
- `Bandler2022/interactive_atlas/barcode_recovery/metadata/{sample_recovery_summary,sample_class_cluster_composition}.tsv`
- `Bandler2022/interactive_atlas/barcode_recovery/audit/{barcode_recovery_validation,match_status_summary,fingerprint_gene_manifest}.tsv`
- `Bandler2022/interactive_atlas/barcode_recovery/figures/ca301_later_atlas_clusters.{png,pdf}`

The combined summary must report both the number/IDs of published samples and
the number/IDs demonstrably represented in the P0 object. Unresolved object
sample IDs and unlinked published samples remain explicit. The requirements
ledger uses `PASS`, `PARTIAL`, `MISSING`, and `NOT_ASSESSED`; the latter is
expected for annotation dictionaries, author-embedding reproduction, broad
class coverage, exact MGE/age-matched counts, and readiness ranking because
those belong to the reviewed next stage. The present visual audit is limited
to observed label hierarchy and the explicit E15/MGE evidence status above;
it does not substitute for the still-missing Mayer label join, the original
Bandler 2022 21-state barcode assignment, or a final reference-readiness
ranking. The later MIND-atlas Bandler E15 barcode join is now complete.

The implementation must remain object-oriented. `SourceRegistry` owns source
validation, `ImmutableSourceCache` owns atomic read-only materialization,
`PublishedMetadataCollector` owns sample evidence, study-specific
`StudyObjectInspector` subclasses own P0 inspection,
`StudyInspectionWorkflow` owns per-study reconciliation, and
`CheckpointPublisher` plus `CurationRequirementEvaluator` own combined output
publication, `HierarchyComposition` and its plotter classes own the visual
audit, and `ReferenceCurationReportBuilder` owns the observed-object report
and label dictionaries. `BandlerAuthorObjectRecovery` owns atomic author and
publisher artifact retrieval, `PublishedSampleLedger` and
`PublishedAnnotationTables` own supplement extraction, and
`BandlerEvidencePublisher` coordinates the Bandler evidence package. Do not
move these responsibilities into a monolithic command.
`ShinySockJsSession` owns the intended public app session,
`GhostscriptRenderer` owns deterministic PDF rendering/vector conversion,
`VectorCellCounter` owns circle extraction and reconciliation, and
`AtlasEvidencePublisher` owns the later-atlas evidence package. These classes
must preserve the distinction between publicly rendered study-level evidence
and the independently validated deposited-barcode join. `DepositedMatrixCache`
owns immutable CA301--303 acquisition, `FeatureShinySession` owns intended
public expression downloads, `VectorPlotParser` owns cell/color extraction,
`OrderedFingerprintMatcher` owns exact subsequence and duplicate resolution,
and `BarcodeRecoveryPublisher` owns the barcode-level evidence package.


TASK 0 — ENVIRONMENT AND REPRODUCIBILITY
----------------------------------------

Use current stable R, Seurat, Python, Scanpy, anndata, pandas, scipy, matplotlib and supporting readers.

Do not convert objects destructively.

Keep original downloaded files read-only.

Any converted Seurat/H5AD object should be an explicitly derived artifact.

Record software versions in:
the run package's `tables/software_versions.tsv` and per-study audit records.

Create scripts rather than one-off interactive commands.

Scripts should fail loudly when an expected annotation, embedding, or barcode field is absent.


TASK 1 — BUILD A SOURCE MANIFEST BEFORE DOWNLOADING
---------------------------------------------------

For EACH of the three papers, audit:

- publication URL / DOI
- GEO/SRA/BioProject accessions
- author GitHub repository
- author lab download site
- supplemental files
- processed expression objects
- cell metadata files
- cluster annotation tables
- UMAP / tSNE coordinates
- count matrices
- raw sequencing availability
- spatial files, if relevant
- code needed to recreate published annotation/UMAP

Create:

the run package's `tables/manifest.tsv`

Columns:

paper
dataset
accession
sample
age
region
genotype
technology
resource_type
filename
format
url
size_if_known
download_priority
downloaded
reason
notes

download_priority should be:

P0 = needed now to reproduce published annotated UMAP / inspect labels
P1 = potentially useful processed metadata/object, small enough to justify
P2 = useful later but DO NOT download now
RAW = raw FASTQ/SRA availability only; VERIFY ACCESS BUT DO NOT DOWNLOAD


TASK 2 — LA MANNO 2021
----------------------

Known starting point:

Author repository:
https://github.com/linnarsson-lab/development-mouse

The repository advertises:

developing_mouse_nervous_system.h5ad
approximately 31,053 genes x 292,495 cells

and:

dev_all.agg.loom
approximately 31,053 genes x 798 clusters

Raw sequencing is available under:
PRJNA637987

The Nature paper states that the study annotated hundreds of developmental cell states and integrated spatial information.

WHAT TO DO:

A. Download the MINIMAL cell-level annotated processed object needed for the audit.
   Prefer:
   developing_mouse_nervous_system.h5ad

B. Do NOT download PRJNA637987 FASTQs.

C. Inspect the H5AD structurally.

Record:

adata.shape
adata.raw status
adata.X type
adata.layers.keys()
adata.obs.columns.tolist()
adata.obsm.keys()
adata.uns.keys()
adata.var columns
sparse/dense representation

Export:

LaManno2021/audit/h5ad_structure.txt
LaManno2021/audit/obs_columns.tsv
LaManno2021/audit/obs_value_counts/

For every plausible annotation column, export unique labels and cell counts.

Potential classes/metadata to search for include terms related to:

cluster
class
subclass
cell_type
celltype
type
subtype
identity
annotation
taxonomy
development
age
embryonic_day
tissue
region
forebrain
telencephalon
ganglionic
MGE
GE
progenitor
radial
neuroblast
microglia

Do not rely only on column names. Inspect values.

D. Determine whether the exact author UMAP coordinates are present in .obsm.

If present:
- identify the correct embedding
- DO NOT recompute UMAP
- reproduce the global paper-like UMAP using author coordinates
- color separately by EACH annotation level

Save:
LaManno2021/figures/01_global_author_umap_<annotation>.pdf/png

E. Identify the subset relevant to MGE / ventral telencephalon.

This must be evidence-based.

Use author metadata, taxonomy and/or paper/code definitions.

Do not simply grep "MGE" and stop.

Determine:
- which cells correspond to MGE directly
- whether MGE is represented as an anatomical region, developmental lineage, cluster family, spatial domain, or combination
- ages represented in those MGE-relevant cells
- number of E14/E15/E16-ish cells if age metadata allows
- whether exact E15 exists

Create:
LaManno2021/metadata/mge_relevant_cells.tsv.gz
or a compact barcode + metadata table if exporting all expression is unnecessary.

F. Make MGE-relevant UMAPs using AUTHOR coordinates and AUTHOR annotations.

For each annotation hierarchy available:
- broad class
- subclass
- cluster
- fine developmental state
- anatomical/tissue label

Save separate plots.

G. Explicitly search the published annotations for the classes I care about:

radial glia / apical progenitors
cycling progenitors
intermediate/basal progenitors
neuroblasts
immature neurons
MGE-derived GABAergic precursors
microglia
OPC / oligodendrocyte lineage
astroglia
endothelial
pericytes / vascular
other contaminating or bona fide non-neural MGE populations

For each class, report:
- exact author label
- exact annotation column
- total cells
- MGE-relevant cells
- age distribution
- marker genes only as supporting context, NOT as a substitute for author annotation

H. Determine whether there is an author-provided hierarchy connecting the 798 states to broader classes.

If there is a hierarchy/table/code mapping:
DOWNLOAD IT if small.
Preserve it exactly.

Do not collapse labels yourself during this audit.


TASK 3 — BANDLER 2022
---------------------

This dataset is biologically especially important because it contains a WT MGE sample collected at E15.5 (GEO shorthand E15).

GEO:
GSE188528

Critical sample:
GSM5684876
CA301
WT MGE E15

Known processed file:
GSM5684876_CA301_filtered_RNA_counts.RDS.gz

GEO describes this as:
"Processed_Transcript: Filtered count matrix"

IMPORTANT:
Do NOT assume this is a Seurat object merely because it is RDS.

A. Download GSM5684876_CA301_filtered_RNA_counts.RDS.gz.

This is P0 because the sample is essentially matched to my experiment.

B. Inspect in R:

obj <- readRDS(...)
class(obj)
str(obj, max.level=2)
dim(obj)

Determine whether it is:
- matrix
- sparse Matrix
- data.frame
- SingleCellExperiment
- Seurat
- other

If Seurat:
inspect:
colnames(obj@meta.data)
Reductions(obj)
Assays(obj)

If matrix:
document exact row/column structure and cell barcode format.

C. Find the published CELL-LEVEL annotation mapping.

Search all of:

1. GSE188528 supplementary files
2. individual GSM records
3. Nature Supplementary Data
4. https://github.com/mayer-lab/Bandler-et-al_lineage
5. scripts and serialized intermediate objects referenced by the code
6. Zenodo/Figshare/lab repositories if linked by authors
7. metadata tables containing cell IDs/barcodes
8. code that assigns cluster identities
9. code that creates the embryonic UMAP
10. code or files mapping CA301 cells to integrated embryonic clusters

The key question is:

CAN THE BARCODE/CELL IDs FROM CA301 WT E15.5 MGE BE MAPPED DIRECTLY TO THE PUBLISHED EMBRYONIC ANNOTATION?

Answer YES/NO/PARTIAL and prove it.

Current proven answer: **YES for the later MIND-atlas annotation, NO for the
original Bandler 2022 21-cluster per-cell taxonomy**. The recovered author
Seurat object is postnatal STICR and contains no CA301 IDs; Supplementary Data
4 defines the 21 original embryonic clusters but has no cell IDs. Separately,
24 intended public MIND RNA-expression vectors and preserved plot order form
an exact fingerprint subsequence of the deposited E15 matrices. They
definitively map 4,481 CA301 barcodes to later-atlas class/cluster labels and
public vector UMAP positions; 35 of 4,516 deposited CA301 cells were not
retained. The same validation maps 2,937 CA302 and 2 CA303 barcodes, totaling
all 7,420 E15 atlas cells. Do not relabel the later-atlas clusters as the
original 2022 assignments.

D. Recover the exact published embryonic annotation vocabulary.

The paper reports embryonic clusters/states including examples such as:

m_Fabp7
m_Top2a
m_Abracl

and immature/postmitotic states such as:

i_Nkx2-1
i_Lhx6/Npy
i_Nnat
i_Sox4
i_Six3
i_Gucy1a3

and differentiated inhibitory trajectories including:

i_Six3/Gucy1a3
i_Ebf1/Isl1
i_Phlda1/Isl1
i_Nr2f2
i_Nxph1

DO NOT take this list as complete.
Extract the COMPLETE label list from the actual published metadata, code, supplement, or object.

Also recover broader cell-class annotations if available.

The Nature supplement includes tables of:
- cell class markers
- cell type cluster markers
- embryonic cluster markers

Use these only to document and cross-check the annotation structure.
Do not infer cell labels from marker tables if a cell-level mapping exists elsewhere.

E. If the integrated annotated object is downloadable:
download the minimal object and inspect it fully.

If it is not directly downloadable but author code reconstructs it from deposited matrices:
DO NOT immediately rerun the whole study.

Instead:
1. identify exactly which small metadata/intermediate artifacts are needed;
2. determine whether the original CA301 barcodes survive integration;
3. determine whether published labels can be reconstructed deterministically;
4. document the shortest path.

If reconstructing the mapping is simple and bounded, do it.

F. Reproduce the published embryonic UMAP using AUTHOR coordinates if available.

Highest priority plots:

1. all embryonic cells colored by published broad class
2. all embryonic cells colored by published embryonic cluster
3. highlight CA301 WT E15.5 MGE
4. CA301-only plot using the same public later-atlas geometry, colored by later-atlas annotation
5. CA301 annotation composition bar plot/table

If exact author UMAP coordinates are unavailable but the repository contains exact code and preprocessed integrated input:
document that distinction before recomputing anything.

G. Explicitly determine whether Bandler can annotate broad MGE dissection classes including:

radial glia
cycling progenitors
IPC/basal progenitors
neuroblasts
immature neurons
glia
microglia
endothelial/vascular

It is possible this study is excellent for neuronal lineage states but less suitable as a census of non-neuronal MGE cells.

Prove which classes are actually represented and annotated.

Do not force absent classes into the taxonomy.


TASK 4 — MAYER 2018
-------------------

Key accessions:

GSE104158 = SuperSeries

GSE104156 = 10x developmental interneuron dataset

GSM2790898 = MGE_E13.5_Lhx6pos

GSE104157 = Smart-seq2 MGE dataset

GSE104156 provides:
GSE104156_digital_expression.csv.gz

Raw sequencing is available in SRA.

A. Determine what processed author-level objects and metadata exist beyond the GEO expression CSV.

Search:
- GEO supplementary files
- paper supplementary material
- Satija lab / Mayer code
- GitHub
- archived analysis scripts
- downloadable Seurat objects / RDS
- cell annotation tables
- cluster tables
- cell barcode metadata
- UMAP/tSNE coordinates

B. Download only the minimum processed data required to:
- identify MGE E13.5 cells
- recover author cell annotations
- reproduce the relevant published embedding if possible

Do not download raw FASTQs.

C. Determine whether the MGE 10x sample is preselected:
GSM2790898 = MGE_E13.5_Lhx6pos

This matters because my E15 sample is an unbiased/dissected MGE sample.

Document the selection strategy and the consequences for reference suitability.

D. Inspect annotation granularity.

Specifically ask whether author labels distinguish:

radial glia / progenitors
cycling cells
intermediate progenitors
neuroblasts
early postmitotic cells
MGE-derived interneuron precursors
emerging PVALB/SST-related trajectories
other classes

Also determine whether microglia/vascular/glial populations were excluded by experimental design or filtering.

E. Reproduce author UMAP/tSNE with original annotations if cell-level coordinates are available.

If the paper used tSNE rather than UMAP, preserve that.
Do not replace the author's embedding merely to standardize appearance.

F. Determine whether the Smart-seq2 MGE data add annotation information not available in the 10x object.

Do not download all Smart-seq2 raw reads.

Only download small processed metadata/expression files if they materially help reconstruct published MGE progenitor annotations.


TASK 5 — RAW DATA ACCESS AUDIT, BUT NO RAW DOWNLOAD
---------------------------------------------------

For all three papers, verify whether complete raw data are publicly accessible.

Do NOT download FASTQ/BAM/SRA data.

Create:

references/developing_mouse_mge/raw_access_manifest.tsv

Columns:

paper
accession
sample
raw_repository
raw_accession
raw_available
assay
reference_genome_if_known
library_chemistry_if_known
estimated_size_if_known
can_reprocess_from_raw
notes

At minimum verify:

La Manno:
PRJNA637987

Bandler:
GSE188528 / associated SRA BioProject

Mayer:
GSE104158 / GSE104156 / GSE104157 / associated SRA

The purpose is simply to establish:
"If I choose this reference later, can I reproduce processing from raw reads?"

Do not spend compute/storage on that now.


TASK 6 — ANNOTATION AUDIT TABLE
-------------------------------

This is one of the MOST IMPORTANT outputs.

Create:

references/developing_mouse_mge/audit_summary.tsv

One row per paper.

Columns:

paper
exact_MGE_age
E15_match
whole_MGE_dissection
selection_bias
processed_cell_object_available
object_format
seurat_object_available
h5ad_available
loom_available
author_embedding_available
author_embedding_type
cell_level_annotations_embedded
separate_cell_metadata_available
barcodes_linkable_to_annotations
broad_annotation_available
fine_annotation_available
hierarchical_annotations_available
radial_glia_label
cycling_progenitor_label
IPC_label
neuroblast_label
immature_neuron_label
MGE_GABA_precursor_labels
microglia_label
OPC_label
astroglia_label
endothelial_label
vascular_pericyte_label
number_MGE_cells
number_age_matched_MGE_cells
raw_data_available
raw_accession
published_umap_reproducible_without_reprocessing
major_strength
major_limitation
reference_readiness

reference_readiness must be one of:

READY_NOW
READY_AFTER_METADATA_JOIN
READY_AFTER_LIGHT_RECONSTRUCTION
REQUIRES_REANALYSIS
NOT_SUITABLE_AS_PRIMARY_REFERENCE


TASK 7 — ANNOTATION DICTIONARIES
--------------------------------

For EACH dataset create a TSV containing every relevant published label.

Example:

LaManno2021/metadata/annotation_dictionary.tsv
Bandler2022/metadata/annotation_dictionary.tsv
Mayer2018/metadata/annotation_dictionary.tsv

Columns:

annotation_column
published_label
parent_label_if_known
annotation_level
n_cells_global
n_cells_MGE
n_cells_age_matched_MGE
paper_figure_or_source
notes

annotation_level should be descriptive, e.g.:

broad_class
subclass
developmental_state
cluster
anatomical_region
developmental_age
lineage
fine_cell_type

If the hierarchy is not explicitly defined by the authors, write:
parent_label_if_known = NA

Do NOT invent parent-child relationships.


TASK 8 — IMMEDIATE FIGURES
--------------------------

The visual priority is to see what we would ACTUALLY get if we selected each paper as a reference.

For each dataset, create the maximum possible subset of these plots using original author embeddings:

A. global_author_embedding_broad_labels
B. global_author_embedding_fine_labels
C. MGE_highlighted_on_global_embedding
D. MGE_only_author_embedding_broad_labels
E. MGE_only_author_embedding_fine_labels
F. E15_or_nearest_age_MGE_embedding
G. annotation_composition_MGE
H. annotation_composition_age_matched_MGE

Use rasterization for large scatter plots when needed.

Do not beautify aggressively.
These are diagnostic reference-audit figures.

Keep author labels unchanged.

Use legends that remain readable.

If an author embedding does not exist, make a text note in place of silently computing a new one.


TASK 9 — VALIDATE THAT A PUBLISHED UMAP IS ACTUALLY REPRODUCED
--------------------------------------------------------------

When claiming that a UMAP/embedding corresponds to the paper:

Provide evidence:

- which object/file supplied coordinates
- exact coordinate column/key
- whether coordinates were author-provided or recomputed
- number of cells plotted
- annotation column used
- paper figure it most closely corresponds to

Write this to a companion TXT/TSV.

Never label a newly calculated UMAP as "published UMAP".


TASK 10 — FINAL HUMAN-READABLE REPORT
-------------------------------------

Create:

references/developing_mouse_mge/REFERENCE_CURATION_REPORT.md

Keep it concise but evidence-driven.

For each paper report:

1. What was downloaded?
2. What is the actual object type?
3. Does the object contain counts?
4. Does it contain normalized data?
5. Does it contain author dimensional reductions?
6. What metadata columns exist?
7. What are the exact annotation levels?
8. What are the exact labels relevant to an E15 MGE dissection?
9. Are broad RG/IPC/neuroblast/microglia/etc. classes present?
10. Are finer MGE developmental states present?
11. Can MGE cells be subset directly?
12. Can E15 or nearest-age cells be subset directly?
13. Is the MGE sample unbiased or experimentally enriched?
14. Can the paper UMAP be reproduced immediately?
15. Are raw data available for later complete reprocessing?
16. What is missing?
17. How much reconstruction would be required before label transfer?

Finish with a side-by-side recommendation table.

BUT DO NOT choose the final reference solely from citation count or journal prestige.

Rank by:

1. MATCH TO E15 DISSECTED MGE
2. EXISTING CELL-LEVEL ANNOTATION QUALITY
3. ABILITY TO RECOVER BROAD CELL CLASSES
4. MGE-SPECIFIC SUBANNOTATION DEPTH
5. ABILITY TO IMMEDIATELY VISUALIZE AUTHOR-ANNOTATED CELLS
6. CELL BARCODE / METADATA TRACEABILITY
7. PROCESSED OBJECT QUALITY
8. RAW REPRODUCIBILITY FOR LATER
9. EXPERIMENTAL SELECTION BIAS


EXPECTED BIOLOGICAL TRADEOFF TO TEST
------------------------------------

Do not assume this is the answer, but explicitly test the following likely tradeoff:

La Manno 2021:
Potentially strongest broad developmental atlas and best coverage of non-neuronal cell types, with a directly accessible H5AD and extensive developmental taxonomy.

Bandler 2022:
Potentially strongest age/anatomy match because CA301 is WT MGE collected at E15.5 (GEO shorthand E15), and potentially excellent MGE neuronal developmental-state annotations, but deposited CA301 RDS is a filtered count matrix and requires a barcode-preserving join to recover published cell labels.

Mayer 2018:
Potentially excellent canonical MGE/interneuron developmental biology and precursor-state resolution, but E13.5 rather than E15 and Lhx6+ enrichment creates important sampling bias relative to my whole E15 MGE dissection.

The audit must PROVE or REFUTE each of these statements.


STOP CONDITIONS / DO NOT DO YET
-------------------------------

At this stage:

DO NOT map my organoid data.
DO NOT map my own E15 mouse data.
DO NOT run scVI.
DO NOT run Seurat anchors.
DO NOT run Harmony.
DO NOT integrate references together.
DO NOT generate de novo cell annotations.
DO NOT download raw FASTQs.
DO NOT normalize/recluster giant datasets just because it is possible.
DO NOT overwrite author metadata.
DO NOT collapse fine labels prematurely.

The deliverable is a CURATED, AUDITED, IMMEDIATELY VIEWABLE set of candidate references.


MINIMUM SUCCESS CRITERIA
------------------------

I should be able to open the output directory and answer, for each of the three papers:

1. What processed data can I actually access?
2. What object/file type is it?
3. What published cell-level labels are actually available?
4. Are those labels embedded in the object or recoverable by barcode?
5. What are ALL annotation levels?
6. Does it distinguish RG, cycling progenitors, IPCs, neuroblasts, immature neurons, microglia, etc.?
7. What finer MGE states exist?
8. How many MGE cells are there?
9. What developmental ages are represented?
10. How close is it to E15?
11. Can I immediately display the authors' annotated UMAP/embedding?
12. Can I later reproduce the dataset from raw reads?
13. What is the shortest path to use this dataset for label transfer?
14. Which of the three is most reference-ready for an E15 whole-MGE sample?

FIRST EXECUTION STEP
--------------------

Begin with the source/metadata audit and smallest P0 downloads.

Do not launch large jobs.

I want an early checkpoint after the three actual processed objects/resources have been inspected, BEFORE any reconstruction or heavy analysis.

At that checkpoint, print a compact table with:

paper
P0 file
actual object type
dimensions
author embedding present?
annotation columns present?
MGE selectable?
age selectable?
cell-level labels immediately usable?
next minimal action

Then proceed with the UMAP/annotation audit.

The implemented first-checkpoint submission entry point is:

    ./paper3_pcdh19/bin/submit_developing_mouse_mge_reference_curation.sh --dry-run
    ./paper3_pcdh19/bin/submit_developing_mouse_mge_reference_curation.sh

The default command creates a new versioned run. To intentionally regenerate
one existing inactive run within this step:

    ./paper3_pcdh19/bin/submit_developing_mouse_mge_reference_curation.sh \
      --replace-run 00_developing_mouse_mge_reference_curation_<timestamp>_<commit>

The submitted chain is source audit, a three-study P0 inspection array, and a
dependency-gated checkpoint/report/visual aggregator. The array also retrieves only the
small GEO/author metadata registries needed to publish the sample/library
inventories above. It stops after the early object/sample-metadata checkpoint;
missing Bandler/Mayer label reconstruction and any final reference selection
must be separately reviewed later stages.

The rerunnable later-atlas follow-up for an existing completed run is:

    ./paper3_pcdh19/bin/submit_mind_public_atlas_capture.sh --dry-run \
      /absolute/path/to/00_developing_mouse_mge_reference_curation_<run>
    ./paper3_pcdh19/bin/submit_mind_public_atlas_capture.sh \
      /absolute/path/to/00_developing_mouse_mge_reference_curation_<run>

It freezes exact executable copies into that run's `code/` directory, uses
that frozen copy for public acquisition on the submission host, and has SLURM
independently validate and parse the captured files. A repeat submission
atomically overwrites only the derived
`Bandler2022/interactive_atlas/` package and regenerates the main report; it
does not replace cached scientific inputs, the parent checkpoint, or any
other PCDH19 step.

The rerunnable deposited E15 barcode follow-up is:

    ./paper3_pcdh19/bin/submit_bandler_e15_barcode_recovery.sh --dry-run \
      /absolute/path/to/00_developing_mouse_mge_reference_curation_<run>
    ./paper3_pcdh19/bin/submit_bandler_e15_barcode_recovery.sh \
      /absolute/path/to/00_developing_mouse_mge_reference_curation_<run>

It freezes its Python, R, shell, report, and SLURM files under `code/`, reuses
immutable deposited matrices, captures public feature PDFs on the submission
host, and overwrites only
`Bandler2022/interactive_atlas/barcode_recovery/` before rebuilding the main
report. Failed jobs remain in `provenance/job_ids.tsv`.

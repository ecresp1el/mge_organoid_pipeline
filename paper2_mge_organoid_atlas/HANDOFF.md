# Paper 2 MGE organoid atlas: operational handoff

Last updated: 2026-08-27

## Scientific objective

Build a cross-study in-vitro MGE organoid single-cell atlas from exactly six
processed datasets:

- Varela DIV30
- Varela DIV90
- Walsh
- Bershteyn 2025
- Bershteyn 2023
- Siebert 2026

The initial workflow is:

```text
six processed organoid objects
  -> audit expression layers and metadata
  -> create paired clean Seurat/AnnData canonical inputs
  -> harmonize genes and metadata
  -> freeze one pre-integration master input
  -> compare integration methods on that same input
  -> evaluate batch removal and biological preservation
  -> select a supported integrated representation
  -> freeze the in-vitro MGE organoid atlas
```

Primary fetal datasets are not members of the initial integration. They will
be used later as independent references for mapping and biological evaluation.

## Non-negotiable scope decisions

1. Do not repeat the six studies' original QC by default.
2. Do not silently remove cells, genes, assays, layers, or metadata during the
   input audit.
3. Do not include fetal cells in the frozen pre-integration organoid object.
4. Do not let individual integration methods receive different starting cells,
   genes, or harmonized covariates unless a method-specific requirement is
   explicitly documented and evaluated.
5. Do not overwrite a frozen input or selected atlas. Active working packages
   may use the explicit replacement mode described below; frozen milestones
   require a new versioned package documenting what they supersede.

## Canonical locations

```text
Code repository:
/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline

Paper 2 code/configuration:
/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline/paper2_mge_organoid_atlas

Existing MGE project/data root:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder

Paper 2 Turbo root:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper2_mge_organoid_atlas
```

The Paper 2 Turbo root owns its own `inputs/`, `results/`, `logs/`, `jobs/`,
and `final_figures/`. Existing processed objects remain read-only upstream
provenance; the six cleaned copies under `inputs/canonical/` are the only
objects downstream Paper 2 workflows are allowed to read.

## Registered input objects

The authoritative candidate-input registry is
[`config/input_objects.tsv`](config/input_objects.tsv). It currently selects
the annotated Walsh final object because that is the object used by the
existing cross-study workflows; the registry documents its relationship to
the post-QC unannotated object.

At this stage, "registered" does not mean "frozen." The six objects become a
frozen pre-integration input only after the following are reviewed:

- file identity and SHA-256 checksum;
- object class and dimensions;
- assays and available counts/data/scale-data layers;
- feature namespaces and cross-study gene overlap;
- metadata columns, types, missingness, and candidate harmonization fields;
- reductions already present, recorded for provenance but not assumed to be
  valid shared integration inputs;
- per-study cell identifiers and a collision-safe atlas cell ID plan.

### Step 00 UMAP visualization contract

Step 00 does not recompute UMAP embeddings. It reads each processed object's
saved UMAP coordinates and regenerates publication-quality cluster inventory
plots in mandatory PNG/PDF and optional editable-text SVG. Rasterized layers
in every enabled format are fixed at 300 dpi.

- Every registered cell with finite saved UMAP coordinates is shown.
- No study-level QC is rerun.
- DIV90 current clusters 6/7 are retained in the input-audit plots even though
  they were excluded from some prior publication-style final figures.
- DIV90 receives the established plotting-only transform
  `UMAP1_plot = UMAP1_original`, `UMAP2_plot = -UMAP2_original`.
- Original coordinates are retained unchanged in the per-cell audit tables.
- DIV30 uses the documented seven-raw-cluster to five-paper/manual-class
  mapping used by the existing cross-study final figures. Both raw and
  displayed cluster counts are reported.
- DIV90 uses `cluster_number_name` from the object when present, while also
  retaining `seurat_clusters` as the raw ID.
- Other studies use the first available documented cluster-name metadata
  field; if none exists, the plot and tables state that only raw IDs were
  available.

Required UMAP outputs include a six-study grid, a dedicated DIV30/DIV90 grid,
one figure per study, cluster-count and label-source tables, a figure manifest,
and the per-cell original-coordinate tables.

## Independent Paper 2 step numbering

Step numbers are local to `paper2_mge_organoid_atlas/`; they do not continue
the historical top-level pipeline numbering.

| Step | Status | Purpose |
| --- | --- | --- |
| `00_input_audit` | Completed, replaced in place, and visually checked; current job `58956196` | Registered and checksummed the six processed objects; inventoried assays, layers, features, metadata, reductions, and saved UMAP/cluster labels; regenerated all-cell UMAP inventory figures. Job `58956196` replaced the original job `58955368` package in the same run directory. |
| `01_canonical_inputs` | Completed and frozen; array job `58958446`, finalizer `58958447` | Created minimal current Seurat RDS and AnnData H5AD pairs for all six audited studies, with exact expression/ID/metadata equivalence and no integration or harmonization. |
| `02_harmonize_genes` | Report-only audit completed; job `58978281`; STOP FOR REVIEW | Mapped feature identifiers to versionless GENCODE 50 Ensembl gene IDs and reported namespaces, ambiguities, unresolved features, duplicate mappings, pairwise overlaps, and two six-way-intersection definitions. No matrix was changed or created. |
| `02b_legacy_gene_id_recovery` | Report-only audit completed; job `58986448`; STOP FOR REVIEW remains active | Used confirmed original Varela Cell Ranger feature tables and otherwise labeled mappings from GENCODE 27/32/35/44/50 as historical-reference candidates. It did not rewrite canonical inputs or Step 02. |
| `03_harmonize_metadata` | Planned; blocked pending Step 02 review | Create a documented common schema for study, dataset, sample, replicate, age/time point, cell labels, and QC provenance. |
| `04_freeze_preintegration` | Planned | Produce and validate the immutable six-study pre-integration master object/package. |
| `10_scvi` | Planned | Run scVI from the frozen input. |
| `11_liger` | Planned | Run LIGER from the frozen input. |
| `12_scpoli` | Planned | Run scPoli from the frozen input. |
| `13_seurat_rpca` | Planned | Run Seurat RPCA from the frozen input. |
| `20_integration_benchmark` | Planned | Compare batch removal and biological preservation using common metrics and panels. |
| `30_freeze_atlas` | Planned | Select, justify, validate, and package the supported organoid atlas representation. |
| `40_fetal_reference_mapping` | Deferred | Map or benchmark the frozen organoid atlas against independent fetal GE references. |

New steps should receive a new number and descriptive name. Existing step
numbers should not be repurposed after they have generated a reviewed output.

## Major-output reproducibility contract

Every major `results/` package and every `final_figures/` package must contain:

```text
README.md                 scientific purpose, inputs, methods, outputs, status
code/                     exact submitted scripts/wrappers used for this run
config/                   submitted and resolved configuration/input registry
figures/png/              high-resolution raster exports when applicable
figures/pdf/              publication/vector exports when applicable
figures/svg/              optional editable SVG exports when enabled
tables/                   machine-readable results and validation summaries
logs/                     SLURM stdout/stderr and pipeline logs
provenance/               command, job ID, Git state, environment, checksums
SUCCESS.txt or FAILED.txt unambiguous terminal state
```

For expensive analyses, figure formatting must use cached, validated tables or
objects through a plot-only step. A cosmetic change must not silently rerun or
replace the integration.

There are two explicit output modes:

- Default versioned mode creates a new timestamped run directory and is used
  for frozen milestones.
- Working replacement mode uses `--replace-run RUN_ID`, clears and regenerates
  that exact completed run directory, and does not create another results
  folder. It is refused for targets outside the step directory and while the
  prior job is pending or running. The replacement records the prior job ID.

PNG and PDF are mandatory for UMAP renders. SVG is controlled by
`--svg true|false` and defaults to false to reduce asset accumulation. When
enabled, SVG rasterized point layers are fixed at 300 dpi; SVG text remains
editable. PNG and PDF rasterized layers are also fixed at 300 dpi.

## SLURM submission and monitoring

Initialize the Turbo directories idempotently:

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline/paper2_mge_organoid_atlas
./bin/initialize_turbo.sh
```

Inspect the first submission without creating or submitting a run:

```bash
./bin/submit.sh 00_input_audit --dry-run
```

Submit:

```bash
./bin/submit.sh 00_input_audit
```

Replace the existing completed working audit instead of making a new folder:

```bash
./bin/submit.sh 00_input_audit \
  --replace-run 00_input_audit_20260827_130337_6ce39f6 \
  --svg false
```

Enable SVG when it is specifically needed:

```bash
./bin/submit.sh 00_input_audit \
  --replace-run 00_input_audit_20260827_130337_6ce39f6 \
  --svg true
```

Versioned submission creates a new timestamped output package; replacement
submission reuses only the explicitly named completed working directory. Both
snapshot code and configuration, copy the exact submitted `.sbatch` into the
Paper 2 `jobs/` directory, record the current Git commit and working-tree
status, and capture the parsable SLURM job ID. The printed run directory is the
durable handle for monitoring:

```bash
./bin/status.sh <run_directory>
```

`status.sh` reports `sacct`, terminal markers, and recent stdout/stderr. A job
is complete only when SLURM reports successful completion, `SUCCESS.txt`
exists, required outputs exist and are nonempty, and the output README and
provenance are sufficient to reconstruct the run.

## Step 00 expected outputs

`00_input_audit` is read-only with respect to all six source objects. Its
versioned package will include:

- audited input registry with file size, modification time, and SHA-256;
- object and default-assay summaries;
- assay/layer inventory;
- metadata-column type, missingness, cardinality, and example-value inventory;
- reduction inventory;
- all-cell UMAP source/cluster-label inventory and raw/display cluster counts;
- six-study, Varela-paired, and per-study cluster UMAP figures in mandatory
  PNG/PDF and optional SVG;
- per-cell saved UMAP coordinates with original IDs, collision-safe candidate
  atlas IDs, sample provenance, raw cluster IDs, displayed names, and mapping
  sources;
- per-study default-assay feature lists;
- pairwise and six-way feature-overlap summaries;
- R/Seurat session information;
- exact input registry, R script, SLURM script, submit command, job ID, Git
  commit/status, package checksum manifest, and logs.

The audit does not decide which expression layer should be used. That decision
belongs to `02_harmonize_genes` and must cite the audit and canonical-input
evidence.

## Step 00 completed evidence and findings

The current completed audit package is:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper2_mge_organoid_atlas/results/00_input_audit/00_input_audit_20260827_130337_6ce39f6
```

The replacement run loaded all six registered objects and found 404,789 cells
with saved UMAP coordinates:

| Study | Cells | Default assay | Raw cluster IDs | Display groups in audit |
| --- | ---: | --- | ---: | ---: |
| Varela DIV30 | 90,631 | RNA | 7 | 5 paper/manual classes |
| Varela DIV90 | 22,338 | RNA | 13 | 13 `cluster_number_name` classes |
| Walsh | 4,519 | RNA | 24 | 156 observed raw-cluster/`walsh_group` combinations |
| Bershteyn 2025 | 124,583 | RNA | 9 | 9 raw IDs; no candidate name column was present |
| Bershteyn 2023 | 98,042 | RNA | 6 | 9 observed raw-cluster/`celltype` combinations |
| Siebert 2026 | 64,676 | SCT | 28 | 28 raw IDs; no candidate name column was present |

All registered cells with finite saved UMAP coordinates were retained in the
audit plots. In particular, the DIV90 audit includes current clusters 6 and 7;
their exclusion from some older published-style figures is not an input-audit
filter.

The default-assay feature union contains 54,746 genes, of which 13,965 occur
in all six objects. The two Varela objects have sample-split Seurat v5 RNA
count/data layers. Siebert defaults to SCT while the other five objects default
to RNA. These are completed audit observations, not harmonization decisions;
the expression layer, feature-identity policy, and method-compatible matrix
must be decided in Step 01.

### Varela DIV30 raw-cluster gap

The DIV30 object contains 90,631 cells and exactly seven observed raw Seurat
cluster IDs: `0, 1, 2, 3, 4, 6, 7`. Both `seurat_clusters` and
`RNA_snn_res.0.2` contain those same values and agree for all cells. There are
no missing cluster assignments. Raw cluster ID `5` is absent from the source
object and from the derived AnnData categorical levels; it was not removed by
the Paper 2 plotting code.

Do not confuse the absent raw cluster ID `5` with publication/display class 5.
Publication class `5 - MGE subpallial neurons` is present and contains 17,287
cells from raw Seurat cluster `2`. The full DIV30 mapping is:

| Raw Seurat cluster | Cells | Paper/manual display class |
| ---: | ---: | --- |
| 0 | 19,148 | 1 - Radial glia |
| 1 | 18,907 | 3 - SST+ cIN |
| 2 | 17,287 | 5 - MGE subpallial neurons |
| 3 | 15,932 | 1 - Radial glia |
| 4 | 13,750 | 4 - PV neuron precursor |
| 5 | 0 / absent | not present in the registered source object |
| 6 | 3,696 | 2 - Inhibitory progenitors |
| 7 | 1,911 | 1 - Radial glia |

The nonconsecutive raw IDs are consistent with an upstream removal, merge, or
manual relabeling before the legacy `Day30.rds` was saved, but the current
repository and Turbo legacy directory do not contain a pre-change object or
source script that establishes which event occurred or why. This unresolved
history must remain an explicit provenance limitation. Do not manufacture a
cluster 5 or infer its biology without locating an earlier object or original
analysis record.

## Step 01 canonical inputs: completed evidence

The permanent canonical input layer was successfully published at:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper2_mge_organoid_atlas/inputs/canonical
```

Downstream Paper 2 workflows must read only from this directory, not directly
from the historical source objects. For each study, it contains a minimal
current Seurat RDS and an equivalent AnnData H5AD. The H5AD `X` and Seurat RNA
`counts` layers contain the source RNA counts. Existing source RNA normalized
data are preserved as H5AD `layers["lognorm"]` and the Seurat RNA `data` layer;
no normalization was calculated. Siebert 2026 had no normalized RNA data layer
to preserve and therefore contains counts only.

| Study | Cells | Genes | Preserved expression |
| --- | ---: | ---: | --- |
| Varela DIV30 | 90,631 | 18,082 | counts + source RNA normalized data |
| Varela DIV90 | 22,338 | 18,082 | counts + source RNA normalized data |
| Walsh | 4,519 | 20,194 | counts + source RNA normalized data |
| Bershteyn 2025 | 124,583 | 45,068 | counts + source RNA normalized data |
| Bershteyn 2023 | 98,042 | 45,068 | counts + source RNA normalized data |
| Siebert 2026 | 64,676 | 32,131 | counts only |

Every pair passed exact checks for ordered cell IDs, ordered gene IDs, sparse
counts, optional normalized expression, selected cell metadata, selected gene
metadata, and absence of analysis embeddings/pairwise arrays. The Seurat files
contain only the RNA assay and have no reductions, graphs, neighbors, command
history, tools, images, scaled matrices, SCT assay, or integrated assay. No
cell/gene subsetting, normalization, gene harmonization, or cell-type
harmonization was performed.

The recorded environments are R 4.4.1, Seurat 5.1.0, SeuratObject 5.0.2,
Python 3.11.15, Scanpy 1.11.5, and AnnData 0.12.14. Exact source and canonical
file SHA-256 checksums, code/config snapshots, SLURM logs, session information,
and validation reports are stored within `inputs/canonical/`.

The successful jobs were array job `58958446` and dependent finalizer job
`58958447`; all tasks and the finalizer completed with exit `0:0`. Two earlier
tasks (`58958104` and `58958105`) were canceled before publication after the
RDS representation was clarified to require actual Seurat objects. Their
incomplete staging directory was removed; they did not publish or alter the
canonical layer.

`inputs/canonical/` is logically frozen: the build launcher refuses to
overwrite it, the `FROZEN.txt` marker records the terminal state, and checksum
manifests detect changes. Turbo's NFS export preserves project-group write
bits even after a successful `chmod`, so POSIX mode bits alone must not be
interpreted as the freeze guarantee.

## Step 02 gene-identifier audit: completed, review required

The completed report package is:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper2_mge_organoid_atlas/results/02_harmonize_genes/02_harmonize_genes_20260827_171616_2483850
```

SLURM job `58978281` completed with exit `0:0` in 50 seconds and peak batch RSS
of approximately 2.2 GiB. It read only H5AD `/var` feature metadata after
verifying the six canonical H5AD checksums. It did not load expression
matrices, modify canonical inputs, concatenate datasets, normalize, select
HVGs, create a cross-study object, or integrate anything.

The common human identity is a versionless Ensembl gene ID from GENCODE release
50 (GRCh38.p14). Mapping used exact Ensembl IDs and exact GENCODE symbols first,
then uniquely resolvable HGNC approved, previous, and alias symbols. No
case-insensitive or heuristic mapping was allowed. The package freezes the
exact GENCODE GTF and HGNC complete-set snapshot with SHA-256 checksums:

| Reference | SHA-256 |
| --- | --- |
| GENCODE 50 comprehensive GTF | `83fba3e9b03f0b8c958f3595c6c350adc55f468abf8b0e47b6d5284cfe13a453` |
| HGNC complete set retrieved 2026-08-27 | `0615a070f1628e6727953f67ad9248dd0f0ddbb16d41a7b40e06aa852fc3f448` |

### Observed feature namespaces

| Study | Input features | Namespace |
| --- | ---: | --- |
| Varela DIV30 | 18,082 | gene symbols |
| Varela DIV90 | 18,082 | gene symbols |
| Walsh | 20,194 | gene symbols |
| Bershteyn 2025 | 45,068 | 45,067 symbol-like IDs plus `THRA1/BTR` |
| Bershteyn 2023 | 45,068 | 45,067 symbol-like IDs plus `THRA1/BTR` |
| Siebert 2026 | 32,131 | mixed: 22,413 symbol-like IDs and 9,718 Ensembl gene IDs |

### Mapping outcomes

| Study | Mapped features | Unique common IDs | Ambiguous | Unresolved | Duplicate common IDs | Strict one-to-one IDs |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Varela DIV30 | 17,976 | 17,975 | 37 | 69 | 1 | 17,974 |
| Varela DIV90 | 17,976 | 17,975 | 37 | 69 | 1 | 17,974 |
| Walsh | 16,702 | 16,680 | 37 | 3,455 | 21 | 16,659 |
| Bershteyn 2025 | 24,318 | 23,577 | 80 | 20,670 | 729 | 22,848 |
| Bershteyn 2023 | 24,318 | 23,577 | 80 | 20,670 | 729 | 22,848 |
| Siebert 2026 | 31,195 | 31,113 | 194 | 742 | 78 | 31,035 |

Many unresolved Bershteyn/Walsh features are historical `AL*`/`AC*`
lncRNA-style symbols absent from current GENCODE 50 and not uniquely rescued by
the frozen HGNC aliases. They remain unresolved rather than being guessed.
Duplicate groups commonly pair a current approved symbol with a previous
symbol retained as a separate source feature, such as `MCUB|CCDC109B` mapping
to the same Ensembl identity.

### Six-way overlap reports

| Definition | Genes |
| --- | ---: |
| Exact raw source-feature-ID intersection | 14,152 |
| Mapped common identity present in all six | 14,483 |
| Strict one-to-one common identity present in all six | 14,112 |

The identity-level result includes any mapped Ensembl identity present in all
six even if one dataset has multiple source features mapping to it. The strict
result excludes such identities in every dataset. There are 371 six-way common
identities present in the identity-level report but excluded by the strict
duplicate rule. No matrix has been collapsed or subset under either policy.

### Plot-only overlap visualization extension

Plot-only job `58981144` completed with exit `0:0` in 19 seconds and peak batch
RSS below 1 GiB. It used only the cached Step 02 mapping/presence tables and did
not reopen canonical objects or read expression matrices. It added four
mandatory PNG/PDF figure pairs at 300 dpi; SVG was disabled:

- `gene_overlap_upset_identity_level`: exact study-presence patterns for all
  mapped common identities;
- `gene_overlap_upset_strict_one_to_one`: the same view after duplicate-to-one
  exclusions;
- `gene_overlap_pairwise_heatmaps`: mapped and strict pairwise intersection
  counts plus Jaccard overlap;
- `gene_mapping_coverage_and_nonoverlap`: mapping coverage, strict-set losses,
  and the distribution of identities across one through six studies.

The figures make the main catalog relationships explicit:

- 14,483 mapped identities occur in all six studies.
- 10,569 mapped identities occur only in Siebert 2026.
- 2,154 occur in both Varela datasets, both Bershteyn datasets, and Siebert,
  but not Walsh.
- 2,068 occur in Walsh, both Bershteyn datasets, and Siebert, but not either
  Varela dataset.
- 1,823 occur in both Bershteyn datasets and Siebert only.
- 1,665 occur only in the two Bershteyn datasets.
- 1,255 occur in both Varela and both Bershteyn datasets, but not Walsh or
  Siebert.

The Varela DIV30/DIV90 mapped feature sets are identical to each other
(Jaccard 1.000), as are the Bershteyn 2023/2025 sets. Among non-identical
pairs, mapped Jaccard overlap is lowest for Varela–Siebert (0.513), followed by
Walsh–Siebert (0.530) and Bershteyn–Siebert (0.601). Varela–Bershteyn is 0.760,
Varela–Walsh is 0.722, and Walsh–Bershteyn is 0.707. These values describe
reference feature catalogs after identifier mapping—not expressed-gene overlap
or biological similarity between cell populations.

Exact plotted pattern counts are stored in
`intersection_pattern_counts_identity_level.tsv` and
`intersection_pattern_counts_strict_one_to_one.tsv`. The figure manifest and
all refreshed package checksums pass. `FIGURES_SUCCESS.txt` retains
`review_stop=YES`.

The full per-feature map, ambiguity list, unresolved list, duplicate groups,
gene-presence matrix, pairwise overlaps, and both six-way intersection tables
are under `tables/`. Every package checksum passes. The run's `SUCCESS.txt`
sets `review_stop=YES`.

### Required review decision

Before any common matrix is created, choose whether downstream work should use
the 14,112 strict one-to-one intersection or define and validate an explicit
duplicate-collapse/resolution policy for the larger 14,483 identity-level
intersection. Also decide whether historical unresolved lncRNA symbols require
a separately versioned legacy-reference mapping pass. Do not advance to Step
03 until these points are reviewed.

### Step 02b legacy-ID recovery decision and evidence rules

The legacy-reference mapping pass has now been explicitly authorized, but it
is an extension of the Step 02 review—not approval to change the data. It is
implemented as `02b_legacy_gene_id_recovery` with its own versioned result
package, exact evidence registry, checksums, code/config snapshot, SLURM logs,
candidate maps, and README.

Original Cell Ranger `features.tsv.gz` files were found for Varela DIV30 and
DIV90. Both contain 37,143 rows, have the same SHA-256 checksum
`51f263c7ce685617100bb30d201411372c4494e73bd2d90f2e74942dd132c164`,
and retain stable Ensembl IDs beside symbols. Unique exact matches to these
tables can therefore be labeled `confirmed_source_feature_table`.

No exact original feature table/reference bundle has yet been found locally
for Walsh, Bershteyn 2025, Bershteyn 2023, or Siebert 2026. Available GEO/NeMO
metadata establish broad Cell Ranger/GRCh38 provenance but not the exact gene
annotation. For these datasets, exact-symbol comparisons across GENCODE 27,
32, 35, 44, and 50 are candidate recovery evidence only:

- one identical Ensembl ID across at least two releases is
  `historical_consensus_unique`;
- a unique match in only one tested release is lower-confidence
  `historical_single_release_unique`;
- multiple IDs remain `historical_reference_ambiguous`;
- names absent from all tested references remain unresolved.

Every proposed mapping and any recalculated six-way intersection is marked
`REPORT_ONLY_NOT_APPLIED`. Step 02b must not open expression matrices, modify
the frozen canonical inputs, overwrite Step 02, concatenate studies, normalize,
select HVGs, or integrate. Step 03 remains blocked until the resulting evidence
tiers and duplicate/intersection policy are reviewed.

### Step 02b completed evidence and results

The completed report package is:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper2_mge_organoid_atlas/results/02b_legacy_gene_id_recovery/02b_legacy_gene_id_recovery_20260827_181028_fd488f5
```

SLURM job `58986448` completed with exit `0:0` in 59 seconds and peak batch
RSS of approximately 1.0 GiB. All 34 package checksums pass. The job did not
open canonical objects or expression matrices and did not modify Step 02.

The audit proposed candidate IDs for 22,947 additional study-feature rows:

| Study | Newly recovered candidates | Confirmed source-table assignments | Historical consensus | Historical single-release | Remaining ambiguous | Remaining unresolved |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Varela DIV30 | 102 | 18,073 total mapped from source table | 0 | 0 | 1 | 3 |
| Varela DIV90 | 102 | 18,073 total mapped from source table | 0 | 0 | 1 | 3 |
| Walsh | 3,107 | 0 | 2,709 | 398 | 343 | 42 |
| Bershteyn 2025 | 9,455 | 0 | 8,040 | 1,415 | 866 | 10,429 |
| Bershteyn 2023 | 9,455 | 0 | 8,040 | 1,415 | 866 | 10,429 |
| Siebert 2026 | 726 | 0 | 721 | 5 | 194 | 16 |

The two Varela source tables exactly agree. Of the 18,082 canonical features
in each Varela dataset, 18,073 have a unique exact symbol-to-Ensembl match in
the source feature table. Five features retain a Step 02 current-reference
fallback because their source symbols map to more than one original Ensembl
ID; `GOLGA8M` remains source-ambiguous; and the suffixed feature names
`TBCE.1`, `HSPA14.1`, and `TMSB15B.1` remain unresolved because the symbol-only
Seurat representation no longer retains which duplicate source row they came
from.

Seventeen mapped Varela feature names per dataset point to a different
original Ensembl ID than the same symbol does in GENCODE 50. These are reported
as `reassigned_by_confirmed_source_table`, not silently applied. This directly
demonstrates why current-symbol matching alone cannot always reconstruct an
older feature identity.

The proposed six-way overlap changes are modest despite the large number of
recovered study-feature rows:

| Definition | Step 02 | Step 02b proposed | Net change |
| --- | ---: | ---: | ---: |
| Identity-level intersection | 14,483 | 14,496 | +13 |
| Strict one-to-one intersection | 14,112 | 14,122 | +10 |

At identity level, 21 Ensembl IDs enter and 8 leave because confirmed Varela
source identities supersede some current-symbol assignments, for a net gain of
13. Under the strict rule, 20 enter and 10 leave, for a net gain of 10. Most of
the 22,947 recovered rows are not present in every study, so they improve
per-study identity recovery without materially enlarging the six-way overlap.
Neither proposed intersection has been adopted.

## Resume point

Current state as of 2026-08-27:

- Paper 2 scope is defined.
- Six candidate processed objects are registered and exist on Turbo.
- The code and Turbo directory scaffold are established.
- Step `00_input_audit` was successfully replaced in place from the six source
  objects as SLURM job `58956196` (`COMPLETED`, exit `0:0`, elapsed
  `00:05:44`, peak batch RSS approximately 60.3 GiB).
- Its current completed run directory is
  `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper2_mge_organoid_atlas/results/00_input_audit/00_input_audit_20260827_130337_6ce39f6`.
- The replacement configuration records `OUTPUT_MODE=replace`,
  `REPLACED_PREVIOUS_JOB_ID=58955368`, `PNG_DPI=300`, `PDF_DPI=300`,
  `SVG_DPI=300`, and `MAKE_SVG=false`.
- The run has a `SUCCESS.txt` marker, complete audit tables, 8 PNG files,
  8 PDF files, no SVG files, and 49 verified package checksum entries.
- All eight PNG files report embedded resolution of approximately
  299.9994 dpi, and every figure-manifest row records PNG/PDF at 300 dpi.
- Exactly one Step 00 results directory exists; the replacement did not create
  another timestamped output directory.
- The Varela paired UMAP was visually checked: DIV30 contains 90,631 cells,
  seven raw clusters mapped to five paper/manual classes; DIV90 contains
  22,338 cells and all 13 current `cluster_number_name` clusters, including
  clusters 6/7, with the documented plotting-only vertical flip.
- Step `01_canonical_inputs` completed successfully for all six studies.
- Twelve permanent analysis-input files (six minimal Seurat RDS and six
  equivalent AnnData H5AD) are published and frozen under
  `inputs/canonical/` with exact pairwise equivalence marked `PASS`.
- Step `02_harmonize_genes` report job `58978281` completed successfully and
  produced verified namespace, mapping, ambiguity, unresolved-feature,
  duplicate, and six-way-overlap reports without touching expression data.
- The reported six-way intersections are 14,483 mapped identity-level genes
  and 14,112 strict one-to-one genes; no choice between them has been applied.
- Step `02b_legacy_gene_id_recovery` job `58986448` completed successfully and
  produced a fully checksummed, report-only candidate map from confirmed
  Varela feature tables and five historical GENCODE releases.
- Step 02b proposes 14,496 identity-level and 14,122 strict six-way genes, but
  explicitly marks both sets `REPORT_ONLY_NOT_APPLIED`; the canonical inputs
  and parent Step 02 mapping remain unchanged.
- No cross-study master dataset has yet been created or frozen.
- No integration method has been run for Paper 2.

The next safe action is review of the completed Step 02 and Step 02b reports.
The workflow must remain stopped until the confirmed-source versus historical-
candidate evidence tiers and the strict-vs-duplicate-resolution intersection
policy are explicitly chosen. Do not create a common matrix or begin Step
`03_harmonize_metadata` yet.

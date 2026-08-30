# GSE94641 rapid MGE reference mapping

## Current status: E15.5-focused label transfer complete

This is a small, independent Step 00 preprocessing workstream for broad
developmental MGE annotation. It does not consume or modify the frozen PCDH19
genotype-classification Steps 03–07. GSE94641 has been downloaded and
validated, and its published labels have now been transferred to all 450,788
Paper 3 cells. The primary calls use E15.5 reference neighbors; the full
E11.5–E17.5 reference, all raw label fields, and an all-age contextual transfer
are retained alongside them.

This completed transfer is not a primary processed or integrated Paper 3
analysis. The next true cleanup is the standalone `primary_processing`
workflow with its own Steps 00–09, starting from the original per-sample Cell
Ranger outputs rather than this mapping package. These GSE94641 calls may be
joined only later as annotation evidence. Its input/output and
decision-to-script plan is in
[`primary_processing/README.md`](primary_processing/README.md); it has not yet
been implemented or run.

GSE94641 was selected because it contains published single-cell RNA-seq from
actual mouse embryonic medial ganglionic eminence at E11.5, E13.5, E15.5, and
E17.5. Its E15.5 cells approximately match the Paper 3 developmental stage,
and its author-supplied fields allow a rapid reference-based annotation without
reclustering or reinventing identities.

A separate sibling workstream completed its first checkpoint for curating La
Manno 2021, Bandler 2022, and Mayer 2018 on 2026-08-29 before selecting one
candidate for possible future reference use. See
[`PCDH19_DEVELOPING_MOUSE_MGE_REFERENCE_CURATION_HANDOFF.md`](PCDH19_DEVELOPING_MOUSE_MGE_REFERENCE_CURATION_HANDOFF.md).
That audit does not replace, reopen, or modify this completed GSE94641 package,
and its first checkpoint separately inventories published versus P0-linked
sample IDs, ages, tissues, library/design metadata, and remaining curation
requirements; it stops before loading or mapping the Paper 3 query.

## Validation result

| Property | Validated result |
| --- | --- |
| Reference cells | 225 |
| Ages | E11.5: 96; E13.5: 48; E15.5: 63; E17.5: 18 |
| Expression dimensions | 36,528 gene rows × 225 cells |
| Processed measurements | size-factor-adjusted nRPKM and uniquely aligning read count |
| Gene IDs | numeric IDs; GEO processing metadata explicitly identifies these as Entrez Gene IDs |
| Symbols directly supplied | No. GPL17021 is a virtual platform and this deposit has no platform annotation table. |
| Query cells loaded | 450,788 across all 12 samples |
| Label transfer run | Yes: E15.5 primary plus all-age context |

The 225 processed files have identical, unique gene-ID rows in identical order.
Every metadata GSM accession matches exactly one processed file. All locked
source identities and all validation checks pass. FASTQs were not downloaded.

## Published annotation fields and actual categories

GEO supplies six `Sample_characteristics_ch1` fields for every cell: `tissue`,
`age`, `extendedphenotype`, `subtype`, `nonneuralcelltype`, and `date`. The
complete per-cell table also preserves every other published `Sample_*` field,
including author cell ID/title, GEO accession, protocols, platform, BioSample,
SRA, and processed-file URL. `date` is an acquisition/preparation date and
should be retained as technical metadata, not treated as a developmental
identity.

### Broad and fine neural labels

`extendedphenotype` is the authors' broad neural developmental program:
`proliferating` corresponds to proliferating neural progenitors and `maturing`
corresponds to post-mitotic immature neurons, not mature neurons. The paper
further characterizes VZ/SVZ progenitors and MGE-derived, mixed MGE/LGE, and
LGE-derived immature-neuron populations. `subtype` is a finer published
partition nested perfectly within that field: `pro1`–`pro5` occur in the
proliferating group and `mat1`–`mat5` in the maturing group. These numbered
subtypes are published transcriptional groups, not standard anatomical cell
type names; they are preserved verbatim rather than renamed. They must not be
interpreted as an ordered maturation scale. The authors did not resolve future
PV versus SST interneuron identity.

| Field | Published label | All cells | E15.5 cells |
| --- | --- | ---: | ---: |
| extendedphenotype | maturing | 133 | 39 |
| extendedphenotype | proliferating | 74 | 12 |
| extendedphenotype | NA | 18 | 12 |
| subtype | mat1 | 16 | 2 |
| subtype | mat2 | 25 | 10 |
| subtype | mat3 | 24 | 3 |
| subtype | mat4 | 8 | 1 |
| subtype | mat5 | 60 | 23 |
| subtype | pro1 | 1 | 0 |
| subtype | pro2 | 7 | 2 |
| subtype | pro3 | 17 | 6 |
| subtype | pro4 | 40 | 4 |
| subtype | pro5 | 9 | 0 |
| subtype | NA | 18 | 12 |

### Non-neural labels

`nonneuralcelltype` is a separate author-supplied identity field that flags
small non-neural or glial/ependymal populations. It is not a clean inverse of
`extendedphenotype`: some cells with one of these labels also retain a neural
extended-phenotype assignment. `NA` therefore means no non-neural label was
assigned, not a new biological class. The mixed label `astrocytes Epend` is
preserved exactly as deposited.

| Published label | All cells | E15.5 cells |
| --- | ---: | ---: |
| NA | 210 | 59 |
| Epend | 3 | 0 |
| Vsmc | 3 | 1 |
| astrocytes | 2 | 1 |
| astrocytes Epend | 2 | 2 |
| endothelial | 1 | 0 |
| microglia | 4 | 0 |

Machine-readable totals for all six characteristics are in
`GSE94641_annotation_counts.tsv`; the complete field-by-label-by-age result is
in `GSE94641_annotation_by_age.tsv`.

## Mapping design

All 225 reference cells and every published metadata field remain in the
reference package. Current NCBI mouse gene information maps deposited Entrez
Gene IDs to query Ensembl IDs; 18,851 unique one-to-one genes overlap. Reference
uniquely aligning read counts and query Flex UMIs are independently
library-size normalized to 10,000 and transformed with `log1p`. The 2,000 most
variable shared reference genes define a 20-component PCA fitted to all four
reference ages. Query cells are projected into that fixed space.

The primary neighbor pool is restricted to the 63 E15.5 cells and uses fixed
`k=5`, selected because `k=20` is too large for this small stage-matched pool.
The exact composite of `extendedphenotype`, `subtype`, and
`nonneuralcelltype`, as well as each constituent field, is transferred by
majority vote. Confidence is the fraction of neighbors supporting the winning
label. A second contextual result uses all 225 cells with fixed `k=20` and also
transfers neighbor age. No k optimization, integration, clustering, or manual
query annotation is performed.

The raw per-field votes and exact composite-tuple vote are intentionally
separate. Because a tuple majority need not equal each marginal majority, they
disagree for 4.68% of cells. Use `GSE94641_broad_state` for the broad neural
state, the individual raw fields for field-specific work, and `GSE94641_label`
only when the exact published-label tuple is required.

## Actual transferred annotations

All 450,788 Cell Ranger filtered cells were mapped and joined to their existing
per-sample vendor UMAP by exact barcode. The E15.5-primary broad result is:

| Sample | Post-mitotic immature neuron | Proliferating neural progenitor | Not assigned neural state |
| --- | ---: | ---: | ---: |
| JZ-1 | 48.02% | 44.95% | 7.03% |
| JZ-2 | 48.69% | 39.85% | 11.47% |
| JZ-3 | 62.08% | 26.11% | 11.80% |
| JZ-4 | 40.90% | 51.58% | 7.51% |
| JZ-5 | 53.26% | 42.30% | 4.44% |
| JZ-6 | 47.62% | 48.17% | 4.21% |
| JZ-7 | 38.82% | 53.74% | 7.44% |
| JZ-8 | 46.41% | 47.40% | 6.19% |
| JZ-9 | 44.44% | 50.89% | 4.67% |
| JZ-10 | 61.38% | 32.78% | 5.84% |
| JZ-11 | 53.08% | 39.10% | 7.82% |
| JZ-12 | 38.05% | 55.28% | 6.67% |
| **Pooled** | **49.91% (225,001)** | **43.38% (195,564)** | **6.70% (30,223)** |

The E15.5 subtype calls are dominated by `mat5` (222,233; 49.30%) and `pro3`
(190,665; 42.30%), with `NA` at 34,071 (7.56%). Rare E15.5 reference labels
are rarely or never majority calls, so these fine subtypes should be treated as
exploratory raw published codes. The median exact-composite vote confidence is
0.60 (mean 0.607); the median broad-state confidence is 0.60 (mean 0.677).

An E15.5 reference leave-one-cell-out diagnostic using the same `k=5` recovered
the published broad state for 56/63 cells (88.89%), subtype for 48/63 (76.19%),
non-neural field for 57/63 (90.48%), and exact tuple for 45/63 (71.43%). This is
an internal small-reference sanity check, not independent biological
validation.

The all-age contextual broad call agrees with the E15.5-primary broad call for
344,309/450,788 cells (76.38%). Its modal neighbor age is E11.5 for 54.64% of
query cells, E15.5 for 37.38%, E13.5 for 7.97%, and E17.5 for 0.007%. This age
shift is why E15.5 remains the declared primary result while the complete
four-age result is retained for context.

## Files and provenance

Downloaded assets are cached under
`paper3_pcdh19/references/GSE94641/source_files/`; the 225 extracted processed
files are under `processed_expression/` in the same reference directory.

| Downloaded file | Purpose | SHA-256 |
| --- | --- | --- |
| `GSE94641_RAW.tar` | 225 processed per-cell expression files | `f237a6c3e12814d9cb8351d145f7b5df7fb0115b03f23759d246e6c45e03bdad` |
| `filelist.txt` | GEO supplementary inventory and byte sizes | `bebbeafd3390032571c4187ec52745dc4f5740659c5dc4c064a0dd710232d6be` |
| `GSE94641_series_matrix.txt.gz` | series metadata and all per-cell annotations | `61cd11ee985b7af87046469ebea29956cb0fd3fd0ddb84b03192e9e18e0529e6` |
| `GSE94641_family.soft.gz` | family-level series metadata | `6a1fb011296442c7f45df854eca92e26bbab221fbae7a37e062b407354cf70f0` |
| `GPL17021_family.soft.gz` | platform metadata and direct-symbol availability check | `8c0c082d2ee55cfd6ce395502e64a944db664e2ca2d57248b5a0ab2f785f9074` |

Validation outputs are atomically published under:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/
  results/mge_reference_mapping_gse94641/reference_validation/
    GSE94641_cell_metadata.tsv
    GSE94641_published_cell_fields.tsv
    GSE94641_annotation_counts.tsv
    GSE94641_annotation_by_age.tsv
    GSE94641_expression_file_manifest.tsv
    GSE94641_gene_identifiers.tsv
    GSE94641_reference_summary.tsv
    source_manifest.tsv
    validation_checks.tsv
    output_manifest.tsv
```

The output-manifest SHA-256 is
`41cf973380a7d9349092a6629e68493058b77bb3d4a6d5403328be010e7ec9e00`.
The checksum lock is
[`config/gse94641_reference_validation.lock.json`](config/gse94641_reference_validation.lock.json),
and the modular implementation is
[`scripts/reference_mapping/01_download_validate_GSE94641.py`](scripts/reference_mapping/01_download_validate_GSE94641.py).

Run or re-verify locally with:

```bash
./paper3_pcdh19/bin/run_step_00_gse94641_reference_validation.sh
```

The corresponding SLURM entry point is:

```bash
sbatch paper3_pcdh19/slurm/step_00_gse94641_reference_validation.sbatch
```

## Mapping outputs

The atomically published mapping package is:

```text
results/mge_reference_mapping_gse94641/label_transfer_e15_5/
  gse94641_query_cell_label_transfer.tsv
  gse94641_reference_cells_with_pca.tsv
  gse94641_shared_selected_genes.tsv
  gse94641_pca_variance.tsv
  gse94641_transferred_label_summary_by_sample.tsv
  mapping_configuration.tsv
  query_input_manifest.tsv
  software_environment.tsv
  figures/
    gse94641_reference_pca_published_labels.png
    query_vendor_umap_gse94641_e15_5_broad_state.png
    query_vendor_umap_gse94641_e15_5_confidence.png
    gse94641_e15_5_transferred_broad_state_proportions_by_sample.png
  output_manifest.tsv
```

The mapping output-manifest SHA-256 is
`ca214c2c65d6aaac6235b7b8e66b2388fbfb2e9f1ac0369de0e9bb0e53280b11`.
The 197 MB per-cell table preserves sample/barcode, sample design metadata,
existing vendor UMAP, query library UMIs, E15.5-primary raw labels and
confidence, broad biological interpretation, distances, and all-age contextual
labels/age/confidence.

Run or checksum-verify the mapping with:

```bash
./paper3_pcdh19/bin/run_step_00_gse94641_e15_5_label_transfer.sh
sbatch paper3_pcdh19/slurm/step_00_gse94641_e15_5_label_transfer.sbatch
```

## Architecture and scope boundary

Reference validation separates `LockedReferenceDownloader`,
`GeoMetadataParser`, `ProcessedExpressionValidator`,
`PlatformAnnotationInspector`, and `ReferenceValidationPublisher`. Mapping adds
`LockedGeneAnnotation`, `ReferencePreprocessor`, `QuerySampleLoader`,
`KNNLabelTransfer`, `MappingPlotter`, and `OutputPublisher` under the same
workflow rather than creating a parallel architecture.

This remains rapid reference-based annotation, not de novo cell-type
discovery. The reference uses Fluidigm C1 read counts whereas the query uses
10x Flex UMIs, and no batch integration was requested; absolute distances and
fine subtype assignments therefore require caution. `not_assigned_neural_state`
is absence of an `extendedphenotype` majority, not a novel cell type. Rare
non-neural labels are poorly supported by only 63 E15.5 reference cells.

PCDH19 WT-like/KO-like classifications were not opened or used as mapping
features. The mapping output is an independent annotation layer joined only by
sample and exact barcode and cannot retroactively modify the frozen genotype
classifier.

## Downstream Step 08 use

Formal Step 08 consumes `GSE94641_broad_state` as an immutable annotation and
joins it by exact sample/barcode to the frozen Step 02a A/B/C probe audit. It
preserves every primary and all-age mapping field in a compressed joined
per-cell table but performs biological-sample comparisons only at the three
broad-state levels. See
[`PCDH19_DEVELOPMENTAL_STATE_PROBE_DETECTABILITY.md`](PCDH19_DEVELOPMENTAL_STATE_PROBE_DETECTABILITY.md).

This downstream use does not revise the mapping, promote cryptic `mat*`/`pro*`
codes to cell-type names, or alter the frozen PCDH19 classifier.

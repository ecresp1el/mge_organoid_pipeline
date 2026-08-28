# Paper 3 Ziobro PCDH19 MGE single-cell RNA-seq: operational handoff

Last updated: 2026-08-27

## Correct allocation and source

The current workspace account belongs to these relevant Unix groups:

```text
umms-parent-turbo
umms-ziobroj-turbo
```

Paper 3 belongs on the second allocation:

```text
Ziobro Turbo root:
/nfs/turbo/umms-ziobroj

Ziobro lab directory:
/nfs/turbo/umms-ziobroj/Ziobro Lab

PCDH19/MGE source directory:
/nfs/turbo/umms-ziobroj/Ziobro Lab/MGE scRNA-seq
```

The current session can read and traverse these paths. The source directories
are group-owned by `umms-ziobroj-turbo`. No source file has been modified,
copied, renamed, or checksummed during discovery.

Paper 3 reads this Ziobro source in place but writes all workflow products to
the existing `umms-parent` MGE project, using the same `parent0` compute account
and output conventions as Paper 2:

```text
Great Lakes output root:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19

Workstation mount equivalent:
/Volumes/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19
```

## Rejected false lead

The following directories are Miranda bulk RNA-seq experiments and were
explicitly rejected by the user as Julia Ziobro Paper 3 data:

```text
/nfs/turbo/umms-parent/Miranda/Miranda RNAseq/01.30.25_PCDH19 D20 FACS sorted
/nfs/turbo/umms-parent/Miranda/Miranda RNAseq/PCDH19 bulk RNAseq_D20
```

They must not be reintroduced into Paper 3 unless the user later changes the
scope explicitly.

## Dataset identity

The correct source is AGC service request `15662-JZ`. Its README identifies it
as 10x Flex v2 data generated on the Ultima Genomics UG 100 using sequencing
recipe version 1.9.0.2.

The Cell Ranger multi configuration records:

- organism/reference: mouse, GRCm39;
- probe set: Chromium Mouse Transcriptome Probe Set v2.0.0,
  `GRCm39-2024-A_EGFP`;
- library FASTQ ID: `15662-JZ-P01-CNVT01`;
- feature type: Gene Expression;
- 12 samples mapped to probe barcodes `B-A03` through `B-H03` and `B-A04`
  through `B-D04`.

The source contains two pooled FASTQs, approximately 212 GB and 237 GB, plus
pooled Cell Ranger products and 12 demultiplexed sample-level filtered
feature-barcode matrices. The delivered sample matrices each contain 19,071
features.

## Loupe inspection

There are two pooled `raw_cloupe.cloupe` copies of approximately 5.3 GB each
and 12 distinct per-sample `sample_cloupe.cloupe` files. The per-sample files
also have byte-size-identical copies in the AGC Dropbox delivery.

The files use Loupe format 9.0. Direct inspection of the embedded index and run
metadata found:

- run name `Sample_15662-JZ-P01`;
- an empty run `Description`;
- GRCm39 and Chromium Mouse Transcriptome Probe Set v2.0.0;
- 450,788 barcodes in the pooled matrix;
- vendor UMAP and t-SNE projections;
- graph-based and k-means cluster assignments;
- no custom cell-data tables; and
- no strings or metadata fields identifying PCDH19 genotype, condition,
  replicate, or other biological sample identities.

Therefore the `.cloupe` files are useful for viewing the vendor expression and
clustering results, but they were not the source of the biological sample key.
The key was subsequently supplied by the user and is registered separately in
[`config/sample_key.csv`](config/sample_key.csv).

## Preliminary technical metrics

The Cell Ranger summary reports 450,788 called cells in the physical library
and 5,387,237,349 analyzed reads. Sample-level values transcribed into
[`config/sample_manifest_draft.tsv`](config/sample_manifest_draft.tsv) include:

| Sample | Probe barcode | Called cells | Median UMIs/cell | Median genes/cell |
| --- | --- | ---: | ---: | ---: |
| `15662-JZ-1` | `B-A03` | 51,229 | 4,076 | 2,622 |
| `15662-JZ-2` | `B-B03` | 37,553 | 3,529 | 2,373 |
| `15662-JZ-3` | `B-C03` | 25,354 | 2,978 | 2,042 |
| `15662-JZ-4` | `B-D03` | 21,440 | 3,895 | 2,542 |
| `15662-JZ-5` | `B-E03` | 41,878 | 4,050 | 2,611 |
| `15662-JZ-6` | `B-F03` | 56,099 | 4,227 | 2,696 |
| `15662-JZ-7` | `B-G03` | 19,623 | 3,828 | 2,486 |
| `15662-JZ-8` | `B-H03` | 20,799 | 3,686 | 2,432 |
| `15662-JZ-9` | `B-A04` | 60,680 | 4,591 | 2,867 |
| `15662-JZ-10` | `B-B04` | 52,753 | 4,050 | 2,637 |
| `15662-JZ-11` | `B-C04` | 52,295 | 4,366 | 2,782 |
| `15662-JZ-12` | `B-D04` | 11,085 | 4,453 | 2,804 |

These are vendor pipeline calls, not final analysis-QC decisions. The wide
cell-count range, especially sample 12, must be evaluated against the
registered biological design during formal QC.

## Registered biological sample key

The inspected AGC configuration and QC files map technical sample IDs to probe
barcodes but leave `Sample description` blank. On 2026-08-27, the user supplied
the biological mapping below. Its canonical machine-readable registration is
[`config/sample_key.csv`](config/sample_key.csv); its provenance is recorded
there as `user_provided_in_codex_thread`, not as metadata recovered from the
AGC delivery.

| Technical sample | Submitted name | Target cells | Organism | Tissue/region | Genotype | Sex | Design group |
| --- | --- | ---: | --- | --- | --- | --- | --- |
| `15662-JZ-1` | `Sample 1.1` | 20,000 | Mouse | Embryonic brain, MGE | WT | M | `WT_M` |
| `15662-JZ-2` | `Sample 1.7` | 20,000 | Mouse | Embryonic brain, MGE | WT | M | `WT_M` |
| `15662-JZ-3` | `Sample 2.1` | 20,000 | Mouse | Embryonic brain, MGE | WT | M | `WT_M` |
| `15662-JZ-4` | `Sample 1.4` | 20,000 | Mouse | Embryonic brain, MGE | WT | F | `WT_F` |
| `15662-JZ-5` | `Sample 1.5` | 20,000 | Mouse | Embryonic brain, MGE | WT | F | `WT_F` |
| `15662-JZ-6` | `Sample 2.5` | 20,000 | Mouse | Embryonic brain, MGE | WT | F | `WT_F` |
| `15662-JZ-7` | `Sample 1.3` | 20,000 | Mouse | Embryonic brain, MGE | HET | F | `HET_F` |
| `15662-JZ-8` | `Sample 1.6` | 20,000 | Mouse | Embryonic brain, MGE | HET | F | `HET_F` |
| `15662-JZ-9` | `Sample 2.4` | 20,000 | Mouse | Embryonic brain, MGE | HET | F | `HET_F` |
| `15662-JZ-10` | `Sample 1.2` | 20,000 | Mouse | Embryonic brain, MGE | KO | M | `KO_M` |
| `15662-JZ-11` | `Sample 1.8` | 20,000 | Mouse | Embryonic brain, MGE | KO | M | `KO_M` |
| `15662-JZ-12` | `Sample 2.2` | 20,000 | Mouse | Embryonic brain, MGE | KO | M | `KO_M` |

The group-level allocation is three WT males, three WT females, three HET
females, and three KO males. This gives two direct sex-matched genotype
contrasts: HET female versus WT female, and KO male versus WT male. It is not a
complete genotype-by-sex factorial because HET males and KO females are absent.
Do not yet assume that the three submitted samples per group are independent
biological replicates: donor/embryo/litter, cell line, differentiation batch,
age/time point, mosaic fraction, and the meaning of the `1.x`/`2.x` submitted
names still require confirmation.

## Authoritative configuration

- [`config/input_candidates.tsv`](config/input_candidates.tsv) records the
  correct source root and technical format.
- [`config/sample_manifest_draft.tsv`](config/sample_manifest_draft.tsv)
  contains only verified technical IDs, probe barcodes, and vendor QC metrics.
- [`config/sample_key.csv`](config/sample_key.csv) is the authoritative
  user-provided mapping from technical IDs to submitted names, target cells,
  organism, tissue/region, genotype, sex, and four design groups. Its current
  SHA-256 is
  `5b20e8596c28f95b6adedcb0fe17019d2a6db46384f18f02d4bfdc58803f1dc8`.
- [`config/greatlakes.env`](config/greatlakes.env) keeps the source on the
  Ziobro allocation and points Paper 3 outputs to the existing `umms-parent`
  MGE project.
- [`PIPELINE_IO_AND_BIOLOGICAL_SCOPE.md`](PIPELINE_IO_AND_BIOLOGICAL_SCOPE.md)
  is the authoritative file-by-file contract for every executable, scientific
  input, output column, validation, and permitted biological use in the locked
  probe audit.
- [`XGFP_PROBE_COMPATIBILITY_AUDIT.md`](XGFP_PROBE_COMPATIBILITY_AUDIT.md)
  records the separate construct-level forensic audit of the custom Flex EGFP
  probes against the Nagy/Kalantry D4/XEGFP reporter.

## Locked X-GFP/Flex sequence-compatibility result

Before GFP counts or biological interpretation, Step
`02b_xgfp_probe_compatibility_audit` established whether the exact custom Flex
panel could recognize Julie's Nagy/Kalantry X-GFP reporter by sequence.

The reporter is the Kalantry laboratory X-GFP line identified in their methods
as D4/XEGFP: `Tg(CAG-EGFP)D4Nagy`, MGI `3055027`, JAX stock `003116`. The
original line derives from pCX-EGFP-expressing R1 ES cells. That construct uses
Clontech enhanced GFP cDNA under the CMV immediate-early enhancer/chicken
beta-actin promoter with rabbit beta-globin polyadenylation sequence.

The delivered `probe_set.csv` contains exactly three included EGFP probes. The
full delivered IDs, sequences, and exact construct-reference alignments are:

| EGFP 5-prime order | Full probe ID | Probe sequence, 5-prime to 3-prime | Clontech EGFP CDS interval | Alignment |
| ---: | --- | --- | ---: | --- |
| 1 | `EGFP|EGFP|probe03` | `ATGGTGCGCTCCTGGACGTAGCCTTCGGGCATGGCGGACTTGAAGAAGTC` | 247--296 | Unique reverse-complement match; 50/50 bases, 0 mismatches, 0 gaps |
| 2 | `EGFP|EGFP|probe02` | `AGGGTGTCGCCCTCGAACTTCACCTCGGCGCGGGTCTTGTAGTTGCCGTC` | 310--359 | Unique reverse-complement match; 50/50 bases, 0 mismatches, 0 gaps |
| 3 | `EGFP|EGFP|probe01` | `GGTAGTGGTCGGCGAGCTGCACGCTGCCGTCCTCGATGTTGTGGCGGATC` | 501--550 | Unique reverse-complement match; 50/50 bases, 0 mismatches, 0 gaps |

All three delivered rows have `included=TRUE` and the 10x annotation
`region=unspliced`. The exact reference is the Clontech-submitted EGFP CDS in
NCBI GenBank `U55762.1`, bases 679--1398. The 720-base CDS SHA-256 is
`e8d734949c2b74ac68d9ecd5aa1016aa75d496d292b124db7d4993f2e5409449`.

The defensible conclusion is that sequence compatibility passes: the custom
Flex panel was theoretically capable of detecting RNA from this reporter.
Sequence incompatibility cannot explain absent GFP counts. This does not prove
reporter transcription, RNA preservation, probe efficiency, or an observed
GFP UMI in any sample or cell.

The exact full D4/XEGFP integrated concatemer and mouse/transgene junction
sequence is not deposited in the cited line records; MGI lists the precise
insertion location as unknown. The exact-match claim is therefore at the
original construct reporter-CDS level, not colony-specific resequencing of
every integrated copy.

Reproduce or verify with:

```bash
./paper3_pcdh19/bin/run_xgfp_probe_audit.sh
```

Validated machine-readable results are under:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/xgfp_probe_audit
```

The output manifest SHA-256 is
`3d4a3affc77aa844cc779ce947ba8805a06261360b0ec4513435833d5d540912`.
No GFP counts, barcodes, UMAP, or cell classification were analyzed in this
step.

## Locked Pcdh19 probe identities and genomic targets

These are the three and only three Pcdh19 probes present in the delivered
Chromium Mouse Transcriptome Probe Set v2.0.0, GRCm39, reference 2024-A. The
full ID—not only the seven-character suffix—is the identity used for every
join and validation. Coordinates below are GRCm39 BED convention: 0-based,
half-open. The corresponding 1-based inclusive interval is provided to avoid
coordinate-system ambiguity.

### A: exon-1 probe `a3f4e22`

- Full ID: `ENSMUSG00000051323|Pcdh19|a3f4e22`
- Probe sequence: `TGTAGGTCCGGATCTCTTTGTTGTCACGCTTGCATTTGATTGCCACGAAG`
- 10x annotation: `region=unspliced`; this label is retained verbatim and is
  not interpreted as meaning intronic.
- BED target: `chrX:132585873-132585923`, `+` BED strand; equivalently
  `chrX:132585874-132585923` in 1-based inclusive coordinates.
- BED12 blocks: one 50-bp block; `blockCount=1`, `blockSizes=50`,
  `blockStarts=0`.
- `transcript_id_set`: `ENSMUST00000060309`, `ENSMUST00000149154`, and
  `ENSMUST00000167944`.
- GTF intersection: the full 50-bp target lies within exon 1 of every targeted
  transcript. It is entirely exonic.
- Position along Pcdh19: first of the three probes in 5-prime-to-3-prime gene
  order.

### B: downstream probe `8215225`

- Full ID: `ENSMUSG00000051323|Pcdh19|8215225`
- Probe sequence: `AAGTAAGGGAGGAGCAACTGACAACATTCATCTTGTCTGTCTCCTCCACG`
- 10x annotation: `region=unspliced`; this does not mean intronic.
- BED target: `chrX:132582020-132582070`, `+` BED strand; equivalently
  `chrX:132582021-132582070` in 1-based inclusive coordinates.
- BED12 blocks: one 50-bp block; `blockCount=1`, `blockSizes=50`,
  `blockStarts=0`.
- `transcript_id_set`: `ENSMUST00000060309`, `ENSMUST00000149154`, and
  `ENSMUST00000167944`.
- GTF intersection: the full target lies within exon 2 of
  `ENSMUST00000060309`, exon 3 of `ENSMUST00000149154`, and exon 2 of
  `ENSMUST00000167944`. It is entirely exonic for every targeted transcript.
- Additional annotation detail: the interval has a 2-bp overlap with exon 3
  of non-targeted Pcdh19 transcript `ENSMUST00000193485`; that transcript is
  not in the matching 10x `transcript_id_set` and does not change the target
  assignment above.
- Position along Pcdh19: second of the three probes in 5-prime-to-3-prime gene
  order.

### C: downstream probe `d013e0b`

- Full ID: `ENSMUSG00000051323|Pcdh19|d013e0b`
- Probe sequence: `AAGTCACACTGGTGTTCAGGACATCATTCACAGCAGTATCGCAGTACAGG`
- 10x annotation: `region=unspliced`; this does not mean intronic.
- BED target: `chrX:132526044-132526094`, `+` BED strand; equivalently
  `chrX:132526045-132526094` in 1-based inclusive coordinates.
- BED12 blocks: one 50-bp block; `blockCount=1`, `blockSizes=50`,
  `blockStarts=0`.
- `transcript_id_set`: `ENSMUST00000060309`, `ENSMUST00000149154`, and
  `ENSMUST00000167944`.
- GTF intersection: the full target lies within exon 4 of
  `ENSMUST00000060309`, exon 5 of `ENSMUST00000149154`, and exon 4 of
  `ENSMUST00000167944`. It is entirely exonic for every targeted transcript.
- Position along Pcdh19: third of the three probes in 5-prime-to-3-prime gene
  order.

Pcdh19 is annotated on the `-` gene strand in the GRCm39-2024-A GTF, whereas
the 10x BED records `+` in its strand column for each probe. Both values are
preserved as supplied. Because the gene is negative-strand, decreasing GRCm39
coordinate gives the validated 5-prime-to-3-prime order `A -> B -> C`.

### How the target assignments were established

1. The delivered `probe_set.csv` header was required to match panel name
   `Chromium Mouse Transcriptome Probe Set v2.0.0`, genome `GRCm39`, and
   reference version `2024-A`.
2. Rows with `gene_id=ENSMUSG00000051323` were extracted. The program required
   the resulting full-ID set to equal exactly the three IDs above—no missing,
   additional, or duplicate Pcdh19 probes.
3. The exact matching 10x v2.0.0/GRCm39/2024-A BED and metadata TSV were
   checksum-validated. Their Pcdh19 full-ID sets were independently required
   to equal the delivered three-probe set exactly.
4. Sequences and `region` values were required to agree between the delivered
   panel and 10x metadata; coordinates and BED12 blocks came only from the
   matching BED; `transcript_id_set` came only from the matching metadata.
5. Each 50-bp BED interval was intersected with Pcdh19 exons in the exact
   GRCm39-2024-A GTF. For every transcript named by 10x, the workflow required
   the entire probe interval to fall inside exactly one exon and required the
   observed exon number to match the frozen lock.
6. The 5-prime order was derived from genomic coordinates plus the GTF Pcdh19
   `-` strand and was required to equal the frozen `A`, `B`, `C` order.

The exact reference inputs and SHA-256 values are:

| Input | Resolved path | SHA-256 |
| --- | --- | --- |
| Delivered `probe_set.csv` | `/nfs/turbo/umms-ziobroj/Ziobro Lab/MGE scRNA-seq/10x_analysis_15662-JZ/Sample_15662-JZ-P01/probe_set.csv` | `70368b977ffd3d7ce183ebf82581d4bf0357c3c4c4bb18d8a99b013f9c77b8f1` |
| Matching 10x BED | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/inputs/pcdh19_probe_audit/references/Chromium_Mouse_Transcriptome_Probe_Set_v2.0.0_GRCm39-2024-A.bed` | `2e6ac3b6029ae2aa7478f9da5cb44e27611bb9e9fde8f2bf307d07a8e6302440` |
| Matching 10x metadata | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/inputs/pcdh19_probe_audit/references/Chromium_Mouse_Transcriptome_Probe_Set_v2.0.0_GRCm39-2024-A.probe_metadata.tsv` | `68d720efa907394ec3dddc5a8b503cf5917d720b911225fefaf28af387c09c54` |
| GRCm39 GTF | `/nfs/turbo/agc-data/refs/Mus_musculus/GRCm39/refdata-gex-GRCm39-2024-A/genes/genes.gtf.gz` | `dfc17ee4c85f3d7665f48fe758f9a1109c3cc4c38a4185943e9e9678fae45fe6` |
| GRCm39 `reference.json` | `/nfs/turbo/agc-data/refs/Mus_musculus/GRCm39/refdata-gex-GRCm39-2024-A/reference.json` | `300f0ad96009d483c66c09841d3816375b077cba96a68ca36bbb2b14b6c0c614` |

The canonical machine-readable results are
`results/pcdh19_probe_audit/references/pcdh19_probe_coordinates.tsv` and
`reference_manifest.tsv` under the Paper 3 Turbo output root. The coordinate
table SHA-256 is
`05c40f1f9f04fd9a6252e09c6efb868004f4b3f5171ec1d73b2050d4d5bdbe1c`.

## Independent Paper 3 steps

### Current workstream and step numbering

`00_source_discovery` is the completed read-only discovery step that found the
correct Ziobro allocation and the `15662-JZ` delivery. It is not the Pcdh19
probe audit. Step `01_sample_key` is now complete from the user-provided
mapping. The locked cross-sample probe work is the independent technical
substep `02a_pcdh19_probe_audit`; the construct-level X-GFP compatibility gate
is completed as Step `02b_xgfp_probe_compatibility_audit`. The broader Step
`02_input_audit` has not yet been run. The completed technical audits remain
unchanged; biological labels enter only through the separate sample-key join
in later work.

| Step | Status | Purpose |
| --- | --- | --- |
| `00_source_discovery` | Completed | Locate the correct Ziobro Turbo allocation and inventory the `15662-JZ` delivery read-only. |
| `01_sample_key` | Completed; experimental-unit details remain | Registered submitted name, target cells, organism, tissue/region, genotype, sex, and design group for all 12 samples. Donor/embryo/litter/batch structure still requires confirmation. |
| `02_input_audit` | Not started | Verify MD5s, choose one authoritative matrix location, audit features/barcodes/QC, and record exact inputs. |
| `02a_pcdh19_probe_audit` | Completed | Checksum-lock the v2.0.0/GRCm39-2024-A probe references and reproduce raw three-probe Pcdh19 counts and binary patterns for all 12 technical samples without biological labels. |
| `02b_xgfp_probe_compatibility_audit` | Completed | Validate the exact three custom Flex EGFP probes against the original Nagy/Kalantry D4/XEGFP construct-level reporter sequence before any GFP count interpretation. |
| `03_canonical_inputs` | Not started | Create and validate a minimal analysis-ready object without altering source files. |
| `04_qc_and_filtering` | Not designed | Define sample-aware cell/gene QC after the biological design and expected cell types are known. |
| `10_primary_analysis` | Not designed | Normalize, integrate only if justified, cluster, annotate, and test approved comparisons. |
| `20_validation` | Not designed | Evaluate robustness, replicate structure, PCDH19 biology, and reference mappings. |
| `30_final_figures` | Not designed | Create provenance-complete, versioned Paper 3 figure packages. |

## Cross-sample Pcdh19 probe-pattern snapshot

The completed `02a_pcdh19_probe_audit` gives the following technical
measurements. `Design group` is shown by a descriptive join to the separately
registered sample key; it was not an input to the frozen audit. `A` is exon-1
probe `a3f4e22`; `B` and `C` are downstream probes `8215225` and `d013e0b`.
`B+C / A-negative` is the observed `B+C` binary pattern, not a mutant-cell
label. The final column uses all cells with any downstream detection (`B or C`)
as its denominator.

| Technical sample | Design group | Filtered cells | Any Pcdh19, % | A+B+C, n | B+C / A-negative, n | Downstream A-negative, % |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `15662-JZ-1` | `WT_M` | 51,229 | 18.595 | 207 | 291 | 74.282 |
| `15662-JZ-2` | `WT_M` | 37,553 | 15.961 | 108 | 182 | 75.928 |
| `15662-JZ-3` | `WT_M` | 25,354 | 12.357 | 55 | 75 | 75.679 |
| `15662-JZ-4` | `WT_F` | 21,440 | 17.397 | 73 | 116 | 75.176 |
| `15662-JZ-5` | `WT_F` | 41,878 | 17.773 | 139 | 230 | 75.800 |
| `15662-JZ-6` | `WT_F` | 56,099 | 19.403 | 220 | 330 | 74.943 |
| `15662-JZ-7` | `HET_F` | 19,623 | 18.738 | 32 | 319 | 94.908 |
| `15662-JZ-8` | `HET_F` | 20,799 | 18.770 | 52 | 219 | 90.386 |
| `15662-JZ-9` | `HET_F` | 60,680 | 20.840 | 135 | 1,193 | 92.300 |
| `15662-JZ-10` | `KO_M` | 52,753 | 20.932 | 2 | 1,319 | 99.946 |
| `15662-JZ-11` | `KO_M` | 52,295 | 21.759 | 2 | 1,505 | 99.903 |
| `15662-JZ-12` | `KO_M` | 11,085 | 19.684 | 0 | 277 | 99.908 |

The notable technical result is a sharp shift in probe combination rather
than a loss of Pcdh19 detection overall. Samples 1--6 have 24--26% exon-1
positivity among downstream-positive cells; samples 7--9 have only 5--10%;
and samples 10--12 have approximately 0.1% or less. Triple-probe detections
are 55--220 cells in samples 1--6 but only 2, 2, and 0 cells in samples
10--12. Meanwhile, any-Pcdh19 detection remains 18.7--21.8% in samples
7--12, so there is no corresponding global disappearance of Pcdh19 signal.

The registered sample key shows that this technical separation follows the
four supplied design groups: both WT male and WT female samples have 74--76%
downstream-positive/A-negative cells (unweighted sample means 75.296% and
75.306%, respectively), HET female samples have 90--95% (mean 92.532%), and KO
male samples have approximately 99.9% (mean 99.919%). The appropriate future
genotype contrasts are HET female versus WT female and KO male versus WT male. This is
a strong descriptive correspondence, not yet a statistical result, and it
does not make an A-negative downstream-positive barcode a mutant cell. Probe
non-detection remains subject to sampling and assay efficiency.

## Next action

Preserve [`config/sample_key.csv`](config/sample_key.csv) as the biological
annotation layer and confirm the experimental-unit structure behind the
submitted names: donor/embryo/litter, batch, age, and whether the three samples
per group are independent biological replicates. Then complete the broader
Step `02_input_audit` before canonical-object construction and formal
sample-level comparisons. Step `00_source_discovery` remains closed unless the
delivered source changes; do not rename the completed probe audit as Step 00.
The X-GFP sequence gate is complete, but GFP UMI extraction remains a separate
future action and must not be represented as part of Step `02b`.

The independent technical Pcdh19 probe audit does not ingest the sample key and
preserves only `15662-JZ-1` through `15662-JZ-12`. Its single local
entry point is `paper3_pcdh19/bin/run_pcdh19_probe_audit_all.sh`; the matching
batch entry is `paper3_pcdh19/slurm/pcdh19_probe_audit_all.sbatch`. The runner
must reproduce the frozen JZ-1 barcode table SHA-256 before it can advance to
samples 2 through 12.

The locked runner completed all 12 technical samples on 2026-08-27. All
per-sample and combined validation rows are `PASS`; the JZ-1 barcode table is
byte-equivalent to the validated prototype; and the final output manifest
covers every published reference, per-sample, combined, and environment file.
The completed result root is:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/pcdh19_probe_audit
```

The production package records Python implementation SHA-256
`731edca099bcbca8fb5e968ed621c2e5d8e31e3f25b6209b5008d2eb082d7b43`.
The computational file remains frozen at that checksum; expanded explanatory
documentation is maintained separately so the existing provenance remains
truthful.

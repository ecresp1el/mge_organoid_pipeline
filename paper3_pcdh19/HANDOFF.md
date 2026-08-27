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
clustering results, but they do not replace the missing biological sample key.

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
cell-count range, especially sample 12, must be evaluated after the biological
design is known.

## Missing biological metadata

The inspected AGC configuration and QC files map technical sample IDs to probe
barcodes but leave `Sample description` blank. They do not establish:

- PCDH19 genotype or mosaic fraction;
- control versus experimental group;
- cell line, clone, donor, sex, or reporter;
- dorsal/ventral or other regional identity;
- organoid age/time point;
- differentiation batch;
- biological replicate structure; or
- intended primary comparisons.

Do not infer these fields from sample order. Obtain the experiment submission
sheet, lab sample key, or confirmation from Julia Ziobro before building the
analysis object.

## Authoritative configuration

- [`config/input_candidates.tsv`](config/input_candidates.tsv) records the
  correct source root and technical format.
- [`config/sample_manifest_draft.tsv`](config/sample_manifest_draft.tsv)
  contains only verified technical IDs, probe barcodes, and vendor QC metrics.
- [`config/greatlakes.env`](config/greatlakes.env) keeps the source on the
  Ziobro allocation and points Paper 3 outputs to the existing `umms-parent`
  MGE project.

## Independent Paper 3 steps

| Step | Status | Purpose |
| --- | --- | --- |
| `00_source_discovery` | Completed | Locate the correct Ziobro Turbo allocation and inventory the `15662-JZ` delivery read-only. |
| `01_sample_key` | Waiting for metadata | Register biological identities, experimental units, and intended comparisons for all 12 samples. |
| `02_input_audit` | Not started | Verify MD5s, choose one authoritative matrix location, audit features/barcodes/QC, and record exact inputs. |
| `02a_pcdh19_probe_audit` | Completed | Checksum-lock the v2.0.0/GRCm39-2024-A probe references and reproduce raw three-probe Pcdh19 counts and binary patterns for all 12 technical samples without biological labels. |
| `03_canonical_inputs` | Not started | Create and validate a minimal analysis-ready object without altering source files. |
| `04_qc_and_filtering` | Not designed | Define sample-aware cell/gene QC after the biological design and expected cell types are known. |
| `10_primary_analysis` | Not designed | Normalize, integrate only if justified, cluster, annotate, and test approved comparisons. |
| `20_validation` | Not designed | Evaluate robustness, replicate structure, PCDH19 biology, and reference mappings. |
| `30_final_figures` | Not designed | Create provenance-complete, versioned Paper 3 figure packages. |

## Next action

Find the biological key for samples `15662-JZ-1` through `15662-JZ-12`.
Likely places to check are the AGC submission records, a lab spreadsheet, email
handoff, or another directory under `Ziobro Lab`. Until that key is recovered,
the correct next task is metadata discovery rather than Seurat processing.

The independent technical Pcdh19 probe audit may be run without the sample key
because it preserves only `15662-JZ-1` through `15662-JZ-12`. Its single local
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

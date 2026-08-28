# X-GFP/Flex probe compatibility audit

Status: completed and validated on 2026-08-27 as Paper 3 Step
`02b_xgfp_probe_compatibility_audit`.

## Question

Was the exact custom Flex panel delivered for `15662-JZ` sequence-compatible
with the Nagy/Kalantry X-GFP reporter used by Julie?

## Answer

Yes at the construct-reference level. The delivered panel contains exactly
three included EGFP probes. All three are unique, exact 50-of-50-base
reverse-complement matches to the Clontech EGFP coding sequence specified for
the original pCX-EGFP construct used to generate D4/XEGFP. There are zero
mismatches and zero gaps.

Therefore the delivered Flex panel was capable, by target sequence, of
detecting transcripts from this X-GFP reporter. This establishes theoretical
sequence compatibility. It does not establish that the reporter was
transcribed in a particular sample or cell, that its RNA was preserved, that
all probes hybridized with equal efficiency, or that any UMI was observed.

## Reporter identity and evidence chain

1. Kalantry laboratory supplementary methods identify their X-GFP as the
   previously described D4/XEGFP line.
2. MGI identifies D4/XEGFP as `Tg(CAG-EGFP)D4Nagy`, MGI `3055027`; JAX lists
   the same line as stock `003116`.
3. The original line reports trace D4/XEGFP to R1 ES cells transfected with
   pCX-EGFP.
4. The original construct report describes pCX-EGFP as Clontech enhanced GFP
   cDNA driven by the CMV immediate-early enhancer/chicken beta-actin promoter
   with rabbit beta-globin polyadenylation sequence.
5. NCBI GenBank `U55762.1` is the 1996 Clontech-submitted pEGFP-N1 record. Its
   EGFP CDS is bases 679--1398, 720 bases including the stop codon. The CDS
   SHA-256 used here is
   `e8d734949c2b74ac68d9ecd5aa1016aa75d496d292b124db7d4993f2e5409449`.

Primary and authoritative records:

- [MGI `Tg(CAG-EGFP)D4Nagy` line record](https://www.informatics.jax.org/allele/MGI:3055027)
- [JAX D4/XEGFP stock 003116](https://www.jax.org/strain/003116)
- [Original X-linked line report, Hadjantonakis et al. 1998](https://doi.org/10.1038/893)
- [Original pCX-EGFP/D4 derivation report, Hadjantonakis et al. 1998](https://doi.org/10.1016/S0925-4773(98)00093-8)
- [pCX-EGFP construct report, Okabe et al. 1997](https://doi.org/10.1016/S0014-5793(97)00313-X)
- [Kalantry laboratory X-GFP supplementary methods](https://websites.umich.edu/~kalantry/Kalantry_Lab/Publications_files/Supplement.pdf)
- [NCBI GenBank `U55762.1`](https://www.ncbi.nlm.nih.gov/nuccore/U55762.1)

## Exact delivered Flex probes and alignments

The accessible delivered panel is:

```text
/nfs/turbo/umms-ziobroj/Ziobro Lab/MGE scRNA-seq/10x_analysis_15662-JZ/Sample_15662-JZ-P01/probe_set.csv
SHA-256 70368b977ffd3d7ce183ebf82581d4bf0357c3c4c4bb18d8a99b013f9c77b8f1
```

Cell Ranger `config.csv` records the AGC custom source name as
`Chromium_Mouse_Transcriptome_Probe_Set_v2.0.0_GRCm39-2024-A_EGFP.csv`.
The current workspace cannot traverse that AGC custom-reference directory, so
the audit uses the exact delivered `probe_set.csv` emitted with the run. The
delivered header is Chromium Mouse Transcriptome Probe Set v2.0.0, GRCm39,
reference 2024-A.

| EGFP 5-prime order | Full probe ID | Probe sequence, 5-prime to 3-prime | EGFP CDS interval, 1-based inclusive | Orientation | Identity |
| ---: | --- | --- | ---: | --- | ---: |
| 1 | `EGFP|EGFP|probe03` | `ATGGTGCGCTCCTGGACGTAGCCTTCGGGCATGGCGGACTTGAAGAAGTC` | 247--296 | Reverse complement | 50/50; 100% |
| 2 | `EGFP|EGFP|probe02` | `AGGGTGTCGCCCTCGAACTTCACCTCGGCGCGGGTCTTGTAGTTGCCGTC` | 310--359 | Reverse complement | 50/50; 100% |
| 3 | `EGFP|EGFP|probe01` | `GGTAGTGGTCGGCGAGCTGCACGCTGCCGTCCTCGATGTTGTGGCGGATC` | 501--550 | Reverse complement | 50/50; 100% |

All three delivered rows have `included=TRUE` and the 10x annotation
`region=unspliced`. The alignment tests both orientations and requires one
unique exact match in the EGFP CDS. Each probe matches only as the reverse
complement, as expected for the reported probe sequence orientation.

## Important sequence limitation

The cited MGI line record reports the X-chromosomal insertion but lists its
precise location as unknown. The cited line records do not deposit the full
integrated concatemer and mouse/transgene junction sequence. Consequently,
this audit makes an exact construct-level statement about the EGFP reporter
CDS specified by the original pCX-EGFP references; it does not claim to have
resequenced every integrated D4/XEGFP copy in Julie's colony.

If colony-specific sequencing later reveals a mutation in any of the three
target intervals, the alignment must be rerun against that sequence. In the
absence of such evidence, all three probes exactly match the authoritative
construct-level EGFP reference.

## Reproducible implementation and outputs

Run:

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
./paper3_pcdh19/bin/run_xgfp_probe_audit.sh
```

Version-controlled inputs:

- `config/xgfp_probe_audit.lock.json`
- `scripts/xgfp_probe_compatibility_audit.py`
- `bin/run_xgfp_probe_audit.sh`

Validated production output:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/xgfp_probe_audit/
  gfp_probe_alignment.tsv
  conclusion.tsv
  validation.tsv
  software_environment.tsv
  output_manifest.tsv
  references/
    U55762.1_EGFP_CDS.fasta
    source_manifest.tsv
```

The runner checksum-validates the delivered panel and NCBI reference,
requires the exact three-probe ID/sequence set, validates all alignments before
publication, writes through a temporary directory, refuses different existing
outputs, and checksum-verifies an existing completed audit on restart.

## Stop point

This step does not inspect GFP UMI counts, sample-level GFP detection,
barcodes, clusters, or UMAP. Those are separate downstream questions. The only
conclusion here is that sequence incompatibility does not explain failure to
detect the Nagy/Kalantry reporter: the delivered Flex panel contains three
perfectly matching EGFP probes.

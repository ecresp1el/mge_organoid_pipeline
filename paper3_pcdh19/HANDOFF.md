# Paper 3 Ziobro PCDH19 MGE single-cell RNA-seq: operational handoff

Last updated: 2026-08-29

> **Classification freeze:** Steps 03–07 are complete and immutable. The
> authoritative resume document is
> [`PCDH19_GENOTYPE_CLASSIFICATION_HANDOFF.md`](PCDH19_GENOTYPE_CLASSIFICATION_HANDOFF.md).
> Do not change any existing model, threshold, classification, or result
> package. Post-freeze analyses must remain separate and downstream.

## Independent developing-mouse MGE reference curation (first checkpoint complete)

The authoritative instructions for auditing La Manno 2021, Bandler 2022, and
Mayer 2018 as candidate WT developing-mouse MGE references are now registered
in
[`PCDH19_DEVELOPING_MOUSE_MGE_REFERENCE_CURATION_HANDOFF.md`](PCDH19_DEVELOPING_MOUSE_MGE_REFERENCE_CURATION_HANDOFF.md).
This is an isolated reference-curation and verification module. Its first
checkpoint is limited to source/metadata audit and inspection of the smallest
P0 processed objects. It must not download raw FASTQs, integrate or map the
Paper 3 cells, invent annotations, recluster giant datasets, or select a final
reference before the evidence audit is complete.

The first checkpoint now creates a curation-ready sample and library ledger,
not just an object-structure row. It preserves published/author sample IDs,
sample counts, ages, tissues/regions, genotype/condition, sex/strain,
pool/replicate/QC, single-cell modality, capture technology, library and
instrument fields, chemistry/protocol, experimental selection, raw/BioSample
accessions, and P0-object membership. Published metadata and proven P0 cell
membership remain separate claims. The combined outputs are
`tables/all_candidate_reference_samples.tsv`, `tables/study_sample_summary.tsv`,
`tables/sample_metadata_data_dictionary.tsv`, and
`tables/reference_curation_requirements_ledger.tsv`; unresolved or future-stage
requirements are marked explicitly.
After the three inspections pass, the checkpoint job also generates
`REFERENCE_CURATION_REPORT.md`, `tables/reference_annotation_availability.tsv`,
and per-study observed-label `annotation_dictionary.tsv` files from the exact
report code copied into that run.

This auxiliary implementation follows the same modular object-oriented rule
as the classification code: registry, immutable cache, metadata collection,
study-specific object inspection, per-study orchestration, readiness
evaluation, and checkpoint publication are separate objects.

This new three-candidate workstream is a sibling of the completed GSE94641
Step 00 mapping work, not a replacement for it and not an input to the frozen
Steps 03–07 classifier. Its first checkpoint completed on 2026-08-29 in run
`00_developing_mouse_mge_reference_curation_20260829_141944_3b2ad52`.
La Manno was directly annotation/embedding-ready; the Bandler and Mayer P0
files were counts-only and require published barcode metadata before their
cell types can be displayed. See that run's `REFERENCE_CURATION_REPORT.md`.
For Bandler, a historical author link has now yielded the separate 65,700-cell
postnatal STICR Seurat object with author labels and UMAP. It contains no
CA301/CA298--303 samples and therefore does not annotate the exact WT E15.5 MGE
matrix. The run-scoped Bandler recovery report preserves that negative join,
the official seven-sample embryonic design, all 21 published embryonic labels,
and the later 2025 interactive-atlas evidence without conflating these artifacts.
The live atlas has now been queried through its intended public Shiny session
and vector-plot download handlers. It contains 18,424 Bandler cells in the
later filtered embryonic atlas: 11,004 E13 and 7,420 E15, with 12 later-atlas
clusters, 2,877 mitotic cells, and 15,547 inhibitory neuron precursors. A
24-gene expression-fingerprint and preserved-order follow-up resolves every
E15 atlas row to a deposited barcode: 4,481 CA301 MGE, 2,937 CA302 CGE, and 2
CA303 LGE. All 7,420 assignments are definitive under the order-constrained
validation. These are later MIND-atlas labels, not the original Bandler 2022
21-cluster per-cell assignments.
The E15/MGE evidence is not equivalent across candidates: Bandler directly
provides 4,516 WT E15.5 MGE cells (GEO shorthand E15); Mayer directly provides 6,515 MGE E13.5
Lhx6-positive cells; and La Manno has exact E15 and ventral-forebrain cells but
no inspected author label that proves exact E15 MGE membership. The run now
also contains an author-label hierarchy/composition plot and a three-study
E15/MGE evidence matrix, with underlying TSVs and exact executed plot code.

The first-checkpoint scaffold is implemented through
`bin/submit_developing_mouse_mge_reference_curation.sh`. Default submission
creates a versioned package beneath
`results/00_developing_mouse_mge_reference_curation/`. Before submission it
copies the exact Python, R, submitter, configuration, and SLURM files into the
run package; the jobs execute those copies. P0 processed objects are cached
once beneath `inputs/developing_mouse_mge/` and are not duplicated per run.
The dependency-gated checkpoint invokes the object-oriented hierarchy plotter
before rebuilding the Markdown report, so a rerun or guarded same-step
replacement regenerates the visual outputs rather than relying on a manual
figure.
The standalone Bandler recovery job likewise executes only frozen run copies
of the downloader, Seurat inspector, supplement publisher, report builder, and
SLURM script. Rerunning it intentionally overwrites derived tables, figures,
and reports inside the same step while reusing validated source-cache
artifacts; the copied executable code remains the primary record of what ran.
The standalone `bin/submit_mind_public_atlas_capture.sh RUN_DIRECTORY`
follow-up follows the same rule: it freezes the exact Python/report/SLURM code
under the run's `code/`, captures the public PDFs from the submission host,
and has SLURM validate/parse those exact files. It intentionally overwrites
only `Bandler2022/interactive_atlas/` and refreshes the main report when
rerun. Job 59175410 completed this package; failed compute-node network
attempts 59175246 and 59175304 remain visible in provenance.
The standalone `bin/submit_bandler_e15_barcode_recovery.sh RUN_DIRECTORY`
follow-up applies the same frozen-code rule, overwrites only
`Bandler2022/interactive_atlas/barcode_recovery/`, reuses immutable deposited
matrices, captures intended public expression PDFs on the submission host,
and performs the fingerprint/order join under SLURM.

Use `--replace-run RUN_ID` only when intentionally regenerating an existing inactive run
within this exact step. The guarded replacement rejects unsafe names and active
jobs, clears only the selected run package, and preserves the shared input
cache. This replace mode does not authorize overwriting the frozen Steps 03–07
packages.

## Current next analysis boundary

Step `03_pcdh19_genotype_classification_setup` is implemented and complete.
It checksum-validates the frozen Step 02a per-cell probe tables, joins the
registered biological sample key, and publishes a 230,269-row WT-male/KO-male
ground-truth candidate table beneath
`results/step_03_pcdh19_genotype_classification_setup/`. The table retains raw
A/B/C counts and patterns and adds direct detection indicators and explicit
WT=`0`/KO=`1` targets. It fits no classifier and does not score HET cells.

All new Python analysis steps from Step 03 forward should retain the modular,
object-oriented pattern established by
`scripts/Step_03_PCDH19_Genotype_Classification_Setup.py`: keep loading,
encoding, splitting, classification, evaluation, and HET application in small
separable components rather than extending one monolithic script.

Step `04_pcdh19_empirical_pattern_classifier` is also implemented and
complete. It consumes only the manifested Step 03 table and publishes the
eight-state empirical WT/KO probability/likelihood-ratio model plus three
diagnostic plots beneath
`results/step_04_pcdh19_empirical_pattern_classifier/`. Every pattern,
including `000`, remains uncalled; Step 04 does not load HET cells or report
held-out performance. The current architecture and future-only module plan are
registered in `PCDH19_CLASSIFICATION_FRAMEWORK.md`.

Step `05_pcdh19_logistic_regression_baseline` is implemented and complete. It
reuses the locked Step 04 encoder, ground-truth loader, empirical estimator,
and probability interface, then adds an unpenalized
`A_detected+B_detected+C_detected` logistic model with WT=`0`, KO=`1`. The
`000` intercept gives P(WT)=0.5099 and P(KO)=0.4901, and remains uncalled. Step
05 reports coefficients, odds ratios, all eight probabilities, Step 04
comparisons, and plots beneath
`results/step_05_pcdh19_logistic_regression_baseline/`.

Step 05 also now owns a leave-one-registered-sample-out validation subpackage
at `sample_level_held_out_validation/`. Each of the six WT-male/KO-male samples
is held out once; the same three-feature model is fit on the other five. Pattern
`000` is always uncalled, and informative patterns use a fixed, untuned 0.5
probability rule. Overall, 43,256/230,269 cells (18.785%) are called; called-
cell accuracy is 80.729%, KO sensitivity is 99.711%, and KO specificity is
55.691%. This asymmetry and the low call percentage must remain explicit in
downstream interpretation. The registered `technical_sample_id` is the
available holdout key; donor/embryo/litter independence is not established.

The male-only package is preserved byte-for-byte. A second Step 05 subpackage,
`wt_male_female_vs_ko_male_validation/`, defines WT as JZ-1--6 WT males plus
WT females, KO as JZ-10--12 KO males, and completely excludes JZ-7--9 HET
females. It holds out each of nine samples once and applies the same unweighted
A/B/C model and fixed rule. Of 349,686 cells, 65,314 (18.678%) are called;
called-cell accuracy is 41.960%, KO sensitivity 20.302%, and WT specificity
55.049%.

WT-female inclusion does not change WT-male specificity (55.691%) or its 8,265
false-KO cells. Those errors are exactly `001`=4,586, `010`=3,131, and
`011`=548; all A-detected states contribute zero. WT-female specificity is
54.506%. Holding out either large KO sample causes `001` and `010` to flip to
WT because the unchanged model is unweighted and the remaining training cohort
is WT-heavy; this explains the aggregate KO-sensitivity collapse. Preserve
this diagnostic result. Do not silently add weighting or optimize the threshold
inside the baseline.

Step 05 still performs no final model selection, threshold optimization, or
HET inference. Step 06 now owns inference-only application of both established
full-fit models; those responsibilities must not be moved back into Step 05.

Step 05 now also owns the separately manifested descriptive module
`a_negative_raw_bc_umi_evidence/`. It restricts the registered WT-M+F and KO-M
controls to exact `A_UMI == 0`, preserves raw integer Probe B and C UMI counts,
and fits no classifier. The output contains 327,204 cells: 103,748 WT-M,
107,394 WT-F, and 116,062 KO-M; no HET-female table is opened.

Increasing B/C evidence is monotonically KO-enriched. For raw B+C total,
P(KO) is 54.255% with one UMI, 67.640% with two, and 80.216% with three or
more; KO:WT likelihood ratios are 2.158, 3.803, and 7.376. B alone and C alone
show the same trend. The 3+ bins are small (158 B, 499 C, and 1,663 B+C
cells), and WT cells remain present, so this is a descriptive result rather
than a count-based genotype rule or classifier specificity. The ten-file
package manifest is
`0a27a9e5c1e983c5b1a6772435b299f72a974bac405d9a39932e0c8146dc0631`.

Step 05 additionally owns `count_informed_vs_binary_validation/`. This package
compares the immutable expanded binary `A_detected+B_detected+C_detected`
model against `A_detected+B_UMI+C_UMI` on the exact same 349,686 cells and
nine registered-sample holdouts. B/C counts are raw probe-level UMI/ligation
evidence, not transcript numbers. Both models leave exact `000` uncalled and
use the same fixed, untuned 0.5 rule. No normalization, weighting, interaction,
nonlinearity, HET inference, or threshold change is present.

The count-informed model raises called-cell accuracy from 41.960% to 44.695%
and KO sensitivity from 20.302% to 27.570%; WT specificity is essentially
unchanged (55.049% to 55.044%). It corrects 2,373 binary errors and introduces
587 regressions. This gain is not sample-consistent: JZ-10 and JZ-11 improve by
10.29 and 10.87 points, JZ-12 declines by 26.81 points, and pooled held-out AUC
drops from 0.254 to 0.201. Treat this as modest, sample-dependent aggregate
improvement; equally weighted mean sample accuracy decreases from 50.613% to
49.983%. This is not proof that the count model generalizes robustly or should
be selected for HET inference. The 15-file package manifest is
`cf68ee398a2d3a9e19e3de05a59be50d9aed2bfe1084b3445f63ab47e556130f`.

The read-only `sample_level_probe_evidence_diagnostics/` package reconciles all
349,686 WT-M/WT-F/KO-M cells to manifested raw probe rows; no eligible cell is
dropped and no HET table is opened. JZ-12 is modestly weaker than JZ-10/JZ-11
in A-negative B+C evidence (mean 0.2536 versus 0.2700/0.2867; zero fraction
80.40% versus 79.11%/78.29%), but remains biologically KO-like. The 2,373
count-model corrections are B+C=2 (2,070) or 3+ (303); 585/587 regressions are
B+C=1 and occur in JZ-12. This identifies a fold-specific one-UMI decision
effect rather than a categorically anomalous JZ-12 probe distribution.

The AUC audit confirms KO=`1`, P(KO) scoring, and consistent orientation.
Reversing to P(WT) yields AUC above 0.5 but is only a diagnostic; do not replace
the immutable stored AUCs. The 29-file diagnostic package passed 46 checks and
has manifest SHA-256
`ad4bb28c6fce5adea9b66916378e43c2407ef6d74d67d47d65d2f10f77ecf7fe`.

Step `06_pcdh19_het_female_inference` is implemented and complete. It verifies
the exact Step 05 expanded binary and count-informed full-fit models, loads
only JZ-7--9 HET-female probe rows for inference, and applies both models
unchanged. The 101,102-row per-cell output preserves sample, barcode, raw
A/B/C probe-level UMI/ligation counts, detection states, pattern, and both
models' WT/KO probabilities. HET rows never enter fitting, coefficient
estimation, weighting, model selection, or threshold optimization. No cell
type is loaded and no HET genotype call is created.

Exact `000` accounts for 81.26%, 81.23%, and 79.16% of JZ-7--9 and is labeled
`uncalled_000`. A+ fractions are 2.40%, 5.16%, and 3.99%. Among A-negative
cells, B+C=2/3+ totals 3.45%, 2.35%, and 3.68%, intermediate between WT-F
(1.23%) and KO-M (4.90%). Non-`000` mean count-informed P(KO) values are
0.492, 0.403, and 0.458, versus 0.250 for WT-F and 0.570 for KO-M. Each HET
sample shows low-P(KO) A-detected and higher-P(KO) A-negative/B/C evidence.
That is consistent with mixed WT-like and KO-like probe evidence, but control
overlap and Step 05's sample-level limitations prevent interpreting these
modes as validated cellular genotypes.

The Step 06 package is at
`results/step_06_pcdh19_het_female_inference/`; all checks and immediate
idempotent verification pass. Its output-manifest SHA-256 is
`8d5bf77448fee3fc67991094960fa355a826def67de0e53777d712a8079d9a2a`.

Step `07_pcdh19_het_female_wt_ko_like_classification` is implemented and
complete. It selects thresholds before opening HET data, using only the exact
Step 05 count-informed leave-one-sample-out WT-M/WT-F/KO-M probabilities. The
prespecified target is at least 95% held-out precision with maximum coverage;
if that is unavailable, the fallback is maximum precision then maximum
coverage.

The frozen rule is WT-like at P(KO)<=0.301037619832, KO-like at
P(KO)>=0.911554020713, uncertain between, and `Uncalled_000` for exact `000`.
WT-like achieves 99.684% precision and 55.042% informative-WT sensitivity.
The KO-like 95% target is not achievable: the retained best KO-enriched tail
has 76.316% precision, 0.471% informative-KO sensitivity, and 152 held-out
calls. Keep this KO-like category, but always state that its validation
confidence is lower than WT-like.

The resulting HET counts are JZ-7: 471 WT-like, 26 KO-like, 3,180 uncertain,
15,946 `000`; JZ-8: 1,074, 10, 2,820, 16,895; and JZ-9: 2,419, 104, 10,123,
48,034. Pooled informative proportions are 19.598% WT-like, 0.692% KO-like,
and 79.710% uncertain. These are inferred PCDH19 probe-evidence states, not
independently observed DNA genotypes.

HET rows influenced no model coefficient, feature, calibration, threshold, or
performance estimate, and no cell type was loaded. The complete package is at
`results/step_07_pcdh19_het_female_wt_ko_like_classification/`; its manifest
SHA-256 is
`aec9a9cf7c8575ca453fbc61172fdd6d5c28e6c9be177d662a3c7a75cb40af9d`.

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
is completed as Step `02b_xgfp_probe_compatibility_audit`; and the raw EGFP
count audit is completed as Step `02c_egfp_probe_audit`. The broader Step
`02_input_audit` has not yet been run. Step `02c` uses the registered sample
key only to create a clearly labeled descriptive sample-design table; it does
not use those fields to alter counts or classify cells.

| Step | Status | Purpose |
| --- | --- | --- |
| `00_source_discovery` | Completed | Locate the correct Ziobro Turbo allocation and inventory the `15662-JZ` delivery read-only. |
| `01_sample_key` | Completed; experimental-unit details remain | Registered submitted name, target cells, organism, tissue/region, genotype, sex, and design group for all 12 samples. Donor/embryo/litter/batch structure still requires confirmation. |
| `02_input_audit` | Not started | Verify MD5s, choose one authoritative matrix location, audit features/barcodes/QC, and record exact inputs. |
| `02a_pcdh19_probe_audit` | Completed | Checksum-lock the v2.0.0/GRCm39-2024-A probe references and reproduce raw three-probe Pcdh19 counts and binary patterns for all 12 technical samples without biological labels. |
| `02b_xgfp_probe_compatibility_audit` | Completed | Validate the exact three custom Flex EGFP probes against the original Nagy/Kalantry D4/XEGFP construct-level reporter sequence before any GFP count interpretation. |
| `02c_egfp_probe_audit` | Completed | Extract raw UMI counts for the same three exact EGFP probes from every vendor-filtered barcode, reproduce Cell Ranger's EGFP row, and report all eight detection patterns across all 12 samples. |
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

## Cross-sample EGFP probe-count result

Step `02c_egfp_probe_audit` used the exact Step `02b`-validated probes:
`E01=EGFP|EGFP|probe01`, `E02=EGFP|EGFP|probe02`, and
`E03=EGFP|EGFP|probe03`. Each per-barcode count is a raw integer UMI extracted
from `sample_raw_probe_bc_matrix.h5`. The denominator is Cell Ranger's existing
vendor-filtered barcode set; the pipeline did not make new cell calls. For
every barcode in all 12 samples, `E01 + E02 + E03` exactly equals the filtered
Cell Ranger EGFP feature count.

| Technical sample | Design group | Filtered cells | EGFP UMIs | EGFP-positive barcodes | Detected, % | E01 | E02 | E03 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `15662-JZ-1` | `WT_M` | 51,229 | 3 | 3 | 0.005856 | 3 | 0 | 0 |
| `15662-JZ-2` | `WT_M` | 37,553 | 0 | 0 | 0 | 0 | 0 | 0 |
| `15662-JZ-3` | `WT_M` | 25,354 | 0 | 0 | 0 | 0 | 0 | 0 |
| `15662-JZ-4` | `WT_F` | 21,440 | 2 | 2 | 0.009328 | 2 | 0 | 0 |
| `15662-JZ-5` | `WT_F` | 41,878 | 2 | 2 | 0.004776 | 1 | 0 | 1 |
| `15662-JZ-6` | `WT_F` | 56,099 | 0 | 0 | 0 | 0 | 0 | 0 |
| `15662-JZ-7` | `HET_F` | 19,623 | 2 | 2 | 0.010192 | 2 | 0 | 0 |
| `15662-JZ-8` | `HET_F` | 20,799 | 0 | 0 | 0 | 0 | 0 | 0 |
| `15662-JZ-9` | `HET_F` | 60,680 | 3 | 3 | 0.004944 | 2 | 1 | 0 |
| `15662-JZ-10` | `KO_M` | 52,753 | 0 | 0 | 0 | 0 | 0 | 0 |
| `15662-JZ-11` | `KO_M` | 52,295 | 0 | 0 | 0 | 0 | 0 | 0 |
| `15662-JZ-12` | `KO_M` | 11,085 | 0 | 0 | 0 | 0 | 0 | 0 |
| **Total** |  | **450,788** | **12** | **12** | **0.002662** | **10** | **1** | **1** |

All 12 detected barcodes have exactly one EGFP UMI. No barcode has two or
three EGFP probes detected: the only nonzero patterns are ten `E01 only`, one
`E02 only`, and one `E03 only`. Thus the apparent lack of GFP is not a
gene-visualization artifact hiding a moderate joint-probe signal. The audit
finds an extremely sparse, non-concordant trace signal. Five of the 12 events
occur in HET-female samples (two in JZ-7 and three in JZ-9; JZ-8 has zero), but
that is not enough evidence to use EGFP as a reporter-based cell classifier or
to test concordance with Pcdh19 probe patterns. No EGFP/Pcdh19 cell label was
created.

The production result is
`/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/egfp_probe_audit`.
It contains 221 validation rows, all `PASS`. The output-manifest SHA-256 is
`ee81480aed00d8f9c86bba1044cc8d2ae66bd5a95a12deb269151ae2f9849f1d`.
The JZ-1 barcode table was recreated from source and matched the independently
frozen SHA-256
`a6e773d19d6ea7f796848dad770c72f2f540e2205413f55f7a2dc25c0722fa4c`.
An immediate full rerun verified every existing checksum and published no
replacement files.

## Next action

Preserve [`config/sample_key.csv`](config/sample_key.csv) as the biological
annotation layer and confirm the experimental-unit structure behind the
submitted names: donor/embryo/litter, batch, age, and whether the three samples
per group are independent biological replicates. Then complete the broader
Step `02_input_audit` before canonical-object construction and formal
sample-level comparisons. Step `00_source_discovery` remains closed unless the
delivered source changes; do not rename the completed probe audit as Step 00.
The X-GFP sequence gate and the separate EGFP raw-count audit are both
complete. Do not merge Step `02c` into Step `02b`: Step `02b` establishes
sequence compatibility, whereas Step `02c` measures observed UMIs. Because
Step `02c` found only 12 isolated single-probe events, defer EGFP/Pcdh19
cell-level concordance and any reporter-based cell recovery unless an
independent reason and explicit low-count interpretation plan are established.

The independent technical Pcdh19 probe audit does not ingest the sample key and
preserves only `15662-JZ-1` through `15662-JZ-12`. Its single local
entry point is `paper3_pcdh19/bin/run_pcdh19_probe_audit_all.sh`; the matching
batch entry is `paper3_pcdh19/slurm/pcdh19_probe_audit_all.sbatch`. The runner
must reproduce the frozen JZ-1 barcode table SHA-256 before it can advance to
samples 2 through 12.

The analogous EGFP entry point is
`paper3_pcdh19/bin/run_egfp_probe_audit_all.sh`; its scheduler wrapper is
`paper3_pcdh19/slurm/egfp_probe_audit_all.sbatch`. It first checksum-validates
the delivered panel and completed Step `02b`, then requires the independent
JZ-1 prototype before advancing. Its descriptive design summary is the only
Step `02c` output that joins `config/sample_key.csv`; all per-barcode tables
retain only technical sample IDs and raw counts.

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
That hash identifies the source at repository commit `22a342e`. The current
working source contains expanded docstrings and therefore has different bytes;
it must not be described as having generated the historical package. The
verified distinction between analyzed inputs, generated assets, run logs, and
current documented source is maintained in
[`ASSET_AND_RUN_INVENTORY.md`](ASSET_AND_RUN_INVENTORY.md).

# Paper 3: Ziobro PCDH19 MGE single-cell RNA-seq

This workstream is reserved for the Ziobro PCDH19 MGE single-cell paper. It is
separate from the Paper 2 cross-study atlas and uses the Ziobro Turbo
allocation, independent step numbering, and independent outputs.

## Critical current query-data state: no joint UMAP yet

The 12 Paper 3 samples have **not** yet undergone a unified Scanpy or Seurat
processing/integration workflow. The source delivery provides a filtered
feature-barcode matrix, graph-cluster assignments, and UMAP coordinates from
Cell Ranger for each sample. Each sample's graph and UMAP were computed
independently; cluster numbers and distances are therefore not comparable
between samples.

Step 01 subsequently created a derived Seurat mapping container with the
original counts, metadata, per-sample Cell Ranger coordinates, Bandler/MIND
predictions, and MapMyCells predictions. Seurat normalization and reference
PCA projection were used internally for label transfer, but this object is
**not** a completed processed or integrated analysis object. No joint batch
correction/integration, neighbor graph, clustering, or UMAP has been run in
Scanpy or Seurat. All current query UMAP plots must remain faceted by sample;
placing the 12 coordinate sets in one panel would imply a geometry that does
not exist.

The **true primary bioinformatics workflow** is now a separate namespace with
its own Steps 00–07. Start with
[`PCDH19_PRIMARY_PROCESSING_HANDOFF.md`](PCDH19_PRIMARY_PROCESSING_HANDOFF.md).
Its expression inputs are the original 12 per-sample Cell Ranger filtered
feature-barcode matrices; the registered sample key and technical manifest are
joined only as metadata. Step 00 constructs the canonical raw-count AnnData
without QC filtering, normalization, reduction, clustering, or annotation.
Every computed step remains `IN_REVIEW` until explicitly approved, and only an
approved checkpoint may feed the next step. Existing Bandler/MIND,
MapMyCells, GSE94641, probe, and classification results remain preserved for
later downstream integration; they are not primary-processing inputs.

Primary-processing code is object-oriented Python under
`scripts/primary_processing/`. Submissions create versioned packages under
`PAPER3_ROOT/results/primary_processing/`, freeze exact executable and metadata
copies before submission, run those copies on Great Lakes, and keep scheduler
scripts and logs under `PAPER3_ROOT/jobs/` and
`PAPER3_ROOT/logs/primary_processing/`, respectively.

Step 00 computed run
`00_input_validation_and_canonical_anndata_20260830_113749_d8b6bf7` with
successful Great Lakes job `59279775`. The 6.60-GB canonical H5AD contains
450,788 cells × 19,071 genes and 1,295,361,777 sparse raw-count entries; all
221 validations passed. The user explicitly **APPROVED** this exact Step 00
run on 2026-08-30. Step 01 run
`01_qc_metrics_20260830_115715_2b57907` then completed in Great Lakes job
`59281063` with all 450,788 cells and 19,071 genes retained, 40/40 validation
checks passing, and no thresholds or filtering. Before approval, Step 01a was
authorized to visualize 3/4/5
scaled-MAD candidate boundaries independently within each technical sample,
without filtering or upper count/gene rules. Review-target Step 01a run
`01a_qc_mad_sensitivity_20260830_121931_d5936f9` completed in job `59282437`
with 71/71 validations passing and all 450,788 cells retained. Any-criterion
candidate percentages are 4.765%, 2.024%, and 0.985% at 3/4/5 MAD and are
dominated by the mitochondrial tail. The user explicitly **APPROVED** the
exact Step 01 and review-target Step 01a runs on 2026-08-30 and authorized
Step 02 to apply the per-sample 5-MAD low-count OR low-gene OR high-mt union.
Step 02 must preserve every exclusion reason, retain all genes/raw counts,
report sample/design before-after counts, and stop for review. Step 03 will use
scDblFinder, but cannot run until the proposed one-`GEX_1`-
capture definition is explicitly approved. The detailed results and execution history are in
[`PCDH19_PRIMARY_PROCESSING_HANDOFF.md`](PCDH19_PRIMARY_PROCESSING_HANDOFF.md).
The capture evidence is in
[`PCDH19_STEP03_SCDBLFINDER_CAPTURE_DECISION.md`](PCDH19_STEP03_SCDBLFINDER_CAPTURE_DECISION.md).

Step 02 run `02_qc_filtering_20260830_124611_97e1bb5` completed in Great Lakes
job `59287494` with 38/38 validation checks passing. It excluded 4,439 reviewed
5-MAD candidates and retained 446,349 cells × all 19,071 genes as sparse raw
integer counts. Exact reasons remain available for every original cell. The
checkpoint is **IN_REVIEW**; scDblFinder was not run.

## New auxiliary work: three-reference MGE curation

A standalone module audits three published WT developing-mouse
references—La Manno 2021, Bandler 2022, and Mayer 2018—before one is selected
for future mapping. Its first checkpoint completed on 2026-08-29 in run
`00_developing_mouse_mge_reference_curation_20260829_141944_3b2ad52`. Start with
[`PCDH19_DEVELOPING_MOUSE_MGE_REFERENCE_CURATION_HANDOFF.md`](PCDH19_DEVELOPING_MOUSE_MGE_REFERENCE_CURATION_HANDOFF.md).

The immediate scope is source discovery, minimal processed-object inspection,
author-annotation recovery, and author-embedding verification. Do not download
raw FASTQs, reconstruct large analyses, integrate references, or map the
PCDH19 cells before the required early checkpoint. This is a sibling of the
existing GSE94641 workstream and is not part of the frozen Steps 03–07
classification unit.

The observed P0 files show that La Manno has embedded `Class`, `Subclass`,
`CellType`, `cluster_id`, age/dissection metadata, and author UMAP/tSNE
coordinates. Bandler's exact WT E15.5 MGE object (GEO shorthand E15) and Mayer's six-sample 10x CSV
contain counts but no cell-level annotations or saved embeddings. The complete
checkpoint interpretation is in the run's `REFERENCE_CURATION_REPORT.md`.

A follow-up Bandler recovery found and inspected the authors' 65,700-cell
postnatal STICR Seurat object, including broad/refined annotations and UMAP.
It is not the missing embryonic object and contains no CA301 cells. Official
supplements preserve the complete 21-state embryonic vocabulary but no
barcode assignments. The later Mayer-lab interactive GE atlas uses an
unpublished local integrated RDS with Bandler embryonic cells; treat it as a
later reanalysis, not as a recovered original object. Its intended public
Shiny vector outputs can nevertheless be captured: they contain 18,424
Bandler study cells (11,004 E13; 7,420 E15) with later-atlas class and cluster
labels. A 24-gene public-expression fingerprint join now maps every E15 atlas
row back to a deposited barcode: 4,481 CA301 MGE, 2,937 CA302 CGE, and 2 CA303
LGE cells. These are later-atlas labels, not the original 2022 21-cluster
assignments. See
`Bandler2022/BANDLER_AUTHOR_OBJECT_RECOVERY_REPORT.md` in the curation run.

For the exact age/anatomy question: Bandler provides 4,516 directly identified
WT E15.5 MGE cells (`CA301`; GEO shorthand E15); Mayer provides 6,515 MGE cells at E13.5 with
Lhx6-positive enrichment; La Manno contains E15 cells and ventral-forebrain
dissections, but its inspected author labels do not prove which E15 cells are
MGE. The run includes a La Manno author hierarchy/composition plot and a
three-candidate E15/MGE evidence matrix in its `figures/` directories, with
the underlying counts/evidence retained as TSV tables and manifested outputs.

The early checkpoint additionally produces a standardized sample/library
inventory for the next curation stage. It separately reports published sample
counts and IDs versus samples demonstrably present in each P0 object, with
age, tissue/region, genotype/design, pool/replicate/QC, modality, technology,
library, instrument/platform, chemistry/protocol, enrichment/selection, raw
accessions, cell counts when proven, and explicit metadata gaps. See
`tables/all_candidate_reference_samples.tsv`, `tables/study_sample_summary.tsv`,
and `tables/reference_curation_requirements_ledger.tsv` in a completed run.
The Python design is object-oriented, with distinct registry/cache, metadata,
study-inspector, orchestration, evaluator, and publisher components.
Successful future checkpoint jobs also publish `REFERENCE_CURATION_REPORT.md`
and observed-label annotation dictionaries plus hierarchy/evidence plots
using the plotting and report code copied into the run package before
submission.

Prepare or submit the required first checkpoint with:

```bash
./paper3_pcdh19/bin/submit_developing_mouse_mge_reference_curation.sh --dry-run
./paper3_pcdh19/bin/submit_developing_mouse_mge_reference_curation.sh
```

For the rerunnable public-atlas follow-up inside a completed curation run:

```bash
./paper3_pcdh19/bin/submit_mind_public_atlas_capture.sh --dry-run /absolute/path/to/run
./paper3_pcdh19/bin/submit_mind_public_atlas_capture.sh /absolute/path/to/run
```

That wrapper freezes copies under `code/`, captures the public PDFs from the
submission host, then submits their validation/vector parsing/report rebuild
to SLURM. It overwrites only the derived `Bandler2022/interactive_atlas/`
package. Completed job: `59175410`.

For the rerunnable deposited-barcode follow-up:

```bash
./paper3_pcdh19/bin/submit_bandler_e15_barcode_recovery.sh --dry-run /absolute/path/to/run
./paper3_pcdh19/bin/submit_bandler_e15_barcode_recovery.sh /absolute/path/to/run
```

This freezes the exact Python/R/shell/SLURM code, reuses immutable GEO
matrices, captures intended public RNA-expression vectors on the submission
host, and overwrites only
`Bandler2022/interactive_atlas/barcode_recovery/`.

Every default submission creates a versioned package under
`PAPER3_ROOT/results/00_developing_mouse_mge_reference_curation/`. The exact
code and settings are copied into the package before submission and executed
from there. An existing inactive run can be intentionally regenerated with
`--replace-run RUN_ID`; replacement is constrained to this step, refused for
active jobs, and does not delete the shared P0 cache under
`PAPER3_ROOT/inputs/developing_mouse_mge/`.

## Current analysis: PRELIMINARY pseudobulk differential expression

Formal Step 09 is a rapid, explicitly preliminary sample-level pseudobulk
analysis. It tests WT-F/WT-M, HET-F/WT-F, and KO-M/WT-M across all cells,
progenitors, and immature neurons. Start with
[`PCDH19_PRELIMINARY_PSEUDOBULK_DE.md`](PCDH19_PRELIMINARY_PSEUDOBULK_DE.md).

No HET-F/WT-F or KO-M/WT-M test currently has an FDR≤0.05 gene. WT-F/WT-M
shows many mostly small effects in all cells and progenitors but none in
immature neurons; that signal requires sex-marker, batch, and sample-influence
QC before biological interpretation. All plots and tables remain preliminary.

## Current preprocessing work: GSE94641 reference mapping

A separate Step 00 workstream now provides the entry point for rapid published
MGE reference mapping. Start with
[`MGE_REFERENCE_MAPPING.md`](MGE_REFERENCE_MAPPING.md). The full 225-cell,
four-age reference is validated and preserved. All 450,788 Paper 3 cells now
have an E15.5-primary published-label transfer plus a complete all-age context
transfer, with raw labels, confidence, exact barcodes, and existing vendor UMAP
coordinates retained.

This preprocessing work is independent of the frozen PCDH19 classification
unit below. It does not read PCDH19 classifications or alter Steps 03–07.

## Current biological analysis: developmental-state probe detectability

Formal Step 08 joins the independent GSE94641 broad developmental states to
the frozen per-cell PCDH19 probe observations. It reports sample-level state
composition and `% A+`, `% B+`, `% C+`, and `% any PCDH19 probe+` by state and
design group. Start with
[`PCDH19_DEVELOPMENTAL_STATE_PROBE_DETECTABILITY.md`](PCDH19_DEVELOPMENTAL_STATE_PROBE_DETECTABILITY.md).

The validated Step 08 package contains all 450,788 cells. Step 08 is downstream
biological description: it neither opens nor modifies any classifier,
coefficient, threshold, probability, or Step 07 classification.
Its `output_manifest.tsv` SHA-256 is
`296a00caf2e11cc41ad81dddc154f77fb60d0b7b2887320a9a42ee87ba95b76b`.

## Frozen PCDH19 classification unit

The PCDH19 genotype-classification workflow is complete and frozen after Step
07. Start with
[`PCDH19_GENOTYPE_CLASSIFICATION_HANDOFF.md`](PCDH19_GENOTYPE_CLASSIFICATION_HANDOFF.md)
when auditing the classification methods. Do not modify the Step 03–07 models,
thresholds, classifications, or result packages. Post-freeze biological steps
must remain downstream and cannot feed results back into the classifier.

## Correct Turbo allocation

The current account is a member of `umms-ziobroj-turbo` and can read and
traverse:

```text
/nfs/turbo/umms-ziobroj
```

The source dataset discovered there is:

```text
/nfs/turbo/umms-ziobroj/Ziobro Lab/MGE scRNA-seq
```

This is service request `15662-JZ`: a mouse 10x Flex v2 Gene Expression
experiment with 12 probe-barcoded samples. It was sequenced on the Ultima
Genomics UG 100 platform and processed against the GRCm39 mouse probe set.

The earlier Miranda bulk RNA-seq directories are not Julia Ziobro data and are
not Paper 3 inputs.

## Observed data layout

The source directory contains:

- two raw paired FASTQs for the pooled gene-expression library;
- a pooled Cell Ranger multi output with raw and filtered matrices, QC reports,
  Loupe files, and alignment products;
- 12 demultiplexed sample-level filtered feature-barcode matrices;
- sample-level Cell Ranger web summaries, metrics, and Loupe files;
- a second copy of the 12 sample matrices inside the AGC Dropbox delivery;
- MD5 and demultiplexing records.

The Cell Ranger QC summary reports 450,788 called cells across the 12 samples.
Each delivered sample matrix contains 19,071 features. Sample-level called-cell
counts range from 11,085 to 60,680.

## Locations

- Version-controlled Paper 3 code/configuration: this directory
- Source registry: [`config/input_candidates.tsv`](config/input_candidates.tsv)
- Technical sample inventory:
  [`config/sample_manifest_draft.tsv`](config/sample_manifest_draft.tsv)
- Registered biological sample key:
  [`config/sample_key.csv`](config/sample_key.csv)
- Operational handoff: [`HANDOFF.md`](HANDOFF.md)
- Primary-processing handoff:
  [`PCDH19_PRIMARY_PROCESSING_HANDOFF.md`](PCDH19_PRIMARY_PROCESSING_HANDOFF.md)
- Three-candidate developing-mouse MGE reference-curation handoff:
  [`PCDH19_DEVELOPING_MOUSE_MGE_REFERENCE_CURATION_HANDOFF.md`](PCDH19_DEVELOPING_MOUSE_MGE_REFERENCE_CURATION_HANDOFF.md)
- Complete pipeline I/O and interpretation contract:
  [`PIPELINE_IO_AND_BIOLOGICAL_SCOPE.md`](PIPELINE_IO_AND_BIOLOGICAL_SCOPE.md)
- GSE94641 rapid MGE reference validation, transfer, and limitations:
  [`MGE_REFERENCE_MAPPING.md`](MGE_REFERENCE_MAPPING.md)
- Step 08 developmental-state composition and PCDH19 probe detectability:
  [`PCDH19_DEVELOPMENTAL_STATE_PROBE_DETECTABILITY.md`](PCDH19_DEVELOPMENTAL_STATE_PROBE_DETECTABILITY.md)
- Preliminary Step 09 sample-level pseudobulk differential expression:
  [`PCDH19_PRELIMINARY_PSEUDOBULK_DE.md`](PCDH19_PRELIMINARY_PSEUDOBULK_DE.md)
- Verified generated-asset and run ledger:
  [`ASSET_AND_RUN_INVENTORY.md`](ASSET_AND_RUN_INVENTORY.md)
- X-GFP construct/Flex compatibility audit:
  [`XGFP_PROBE_COMPATIBILITY_AUDIT.md`](XGFP_PROBE_COMPATIBILITY_AUDIT.md)
- Paper 3 output root:
  `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19`

On a workstation where the same Turbo allocation is mounted under `/Volumes`,
the corresponding output path is:

```text
/Volumes/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19
```

The source directory remains read-only. Paper 3 analysis products use the same
`inputs/`, `results/`, `logs/`, `jobs/`, and `final_figures/` layout as Paper 2.

### Production provenance note

The existing Step 02a-02c result packages were generated from commit
`22a342e` and record exact Python-source hashes. The current working sources
contain expanded scientific docstrings, so their file hashes intentionally
differ even though their executable Python AST is unchanged. Do not run the
documentation-edited files in place and describe that as verification by the
historical implementation. Use the recorded production commit/hash for exact
reproduction, or create a newly versioned output package. The complete
old-versus-current identities are in
[`ASSET_AND_RUN_INVENTORY.md`](ASSET_AND_RUN_INVENTORY.md).

## Safe access check

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline/paper3_pcdh19
./bin/check_candidate_access.sh
```

This checks path metadata only. It does not open, checksum, copy, or modify the
large scientific files.

## Loupe files

The delivery includes a pooled 5.3 GB `raw_cloupe.cloupe` and one
`sample_cloupe.cloupe` for each of the 12 technical samples. The per-sample
files are also duplicated in the AGC Dropbox delivery directory.

Inspection of the embedded Loupe v9 metadata found a blank run description,
no custom cell-data tables, and no PCDH19 genotype or condition labels. The
files contain the expression matrices, vendor UMAP/t-SNE, graph and k-means
clusters, and technical Cell Ranger metadata. The biological sample key was
subsequently supplied by the user and is registered separately in
[`config/sample_key.csv`](config/sample_key.csv); it was not recovered from
the Loupe or AGC delivery.

## Registered experimental groups

The user-provided sample key maps the 12 technical samples to four equal-sized
groups:

- `15662-JZ-1` through `-3`: WT male;
- `15662-JZ-4` through `-6`: WT female;
- `15662-JZ-7` through `-9`: HET female; and
- `15662-JZ-10` through `-12`: KO male.

All are recorded as mouse embryonic brain, MGE, with a submitted target of
20,000 cells. Each group contains three submitted samples. The design supports
the sex-matched contrasts HET female versus WT female and KO male versus WT
male, but it is not a complete genotype-by-sex factorial. The submitted names
(`Sample 1.1`, `Sample 1.7`, and so on) must not be interpreted as donor,
litter, batch, or biological-replicate structure without additional metadata.

## Locked Pcdh19 probe audit

The technical Pcdh19 Flex audit is implemented as a biological-metadata-free,
version-locked pipeline step. It checks the delivered panel, exact 10x v2.0.0
BED and metadata, and GRCm39-2024-A GTF by SHA-256 before extracting the three
validated Pcdh19 probes from every technical sample.

Run all 12 samples locally or in an interactive Great Lakes allocation:

```bash
cd /home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
./paper3_pcdh19/bin/run_pcdh19_probe_audit_all.sh
```

Submit the same entry point through SLURM:

```bash
sbatch paper3_pcdh19/slurm/pcdh19_probe_audit_all.sbatch
```

The runner always validates `15662-JZ-1` against the frozen prototype before
processing samples 2 through 12. Existing validated outputs are checksum-
verified and retained; incomplete or different outputs cause a hard failure.
Final products are written to:

```text
results/pcdh19_probe_audit/
  references/
    pcdh19_probe_coordinates.tsv
    reference_manifest.tsv
  per_sample/<technical_sample_id>/
    pcdh19_probe_patterns.tsv
    pcdh19_probe_summary.tsv
    validation.tsv
    checksums.sha256
  combined/
    all_samples_pcdh19_probe_summary.tsv
    all_samples_pcdh19_pattern_summary.tsv
    all_samples_validation.tsv
  software_environment.tsv
  output_manifest.tsv
```

The lock is
[`config/pcdh19_probe_audit.lock.json`](config/pcdh19_probe_audit.lock.json),
and the implementation is
[`scripts/pcdh19_probe_audit.py`](scripts/pcdh19_probe_audit.py). This step
contains technical sample IDs only by design. The registered sample key is a
separate downstream annotation layer; it does not retroactively alter the
frozen probe-audit files or their checksums.

For the exact role of every shell, SLURM, configuration, Python input, output
column, validation, and biological boundary, read
[`PIPELINE_IO_AND_BIOLOGICAL_SCOPE.md`](PIPELINE_IO_AND_BIOLOGICAL_SCOPE.md).
The validated production files record the SHA-256 of the frozen Python
implementation; documentation changes must not be represented as having
generated those existing results.

## Locked X-GFP sequence-compatibility audit

Before any GFP count interpretation, Step `02b` tests whether the custom Flex
EGFP probes could recognize the Nagy/Kalantry D4/XEGFP reporter sequence. The
delivered panel contains three included EGFP probes, and all three are unique
exact 50/50 reverse-complement matches to the Clontech EGFP CDS specified for
the original pCX-EGFP construct. Sequence compatibility therefore passes.

Run or verify the idempotent audit with:

```bash
./paper3_pcdh19/bin/run_xgfp_probe_audit.sh
```

See [`XGFP_PROBE_COMPATIBILITY_AUDIT.md`](XGFP_PROBE_COMPATIBILITY_AUDIT.md)
for the exact probes, sequence coordinates, evidence chain, limitations, and
machine-readable output paths. This result establishes theoretical probe
compatibility only; it does not establish reporter expression or observed GFP
UMIs in any sample or cell.

## Locked EGFP raw-count audit

Step `02c` applies the same frozen per-probe machinery used for Pcdh19 to the
three Step `02b`-validated EGFP probes. It reads raw probe UMIs for every
Cell Ranger filtered barcode, requires the three-probe sum to equal Cell
Ranger's EGFP feature for every barcode, and emits all eight binary patterns.
It does not normalize, create new cell calls, label cells, infer genotype, or
perform an EGFP/Pcdh19 concordance analysis.

Run all 12 samples locally:

```bash
./paper3_pcdh19/bin/run_egfp_probe_audit_all.sh
```

Or submit the identical entry point:

```bash
sbatch paper3_pcdh19/slurm/egfp_probe_audit_all.sbatch
```

JZ-1 is always checked against an independently frozen barcode-table checksum
before the remaining samples can run. Existing validated outputs are verified
and retained; different artifacts are never silently overwritten. Outputs are
written beneath `results/egfp_probe_audit/` as per-sample barcode/summary/
validation files, combined probe and eight-pattern summaries, a descriptive
sample-design join, exact probe-reference provenance, environment versions,
and a complete SHA-256 output manifest.

The completed audit found 12 raw EGFP UMIs in 12 of 450,788 filtered barcodes:
10 from `probe01`, one from `probe02`, and one from `probe03`. Every detected
barcode has one UMI and only one probe; no joint two- or three-probe pattern is
observed. This is an extremely sparse trace signal, not a hidden moderate GFP
population, and it is not used for reporter-positive cell classification.

## Step 03 PCDH19 genotype-classification setup

Step `03_pcdh19_genotype_classification_setup` is the first formal biological
classification-preparation step. It does not change or rewrite any Step 02
artifact. Instead, it checksum-validates the existing Step 02a per-cell probe
tables, joins the separately registered sample key, and selects cells from the
three known WT-male and three known KO-male technical samples as ground-truth
training candidates.

Run the step locally or in an interactive Great Lakes allocation with:

```bash
./paper3_pcdh19/bin/run_step_03_pcdh19_genotype_classification_setup.sh
```

Or submit the same entry point through SLURM:

```bash
sbatch paper3_pcdh19/slurm/step_03_pcdh19_genotype_classification_setup.sbatch
```

The object-oriented implementation is
[`scripts/Step_03_PCDH19_Genotype_Classification_Setup.py`](scripts/Step_03_PCDH19_Genotype_Classification_Setup.py).
It separates configuration, sample annotation, per-cell probe representation,
validated loading, output writing, and orchestration responsibilities. The
main intermediate output is:

```text
results/step_03_pcdh19_genotype_classification_setup/
  pcdh19_wt_ko_male_genotype_classification_ready_cells.tsv
  pcdh19_genotype_classification_ground_truth_summary.tsv
  pcdh19_genotype_classification_setup_validation.tsv
  software_environment.tsv
  output_manifest.tsv
```

The classification-ready table contains 230,269 cells and preserves technical
sample ID, submitted sample name, genotype, sex, design group, cell barcode,
raw A/B/C UMI counts, direct A/B/C detection flags, total Pcdh19 UMI, original
detection pattern, and explicit WT=0/KO=1 ground-truth targets. The label is a
registered sample-genotype target, not a probe-derived prediction. No split,
classifier, model score, confusion matrix, performance claim, or HET-cell
prediction is produced in Step 03.

## Step 04 empirical PCDH19 pattern classifier

Step `04_pcdh19_empirical_pattern_classifier` is the first formal classifier
step. It reads only the manifested Step 03 WT-male/KO-male table and uses only
the binary A/B/C detection state. For all eight patterns it reports WT, KO, and
total cells; `P(WT | pattern)` and `P(KO | pattern)`; class-conditional pattern
frequencies; raw cell ratios; and reciprocal likelihood ratios.

Run locally or interactively with:

```bash
./paper3_pcdh19/bin/run_step_04_pcdh19_empirical_pattern_classifier.sh
```

Or submit the same entry point through SLURM:

```bash
sbatch paper3_pcdh19/slurm/step_04_pcdh19_empirical_pattern_classifier.sbatch
```

Outputs are written beneath:

```text
results/step_04_pcdh19_empirical_pattern_classifier/
  step_04_pcdh19_empirical_pattern_classifier.tsv
  step_04_pcdh19_pattern_distribution.tsv
  step_04_wt_vs_ko_pattern_frequency.png
  step_04_wt_ko_conditional_probability_by_pattern.png
  step_04_pattern_cell_count_and_proportion.png
  step_04_pcdh19_empirical_pattern_classifier_validation.tsv
  software_environment.tsv
  output_manifest.tsv
```

Pattern `000` remains explicitly represented and has no hard call. Step 04
does not split data, compute a confusion matrix, fit logistic regression, load
HET cells, or generate HET predictions. The class/module ownership and planned
extensions are documented in
[`PCDH19_CLASSIFICATION_FRAMEWORK.md`](PCDH19_CLASSIFICATION_FRAMEWORK.md).

## Step 05 PCDH19 logistic-regression baseline

Step `05_pcdh19_logistic_regression_baseline` adds the first parametric model
through the existing Step 04 probability interface. It fits the unpenalized
main-effects equation
`logit(P(KO)) = intercept + beta_A*A + beta_B*B + beta_C*C`, with WT=`0` and
KO=`1`. Only the three binary detection indicators are predictors.

Run locally or interactively with:

```bash
./paper3_pcdh19/bin/run_step_05_pcdh19_logistic_regression_baseline.sh
```

Or submit through SLURM:

```bash
sbatch paper3_pcdh19/slurm/step_05_pcdh19_logistic_regression_baseline.sbatch
```

Outputs are written beneath:

```text
results/step_05_pcdh19_logistic_regression_baseline/
  step_05_pcdh19_logistic_regression_coefficients.tsv
  step_05_pcdh19_logistic_pattern_probabilities.tsv
  step_05_pcdh19_empirical_vs_logistic_comparison.tsv
  step_05_pcdh19_logistic_regression_diagnostics.tsv
  step_05_logistic_predicted_genotype_probability_by_pattern.png
  step_05_empirical_vs_logistic_probability_by_pattern.png
  step_05_logistic_coefficient_odds_ratio.png
  step_05_pcdh19_logistic_regression_validation.tsv
  software_environment.tsv
  output_manifest.tsv
  sample_level_held_out_validation/
    step_05_loso_held_out_cell_predictions.tsv
    step_05_loso_confusion_matrix.tsv
    step_05_loso_per_sample_metrics.tsv
    step_05_loso_overall_metrics.tsv
    step_05_loso_fold_model_coefficients.tsv
    step_05_loso_held_out_confusion_matrix.png
    step_05_loso_per_sample_accuracy_called.png
    step_05_loso_per_sample_percent_called.png
    step_05_loso_validation_checks.tsv
    software_environment.tsv
    output_manifest.tsv
  wt_male_female_vs_ko_male_validation/
    step_05_wt_mf_ko_m_held_out_cell_predictions.tsv
    step_05_wt_mf_ko_m_confusion_matrix.tsv
    step_05_wt_mf_ko_m_per_sample_metrics.tsv
    step_05_wt_mf_ko_m_metrics_by_group.tsv
    step_05_wt_mf_ko_m_overall_metrics.tsv
    step_05_wt_mf_ko_m_fold_model_coefficients.tsv
    step_05_wt_mf_ko_m_full_fit_coefficients.tsv
    step_05_wt_mf_ko_m_full_fit_pattern_probabilities.tsv
    step_05_male_only_wt_false_ko_errors_by_pattern.tsv
    step_05_wt_mf_ko_m_wt_false_ko_errors_by_pattern_and_sex.tsv
    step_05_ground_truth_cohort_validation_comparison.tsv
    five diagnostic PNGs
    step_05_wt_mf_ko_m_validation_checks.tsv
    software_environment.tsv
    output_manifest.tsv
  a_negative_raw_bc_umi_evidence/
    step_05_a_negative_raw_bc_cells.tsv
    step_05_a_negative_raw_bc_cohort_summary.tsv
    step_05_a_negative_raw_umi_distribution.tsv
    step_05_a_negative_umi_bin_enrichment.tsv
    step_05_a_negative_joint_bc_count_combinations.tsv
    step_05_a_negative_raw_b_c_umi_distributions.png
    step_05_a_negative_umi_bin_ko_enrichment.png
    step_05_a_negative_joint_bc_count_combinations.png
    step_05_a_negative_raw_bc_validation_checks.tsv
    software_environment.tsv
    output_manifest.tsv
  count_informed_vs_binary_validation/
    step_05_binary_vs_count_informed_held_out_predictions.tsv
    step_05_count_informed_confusion_matrix.tsv
    step_05_count_informed_per_sample_metrics.tsv
    step_05_count_informed_metrics_by_group.tsv
    step_05_count_informed_overall_metrics.tsv
    step_05_count_informed_fold_model_coefficients.tsv
    step_05_count_informed_full_fit_coefficients.tsv
    step_05_binary_vs_count_informed_overall_comparison.tsv
    step_05_binary_vs_count_informed_per_sample_comparison.tsv
    step_05_binary_vs_count_informed_correctness_transitions.tsv
    three diagnostic PNGs
    step_05_count_informed_validation_checks.tsv
    software_environment.tsv
    output_manifest.tsv
  sample_level_probe_evidence_diagnostics/
    step_05_sample_cell_flow.tsv
    step_05_genotype_cell_flow.tsv
    step_05_per_sample_probe_summary.tsv
    step_05_per_sample_Anegative_probe_summary.tsv
    step_05_per_sample_Apositive_probe_summary.tsv
    step_05_per_sample_probe_umi_distribution.tsv
    three UMI-bin composition tables
    step_05_prediction_performance_by_probe_state.tsv
    step_05_binary_vs_count_error_transitions_by_probe_state.tsv
    step_05_sample_pairwise_probe_comparison.tsv
    step_05_auc_orientation_audit.tsv
    step_05_sample_probe_interpretation.tsv
    12 diagnostic PNGs
    scope, validation, environment, and output-manifest TSVs
```

All `000` cells remain in the fit, making the intercept their baseline log
odds. The fitted `000` probabilities are P(WT)=0.5099 and P(KO)=0.4901, but no
held-out call is assigned. The validation extension leaves out one registered
sample at a time, refits on the other five, and predicts every cell in the
held-out sample. No cell-level random split is used. For non-`000` patterns,
the fixed rule calls KO when P(KO) > 0.5 and WT when P(WT) > 0.5; exact ties
remain uncalled, and the threshold is never optimized on held-out cells.

Across all six held-out folds, 43,256/230,269 cells (18.785%) are called.
Accuracy among called cells is 80.729%; with KO as the positive class,
sensitivity is 99.711% and specificity is 55.691%. These denominators and the
large class asymmetry should accompany any interpretation of the overall
accuracy. The original binary logistic classifier still excludes UMI-count
predictors, interactions, nonlinear terms, transcriptome/cell-type features,
HET cells, model selection, and threshold optimization.

The registered Step 03 `technical_sample_id` is the available holdout unit and
is labeled `biological_sample_id` in validation outputs. This is leakage-safe
at the registered-sample level, but the current metadata do not establish
donor, embryo, or litter independence.

The second Step 05 cohort keeps the male-only package byte-for-byte unchanged
and defines WT as WT male plus WT female (six samples, 233,553 cells), KO as KO
male (three samples, 116,133 cells), and HET female as excluded/not loaded. It
performs nine leave-one-sample-out folds with the same unweighted model and
fixed call rule. Of 349,686 cells, 65,314 (18.678%) are called. Called-cell
accuracy is 41.960%, KO sensitivity is 20.302%, and WT specificity is 55.049%.

Adding WT females does not improve the directly comparable WT-male result:
WT-male specificity remains 55.691% and the same 8,265 WT-male cells are called
KO. Every error is A-absent: `001` accounts for 4,586 (55.49%), `010` for 3,131
(37.88%), and `011` for 548 (6.63%); `100`, `101`, `110`, and `111` account for
zero. The expanded cohort also shows a fold-prevalence limitation: holding out
either large KO sample causes `001` and `010` to become WT-favored, sharply
reducing KO sensitivity. Step 05 records this result without adding class
weights or tuning the 0.5 threshold.

Step 05 also contains a separate descriptive module for the 327,204 controls
with exact `A_UMI == 0`: 211,142 WT cells (WT male plus WT female) and 116,062
KO-male cells. It preserves raw integer Probe B and C molecule counts without
normalization and fits no additional classifier. For B+C total, known-KO
enrichment increases from 54.255% at one UMI to 67.640% at two and 80.216% at
three or more; the corresponding KO:WT likelihood ratios are 2.158, 3.803,
and 7.376. Probe B and Probe C separately show the same monotonic trend. The
3+ B+C bin contains only 1,663 cells, so this is evidence of increasing KO
enrichment, not a validated count-based call rule or a sensitivity/specificity
claim. HET-female tables are not loaded.

The next Step 05 implementation compares the original binary model directly
with `A_detected + B_UMI + C_UMI` on the identical nine expanded-cohort LOSO
folds. B/C values remain raw probe-level UMI/ligation evidence and are not
normalized or described as transcript numbers. Both models leave exact `000`
uncalled and apply the same fixed, untuned 0.5 rule to informative cells.

Among the same 65,314 called cells, count-informed accuracy is 44.695% versus
41.960% for binary A/B/C, and KO sensitivity is 27.570% versus 20.302%. WT
specificity is essentially unchanged at 55.044% versus 55.049%. The count
model corrects 2,373 binary-model errors while regressing 587 previously
correct calls. The improvement is not sample-consistent: JZ-10 and JZ-11 gain
10.29 and 10.87 accuracy points, but JZ-12 loses 26.81 points. Brier score and
log loss improve only slightly, and pooled held-out AUC decreases from 0.254
to 0.201. Equally weighted mean per-sample accuracy also decreases from
50.613% to 49.983%. This supports modest pooled-cell information beyond
detection state, but not robust biological-sample generalization or a final
model choice. No
weighting, threshold optimization, interaction, nonlinear term, or HET
inference is introduced.

The biological-sample diagnostic verifies all 349,686 paired cells against the
raw manifested Step 02a rows and excludes no eligible WT/KO control. All
284,372 exact-zero cells remain present but uncalled; both models call the
same 65,314 cells. JZ-10/JZ-11/JZ-12 have A-negative B+C means of 0.2700,
0.2867, and 0.2536 and exact-zero proportions of 79.11%, 78.29%, and 80.40%,
respectively. JZ-12 is therefore modestly weaker in B+C evidence, not a
qualitatively different KO sample. Count-model corrections occur at B+C=2
(2,070) and 3+ (303); 585/587 regressions occur at B+C=1 and specifically in
JZ-12. The AUC audit confirms KO=`1` and P(KO) are oriented correctly; reversed
P(WT) gives 0.746/0.799 only as a diagnostic, not a corrected replacement.

## Step 06 PCDH19 HET-female inference

Step 06 applies the two immutable Step 05 WT-M+F/KO-M full-fit models to the
three registered HET-female samples. It performs no fitting, coefficient
estimation, weighting, model selection, threshold optimization, genotype
calling, or cell-type stratification.

Run directly or through Great Lakes:

```bash
./paper3_pcdh19/bin/run_step_06_pcdh19_het_female_inference.sh
sbatch paper3_pcdh19/slurm/step_06_pcdh19_het_female_inference.sbatch
```

Outputs are published atomically beneath:

```text
results/step_06_pcdh19_het_female_inference/
  step_06_het_female_cell_probabilities.tsv
  step_06_het_female_sample_summary.tsv
  step_06_het_female_pattern_distribution.tsv
  step_06_het_female_a_negative_bc_evidence_distribution.tsv
  step_06_het_and_control_probability_summary.tsv
  step_06_het_and_control_probability_histograms.tsv
  step_06_frozen_step_05_model_identity.tsv
  four step_06 diagnostic PNGs
  step_06_het_female_inference_validation.tsv
  software_environment.tsv
  output_manifest.tsv
```

The per-cell table contains all 101,102 JZ-7--9 cells, sample/barcode, raw
A/B/C probe-level UMI/ligation counts, detection states, pattern, and both
models' WT/KO probabilities. All 80,875 exact `000` cells retain auditable
probabilities but are labeled `uncalled_000`; every other row is explicitly
`probability_only_no_genotype_call`.

The HET samples are intermediate between control distributions. Their `000`
fractions are 81.26%, 81.23%, and 79.16%; A-positive fractions are 2.40%,
5.16%, and 3.99%. Among A-negative cells, B+C=2 or 3+ comprises 3.45%, 2.35%,
and 3.68%, compared with 1.23% in WT-F and 4.90% in KO-M. Among non-`000`
cells, mean count-informed P(KO) is 0.492, 0.403, and 0.458, versus 0.250 in
WT-F and 0.570 in KO-M. The control-relative figures show both low-P(KO)
A-detected evidence and higher-P(KO) A-negative B/C evidence in every HET
sample. This is consistent with mixed WT-like and KO-like probe evidence, but
it is not by itself proof of two biologically separable cell genotypes because
the control distributions overlap and the Step 05 models have known
sample-level limitations.

## Step 07 HET-female WT-like/KO-like classification

Step 07 uses the frozen Step 05 count-informed model probabilities. Before any
Step 06 HET row is opened, it evaluates all unique informative held-out-control
P(KO) values as WT-like lower-tail or KO-like upper-tail cutoffs. The
prespecified objective is at least 95% empirical precision with maximum
coverage; if the target is unavailable, the rule retains the maximum
achievable precision and maximum coverage among ties.

```bash
./paper3_pcdh19/bin/run_step_07_pcdh19_het_female_wt_ko_like_classification.sh
sbatch paper3_pcdh19/slurm/step_07_pcdh19_het_female_wt_ko_like_classification.sbatch
```

The frozen rule is:

```text
000                                  -> Uncalled_000
non-000 and P(KO) <= 0.301037619832 -> WT_like
non-000 and P(KO) >= 0.911554020713 -> KO_like
between thresholds                   -> Uncertain
```

The WT-like cutoff achieved 99.684% held-out precision and 55.042%
informative-WT sensitivity. The 95% target was not achievable for KO-like
calls. The retained best validated KO-like cutoff has 76.316% precision and
0.471% informative-KO sensitivity; only 152 held-out controls entered that
tail. It is therefore a usable KO-enriched probe-evidence category, but not a
confidence tier equivalent to WT-like.

Across held-out controls, 22,631/349,686 cells are called (6.472% of all cells,
34.650% of informative cells), 42,683 informative cells are uncertain, and
284,372 are `000`. The called-cell confusion counts are WT->WT-like 22,408,
WT->KO-like 36, KO->WT-like 71, and KO->KO-like 116.

| Sample | WT-like | KO-like | Uncertain | Uncalled-000 | WT/KO/uncertain among informative |
| --- | ---: | ---: | ---: | ---: | ---: |
| JZ-7 | 471 | 26 | 3,180 | 15,946 | 12.81% / 0.71% / 86.48% |
| JZ-8 | 1,074 | 10 | 2,820 | 16,895 | 27.51% / 0.26% / 72.23% |
| JZ-9 | 2,419 | 104 | 10,123 | 48,034 | 19.13% / 0.82% / 80.05% |
| Pooled | 3,964 | 140 | 16,123 | 80,875 | 19.60% / 0.69% / 79.71% |

Outputs are in
`results/step_07_pcdh19_het_female_wt_ko_like_classification/`. The per-cell
table is `step_07_het_female_cell_classifications.tsv`; rule, full threshold
tradeoff, overall/per-sample control validation, confusion/category tables,
four figures, validation, environment, and manifest files accompany it.

HET cells contributed zero observations to model fitting, calibration,
threshold selection, or validation. `WT_like` and `KO_like` describe inferred
PCDH19 probe evidence and are not independently observed DNA genotypes.

## Step 08 developmental-state PCDH19 detectability

Step 08 joins all Step 00 transferred labels to the exact frozen Step 02a
per-cell A/B/C evidence. It preserves all 450,788 cells and all primary/all-age
annotation fields, while summarizing composition and probe detection with
biological sample—not cell—as the replicate.

```bash
./paper3_pcdh19/bin/run_step_08_pcdh19_developmental_state_probe_detectability.sh
sbatch paper3_pcdh19/slurm/step_08_pcdh19_developmental_state_probe_detectability.sbatch
```

The direct WT-M versus KO-M comparison uses JZ-1–3 and JZ-10–12 and all 20
exact three-versus-three sample-label assignments. WT-F and HET-F are
descriptive context. Progenitors, rather than immature neurons, have the
highest PCDH19 detectability in every group. Within both assigned states, A is
appreciable in WT, nearly absent in KO, and intermediate in HET; B/C remain
KO-enriched. No reproducible WT-M/KO-M broad-state composition shift is
supported by the six samples.

Outputs are under
`results/step_08_pcdh19_developmental_state_probe_detectability/`; interpretation,
exact effects, and limitations are in
[`PCDH19_DEVELOPMENTAL_STATE_PROBE_DETECTABILITY.md`](PCDH19_DEVELOPMENTAL_STATE_PROBE_DETECTABILITY.md).

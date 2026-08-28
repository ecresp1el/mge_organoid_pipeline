# Paper 3: Ziobro PCDH19 MGE single-cell RNA-seq

This workstream is reserved for the Ziobro PCDH19 MGE single-cell paper. It is
separate from the Paper 2 cross-study atlas and uses the Ziobro Turbo
allocation, independent step numbering, and independent outputs.

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
- Complete pipeline I/O and interpretation contract:
  [`PIPELINE_IO_AND_BIOLOGICAL_SCOPE.md`](PIPELINE_IO_AND_BIOLOGICAL_SCOPE.md)
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
```

All `000` cells remain in the fit, making the intercept their baseline log
odds. The fitted `000` probabilities are P(WT)=0.5099 and P(KO)=0.4901, but no
hard call is assigned. Step 05 does not use UMI counts, interactions,
nonlinear terms, transcriptome/cell-type features, HET cells, thresholds, or
held-out performance claims.

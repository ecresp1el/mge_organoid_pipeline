# PCDH19 genotype-classification handoff

**Classification state: frozen after Step 07 on 2026-08-28.** Steps 03–07 are
one completed methodological unit. Do not refit a model, alter a
threshold, relabel a cell, or overwrite a result package. This document is the
entry point for resuming the PCDH19 classification work.

The next true bioinformatics cleanup is a separate primary-processing
sequence, reserved as Steps 03/04/10 in the overall Paper 3 workflow. It must
start from the original per-sample Cell Ranger outputs, build a new canonical
object, perform sample-aware QC/filtering, assess integration, and compute new
reductions, neighbors, clusters, and a unified UMAP. It must not use this
frozen classification unit or the Step 01 mapping RDS to define cell
retention, geometry, or clusters. Frozen classification and reference labels
may be joined only downstream after cleanup. This primary sequence has not yet
been implemented or run.

Post-freeze Step 08 has since begun a separate biological analysis of broad
developmental-state composition and PCDH19 probe detectability. It does not
open Step 07 as an input or modify this methodological unit. See
[`PCDH19_DEVELOPMENTAL_STATE_PROBE_DETECTABILITY.md`](PCDH19_DEVELOPMENTAL_STATE_PROBE_DETECTABILITY.md).

Preliminary Step 09 subsequently performs ordinary design-group pseudobulk
differential expression without loading HET WT-like/KO-like classifications.
It also remains downstream and cannot modify this frozen unit. See
[`PCDH19_PRELIMINARY_PSEUDOBULK_DE.md`](PCDH19_PRELIMINARY_PSEUDOBULK_DE.md).

An additional auxiliary workstream curates La Manno 2021, Bandler 2022, and
Mayer 2018 as candidate developing-mouse MGE references. Its first checkpoint
completed on 2026-08-29 without loading any Paper 3 cell. Its authoritative
scope and stop conditions are in
[`PCDH19_DEVELOPING_MOUSE_MGE_REFERENCE_CURATION_HANDOFF.md`](PCDH19_DEVELOPING_MOUSE_MGE_REFERENCE_CURATION_HANDOFF.md).
It is reference-only at this stage: it cannot load or map the Paper 3 query,
alter GSE94641 outputs, or feed information into this frozen classification
unit.
Its later public-atlas follow-up recovered exact Bandler study-level E13/E15
class and cluster counts from intended Shiny vector exports. A subsequent
24-gene expression-fingerprint/order join now resolves the 7,420 E15 atlas
cells as 4,481 CA301 MGE, 2,937 CA302 CGE, and 2 CA303 LGE deposited barcodes,
with later-atlas labels and public-vector plot coordinates. This auxiliary
result does not alter any Step 03–07 classification, and it is not the original
Bandler 2022 21-cluster per-cell taxonomy. Barcode-recovery job 59179410
completed successfully; failed environment attempt 59177809 remains recorded.

That auxiliary module uses versioned run packages by default and has an
explicit guarded `--replace-run RUN_ID` mode for intentional inactive-run reruns inside its
own step. Its SLURM jobs execute code copied into each run package before
submission. This does not relax the freeze here: Steps 03–07 must not be
replaced through the auxiliary workflow.

Its early checkpoint also publishes standardized per-sample and per-study
curation tables covering sample counts/IDs, age, tissue/region, design and QC,
single-cell modality, library/instrument/chemistry fields, selection bias, raw
accessions, and proven P0-object membership. It preserves all original
GEO/author metadata long-form and records later annotation/embedding/readiness
requirements as `NOT_ASSESSED` rather than silently dropping them. The module
is object-oriented, publishes an observed-object report only after all three
inspections pass, and remains fully separate from this frozen classifier.

## 1. Biological goal

The goal was to determine whether PCDH19 Flex probe-level evidence from known
WT and KO controls could support inference of WT-like versus KO-like PCDH19
states in HET-female cells. WT male, WT female, and KO male samples supplied
the final known-genotype controls. HET females were withheld from model fitting,
feature selection, calibration, threshold selection, and performance
estimation.

## 2. Probe interpretation

- Probe A is strongly WT-associated. A-detected cells drive the highly precise
  WT-like tier.
- Among A-negative cells, Probe B and Probe C evidence is enriched in KO
  controls.
- Increasing raw B, C, or B+C evidence progressively increases KO enrichment.
  For A-negative controls, P(KO) for B+C was 54.255% at one UMI, 67.640% at two,
  and 80.216% at three or more.
- Flex probe-level UMIs are UMI-deduplicated probe-pair ligation evidence. They
  are **not literal transcript counts**.
- Pattern `000` means that no informative PCDH19 probe event was observed in
  that cell. It does **not** mean biological absence of PCDH19 expression and
  must not be interpreted as KO.

## 3. Step-by-step workflow

### Step 03 — classification-ready ground-truth data

- **Purpose:** establish a clean, modular per-cell table for known WT-male and
  KO-male controls.
- **Input:** manifested Step 02 PCDH19 per-cell A/B/C probe tables plus the
  registered sample key.
- **Output:** 230,269 cells with sample, sex, genotype, barcode, raw A/B/C
  counts, detection indicators, pattern, and WT=`0`/KO=`1` target.
- **Conclusion:** the stable ground-truth table was ready for classification;
  no model was fitted.

### Step 04 — empirical eight-pattern baseline

- **Purpose:** quantify WT and KO evidence for all eight binary A/B/C states.
- **Input:** the manifested Step 03 known-control table.
- **Output:** empirical WT/KO counts, conditional probabilities, and likelihood
  ratios for `000`–`111`, with three pattern diagnostic figures.
- **Conclusion:** probe patterns carry genotype information, while `000`
  remains observed but uninformative and uncalled.

### Step 05 — logistic and probe-evidence framework

- **Purpose:** establish parametric baselines and validate them at the
  registered-sample level.
- **Inputs:** Step 03/04 artifacts, manifested WT-female Step 02 rows, and only
  known WT-M, WT-F, and KO-M controls.
- **Outputs:**
  - preserved binary `A_detected+B_detected+C_detected` logistic baseline;
  - leave-one-registered-sample-out validation, never a random cell split;
  - preserved WT-M/KO-M benchmark and expanded WT-M+F/KO-M validation;
  - descriptive A-negative raw B/C evidence analysis;
  - paired binary-versus-count-informed comparison using
    `A_detected+B_UMI+C_UMI`; and
  - biological-sample probe-evidence diagnostics.
- **Conclusion:** Probe A was strongly WT-associated and increasing raw B/C
  evidence was progressively KO-enriched. The count-informed model changed
  predictions beyond binary detection, but improvement was modest and
  sample-dependent. HET cells remained completely excluded.

### Step 06 — frozen-model HET-female inference

- **Purpose:** describe HET probe evidence before defining calls.
- **Input:** exact Step 05 full-fit binary and count-informed models plus the
  manifested JZ-7–9 HET-female probe tables.
- **Output:** 101,102 HET cells with raw evidence, patterns, and both frozen
  models' WT/KO probabilities; sample summaries and control-relative plots.
- **Conclusion:** all HET samples contained both low-P(KO), A-detected evidence
  and higher-P(KO), A-negative B/C evidence. No HET genotype call or threshold
  was created in this step.

### Step 07 — control-derived WT-like/KO-like classification

- **Purpose:** derive conservative probe-state calls and apply them to HET
  females without HET feedback.
- **Input:** immutable Step 05 count-informed held-out control probabilities;
  only after the rule was frozen, the immutable Step 06 HET probabilities.
- **Output:** the full control threshold tradeoff, frozen rule, overall and
  per-sample control validation, four-state classifications for every HET cell,
  and four concise figures.
- **Conclusion:** WT-like calls form a highly precise tier. A validated
  KO-enriched tier is identifiable but has lower precision and extremely low
  sensitivity. Most informative HET cells remain uncertain.

Implementation ownership and detailed schemas remain in
[PCDH19_CLASSIFICATION_FRAMEWORK.md](PCDH19_CLASSIFICATION_FRAMEWORK.md) and
[PIPELINE_IO_AND_BIOLOGICAL_SCOPE.md](PIPELINE_IO_AND_BIOLOGICAL_SCOPE.md).

## 4. Final frozen classification rule

The primary score is the frozen Step 05 count-informed logistic P(KO):

```text
WT_like:      P(KO) <= 0.301037619832
KO_like:      P(KO) >= 0.911554020713
Uncertain:    0.301037619832 < P(KO) < 0.911554020713
Uncalled_000: probe pattern 000, regardless of its model probability
```

The thresholds were derived **exclusively from known-genotype held-out WT-M,
WT-F, and KO-M controls before HET classification**. The prespecified target
was at least 95% empirical precision with maximum coverage. If unavailable,
the frozen fallback was maximum attainable precision and then maximum coverage
among ties. HET cells contributed zero observations to this process.

## 5. Validation performance

| Held-out-control measure | Frozen result |
| --- | ---: |
| WT-like precision | 99.684% |
| KO-like precision | 76.316% |
| WT sensitivity among informative WT controls | 55.042% |
| KO sensitivity among informative KO controls | 0.471% |
| Informative controls classified | 34.650% |
| Informative controls uncertain | 65.350% |
| All controls with pattern `000` | 81.322% |

WT-like calls are highly precise. KO-like calls represent validated
KO-enriched evidence but are substantially less precise and dramatically less
sensitive. The best KO-like tier was retained rather than discarded when it
missed the 95% target, but its 76.316% precision must always remain explicit.

**WT-like and KO-like counts must not be interpreted as directly comparable
estimates of the underlying WT:KO cellular ratio.** The asymmetric precision,
sensitivity, and coverage make such a ratio invalid.

## 6. Final HET result

The frozen Step 07 result contains 101,102 HET-female cells:

| Sample | Total | WT_like | KO_like | Uncertain | Uncalled_000 | WT / KO / uncertain among informative |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| JZ-7 | 19,623 | 471 | 26 | 3,180 | 15,946 | 12.81% / 0.71% / 86.48% |
| JZ-8 | 20,799 | 1,074 | 10 | 2,820 | 16,895 | 27.51% / 0.26% / 72.23% |
| JZ-9 | 60,680 | 2,419 | 104 | 10,123 | 48,034 | 19.13% / 0.82% / 80.05% |
| **Pooled** | **101,102** | **3,964** | **140** | **16,123** | **80,875** | **19.60% / 0.69% / 79.71%** |

Thus 80,875 cells are `Uncalled_000`, 3,964 are `WT_like`, 140 are
`KO_like`, and 16,123 are `Uncertain`. Among the 20,227 informative non-`000`
cells, 19.60% are WT-like, 0.69% are KO-like, and 79.71% are uncertain.

## 7. What can and cannot currently be concluded

### Can conclude

- PCDH19 Flex probe evidence contains reproducible genotype information.
- A-positive evidence supports highly precise WT-like classification.
- Stronger A-negative B/C evidence progressively enriches for KO.
- A control-validated KO-like tier can be identified, although with limited
  precision and extremely low sensitivity.
- HET females contain both WT-like and KO-enriched PCDH19 probe evidence.

### Cannot conclude

- That all HET cells can be genotyped.
- That `000` cells are KO.
- That WT-like:KO-like counts estimate the biological mosaic ratio.
- That KO-like calls have confidence equivalent to WT-like calls.
- That Flex probe UMIs are transcript counts.
- That WT-like or KO-like is an independently observed DNA genotype.

## 8. Exact stopping point and next analysis when resumed

**The genotype-classification framework is frozen after Step 07.** Step 08 is
not part of this methodological unit. Do not retrain, recalibrate, change the
features, move either threshold, relabel existing cells, or regenerate the
published result packages.

The resumed biological analysis is to map
the frozen HET classifications onto existing cell identities/states and ask
how WT-like, KO-like, Uncertain, and Uncalled cells are distributed across
biologically relevant PCDH19-expressing populations.

Any such analysis must consume the frozen Step 07 classifications as an
immutable input. Cell-type or state results must never be used to retroactively
alter the classifier, probability model, thresholds, or Step 07 labels.

## 9. Reproducibility map

Manifest hashes below are SHA-256 values of each published
`output_manifest.tsv`, verified read-only on 2026-08-28.

| Step | Script | Primary input | Primary output | Frozen manifest SHA-256 | Purpose |
| --- | --- | --- | --- | --- | --- |
| 03 | `scripts/Step_03_PCDH19_Genotype_Classification_Setup.py` | Step 02 per-cell probe tables + `config/sample_key.csv` | `results/step_03_pcdh19_genotype_classification_setup/pcdh19_wt_ko_male_genotype_classification_ready_cells.tsv` | `8d0627c0e93e0fa7ec37ebc0d650896565a99cc3cecb2e5e099072edf1cd0098` | Build classification-ready known controls. |
| 04 | `scripts/Step_04_PCDH19_Empirical_Pattern_Classifier.py` | Manifested Step 03 table | `results/step_04_pcdh19_empirical_pattern_classifier/step_04_pcdh19_empirical_pattern_classifier.tsv` | `2f54c0101f1379438aedbf7a627b31776ae7fde6fae2dc450d0444419553ecc0` | Quantify empirical evidence for all eight patterns. |
| 05 | `scripts/Step_05_PCDH19_Logistic_Regression_Baseline.py` | Steps 03/04, manifested WT-F rows, known WT-M/WT-F/KO-M controls | `results/step_05_pcdh19_logistic_regression_baseline/` and five manifested subpackages | Base: `1507fa4c4b6cacb5f7bae9580ed4b04d4dc2483869ebe4eec065e834c8b9f8ff`; male LOSO: `8b4012572b7bfbbbab8548d4a7b1e006aabf642139e84cd9ce6266540b28c239`; expanded: `bcb22b87e567d3ec8ccee15cc17eefd5e36abfe499f25406753b651fb8deac7d`; raw B/C: `0a27a9e5c1e983c5b1a6772435b299f72a974bac405d9a39932e0c8146dc0631`; count comparison: `cf68ee398a2d3a9e19e3de05a59be50d9aed2bfe1084b3445f63ab47e556130f`; sample diagnostic: `ad4bb28c6fce5adea9b66916378e43c2407ef6d74d67d47d65d2f10f77ecf7fe` | Fit and validate binary/count baselines; characterize probe evidence. |
| 06 | `scripts/Step_06_PCDH19_HET_Female_Inference.py` | Frozen Step 05 models + JZ-7–9 probe tables | `results/step_06_pcdh19_het_female_inference/step_06_het_female_cell_probabilities.tsv` | `8d5bf77448fee3fc67991094960fa355a826def67de0e53777d712a8079d9a2a` | Apply frozen probabilities to HET cells without calls. |
| 07 | `scripts/Step_07_PCDH19_HET_Female_WT_KO_Like_Classification.py` | Step 05 held-out controls + frozen Step 06 HET probabilities | `results/step_07_pcdh19_het_female_wt_ko_like_classification/step_07_het_female_cell_classifications.tsv` | `aec9a9cf7c8575ca453fbc61172fdd6d5c28e6c9be177d662a3c7a75cb40af9d` | Freeze control-derived thresholds and assign four-state HET labels. |

Detailed references:

- [PCDH19_CLASSIFICATION_FRAMEWORK.md](PCDH19_CLASSIFICATION_FRAMEWORK.md):
  module ownership, interfaces, and guardrails.
- [PIPELINE_IO_AND_BIOLOGICAL_SCOPE.md](PIPELINE_IO_AND_BIOLOGICAL_SCOPE.md):
  exact inputs, schemas, outputs, and biological scope.
- [ASSET_AND_RUN_INVENTORY.md](ASSET_AND_RUN_INVENTORY.md): published assets,
  run evidence, and source/result provenance.
- [HANDOFF.md](HANDOFF.md): broader Paper 3 operational context.

Frozen output directories:

```text
results/step_03_pcdh19_genotype_classification_setup/
results/step_04_pcdh19_empirical_pattern_classifier/
results/step_05_pcdh19_logistic_regression_baseline/
results/step_06_pcdh19_het_female_inference/
results/step_07_pcdh19_het_female_wt_ko_like_classification/
```
REFERENCE-MAPPING AUXILIARY STEP (DO NOT CONFUSE WITH GENOTYPE CLASSIFICATION)
-----------------------------------------------------------------------------

Standalone Step 01 run
`01_e15_mge_mapping_diagnostic_20260829_172214_a3f7ca7` completed two cell-type
mapping approaches without reading or changing any PCDH19 classifier output.
Successful jobs 59181267, 59181383, 59182039, 59182287, and 59185031 preserve
original Cell Ranger counts, sample-qualified clusters, and per-sample UMAPs.
The combined Seurat object adds only MIND and MapMyCells mapping fields. The
mapping outputs must remain covariates/annotations separate from the
genotype-classification evidence. See
`PCDH19_E15_MGE_MAPPING_DIAGNOSTIC_HANDOFF.md`.

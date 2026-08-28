# PCDH19 modular classification framework

This document defines ownership boundaries for the Paper 3 PCDH19 genotype-
classification workflow. It is intended to keep later models interoperable and
prevent scripts from duplicating input loading, pattern encoding, evaluation,
plotting, or publication logic.

## Current data flow

```text
Step 03
classification-ready known WT/KO male cells
    ↓
Step 04
empirical A/B/C pattern classifier
    ↓
Step 05
three-feature logistic-regression baseline plus
two leave-one-registered-sample-out cohorts:
  A. WT-M versus KO-M (preserved benchmark)
  B. WT-M+F versus KO-M (female-aware reference)
and a descriptive A-negative raw B/C UMI evidence module
using WT-M+F and KO-M controls (no classifier fit)
and a paired nine-fold model comparison:
  binary A/B/C detection versus
  A detection + raw B/C probe-level UMI evidence
    ↓
future model comparison and confidence logic
    ↓
eventual HET-cell inference
```

Steps 03 through 05 are implemented. Step 05 compares the preserved male-only
benchmark with the expanded WT-M+F/KO-M cohort and now compares the binary
model with one count-informed main-effects implementation. No current step
performs final model selection, class reweighting, threshold optimization, or
HET inference.

## Step ownership

| Step | Owns | Explicitly does not own |
| --- | --- | --- |
| Step 03 | Validated loading of Step 02a per-cell PCDH19 probe observations; registered WT-male/KO-male annotation; the stable classification-ready cell table. | Pattern-frequency estimation, model fitting, splitting, performance evaluation, hard calls, or HET inference. |
| Step 04 | Binary A/B/C encoding; empirical eight-pattern WT/KO counts and probabilities; likelihood-ratio evidence; descriptive pattern diagnostics. | Count-based predictors, logistic regression, train/test splitting, confusion matrices, hard calls, confidence thresholds, or HET inference. |
| Step 05 | The original unpenalized binary A/B/C logistic baseline; immutable WT-M/KO-M and WT-M+F/KO-M packages; leave-one-registered-sample-out validation; descriptive raw B/C UMI evidence; and a paired count-informed model using `A_detected + B_UMI + C_UMI` on the same expanded folds. | Random cell splitting, class weighting, threshold optimization, forced calls for `000`, count normalization, treating probe UMIs as transcript numbers, interactions, nonlinear terms, transcriptome/cell-type features, final model selection, or HET inference. |
| Future model steps | A new model implementation behind a shared probability interface. | Reimplementation of Step 03 loading or Step 04 pattern encoding and publication conventions. |
| Future validation extensions | Confirmed biological-unit validation, classifier comparison, and any expanded performance outputs. | Refitting hidden model logic inside evaluation code or optimizing a rule on held-out observations. |
| Future HET application step | Application of a selected, validated, frozen model to HET cells. | Model selection or performance claims based on HET predictions. |

## Current modules

### Step 03: classification-ready data

Implementation:
`scripts/Step_03_PCDH19_Genotype_Classification_Setup.py`.

- `StepConfiguration` owns the locked Step 03 schema and known class mapping.
  It does not load cell rows or define a model.
- `SampleAnnotation` represents registered sample identity, genotype, sex,
  design group, and ground-truth target. It does not infer any annotation from
  probe counts.
- `ProbeObservation` represents and validates one cell's raw A/B/C counts,
  total, and existing pattern. It does not label genotype.
- `PCDH19ProbeCellLoader` validates the Step 02a manifest and streams registered
  ground-truth cells. It does not transform predictors or fit a classifier.
- `ClassificationReadyTableWriter` owns deterministic Step 03 tables,
  validation, provenance, and atomic publication. It contains no modeling
  logic.
- `PCDH19GenotypeClassificationSetup` only coordinates those components.

Step 03's main table remains the stable male-only input contract. Step 05 does
not rewrite it. The expanded female-aware cohort extends the exact Step 03 male
records with a dedicated manifested reader for the three Step 02a WT-female
tables because those rows are deliberately absent from the historical Step 03
output. It uses the same registered metadata and binary encoder and never opens
the HET-female tables.

### Step 04: empirical pattern model

Implementation:
`scripts/Step_04_PCDH19_Empirical_Pattern_Classifier.py`.

- `Step04Configuration` owns the locked Step 03 input identity, expected
  classes, eight-pattern order, numerical serialization, and plotting scope.
  It contains no scientific calculation.
- `ProbePatternEncoder` is the canonical A/B/C binary encoder. It maps
  `A_detected`, `B_detected`, and `C_detected` to `000` through `111` in that
  exact order and supplies human-readable labels. It does not inspect UMI
  magnitudes or genotype.
- `Step03ClassificationTableReader` verifies the Step 03 manifest/table and
  streams only a pattern code and known WT/KO class to estimation. It rejects
  HET or non-male input and does not calculate probabilities.
- `LabeledPatternObservation` is the minimal fit record: encoded pattern plus
  known class. No raw count enters the estimator.
- `EmpiricalPatternEstimator` counts WT and KO observations for all eight
  states and constructs the model. It does not predict HET cells or select
  thresholds.
- `ProbabilisticClassifier` defines the small shared model interface:
  `predict_proba(features)` returns WT and KO probabilities. Future
  probabilistic models should implement this interface without changing the
  evaluator's model-facing contract.
- `EmpiricalPatternClassifier` is the fitted Step 04 model object. It exposes
  pattern evidence and `predict_proba({"pattern_code": code})`. It returns
  probabilities only and has no hard-call method.
- `ClassificationEvaluator` currently performs descriptive integrity checks
  and creates pattern-distribution data. It explicitly does not calculate
  accuracy, sensitivity, specificity, confusion matrices, or held-out
  performance; those require a future validation design.
- `ClassificationPlotter` renders plots from completed model/distribution rows.
  It does not read cells or recompute the classifier.
- `Step04OutputPublisher` owns TSV serialization, plots, software/input
  provenance, manifests, restart verification, and atomic publication. It does
  not own scientific estimation.
- `PCDH19EmpiricalPatternClassifierStep` only coordinates these components.

### Step 05: logistic-regression baseline

Implementation:
`scripts/Step_05_PCDH19_Logistic_Regression_Baseline.py`.

Step 05 imports the exact checksum-locked Step 04 framework script rather than
creating another encoder, ground-truth loader, or classifier protocol.

- `Step05Configuration` owns the locked predictors, outcome encoding,
  no-interaction/no-penalty constraints, optimizer settings, and upstream
  identities. It does not fit the model.
- `Step04EmpiricalModelReader` validates and loads the serialized Step 04 model
  for probability comparison. It does not refit or alter Step 04.
- `LogisticRegressionEstimator` fits an unpenalized grouped-binomial maximum-
  likelihood model by Newton/IRLS. The eight grouped pattern rows are exactly
  likelihood-equivalent to expanding the 230,269 cells because all cells in a
  pattern share the same predictors. It does not create calls or metrics.
- `LogisticRegressionClassifier` implements the existing
  `ProbabilisticClassifier.predict_proba(features)` interface. It accepts the
  canonical pattern code or the three direct binary features and returns WT
  and KO probabilities. It also exposes the linear predictor and fitted
  coefficients; it has no hard-call method.
- `LogisticModelEvaluator` compares Step 05 probabilities with the manifested
  Step 04 empirical probabilities and reports descriptive probability
  differences. It does not own sample splitting or held-out calls.
- `LogisticRegressionPlotter` renders probability, empirical-comparison, and
  coefficient/odds-ratio plots from completed model rows. It does not fit or
  load the model.
- `Step05OutputPublisher` owns Step 05 serialization, validation/provenance,
  plots, manifest verification, and atomic publication for the original base
  package. It has no fitting responsibility.
- `SampleAwarePatternRecord` adds registered sample identity and cell barcode
  plus sex/design group to the minimal pattern/ground-truth record. It does not
  predict a class.
- `Step03SampleAwarePatternReader` verifies and reads the Step 03 contract
  while retaining the registered `technical_sample_id` used for holdout. It
  does not split cells randomly or claim donor/embryo/litter independence.
- `ManifestedProbePatternCohortReader` validates the Step 03 WT-M/KO-M records,
  registered sample key, Step 02a manifest, and three WT-female per-cell tables.
  It explicitly excludes JZ-7--9 HET females and does not infer genotype from
  probes.
- `HeldOutCallingPolicy` owns the fixed decision rule: `000` is always
  `uncalled`; otherwise KO is called when P(KO) > 0.5, WT when P(WT) > 0.5,
  and an exact tie is `uncalled`. It does not fit or tune the threshold.
- `LeaveOneSampleOutLogisticValidator` holds out each registered sample in
  turn, fits `LogisticRegressionEstimator` on all other cohort samples, and
  produces held-out cell probabilities. It supports both the six-sample and
  nine-sample definitions without changing the model.
- `HeldOutValidationEvaluator` creates called-cell confusion and overall/per-
  sample metrics with KO as the positive class. It does not fit a model or
  change calls.
- `HeldOutValidationPlotter` draws only the held-out confusion matrix, per-
  sample called-cell accuracy, and per-sample call percentage. It does not
  recompute results.
- `ExistingMaleOnlyValidationReader` verifies the historical male-only package
  byte-for-byte and exposes its predictions for comparison; it never rewrites
  that benchmark.
- `PatternErrorAnalyzer` attributes WT-to-KO calls to the seven informative
  A/B/C patterns, separately for the male benchmark and expanded WT sexes.
- `ExpandedCohortEvaluator` owns WT-M, WT-F, and KO-M group metrics and the
  direct male-only-versus-expanded comparison. It does not refit models.
- `ExpandedFullFitBuilder` fits the same unweighted three-feature model on all
  nine expanded ground-truth samples and publishes probabilities without hard
  cell calls.
- `ExpandedValidationPlotter` owns the expanded generalization, pattern-error,
  and cohort-comparison plots. It contains no fitting logic.
- `ExistingStep05BasePackageVerifier` protects the already published Step 05
  base fit and its historical manifest from being rewritten by the extension.
- `ANegativeRawCountRecord` represents one known-genotype A-negative cell and
  preserves sample, sex, design group, barcode, and integer raw A/B/C molecule
  counts. It does not encode a model feature or make a genotype call.
- `ANegativeRawCountCohortReader` reuses the exact Step 03 WT-M/KO-M rows and
  reads only the three manifested WT-female Step 02a tables. It requires
  `A_UMI == 0`, verifies raw count arithmetic, and never opens HET-female
  tables.
- `RawBCUMIEvidenceSummarizer` computes exact B, C, and B+C count
  distributions, `0`/`1`/`2`/`3+` enrichment summaries, and observed joint
  B/C count combinations. It fits no classifier and applies no threshold.
- `RawBCUMIEvidencePlotter` renders distributions, enrichment by molecule
  count, and joint-count diagnostics from completed summary rows. It does not
  load cells or recompute counts.
- `ExistingExpandedGroundTruthPackageVerifier` protects the already published
  expanded-cohort package before the raw-count module runs.
- `RawBCUMIEvidenceOutputPublisher` owns the separately manifested raw-count
  tables, plots, provenance, restart validation, and atomic publication. It
  contains no scientific classification logic.
- `CountInformedCellRecord` represents one comparator cell with raw A/B/C
  probe-level UMI/ligation evidence. It exposes exactly `A_detected`, `B_UMI`,
  and `C_UMI` as model features and does not interpret counts as transcript
  numbers.
- `CountInformedGroundTruthReader` loads all nine manifested WT-M/WT-F/KO-M
  control tables and requires an exact cell, label, sex, group, and binary
  pattern match to the immutable expanded binary predictions. It never opens
  HET-female tables.
- `CountInformedLogisticEstimator` fits the unpenalized main-effects model by
  exact grouped-binomial IRLS. Grouping identical count vectors is
  likelihood-equivalent to expanding every cell; it adds no weighting,
  interaction, transformation, or threshold.
- `CountInformedLogisticClassifier` implements the shared
  `predict_proba(features)` contract for `A_detected`, raw `B_UMI`, and raw
  `C_UMI`. B/C coefficients are per additional probe-level UMI/ligation event.
- `LeaveOneSampleOutCountLogisticValidator` refits this model on the same nine
  registered-sample folds used by the expanded binary comparator.
- `PairedBinaryCountModelEvaluator` matches predictions cell for cell and
  reports called-cell metrics, Brier score, log loss, pooled held-out AUC,
  per-sample accuracy, and correctness transitions. It does not select a
  winner or alter calls.
- `CountModelComparisonPlotter` renders only the overall called metrics,
  paired per-sample accuracy, and full-fit odds ratios.
- `ExistingRawBCEvidencePackageVerifier` protects the descriptive raw-count
  package byte-for-byte after the count-informed implementation changes the
  Step 05 source.
- `CountModelComparisonOutputPublisher` owns the separate paired prediction,
  evaluation, plot, provenance, restart, and manifest contract.
- `HeldOutValidationOutputPublisher` owns serialization, provenance,
  validation, plots, manifest verification, and atomic publication for the
  sample-level validation subpackage.
- `ExpandedCohortOutputPublisher` owns the separate WT-M+F/KO-M package and its
  18-file manifest. It does not modify the base or male-only packages.
- `PCDH19LogisticRegressionBaselineStep` coordinates the reused Step 04
  components, the original Step 05 fit, and the held-out extension.

Step 05 retains Step 03's numeric targets: WT=`0`, KO=`1`. Its model is:

```text
logit(P(KO)) = intercept
             + beta_A * A_detected
             + beta_B * B_detected
             + beta_C * C_detected
```

Thus coefficients and odds ratios describe KO-versus-WT odds. The intercept is
the `000` log odds because all predictors are zero. Pattern `000` is included
in fitting and remains uncalled regardless of which fitted probability is
slightly greater than 0.5.

For held-out evaluation, the available registered sample key is
`technical_sample_id`, exposed as `biological_sample_id` in outputs. The
male-only benchmark holds out one of six samples; the expanded cohort holds out
one of nine WT-M/WT-F/KO-M samples. This prevents cell-level leakage, but the
sample key alone does not establish donor, embryo, or litter independence.

The expanded cohort is intentionally unweighted: cells contribute to the
likelihood exactly as in the original baseline. Adding WT females therefore
changes the observed class prevalence in each fold. It does not change the
WT-male call rule: the same 8,265 WT-male cells remain false KO, entirely from
`001` (4,586), `010` (3,131), and `011` (548). When either large KO sample is
held out, the expanded training fold makes `001` and `010` WT-favored at the
fixed 0.5 threshold; aggregate KO sensitivity consequently falls to 20.302%.
This is a documented baseline result, not a reason to silently add class
weights or tune the threshold.

The raw-count module is a descriptive extension of Step 05, not a change to
the logistic model. It restricts the registered WT-M+F and KO-M controls to
exact `A_UMI == 0`, preserves integer molecule counts without normalization,
and asks whether larger B/C counts carry progressively stronger KO enrichment.
Across 327,204 A-negative controls, P(KO) for B+C total is 54.255% at one UMI,
67.640% at two UMIs, and 80.216% at three or more UMIs; the corresponding
KO:WT likelihood ratios are 2.158, 3.803, and 7.376. Probe B and Probe C show
the same monotonic pattern. These are descriptive cell-level enrichment
estimates, not specificity from a classifier, and the 3+ B+C bin contains only
1,663 cells. No HET cell is loaded.

The count-informed comparison is a new classifier implementation inside Step
05, but it leaves the original binary packages untouched. Both models use the
same 349,686 WT-M+F/KO-M cells, nine leave-one-registered-sample-out folds,
fixed 0.5 rule for informative cells, and uncalled exact `000` state. The new
model is:

```text
logit(P(KO)) = intercept
             + beta_A * A_detected
             + beta_B * raw B probe UMI
             + beta_C * raw C probe UMI
```

Called-cell accuracy increases from 41.960% to 44.695% and KO sensitivity from
20.302% to 27.570%, while WT specificity is essentially unchanged (55.049% to
55.044%). The count model corrects 2,373 binary errors and introduces 587
regressions. However, the change is not consistent across biological samples:
JZ-10 and JZ-11 improve by 10.29 and 10.87 percentage points, whereas JZ-12
declines by 26.81 points. All-cell Brier score and log loss improve only
slightly, while pooled held-out AUC decreases from 0.254 to 0.201. Therefore
raw B/C evidence changes predictions and yields modest aggregate called-cell
improvement, but the equally weighted mean per-sample accuracy decreases from
50.613% to 49.983%. It does not demonstrate robust sample-level
generalization.

## Step 04 empirical model definition

For pattern `p` and known classes WT and KO, Step 04 reports:

```text
P(WT | p) = WT cells with p / all cells with p
P(KO | p) = KO cells with p / all cells with p

P(p | WT) = WT cells with p / all WT cells
P(p | KO) = KO cells with p / all KO cells

WT:KO likelihood ratio = P(p | WT) / P(p | KO)
KO:WT likelihood ratio = P(p | KO) / P(p | WT)
```

The output also reports raw WT:KO and KO:WT cell-count ratios so those are not
confused with class-conditional likelihood ratios. There is no pseudocount or
smoothing. A zero denominator is represented as `Inf` when the numerator is
positive and as a missing ratio when both are zero.

The `hard_call` field is blank for every row. Pattern `000` is retained as an
observed no-probe-detection state; its presence in the empirical model does not
make it an informative or callable genotype state.

## Current output contracts

Step 03 output root:

```text
results/step_03_pcdh19_genotype_classification_setup/
```

Step 04 output root:

```text
results/step_04_pcdh19_empirical_pattern_classifier/
```

Step 05 output root:

```text
results/step_05_pcdh19_logistic_regression_baseline/
```

Its sample-level held-out extension is:

```text
results/step_05_pcdh19_logistic_regression_baseline/
  sample_level_held_out_validation/
```

The preserved benchmark comparison and expanded cohort are:

```text
results/step_05_pcdh19_logistic_regression_baseline/
  wt_male_female_vs_ko_male_validation/
```

The descriptive A-negative raw-count package is:

```text
results/step_05_pcdh19_logistic_regression_baseline/
  a_negative_raw_bc_umi_evidence/
```

The paired binary-versus-count model package is:

```text
results/step_05_pcdh19_logistic_regression_baseline/
  count_informed_vs_binary_validation/
```

The Step 04 classifier TSV is the first serialized model artifact. The
distribution TSV and plots are diagnostics derived from it. Validation,
environment, and output-manifest files protect input identity, computation
scope, code/dependency identity, and published bytes.

## Current baselines and planned modular extensions

The first seven entries are implemented. The remaining entries are planned
extensions only.

1. **Empirical pattern classifier — Step 04 (implemented baseline).** Preserve as
   the simplest interpretable reference model.
2. **Logistic regression baseline — Step 05 (implemented baseline).** Preserve
   as the simplest parametric model using the same binary predictors and shared
   probability interface.
3. **Sample-level validation — Step 05 (implemented first pass).** Leave one
   registered sample out, refit on the other samples, and evaluate held-out
   cells. No random cell split is used. Future work should confirm the true
   biological independence unit before stronger replicate-level claims.
4. **Classifier evaluation/confusion matrix — Step 05 (implemented first
   pass).** Consume held-out probabilities and a separately defined fixed
   calling policy. Report confusion and metrics only among called cells while
   reporting the call denominator separately.
5. **Ground-truth cohort comparison — Step 05 (implemented).** Preserve the
   WT-M/KO-M result and compare it with WT-M+F/KO-M under the same unweighted
   model and fixed call rule. This is a cohort comparison, not a comparison of
   different classifier families.
6. **A-negative raw B/C UMI evidence — Step 05 (implemented descriptive
   module).** Compare unnormalized B, C, B+C, and joint B/C molecule counts in
   known WT-M+F and KO-M controls. Do not treat this analysis as a fitted
   count-based classifier.
7. **Count-informed logistic comparison — Step 05 (implemented).** Compare
   `A_detected + B_UMI + C_UMI` directly with the immutable binary A/B/C model
   on identical expanded-cohort folds, cells, call eligibility, and threshold.
   Preserve probe counts as raw UMI/ligation evidence, not transcript numbers.
8. **Model comparison.** Compare frozen candidate models on the same split and
   evaluation contract rather than allowing each model script to define its
   own denominators.
9. **Confidence/uncertain-call logic.** Keep probability estimation separate
   from decision thresholds. Preserve uninformative/uncertain states rather
   than forcing WT or KO.
10. **Application of a validated classifier to HET cells.** Load HET cells only
   after model selection and validation are complete, record the frozen model
   identity, and keep inference separate from performance estimation.

## Guardrails for future steps

- Reuse the Step 03 table contract for the male benchmark; use the registered,
  manifested Step 02a cohort reader when female ground truth is required.
- Reuse `ProbePatternEncoder` whenever the predictor is the binary A/B/C state.
- Implement new probabilistic classifiers behind `predict_proba(features)` so
  they can share later validation and comparison machinery.
- Keep feature computation, fitting, probability prediction, decision rules,
  evaluation, plotting, and I/O in separate components.
- Do not interpret cell counts as independent biological replicates.
- Do not introduce HET cells into model fitting, model selection, or threshold
  tuning.
- Do not silently turn missing probe detection into a genotype call.
- Register every new formal step in the pipeline I/O/scope document and protect
  upstream/output identities through the existing lock and manifest framework.

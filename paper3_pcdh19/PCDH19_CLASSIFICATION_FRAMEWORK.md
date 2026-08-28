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
leave-one-registered-sample-out validation
    ↓
future model comparison and confidence logic
    ↓
eventual HET-cell inference
```

Steps 03 through 05 are implemented. Step 05 now includes the first sample-
held-out evaluation, but no current step performs model selection, threshold
optimization, or HET inference.

## Step ownership

| Step | Owns | Explicitly does not own |
| --- | --- | --- |
| Step 03 | Validated loading of Step 02a per-cell PCDH19 probe observations; registered WT-male/KO-male annotation; the stable classification-ready cell table. | Pattern-frequency estimation, model fitting, splitting, performance evaluation, hard calls, or HET inference. |
| Step 04 | Binary A/B/C encoding; empirical eight-pattern WT/KO counts and probabilities; likelihood-ratio evidence; descriptive pattern diagnostics. | Count-based predictors, logistic regression, train/test splitting, confusion matrices, hard calls, confidence thresholds, or HET inference. |
| Step 05 | Unpenalized main-effects logistic regression using only binary A/B/C detection; coefficient/odds-ratio reporting; probability comparison with Step 04; leave-one-registered-sample-out validation; a fixed, separately owned calling rule; called-cell metrics and generalization plots. | Random cell splitting, threshold optimization, forced calls for `000`, UMI counts, interactions, nonlinear terms, transcriptome/cell-type features, model selection, or HET inference. |
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

Step 03's main table is the stable input contract for classification work. New
models should consume that contract rather than reopen Cell Ranger matrices or
duplicate Step 02a extraction.

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
  to the minimal pattern/ground-truth record. It does not predict a class.
- `Step03SampleAwarePatternReader` verifies and reads the Step 03 contract
  while retaining the registered `technical_sample_id` used for holdout. It
  does not split cells randomly or claim donor/embryo/litter independence.
- `HeldOutCallingPolicy` owns the fixed decision rule: `000` is always
  `uncalled`; otherwise KO is called when P(KO) > 0.5, WT when P(WT) > 0.5,
  and an exact tie is `uncalled`. It does not fit or tune the threshold.
- `LeaveOneSampleOutLogisticValidator` holds out each registered sample in
  turn, fits `LogisticRegressionEstimator` on the other five samples, and
  produces held-out cell probabilities. It does not calculate metrics or
  draw plots.
- `HeldOutValidationEvaluator` creates called-cell confusion and overall/per-
  sample metrics with KO as the positive class. It does not fit a model or
  change calls.
- `HeldOutValidationPlotter` draws only the held-out confusion matrix, per-
  sample called-cell accuracy, and per-sample call percentage. It does not
  recompute results.
- `ExistingStep05BasePackageVerifier` protects the already published Step 05
  base fit and its historical manifest from being rewritten by the extension.
- `HeldOutValidationOutputPublisher` owns serialization, provenance,
  validation, plots, manifest verification, and atomic publication for the
  sample-level validation subpackage.
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

For held-out evaluation, the only available registered sample key is the Step
03 `technical_sample_id`. Step 05 exposes it as `biological_sample_id` in
validation outputs and holds out one of the six registered WT-male/KO-male
samples at a time. This prevents cell-level leakage, but the sample key alone
does not establish donor, embryo, or litter independence.

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

The Step 04 classifier TSV is the first serialized model artifact. The
distribution TSV and plots are diagnostics derived from it. Validation,
environment, and output-manifest files protect input identity, computation
scope, code/dependency identity, and published bytes.

## Current baselines and planned modular extensions

The first four entries are implemented. The remaining entries are planned
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
5. **Model comparison.** Compare frozen candidate models on the same split and
   evaluation contract rather than allowing each model script to define its
   own denominators.
6. **Confidence/uncertain-call logic.** Keep probability estimation separate
   from decision thresholds. Preserve uninformative/uncertain states rather
   than forcing WT or KO.
7. **Application of a validated classifier to HET cells.** Load HET cells only
   after model selection and validation are complete, record the frozen model
   identity, and keep inference separate from performance estimation.

## Guardrails for future steps

- Reuse the Step 03 table contract for known WT/KO model development.
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

# Paper 3 standalone primary-processing workflow plan

Status: design registered; no processing code, run directory, or SLURM job has
been created or executed.

## Purpose and boundary

This is the planned **true primary bioinformatics cleanup and analysis** of the
12 Paper 3 samples. It is a standalone workflow named `primary_processing`
with its own `00`–`09` numbering. These numbers do not refer to the frozen
PCDH19 classification steps, either reference-curation workflow, or either
reference-mapping workflow.

The expression source is the original read-only Cell Ranger delivery:

`/nfs/turbo/umms-ziobroj/Ziobro Lab/MGE scRNA-seq/10x_analysis_15662-JZ/Sample_15662-JZ-P01`

The registered [`../config/sample_key.csv`](../config/sample_key.csv) may be
joined as sample metadata. The 31.84-GB Step 01 mapping RDS, transferred
labels, Cell Ranger UMAPs, and Cell Ranger graph clusters are not
primary-processing inputs. Reference and frozen classification fields may be
joined only in Step 09, after the cleaned geometry and clusters have been
finalized.

## Decision-control model

No scientific threshold or algorithm choice may be hidden in a script. The
human-readable source of decisions is
[`config/decision_registry.tsv`](config/decision_registry.tsv). Each row has a
stable decision ID, workflow stage, planned configuration key, status, allowed
values, evidence requirement, planned implementation class, and affected
outputs.

The implementation must use four layers:

1. `decision_registry.tsv` exposes every decision and its rationale. A value
   remains `TBD` until it is approved or resolved by a declared diagnostic.
2. A validated, user-editable workflow configuration will hold the executable
   values. Scripts may not provide silent scientific defaults.
3. Each submitted run will copy the exact registry, executable configuration,
   scripts, shell wrappers, and SLURM files into its run package before
   execution. Copies of what ran are required; hashes alone are insufficient.
4. Each stage will publish a resolved decision ledger containing requested
   values, resolved values, source/rationale, and the outputs affected. The
   ledger and exact code copies make the run reproducible.

Changing a decision intentionally creates a new versioned run by default. A
guarded `--replace-run RUN_ID` mode may overwrite only an explicitly selected,
inactive run inside the same standalone stage. Source Cell Ranger files and
upstream result packages remain read-only.

## Planned stage map

All stage roots will live under:

`$PAPER3_ROOT/results/primary_processing/<stage>/<run_id>/`

Every successful stage package will contain `code/`, `config/`, `objects/`,
`tables/`, `figures/`, `logs/`, `provenance/`, a stage report, an output
manifest, validation results, and `SUCCESS.txt`.

| Stage | Purpose | Exact upstream inputs | Required primary outputs | Planned driver and main object |
| --- | --- | --- | --- | --- |
| `00_input_registration` | Resolve and validate the 12 Cell Ranger samples without loading a derived analysis object. | Read-only Cell Ranger root; `paper3_pcdh19/config/sample_key.csv`; path/environment configuration. | `tables/cellranger_input_manifest.tsv`; `tables/sample_input_inventory.tsv`; `tables/input_validation.tsv`; `PRIMARY_00_INPUT_REPORT.md`. | `scripts/00_register_inputs.py`; `CellRangerInputRegistrar`. |
| `01_canonical_counts` | Load Cell Ranger count matrices, create globally unique cell IDs, join sample metadata, and preserve raw integer counts without biological filtering. | Successful Step 00 manifest; exact per-sample Cell Ranger filtered feature-barcode matrices; approved feature-type rule; sample key. | `objects/paper3_canonical_raw_counts.<format>`; `tables/cell_identity.tsv.gz`; `tables/feature_inventory.tsv`; `tables/import_summary.tsv`; `PRIMARY_01_CANONICAL_COUNTS_REPORT.md`. | `scripts/01_build_canonical_counts.py`; `CanonicalCountObjectBuilder`. |
| `02_qc_diagnostics` | Compute and display per-cell and per-sample QC evidence without removing cells. | Manifest-locked Step 01 canonical raw-count object; optional Cell Ranger raw-droplet matrices only for declared ambient-RNA diagnostics. | `objects/paper3_qc_metrics_unfiltered.<format>`; `tables/per_cell_qc_metrics.tsv.gz`; `tables/per_sample_qc_distributions.tsv`; diagnostic figures; `PRIMARY_02_QC_DIAGNOSTIC_REPORT.md`. | `scripts/02_compute_qc_diagnostics.py`; `QCDiagnosticWorkflow`. |
| `03_qc_filtering` | Apply only approved, visible cell/gene filtering rules and preserve every keep/drop reason. | Manifest-locked Step 02 object and diagnostics; approved QC decisions from the registry/configuration. | `objects/paper3_qc_filtered_counts.<format>`; `tables/cell_disposition.tsv.gz`; `tables/gene_disposition.tsv.gz`; `tables/filter_retention_by_sample.tsv`; before/after figures; `PRIMARY_03_QC_FILTER_REPORT.md`. | `scripts/03_apply_qc_filtering.py`; `QCFilterWorkflow`. |
| `04_artifact_assessment` | Evaluate doublets, ambient signal, and other declared artifacts; remove or retain them only under approved rules. | Manifest-locked Step 03 object; optional raw-droplet matrices inventoried in Step 00; declared artifact methods and thresholds. | `objects/paper3_artifact_clean_counts.<format>`; `tables/artifact_scores.tsv.gz`; updated `tables/cell_disposition.tsv.gz`; per-sample artifact summaries and figures; `PRIMARY_04_ARTIFACT_REPORT.md`. | `scripts/04_assess_artifacts.py`; `ArtifactAssessmentWorkflow`. |
| `05_normalization_features` | Preserve raw counts, create the selected normalized representation, and select analysis features reproducibly. | Manifest-locked Step 04 cleaned count object; approved normalization, transformation, covariate, and feature-selection settings. | `objects/paper3_normalized_features.<format>` with raw counts retained; `tables/selected_features.tsv`; normalization/HVG diagnostics; `PRIMARY_05_NORMALIZATION_REPORT.md`. | `scripts/05_normalize_select_features.py`; `NormalizationFeatureWorkflow`. |
| `06_unintegrated_baseline` | Establish the mandatory no-correction baseline before any integration decision. | Manifest-locked Step 05 object; approved PCA/latent, neighbor, clustering-diagnostic, and visualization settings. | `objects/paper3_unintegrated_baseline.<format>`; PCA/variance tables; sample/design mixing diagnostics; sample-faceted and joint diagnostic views; `PRIMARY_06_UNINTEGRATED_REPORT.md`. | `scripts/06_build_unintegrated_baseline.py`; `UnintegratedBaselineWorkflow`. |
| `07_integration_decision` | Compare no correction with explicitly approved candidate corrections and select one representation using declared biological-preservation and batch-mixing criteria. Integration is not automatic. | Step 06 baseline; approved candidate methods; a valid nuisance/batch variable distinct from genotype; evaluation metrics and acceptance rules. | `objects/integration_candidates/`; `tables/integration_benchmark.tsv`; `tables/integration_decision.tsv`; comparison figures; selected representation pointer; `PRIMARY_07_INTEGRATION_REPORT.md`. | `scripts/07_evaluate_integration.py`; `IntegrationBenchmarkWorkflow` with strategy objects. |
| `08_graph_clusters_umap` | Construct the final neighbor graph, evaluate clustering resolution/stability, choose the final clusters, and compute the unified UMAP from the selected representation. | Manifest-locked Step 07 decision and selected representation; approved graph, clustering, stability, UMAP, and random-seed settings. | `objects/paper3_primary_processed.<format>`; `tables/cluster_membership.tsv.gz`; `tables/clustering_stability.tsv`; final QC/sample/design/cluster UMAP figures; `PRIMARY_08_GRAPH_UMAP_REPORT.md`. | `scripts/08_build_graph_clusters_umap.py`; `GraphClusterEmbeddingWorkflow`. |
| `09_annotation_publish` | Characterize finalized clusters and optionally attach independent reference/classification evidence without changing cleanup, graph, clusters, or UMAP. | Manifest-locked Step 08 object; marker settings; explicitly selected external label tables and frozen classification fields joined by validated stable cell ID. | `objects/paper3_primary_annotated.<format>`; optional interoperable export; marker tables; annotation concordance/confidence tables; final decision ledger; complete run report; `PRIMARY_09_FINAL_REPORT.md`. | `scripts/09_annotate_publish.py`; `AnnotationPublicationWorkflow`. |

`<format>` is deliberately unresolved. The analysis backend and canonical
object format are visible engineering decisions in the registry; they must be
approved before implementation rather than being implied by a filename.

## Planned object-oriented package

Stage drivers will orchestrate reusable objects rather than contain scientific
logic. The planned shared package is `primary_processing/scripts/pipeline/`:

| Planned module | Responsibility | Principal planned objects |
| --- | --- | --- |
| `configuration.py` | Load and validate every decision; reject missing/TBD executable settings. | `PrimaryProcessingConfiguration`, `DecisionRegistry`, `ConfigurationValidator`. |
| `contracts.py` | Define stage input/output schemas and stable cell/feature identities. | `StageContract`, `CellIdentityContract`, `FeatureContract`. |
| `cellranger.py` | Resolve and read only manifested Cell Ranger artifacts. | `CellRangerInputRegistrar`, `CellRangerMatrixRepository`. |
| `qc.py` | Compute QC evidence and apply separately configured keep/drop policies. | `QCDiagnosticWorkflow`, `QCMetricCalculator`, `QCFilterPolicy`, `QCFilterWorkflow`. |
| `artifacts.py` | Score artifacts separately from the decision to remove them. | `ArtifactScorer`, `DoubletStrategy`, `AmbientStrategy`, `ArtifactAssessmentWorkflow`. |
| `normalization.py` | Preserve counts and construct the configured normalized/feature representation. | `NormalizationStrategy`, `FeatureSelectionStrategy`, `NormalizationFeatureWorkflow`. |
| `integration.py` | Build and compare no-correction and correction candidates. | `IntegrationStrategy`, `NoCorrectionStrategy`, backend adapters, `IntegrationBenchmarkWorkflow`. |
| `graph.py` | Construct graphs, evaluate resolution/stability, and compute the selected embedding. | `NeighborGraphBuilder`, `ClusteringEvaluator`, `EmbeddingBuilder`, `GraphClusterEmbeddingWorkflow`. |
| `annotation.py` | Compute markers and join external evidence only after Step 08 is frozen. | `MarkerCharacterizer`, `ExternalAnnotationJoiner`, `AnnotationPublicationWorkflow`. |
| `validation.py` | Fail closed on input, schema, cell-identity, count, or decision-contract violations. | `ValidationRecorder`, `ManifestValidator`. |
| `publishing.py` | Atomically publish reports, manifests, exact code/config copies, and final objects. | `RunPublisher`, `DecisionLedgerPublisher`, `EnvironmentRecorder`. |

The strategy-object pattern makes methods tunable without conditional logic
scattered through stage drivers. Configuration validation must fail before
data loading if a required decision is absent, still `TBD`, internally
inconsistent, or prohibited by the workflow contract.

## Required decision gates

- Step 00 must pass before any matrix is loaded.
- Step 02 is diagnostic only. Step 03 filtering settings must be approved from
  its evidence and written to configuration before filtering runs.
- Artifact scores and artifact removal are separate decisions. A method call
  alone never authorizes removal.
- Step 06 no-correction results are mandatory. Step 07 may select an
  integration method only if a defensible nuisance/batch variable exists and
  declared acceptance criteria are met. Genotype, sex, or expected biology
  cannot be used as a batch variable.
- Step 08 graph, resolution, and UMAP settings must be fixed before Step 09
  attaches reference or PCDH19-classification labels. Annotation cannot feed
  backward into cleanup or geometry during this first primary pass.
- Each stage stops with a report and decision ledger. Downstream execution is
  separately submitted so diagnostic gates remain inspectable and tunable.

## Execution contract to implement later

Each stage will eventually have a matched set of files:

- `scripts/<NN>_<stage>.py`: thin stage driver using the shared objects;
- `bin/submit_<NN>_<stage>.sh`: dry-run, versioned-run, and guarded replacement
  interface;
- `slurm/<NN>_<stage>.sbatch`: resource allocation only, with no scientific
  settings;
- `config/<NN>_<stage>.lock.json`: executable decisions approved for that
  stage;
- `config/<NN>_<stage>.requirements.txt`: environment declaration.

The submitter will copy all shared package modules as well as the stage driver,
wrapper, SLURM file, configuration, requirements, and this handoff into the run
package. Jobs must execute the copies in `code/`, preserving the exact source
that ran.

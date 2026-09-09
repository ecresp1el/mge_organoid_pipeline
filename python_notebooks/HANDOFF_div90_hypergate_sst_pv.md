# DIV90 cortical Hypergate — completed first phase and next-chat handoff

Updated: 2026-09-09

Status: Both requested target-definition analyses are complete. The next phase has **not** been run. The user is opening a new chat to continue from these results.

## Read this first

The user wants to build on their recovered **cortical Loupe recluster**, especially the location of SST-expressing cells within it. Preserve the saved cell identities and local UMAP coordinates. The identity recovery work was essential: an existing extractor recovered annotation tracks and barcodes from the saved `.cloupe` file and joined those identities to expression by `cell_id`.

Do not restart cortical-versus-subpallial analysis, redo identity recovery, reconstruct trajectories, or block on the scarcity of PVALB. Do not treat SST-negative cells as established PV fate. The accepted first pass was deliberately crude SST detection versus non-detection, followed by a stricter module comparison. **Both of those passes already ran; “next phase” does not mean that the module pass is still pending.**

The current request was to update documentation and prepare a handoff, not execute additional analyses. Follow the user's new-chat direction for the next phase. The possible follow-up work below is a proposal, not an already accepted experimental design.

## Workspace and entry points

Repository:

```text
/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline
```

Completed outputs are in the external Turbo runtime workspace (`PROJECT_ROOT`), not the Git checkout:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_hypergate_sst_pv/
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_hypergate_sst_pv_results.zip
```

The following paths are relative to `PROJECT_ROOT`, not the repository. Start with:

- `results/div90_hypergate_sst_pv/REPORT.md`
- `results/div90_hypergate_sst_pv/decision_summary.json`
- `results/div90_hypergate_sst_pv/run_parameters.json`
- `results/div90_hypergate_sst_pv/provenance/validation.json`

The full expression matrix and original Loupe sources are on the mounted Turbo workspace, not in the local Git checkout.

## Frozen population and target definitions

Entry criteria: recovered `cortical_only` membership, LHX6 > 0 and ERBB4 > 0 in log1p(CP10K), excluding original stressed clusters `6` and `7` exactly as in the preceding audit.

There are **4,768 entry cells**. All retain their original cell IDs, local recluster coordinates, recovered subcluster labels and sample IDs.

| Definition | Positive class | Other class | Ambiguous |
|---|---:|---:|---:|
| Crude detection | 1,849 `non-SST/PV-candidate` (SST == 0) | 2,919 `SST-like` (SST > 0) | 0 |
| Module extremes | 638 `PV-like` | 615 `SST-like` | 3,515 |

Both target directions were fitted for each definition, with beta = 0.5, 1, 2.

Module scores are means of gene-wise z-scores within the 4,768 cells; population standard deviation (`ddof=0`) was used. High means >=65th percentile; low means <=35th percentile. A label requires its own module high and the opposing module low. Exact centering, scaling and cutoffs are in `run_parameters.json`.

- SST module: `SST, NR2F2, GRIK1`.
- PV-like module: `MEF2C, MAF, MAFB, KCNC1, KCNC2, GAD1, GAD2`.
- ERBB4 is omitted from module scoring because it is detected in every entry cell; it remains an allowed gate feature.
- All module-defining genes are excluded from gate features to avoid direct target-label leakage.
- These are the requested pragmatic programs, not validated fate assignments. Do not silently change their definitions when comparing against phase one.

## Completed search

Actual **R hypergate 0.8.5** was installed and used. This was not a Python approximation of Hypergate.

1. First, all one-, two-, and three-gene combinations of `ERBB4, CXCR4, ACKR3, PLXNA2, NRP1, NRP2` were run on the crude labels: 246 fits.
2. A curated CSPA human surfaceome was intersected with expressed genes. High-confidence transmembrane/GPI proteins were retained; target genes, obvious housekeeping genes and selected extracellular-matrix/isoform ambiguities were excluded. Detection >=5% and variation were required. Final allowed list: **314 genes**.
3. All allowed genes were screened univariately in both directions and at all three beta values. Six additional genes were selected per definition through the documented round-robin ranking; the six core genes remained in each panel.
4. All one-, two-, and three-gene combinations within each resulting **12-gene panel** were searched. This is not an exhaustive search of every triplet among all 314 genes.
5. Positive and negative threshold sweeps benchmarked ERBB4/CXCR4 and the other requested hypotheses.

Total: **3,576 Hypergate fits, 4,448 benchmark evaluations, 2,125 distinct gate rules**. All R confusion counts were independently reproduced by applying the exact saved thresholds in Python. The saved UMAP was never recomputed.

Expanded crude panel additions: `NEGR1, FAT3, ITGA9, EPHA5, PTPRM, PLXNA4`.

Expanded module panel additions: `NEGR1, FAT3, ITGA9, EFNB2, CDH4, CDH10`.

Surface annotation: [CSPA supplementary human surfaceome](https://doi.org/10.1371/journal.pone.0121314.s003). A downloaded copy, its checksum, annotation/filter audits and gene matches are saved in `provenance/` and `run_parameters.json`.

## Main results to preserve

Best shared two-marker candidate rule, inside the existing LHX6+/ERBB4+ population:

```text
FAT3 <= 0.764341 AND PTPRM <= 0
```

These are RNA thresholds in log1p(CP10K), not fluorescence thresholds. PTPRM <=0 means undetected RNA in this nonnegative matrix. The rule captures **3,168 entry cells**.

| Evaluation target | Purity | Yield | F1 | Fold enrichment |
|---|---:|---:|---:|---:|
| Crude non-SST/PV-candidate | 52.24% | 89.51% | 0.65976 | 1.347 |
| Module-defined PV-like | 71.36% | 89.81% | 0.79528 | 1.401 |

The shared rule was chosen by maximizing the lower F1 across the two definitions, then their mean F1. **Module purity is calculated only among the 1,253 unambiguous cells. It does not mean 71% of all 3,168 captured cells have established PV-like labels.** Ambiguous captures and total entry-cell capture counts are separate columns in the results.

Best two-marker rule by F1 for the module definition alone:

```text
FAT3 <= 0 AND ITGA9 <= 0
```

Purity 69.85%, yield 94.04%, F1 0.80160, enrichment 1.372x among the unambiguous module labels.

Best crude three-marker rule by F1:

```text
FAT3 <= 0.764341 AND ITGA9 <= 0 AND PTPRM <= 0
```

Rounded performance: purity 54.6%, yield 88.6%, F1 0.676. See the saved table for full precision. Marker counts in fitted rules are distinct genes; `surface_markers_including_entry` also counts the ERBB4 entry marker.

**ERBB4 + CXCR4 was outperformed by other pairs by F1.** Its best positive-positive threshold sweep gave F1 0.53928 for the crude definition and 0.58115 for the module definition. This does not mean it has no enrichment: the tradeoff is higher purity at lower yield, documented in `benchmark_best_by_hypothesis.tsv`.

## Limitations relevant to phase two

- The fits and gene screening used the full labeled populations. Per-sample evaluations applied the same global gate; they are **not held-out validation** and no per-sample gates were optimized.
- Shared-gate crude purity ranged from approximately 6% to 73% across samples, with different target prevalences. Inspect baseline fractions and sample counts rather than assuming uniform performance.
- Sample `10496-MW-5` has only 37 entry cells and 11 unambiguous module cells, with no module-defined PV-like targets. Its module recall is undefined, not zero.
- RNA-zero gates can reflect detection depth as well as biology. No depth-matched analysis or adjustment for total UMIs/detected genes has yet been done in this Hypergate run.
- A curated surface annotation does not establish protein abundance or live-cell sorting performance. **No human FACS antibody clone/catalog/epitope validation has been completed.**
- SST negativity is a first-pass target, not proven PV identity. The module labels are more selective operational labels, not lineage reconstruction.
- No additional phase-two computation, experiment, ordering, or message to collaborators has been performed.

## Files ready for immediate reuse

All paths in this table are relative to `PROJECT_ROOT/results/div90_hypergate_sst_pv/`.

| File | Purpose |
|---|---|
| `target_labels.tsv.gz` | Cell IDs, sample IDs, original/recovered identities, local UMAP coordinates, SST, both labels and module scores |
| `hypergate_input.tsv.gz` | Barcode-aligned expression of the allowed surface genes plus labels/metadata |
| `allowed_surface_markers.tsv` | Final allowed surface panel, detection fractions and annotation provenance |
| `hypergate_results_all.tsv` | Every fit and benchmark, thresholds, precision, yield, F1/F0.5/F2, enrichment and capture counts |
| `top_gates.tsv` | Top 20 distinct rules per definition and target direction, ranked by F1 |
| `per_sample_validation.tsv` | Fixed global gate performance per sample, with baseline target fractions |
| `shared_pair_comparison.tsv` | Exact same two-marker rules evaluated across both definitions |
| `benchmark_best_by_hypothesis.tsv` | Best F1 threshold sweep per requested hypothesis |
| `best_gate_cells.tsv.gz` | All entry cells with capture flags for the shared and definition-specific best pairs |
| `best_gate_rules.json` / `gate_rules.json` | Exact machine-readable gate rules |
| `surface_univariate_screen.tsv` | Feature-ranking results used to select expanded panels |
| `core_fits.tsv`, `expanded_*_fits_*.tsv` | Completed raw R fit checkpoints |
| `figures/` | Virtual FACS plots, capture overlays, SST/module maps, label maps, precision–yield plots, top 20 and per-sample plots |
| `provenance/validation.json` | Count and feature-leakage verification |

Upstream assets:

```text
# Recluster coordinates and guidance expression including SST (Turbo):
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/final_figures/div90_guidance_with_sst/tables/joined_with_sst.tsv.gz

# Full expression, counts, QC metadata and canonical embeddings (Turbo):
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/varela_div90.h5ad

# Recovered identity and sample mapping (Turbo):
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/final_figures/fig_div90_loupe_recluster_annotations_v1/tables/div90_loupe_recluster_membership.tsv.gz

# Original Loupe sources, exports and recovered tracks (Turbo):
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/data/raw/loupe_exports/div90_recluster_annotations/
```

The full H5AD contains 22,338 cells × 18,082 genes, a `counts` layer, original cluster/sample metadata and QC fields including `nCount_RNA` and `nFeature_RNA`. Only the frozen 4,768-cell cortical subset was analyzed here.

## Code and runtime

All active code is under `python_notebooks/scripts/` in the repository. Shared `div90_hypergate_paths.py` resolves `PROJECT_ROOT` from the environment (default: the Turbo workspace) and rejects runtime roots inside the Git checkout. Tables, plots, caches, logs, checkpoints and ZIPs must be written outside the checkout.

Scripts:

- `run_div90_hypergate_workflow.py`: full crude-first entry point from existing assets.
- `prepare_div90_hypergate.py`: surfaceome filtering, module labels and expanded job tables.
- `run_div90_hypergate.R`: actual R fits, exact rules and confusion counts.
- `run_div90_hypergate_expanded.py`: four R workers per definition; expanded crude precedes module fitting.
- `summarize_div90_hypergate.py`: independent validation, benchmark sweeps, reports, figures and ZIP.

Environment:

```text
/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python
/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/Rscript
R 4.5.3; hypergate 0.8.5
```

Full rerun command, only if needed:

```bash
/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python python_notebooks/scripts/run_div90_hypergate_workflow.py
```

Completed fit checkpoints are reused; incomplete ones cause an explicit failure. Do not rerun phase one simply to resume the conversation. For changed labels, panels, or evaluation design, use a separate output directory and retain the completed phase-one results.

Worktree state: code and documentation changes remain uncommitted. Runtime artifacts created during this conversation have been moved out of the checkout to Turbo, with per-file SHA256 verification before removal of the checkout copies. See `PROJECT_ROOT/results/div90_hypergate_sst_pv/provenance/runtime_migration.json`. Historical logs retain their original paths as provenance; those paths are not current output locations. A pre-existing modification to `render_div90_loupe_pooled_subpopulation_violins.py` adds the SST-compatible violin range option. Preserve these changes; no commit or push was requested.

## Proposed next phase — not executed

Subject to the user's next-chat direction, useful follow-up is to evaluate whether the discovered surface phenotype is robust and experimentally actionable:

1. Begin by reviewing the saved gate and label overlays on the recovered cortical recluster, including both target directions and the purity–yield tradeoff.
2. Evaluate RNA detection-depth and sample-composition sensitivity using the existing cell IDs and H5AD QC/count data. Preserve the original labels and thresholds as a reference.
3. If predictive validation is requested, learn feature selection, score calibration where applicable, and thresholds within training partitions; evaluate frozen rules on held-out samples. Do not call the existing per-sample table cross-validation.
4. Check human extracellular/live-cell flow antibody evidence for shortlisted markers and compare a reagent-constrained panel with the exploratory surfaceome panel. Keep catalog/clone/application evidence separate from RNA performance.
5. If revising module genes or quantile cutoffs, version the labels and report sensitivity against both completed definitions. Do not replace the crude SST first pass or require PVALB positivity.

Suggested separate output directory if the user starts follow-up computation:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_hypergate_sst_pv_phase2/
```

## Trigger phrase for phase two

**Resume DIV90 Hypergate phase 2 — Turbo runtime only.**

This phrase means: read this handoff, use the completed phase-one results on Turbo, preserve the recovered cortical identities and UMAP, and begin the follow-up described below or the user’s more specific new-chat instructions. It does not mean rerun the completed crude or module passes.

Default follow-up: first assess sample/detection-depth robustness of the existing gates, then evaluate human extracellular FACS-antibody feasibility for the shortlisted markers. Keep all new runtime files under the separate phase-two Turbo directory. Do not reconstruct trajectories.

## Prompt for the new chat

> Resume DIV90 Hypergate phase 2 — Turbo runtime only. Read `python_notebooks/HANDOFF_div90_hypergate_sst_pv.md` and the linked completed results. Both the crude SST+/SST− Hypergate pass and the stricter module-label pass are already complete. Preserve the recovered cortical cell identities, original local recluster UMAP, and phase-one outputs. Do not redo trajectories, cortical/subpallial comparisons, or identity recovery. Start with robustness and antibody-feasibility follow-up, subject to my additional instructions. Write runtime artifacts only under `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_hypergate_sst_pv_phase2/`; retain only code/configuration/documentation in the Git checkout. The proposed follow-up has not already been run.

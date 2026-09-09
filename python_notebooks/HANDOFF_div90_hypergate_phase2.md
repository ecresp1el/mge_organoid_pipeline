# DIV90 Hypergate phase 2: independent programs and culture-condition analysis

Updated 2026-09-09. Continue from the phase-two artifacts below; do not repeat phase one or replace its scores, labels, identities or coordinates.

## Runtime and entry points

All scientific runtime outputs, caches, logs and interactive assets are under:

`/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_hypergate_sst_pv_phase2/`

Open `index.html`, `REPORT.md`, `CONDITION_REPORT.md`, and `interactive_pv_sst_landscape.html` first. Four standalone HTML assets work offline. The figures directory contains ten core figures plus supplements in PDF/SVG/600-dpi PNG. `figure_manifest.json` and `figures/FIGURE_CAPTIONS.md` describe them.

The completed phase-one directory remains read-only input. Its SHA256 preservation audit is `provenance/phase1_preservation.json`. All 4,768 cell IDs, recovered annotations, original local Loupe UMAP coordinates, and original phase-one labels/scores are retained. No trajectory, UMAP, identity recovery, or cortical/subpallial analysis was repeated.

## Primary interpretation

PV and SST are independent continuous axes. The primary PV module is `MEF2C, KCNC1, KCNC2, TAC1`; the SST module is `SST, SATB1, NR2F2, CDK14, CACNG3`. Scores are equal-weight means of population gene-wise z-scores of existing log1p(CP10K), bounded per gene to [-3,3]. Generic GAD genes and ERBB4 do not define the primary scores. See [the evidence note](DIV90_HYPERGATE_PHASE2_EVIDENCE.md) for each gene's rationale and limitations.

Pooled median thresholds define operational groups, not fates or demonstrated natural clusters. All cells remain assigned:

| State | N |
|---|---:|
| PV-biased | 1,076 |
| PV/SST hybrid | 1,309 |
| SST-biased | 1,075 |
| Unresolved/immature | 1,308 |

The 1,309 dual-high cells are not all validated biological hybrids. Only 356 detect at least two genes in each module; 157 detect at least two PV genes and two non-SST genes from the SST module. Total dual-high abundance does not exceed the independent-gene null preserving sample and UMI quintile. Within-module correlations are weak; RNA complexity contributes. SST supplies about 25.1% of SST-score variance, rather than dominating it alone. Gene-drop and threshold sensitivities are saved for every cell. Quantiles at the all-undetected score floor assign the entire tie to low, never split cells arbitrarily.

Per-sample Scrublet used all 22,338 available DIV90 cells, with an assumed 5% expected rate and automated thresholds. Five entry cells were flagged, including three dual-high cells. These are exploratory scores, not validated doublet truth; homotypic doublets and prior filtering remain limitations. No primary cells were excluded.

## Culture metadata

The user supplied CV/MW assignments and an approximately twofold glucose difference during this conversation. The new configuration [culture_condition_user_20260909.tsv](../metadata/culture_condition_user_20260909.tsv) supplements the original sample/line mapping without changing its labels. All sample IDs were checked against the existing metadata. The `MW` substring in every run ID is not the culture operator.

- DIV90: CV odd samples 1/3/5, approximately 2x glucose; MW even samples 2/4/6, approximately 1x glucose.
- DIV30 metadata is preserved for future use: CV 1/3/5/7/8; MW 2/4/6/9. No DIV30 analysis was run here.
- H9 comprises DIV90 samples 1/2, 79B 3/4, 2E 5/6. There is one sample per line and condition. CV/MW is confounded with operator and other culture differences; no glucose-specific causal claim is supported.
- CV contributes 1,581 entry cells and MW 3,187. H9 and 79B shift toward higher PV/lower SST scores in CV; 2E reverses that pattern with only 37 CV cells. Inspect line-specific results.
- FAT3/PTPRM expression is lower in CV within the observed states as well as in pooled data. Fixed-gate performance therefore requires condition-specific interpretation.

## Gate conclusions

The old `FAT3 <= 0.764341 AND PTPRM <= 0` retention gate still captures exactly 3,168 cells. Under the new operational model it retains 79.09% of PV-biased but only 47.90% of hybrid cells, losing 682 hybrids. It is not a hybrid-preserving strategy.

The recovery-weighted unconstrained optimum makes almost no depletion (`SLC6A6 > 2.9582818`, removes 22 cells). This is an informative limit of the data, not a successful purification.

A separately labeled analyst-selected constraint requires at least 20% SST-biased removal, with 10/20/30/50% sensitivity. Under it the tested practical pair is:

`REMOVE FGFR2 > 1.2491154 OR PTPRS > 3.4076235; RETAIN the complement.`

Use full-precision rules in `gate_summary.json` for evaluation, not rounded display text. The rule removes 744 and retains 4,024 cells:

| Quantity | Result |
|---|---:|
| PV-biased retained | 954 / 1,076 (88.66%) |
| Hybrid retained | 1,169 / 1,309 (89.30%) |
| PV-biased lost | 122 |
| Hybrid lost | 140 |
| SST-biased removal | 219 / 1,075 (20.37%) |
| Retained SST-biased contamination | 856 / 4,024 (21.27%) |
| Starting SST-biased contamination | 22.55% |

This is modest depletion with substantial costs. It is an RNA-derived experimental hypothesis. Human surface-flow evidence exists for FGFR2 clone 98739 and PTPRS AF3430, but neither the antibody behavior nor a fluorescence threshold is validated in DIV90 neurons. ERBB4-only live entry also does not automatically reproduce the recovered cortical/LHX6+ RNA entry criteria. The report proposes a protein pilot and matched maturation comparison of removed, retained and unsorted fractions.

## Computation and verification

All 314 previously curated expressed surface features were screened after target-gene exclusion. There were 471 actual R Hypergate fits, distinct from the Python depletion sweeps. All R confusion counts were independently reproduced.

The full Python search evaluates 2,885,531 rule boundaries, including 305,661 exact one-marker thresholds, 2,457,050 broad positive AND/OR pairs, constrained/unconstrained dense refinements, and bounded three-marker extensions. These are evaluation counts and can include repeated boundaries. Pair results are grid-based; triplets are not exhaustive. Complete metrics and rule definitions occupy two HDF5 files in `tables/`. Independent verification checks count conservation for every candidate and directly reproduces 700 sampled masks, plus all material gate counts and assignments.

No held-out predictive validation was claimed. Per-sample, cell-line and CV/MW comparisons apply fixed global gates. Unresolved exclusion is an evaluation sensitivity, not a demonstrated additional physical sorting gate.

Important files: `cells.tsv.gz`, `module_genes.tsv`, `state_parameters.json`, `hybrid_validation.json`, `gate_summary.json`, `gate_cell_assignments.tsv.gz`, `tables/gate_module_threshold_sensitivity.tsv`, `tables/gate_depth_quintile_sensitivity.tsv`, `tables/condition_*`, `provenance/final_validation.json`, and `provenance/interactive_validation.json`.

## Code and resumption

Main code is [run_div90_hypergate_phase2.py](scripts/run_div90_hypergate_phase2.py), with separate `div90_hypergate_phase2_{gates,doublets,conditions,visuals,validate,report}.py` modules. Only code, configuration and documentation belong in the Git checkout. Runtime dependencies installed for exploratory doublet scoring and browser verification reside under the phase-two Turbo directory; they do not modify the shared conda environment.

Use the existing environment with `PYTHONDONTWRITEBYTECODE=1`; configure Matplotlib/Numba caches under phase-two Turbo. The main script supports `--stage prepare|doublets|qc|gates|conditions|figures|validate|report|verify|all`. Completed gate search checkpoints are reused. For a changed scientific definition, use a newly versioned output instead of silently replacing this run. Do not rerun the whole analysis merely to inspect or continue this handoff.

No experiment, antibody order, external publication, or collaborator message was performed.

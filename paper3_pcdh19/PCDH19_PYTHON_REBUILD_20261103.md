# PCDH19 mouse MGE: current Python workflow and Methods reference

> **Storage/job update, 4 November:** See [scratch migration and current job pointer](SCRATCH_STORAGE_HANDOFF.md). Migration65352824 → scratch array65352825 (max20) → consensus65352826. Prior job IDs below are historical. Transfer is submitted; verify its completion seal.

**Latest concurrency/status snapshot: 2026-11-04 08:39 EST.** Earlier stage descriptions below retain their explicitly dated historical context. This is the current ordered workflow through
MapMyCells reference preparation. It replaces the chronological notes in
[the execution history](PCDH19_REBUILD_EXECUTION_HISTORY_20261103.md).
Do not use the older 446,349-cell / 34-cluster analysis as the outcome of this
rerun. The frozen approved main figure remains an earlier analysis.

## 1. Current milestones

| Stage | Verified status at this snapshot | Evidence |
| --- | --- | --- |
| Canonical counts and original QC | Completed; source checkpoints preserved | Original Step00/Step02 checkpoints and exclusion tables |
| Additional minimum-gene filter | Completed; 445,254 cells retained | `cell_cycle_comparison_min500/SUCCESS.json` |
| Matched before/after regression PCA/UMAP | Completed on the same 445,254 cells | `cell_cycle_comparison_min500/parameters.json` and embeddings |
| Revised QC/cell-cycle supplement | Completed, PNG/PDF/SVG | `supplementary_QC_cell_cycle/VALIDATION.json` |
| HiCAT normalized input and global cycle model | Completed | `hicat_min500_cycle/PREPARED.json` |
| HiCAT recursive fits / 100-iteration consensus | 32/100 fits completed; 20 running, 48 pending; full initializer completed; consensus pending | `hicat_min500_cycle/JOBS.json`, per-fit progress and completion files |
| Final new cluster count | Not available; requires completed consensus and separation review | Future `hicat_min500_cycle/final/SUCCESS.json` |
| MapMyCells developmental reference | Method selected; metadata/taxonomy validation completed, final reference selection and marker fitting pending | `mapmycells_developmental_reference/REFERENCE_VALIDATION.json` |
| Mapping the new query / URD | Not submitted / not started | No new mapping or URD results are claimed |

The first HiCAT fit completed: 33 initial root graph communities → 11 after
root DE merging → 71 recursive leaves → **67 clusters after cross-branch
merging**. Its pairwise audit flagged 65 unresolved separation pairs among
2,211 evaluated pairs. These are first-fit diagnostics, not a finalized taxonomy.
**67 is one subsample fit's cluster count, not the final consensus count.**
Four further fits and the separate full-population initializer are running.

## 2. Study, inputs and sample identity

This dataset comprises **12 dissected E14.5 mouse medial ganglionic eminence
(MGE) samples**, not organoids. The input was the Cell Ranger **filtered-barcode**
gene-expression delivery: **450,788 cells × 19,071 genes**. This is not a count
of every sequenced droplet. The Cell Ranger/software/reference details should
be taken from the original delivery manifest; they are not inferred here.

[Registered sample key](config/sample_key.csv): JZ-1–3 = WT Male; JZ-4–6 = WT
Female; JZ-7–9 = HET Female; JZ-10–12 = KO Male. Stable cell IDs and
`technical_sample_id` are preserved. Three technical sample IDs per group do
not, by themselves, establish independent animal/litter structure.

Runtime root (`ROOT`):
`/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19`.
Current rerun (`RUN`): `ROOT/results/rebuild_20261103`.

The authoritative original QC count checkpoint is:

```text
ROOT/results/primary_processing/02_qc_filtering/
  02_qc_filtering_20260830_124611_97e1bb5/
  objects/pcdh19_step02_qc_filtered.h5ad
```

It contains sparse original counts for 446,349 cells × 19,071 genes.
SHA256: `fadba4a25a7b6b7320219b21c189b6325687493519ba0fe1bd27efd79606b103`.
The 500-gene filter selects rows from this checkpoint; it does not overwrite
its counts. The raw Step00 checkpoint also remains preserved.

## 3. Cell QC and retention

Original QC was computed separately within each technical sample. Define
scaled MAD as `1.4826 × median(abs(x − median(x)))`. A cell was excluded if
any of the following applied:

- `log1p(total UMI count)` below its sample median minus 5 scaled MADs;
- `log1p(number of detected genes)` below its sample median minus 5 scaled MADs;
- mitochondrial UMI percentage above its sample median plus 5 scaled MADs.

Low-count/gene boundaries are back-transformed into original units. The union
removed **4,439 cells**: 0 low-UMI, 67 low-gene, 4,372 high-mitochondrial; these
reasons did not overlap in the approved result. **446,349 cells** remained.

The user then requested a fixed minimum detected-gene filter. We retained
cells with **at least 500 detected genes** (`n_genes_by_counts >= 500`).
This removed **1,095 additional cells**, leaving **445,254**. Total exclusions
relative to the delivered filtered-barcode cells are **5,534 (1.2276%)**;
retention is **98.7724%**. All **19,071 genes** and all 12 samples remain.
The threshold is an explicit QC choice, not an empirically optimized cutoff.
A minimum of 200 genes would exclude no additional approved Step02 cells.

| Sample | Group | Input cells | Additional <500-gene exclusions after original QC | Final retained cells |
| --- | --- | ---: | ---: | ---: |
| JZ-1 | WT Male | 51,229 | 122 | 50,815 |
| JZ-2 | WT Male | 37,553 | 156 | 36,820 |
| JZ-3 | WT Male | 25,354 | 97 | 24,137 |
| JZ-4 | WT Female | 21,440 | 55 | 21,273 |
| JZ-5 | WT Female | 41,878 | 91 | 41,498 |
| JZ-6 | WT Female | 56,099 | 98 | 55,659 |
| JZ-7 | HET Female | 19,623 | 88 | 19,410 |
| JZ-8 | HET Female | 20,799 | 66 | 20,611 |
| JZ-9 | HET Female | 60,680 | 83 | 60,245 |
| JZ-10 | KO Male | 52,753 | 97 | 52,204 |
| JZ-11 | KO Male | 52,295 | 96 | 51,591 |
| JZ-12 | KO Male | 11,085 | 46 | 10,991 |

No upper UMI/gene cutoff, doublet-based exclusion, ambient RNA correction,
batch integration or new broad cell-type exclusion was added. The previously
rejected scDblFinder calls are not used. Cell-cycle regression removes no cells.

Exact code and outputs:

- [Original MAD calculations](scripts/primary_processing/step01a_outliers.py) and [approved QC filtering](scripts/primary_processing/step02_filtering.py).
- [New filter and comparison, `run()`](scripts/compare_pcdh19_cell_cycle_umaps.py).
- `RUN/supplementary_QC_cell_cycle/per_cell_QC_disposition.tsv.gz`: original flags plus minimum-gene flag and final exclusion decision.
- `RUN/supplementary_QC_cell_cycle/sample_cell_retention.tsv`: exact per-sample totals.

## 4. Data representations: what enters each calculation

| Representation | How generated | Used for |
| --- | --- | --- |
| Original counts | Measured UMI counts, retained unchanged | Query input for future MapMyCells; appropriate future count-based expression analyses |
| Normalized, unregressed expression for UMAP | `ln(1 + counts / cell_total × 10,000)` | Cell-cycle scoring and starting expression for the paired UMAP analysis |
| Corrected UMAP representation | OLS residuals for intercept + S/G2M scores, then SD scaling / upper clip 10 | Corrected PCA, neighbor graph and UMAP |
| Normalized, unregressed HiCAT expression | `ln(1 + counts / cell_total × 1,000,000)` | HiCAT local gene eligibility/HVG selection and expression/detection tests |
| Corrected HiCAT representation | Global S/G2M residuals of HiCAT-normalized expression, SD scaling / upper clip 10 | HiCAT local PCA/graphs and held-out marker-prototype mapping |

**HiCAT does not cluster the UMAP coordinates.** Its corrected high-dimensional
expression representation determines PCA/graph clustering. The corrected UMAP
is the display on which the eventual HiCAT labels will be overlaid.

The UMAP and HiCAT pipelines deliberately retain their respective earlier
normalization targets (10,000 and 1,000,000); do not describe them as an identical
matrix or numerically identical PCA space. Both adjust only S/G2M scores.
Regression/scaling is applied as an analysis representation and never replaces
measured counts. No genotype, sex, sample, total-count or mitochondrial covariate
is regressed.

## 5. Matched before/after regression and UMAP

Both arms are refitted on the identical **445,254 post-QC cells** using the
same **3,000 HVGs**, newly selected from their raw counts. HVG selection uses
Scanpy `seurat_v3`, `n_top_genes=3000`, `batch_key='technical_sample_id'`.
Sample-aware HVG selection is not batch integration.

Cells were normalized to 10,000 total counts and natural-log `log1p`
transformed. Cell-cycle scoring uses Scanpy `score_genes_cell_cycle`, with
**43 matched S-phase and 53 matched G2/M genes**, resolved from the preserved
reference gene lists and aliases. One G2/M reference symbol was unmatched.
The original unregressed scores/phase labels are retained for both embeddings.

| Control | Value |
| --- | --- |
| Regression | Ordinary least squares with intercept, S_score, G2M_score |
| Corrected genes | Selected 3,000 HVGs for the UMAP comparison |
| Regression implementation | Blocked QR projection; block size 64 genes; float32 working residual matrix |
| Scaling in both arms | Divide by sample SD (`zero_center=False`), upper clipping at +10 |
| PCA | Centered ARPACK, 50 PCs |
| Neighbors | 30, cosine distance, explicit `use_rep='X_pca'`, 50 PCs |
| UMAP | `min_dist=0.3`, `spread=1.0` |
| Random seed | 20260903 for scoring, PCA, neighbors and UMAP |

The uncorrected arm scales normalized expression without regression. The
corrected arm regresses and then scales expression. Both arms recompute PCA,
neighbors and UMAP. Layouts can rotate independently; coordinate position
itself does not identify a corresponding cell state.

Core code, abridged from the executed script:

```python
keep = a.obs['n_genes_by_counts'].to_numpy() >= 500
filtered = a[keep].copy()
sc.pp.highly_variable_genes(filtered, flavor='seurat_v3',
                          n_top_genes=3000, batch_key='technical_sample_id')
sc.pp.normalize_total(filtered, target_sum=10000)
sc.pp.log1p(filtered)
sc.tl.score_genes_cell_cycle(filtered, s_genes=s_genes,
                           g2m_genes=g2m_genes, random_state=20260903)
# On the HVG matrix: residual = Y - Q @ (Q.T @ Y),
# where Q spans [intercept, S_score, G2M_score].
```

See [regression and score-gene resolution](scripts/compare_pcdh19_cell_cycle_umaps.py)
and [the numerical tests](tests/test_cell_cycle_comparison.py). Gene lists and
aliases are in [the program configuration](config/hicat_validation_programs.json).
The exact matched/missing audit is saved in `cell_cycle_gene_coverage.json`.

Completed output: `RUN/cell_cycle_comparison_min500/`:

- `before_after_embeddings.npz`: exact cell IDs and both PCA/UMAP coordinate sets;
- `cells_and_cycle_scores.tsv.gz`: sample metadata, original S/G2M scores and phase;
- `selected_HVGs.tsv`, `regression_pca_model.npz`, corrected graph files;
- sample overlays, individual-sample highlights, cycle and lineage-marker overlays;
- `parameters.json`, `SUCCESS.json`, `progress.jsonl`, source/config copies.

Maximum absolute correlations of S/G2M scores with any of the first 50 PCs:
S **0.383973 → 0.028507**; G2/M **0.803819 → 0.038647**. These diagnose the
intended reduction of score-associated variation, not complete removal of
biological cell-cycle effects.

The supplement is `RUN/supplementary_QC_cell_cycle/PCDH19_QC_cell_cycle_supplement.{png,pdf,svg}`.
It combines sample-level pre/post-QC counts/distributions with matched
before/after-regression embeddings. Both UMAPs are **post-QC**; no pre-QC
UMAP is implied. Fits use all retained cells; plots display at most 6,000
cells per sample with fixed sampling and matching colors.

## 6. HiCAT: modified rerun using the previous consensus settings

**Clustering is running: 32/100 fits completed at 2026-11-04 08:39 EST; no final consensus cluster count is available.** This uses the
previous pinned Allen Python fitter and consensus implementation, with a new
input cell universe and explicit S/G2M-adjusted clustering representation.
It is not an unchanged rerun of the old 34-cluster / 98-iteration result.

Frozen configuration: `RUN/hicat_min500_cycle/config.json`. Prior full-data
settings were retained from `RUN/baseline/previous_hicat_run.json`.
The complete historical algorithm trace is in
[the previous HiCAT methods guide](HICAT_PARAMETERS_AND_CONSENSUS.md); the
changes described in this document govern the current rerun.

### 6.1 Expression preparation and cycle correction

Read original Step02 counts, select the same `>=500`-gene cells and verify
exact ID ordering against the completed UMAP/scoring output. Create one shared
sparse float64 natural-log CPM object: `full_log1p_cpm.h5ad`.

Fit an all-gene global OLS model to that expression with design
`[intercept, S_score, G2M_score]`. The scores are those already computed from
unregressed CP10K expression. Save regression coefficients and sample residual
SDs in `cycle_model.npz`; the model is fitted on the full retained population
before drawing 80% fitting subsets. Dense all-gene residuals are not stored.
Each node materializes only its selected corrected genes.

```python
# Actual core algebra in fit_cycle_model() / corrected_block():
beta = np.linalg.solve(design.T @ design, np.asarray(X.T @ design).T)
residual = X_selected - design @ beta_selected
corrected = np.minimum(residual / residual_sample_SD, 10)
```

Source: [cycle preparation, representation and workflow](scripts/hicat/cycle_rebuild.py)
(`fit_cycle_model`, `corrected_block`, `CycleEngine._representation`, `prepare`).
Source counts are never replaced by residuals.

### 6.2 Per-fit parameters

| Scientific control | Actual setting / behavior |
| --- | --- |
| Input | 445,254 cells × 19,071 genes; all 12 samples |
| Subsample iterations | 100 independent fits |
| Sampling | Uniform pooled 80% without replacement; 356,203 fitted / 89,051 held-out cells per iteration |
| Normalization | Natural log(1 + CPM), target 1,000,000 |
| Gene eligibility | Unregressed ln(1+CPM) >1 in at least 4 node cells |
| Local feature selection | Original Allen HVG rule, up to 3,000; recomputed at each recursive node |
| PCA representation | S/G2M-regressed, SD-scaled and upper-clipped expression of selected genes |
| PCA fit | All cells at each node; randomized solver; fit up to 50 PCs subject to node dimensions |
| Retained PCs | Allen explained-variance elbow selection; at most 20 |
| Graph neighbors | Requested k=15; actual k=min(15, retained PC dimensions, node cells−1) |
| Distance / approximate search | Euclidean; Annoy, 50 trees |
| Edge weighting / clustering | Jaccard; taynaud Louvain; resolution=1.0 |
| Recursion | Recompute features/PCA/graph/merging in each branch; branches with <40 cells terminate |
| Depth guard | 20; reaching it is a failure, not convergence |
| Cross-branch merge | Enabled, using corrected PCA on the fit-derived marker union |
| Target cluster count | None |

Local HVG eligibility/ranking still uses unregressed expression, as in the
prior fitter. This can influence which genes are offered to the corrected PCA.
The implementation does not claim complete statistical independence from
cell-cycle biology in every stage.

### 6.3 HiCAT expression-based merging and its interpretation

Candidate graph groups undergo the original Allen empirical-Bayes expression
merging. Neighbor-pair prioritization uses the corrected representation;
expression means, variances and detection fractions use **unregressed ln(1+CPM)**.
The same distinction applies to final DE merging and audit. This is deliberate:
the existing detection/expression thresholds retain their original meaning.
Do not describe the DE tests as cell-cycle-adjusted tests.

| Gene / pair control | Actual rule |
| --- | --- |
| Expression detection | ln(1+CPM) >1 |
| Foreground detection q1 | >0.4 |
| Background detection q2 | None; no independent ceiling |
| qdiff | abs(qA−qB)/max(qA,qB) >0.7; 0 if both fractions are 0 |
| Foreground detected-cell minimum | 20 |
| Adjusted p | <0.05, Holm adjustment in the pinned Python implementation |
| Expression difference | Absolute difference of mean ln(1+CPM) >1; not log2 fold change |
| DE score | Sum of −log10(adjusted p) for qualifying up/down genes; threshold 150 |
| Declared endpoint audit | Score ≥150 and at least 5 qualifying genes |
| Merge-neighbor count | 2 |

The pinned merge loop has a score-only early exit and an inner gene-count
comparison `num > 5`; the endpoint audit uses `num >= 5`. Neighbor merging
also does not compare every possible pair at every pass. Surviving weak pairs
are reported for review. The algorithm does not guarantee every final pair
passes the declared separation rule.

These tests concern cluster separation, not sample-replicated PCDH19 genotype
DE. A final cluster count remains provisional until its separation evidence,
sample composition and biological marker context have been reviewed.

### 6.4 Held-out assignment and consensus

For each 80% fit, compute each sampled cluster's mean prototype in **corrected,
scaled marker expression**. Assign each held-out cell to the highest Pearson
correlation prototype; preserve sampled-cell labels. The markers are the union
of the first 20 entries in each up/down separation-gene list captured from the
fit DE calls, not a hand-selected cell-identity marker panel.

There is no held-out rejection threshold. Undefined correlations become zero;
exact ties select the first sorted numeric cluster ID. One-cluster fits assign
the only possible label without inventing a perfect correlation score.
Every iteration therefore assigns all 445,254 cells once.

The separate full-population initializer uses the same corrected recursive
fitter but is not counted among the 100 subsample partitions. Consensus uses
all 100 completed partitions; failures are not silently omitted.

| Consensus control | Actual value |
| --- | --- |
| Pairwise coassignment | P(i,j)=number of common-cluster assignments /100 |
| Storage | Sparse iteration-cluster × cell membership B; never a dense N×N matrix |
| Large-graph branch gate | First-fit sum(cluster_size²) ≥1e9; smaller branch stops pending qualification |
| Initial partition | Separate full-population HiCAT fit |
| Coassignment merge | Between affinity >0.1 and max(within affinities)−between <0.25 |
| Refinement | Max-affinity reassignment; median confusion threshold 0.6; minimum group size 20 |
| Refinement stopping | Up to 50 iterations per outer pass, <1% proposed changes checked before update, or no increase in correct count |
| Final stage | Corrected-marker-PCA-prioritized expression DE merge, then all-pairs audit |

The current final DE prioritization uses corrected marker PCA. The older
production adapter used marker expression directly at that final stage; this
is another explicitly recorded consequence of applying cycle correction to
the current clustering representation.

Implementation: [current workflow](scripts/hicat/cycle_rebuild.py),
[recursive engine](scripts/hicat/fullscale_engine.py),
[feature selection / graph / merging](scripts/hicat/engine.py), and
[consensus operations](scripts/hicat/consensus_core.py).
The **frozen runtime copies**, not later edits to repository source, are
executed by the submitted jobs.

### 6.5 Seeds, resources, checks and storage

- Iteration 0 sampling/fitting seed: 20260908.
- Iterations 1–99: `np.random.SeedSequence(20260908).spawn(99)`, generating two uint32 values per child (sampling then fitting). Exact values: `hicat_min500_cycle/seeds.tsv`.
- Full-population initializer seed: 20260909. Final DE seed: 20260908.
- Preserve the original positive seeds with the qualified Annoy uint64 bridge; no modulo/sign wrapping. All 100 seeds passed native index/KNN/graph qualification.
- HiCAT environment: Python 3.8, Scanpy 1.9.8, AnnData 0.8.0, NumPy 1.24.4, SciPy 1.10.1, pandas 1.5.3, scikit-learn 1.3.2, Annoy 1.17.3, python-louvain 0.16. Exact `pip freeze`: `environment.txt`.
- Slurm: 8 CPUs / 150 GiB / 48 h per fit; Allen workers=1; at most four array fits concurrently, plus the separate initializer when scheduled.
- Sparse computation / mapping blocks: 4,096 cells. This is a memory bound, not biological downsampling.
- Large recursive graphs/projections/indexes are job-local temporary files. Retain labels, models, membership matrices and compressed feature/PCA/DE evidence archives on NFS.
- Completed checks: OLS/SD against dense least squares; block held-out mapping against the earlier Pearson implementation applied to explicit residuals; real-count reader/QC agreement; actual recursive, cross-branch and final DE smoke tests; exact-seed Annoy qualification. Evidence: `hicat_min500_cycle/qualification/`.

The input H5AD from newer AnnData has `uns` encodings unsupported by AnnData
0.8. The current HiCAT reader uses HDF5 to read only counts and required
metadata, then builds a compatible in-memory object; it does not reinterpret
those unsupported annotations. Raw counts and QC metrics were checked for
agreement on a 100-cell real input slice.

HiCAT source hashes/configuration are frozen in `FROZEN.json`; pinned upstream
Python source hashes are in `REFERENCE.json`. Allen Python commit:
`99154957c74023235763025fda9dc4eb3ca952c6`. The inherited Allen R consensus
control-flow reference commit is `9af2f04cb837a7b53b005dee88a92eae337afead`.
This is a Python adapter, not a claim of numerical identity to unchanged R.

## 7. MapMyCells: selected method, reference preparation not yet mapping

**The user selected the actual MapMyCells approach.** Use the local Allen
`cell_type_mapper` engine, with a custom developmentally appropriate reference.
Do not substitute the earlier Seurat/MIND transfer or a generic nearest-neighbor
label transfer and call it MapMyCells. The official implementation supports
[arbitrary reference taxonomies](https://github.com/AllenInstitute/cell_type_mapper/blob/main/docs/mapping_cells.md).

The local engine is installed at
`ROOT/software/mapmycells-v1.7.4/bin/python`; package version **1.7.4**.
It includes reference-statistic precomputation, marker selection and
hierarchical query assignment. The new custom reference has **not** yet had
its statistics/markers fitted, and the new query has **not** been mapped.

### 7.1 Verified available references and limitations

| Candidate | Verified evidence | Current limitation |
| --- | --- | --- |
| La Manno 2021 developing mouse nervous system | Cached 292,495-cell × 31,053-gene H5AD; X sample contains integer counts; published Class/Subclass/CellType and developmental ages, including E14.5 | Exact MGE membership has not been demonstrated; cannot relabel all ventral forebrain cells as MGE |
| Bandler CA301 | 4,481 definitively barcode-joined WT E15.5 MGE cells; later-MIND two-class / 12-state labels | These are later-MIND labels, not the unrecovered original 2022 taxonomy; one leaf has only one reference cell; incomplete whole-dissection/non-neural coverage |
| Allen cached WMB-10X | Official adult whole-brain reference, version 20240831 | Adult comparator, not selected as definitive embryonic MGE reference |
| Allen developing visual cortex | Official MapMyCells developmental taxonomy, CCN20260131, E11.5–P56 | Visual-cortex coverage does not establish embryonic MGE progenitor coverage |

La Manno and Bandler label hierarchies were tested with the **actual installed
MapMyCells `TaxonomyTree` API**. La Manno repeats generic labels such as
`Undefined` beneath different parents. Its original strings cannot directly
serve as globally unique tree-node IDs. Path-qualified JSON machine IDs now
represent the existing published parent path, while preserving original labels
verbatim in separate columns; this creates no new biological assignments.
The path-qualified tree and Bandler tree both pass construction.

Code: [developmental reference validator](scripts/reference_mapping/validate_developmental_mapmycells_reference.py).
Output: `RUN/mapmycells_developmental_reference/`, including
`REFERENCE_VALIDATION.json`, age/class composition, reference-state counts,
original labels / machine IDs and validated taxonomy JSON files.

### 7.2 Remaining mapping work and input contract

1. Finalize age/region coverage, hierarchy support and handling of reference states with too few cells. Preserve broad non-neural/progenitor classes instead of forcing every query cell into an MGE neuronal state.
2. Use MapMyCells to compute the selected reference's summary statistics and reference/query marker assets. Preserve exact gene IDs, normalization settings and provenance.
3. Map **original raw counts for the same 445,254 cells**. Do not use negative regression residuals or scaled expression as MapMyCells input. Query gene IDs are mouse Ensembl IDs.
4. Preserve every returned hierarchy level, assignment score, correlation, bootstrap support and runner-up field; join by stable cell ID to the new HiCAT labels and corrected UMAP.
5. Review reference support, ambiguous calls and marker consistency before using annotations in URD.

MapMyCells assigns individual cells and does not require the final HiCAT
cluster number. Its labels can be joined to clusters afterward. No new
bootstrap/marker parameters or acceptance cutoffs are reported as executed.
The previous adult comparator used 100 bootstraps, factor0.5 and seed20260829;
those are historical settings, **not yet a frozen custom-reference run**.

The existing [MapMyCells query runner](scripts/reference_mapping/run_mapmycells_comparator.py)
is a code reference for preserving hierarchy/confidence output. It is not
proof that the current custom-reference query mapping has run.

Official resources consulted:

- [MapMyCells custom-taxonomy workflow](https://github.com/AllenInstitute/cell_type_mapper/blob/main/docs/mapping_cells.md).
- [Current Allen MapMyCells references](https://brain-map.org/bkp/analyze/mapmycells).
- [La Manno et al. 2021](https://www.nature.com/articles/s41586-021-03775-x).
- [Bandler et al. 2022](https://www.nature.com/articles/s41586-021-04237-0).

## 8. Reproducibility, job lineage and navigation

### Code portions used in this workflow

| Portion | Entry points / functions | Exact executed source |
| --- | --- | --- |
| Filter, scoring, paired PCA/UMAP | `run`, `resolve_cycle_genes`, `regress_blocks` | [Frozen comparison source](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/rebuild_20261103/min500_code/scripts/compare_pcdh19_cell_cycle_umaps.py) |
| QC supplement and cell disposition | `render` | [Frozen figure source](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/rebuild_20261103/min500_code/scripts/render_pcdh19_qc_cycle_supplement.py) |
| HiCAT global regression and scaling | `fit_cycle_model`, `corrected_block` | [Frozen current HiCAT adapter](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/rebuild_20261103/hicat_min500_cycle/code/hicat/cycle_rebuild.py) |
| HiCAT corrected local representation | `CycleEngine._representation` | [Frozen current HiCAT adapter](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/rebuild_20261103/hicat_min500_cycle/code/hicat/cycle_rebuild.py) |
| HiCAT fit, held-out assignment, final workflow | `fit`, `map_corrected`, `aggregate` | [Frozen current HiCAT adapter](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/rebuild_20261103/hicat_min500_cycle/code/hicat/cycle_rebuild.py) |
| Recursive graph/DE fitter | `FullScaleEngine._visit`, `PilotEngine._merge`, `PilotEngine.run` | [Frozen full-scale engine](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/rebuild_20261103/hicat_min500_cycle/code/hicat/fullscale_engine.py), [base fitter](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/rebuild_20261103/hicat_min500_cycle/code/hicat/engine.py) |
| Coassignment merging/refinement | `Membership.merge_by_co`, `Membership.refine` | [Frozen consensus core](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/rebuild_20261103/hicat_min500_cycle/code/hicat/consensus_core.py) |
| Developmental label-tree audit | `main` | [Executed validator snapshot](/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/paper3_pcdh19/results/rebuild_20261103/mapmycells_developmental_reference/validate_developmental_mapmycells_reference.py) |

The code ledger includes hashes for repository and executed snapshots. Three
inherited HiCAT files have later explanatory docstrings in the repository;
their executable abstract syntax trees were checked and are identical to the
frozen files. Current regression, figure and reference-validation entry points
are byte-identical to their executed snapshots. No running code was edited
while preparing this Methods document.


| Artifact | Current location relative to RUN |
| --- | --- |
| Current job registry | `job_ids.json` |
| Filtered comparison source/config snapshot | `min500_code/` and its `MANIFEST.json` |
| HiCAT scientific controls / exact seeds | `hicat_min500_cycle/config.json`, `seeds.tsv` |
| HiCAT submitted dependency chain | `hicat_min500_cycle/JOBS.json` |
| HiCAT frozen adapter / reference hashes | `hicat_min500_cycle/FROZEN.json`, `REFERENCE.json` |
| HiCAT software environment | `hicat_min500_cycle/environment.txt` |
| Prepared HiCAT input identity | `hicat_min500_cycle/PREPARED.json` |
| Per-iteration outputs | `hicat_min500_cycle/iterations/000` through `099` |
| Full-population initializer | `hicat_min500_cycle/fullfit/` |
| Consensus endpoint and pairwise audit | `hicat_min500_cycle/final/` when completed |
| Reference audit and code snapshot | `mapmycells_developmental_reference/` |
| Code/function Methods ledger | `METHODS_CODE_LEDGER.tsv` |
| Revised QC figure | `supplementary_QC_cell_cycle/` |

Some inherited `fitting` JSON fields, including the old benchmark scope and
`executed_subsample_iterations=1`, are historical metadata. They are not
current completion counts. The current runner uses top-level `iterations=100`;
per-iteration completion files and `final/SUCCESS.json` establish progress.

Submitted jobs: comparison **65323559** (completed); HiCAT preparation
**65323905** (completed); revised supplement **65323910** (completed);
first fit **65323906** (completed, 01:04:31); full initializer **65323907**
(running); iteration array **65323908** (tasks1–4 running, tasks5–99 pending
at the four-task concurrency limit); final consensus/DE **65323909** (waiting
on those dependencies). There is no MapMyCells or URD job for this
new cell universe. Job IDs alone do not establish completion.

The first-fit hold during discussion of unregressed HiCAT DE was subsequently
resolved after the user clarified the representation and requested proceeding.
The first fit was confirmed running. `CLUSTERING_HOLD.json` records that
resolution; it is not an outstanding instruction to stop the current run.

Software for UMAP/QC: Python3.11, Scanpy1.11.5, AnnData0.12.14, NumPy2.4.4,
SciPy1.17.1, umap-learn0.5.12, scikit-learn1.8.0. Actual versions are captured
in the completed comparison's `parameters.json`. HiCAT uses its separately
pinned Python3.8 environment described above.

Old replaceable analyses were removed under user authorization; preserve the
frozen final figure `ROOT/final_figures/mge_activity_IEG_FINAL_20261101/`,
raw Step00, approved Step02 and reference/software assets. Cleanup records
are in `cleanup_status.json` / `cleanup_manifest.json`. The old main figure
is not silently relabeled as a result of this rerun. No git commit is implied
by the existence of a runtime code snapshot.

## 9. Methods-ready prose: completed preprocessing only

Single-cell RNA-sequencing data from 12 dissected E14.5 mouse MGE samples were
processed in Python. The Cell Ranger filtered-barcode delivery contained
450,788 cells and 19,071 genes. Cells were filtered separately within each
sample using a five-scaled-MAD union rule for low log-transformed UMI counts,
low log-transformed detected-gene counts, or high mitochondrial UMI percentage,
removing 4,439 cells. An additional requirement of at least 500 detected genes
removed 1,095 cells, yielding 445,254 retained cells. No doublet-based exclusion,
ambient correction or batch integration was applied in this workflow. For the
UMAP comparison, 3,000 variable genes were selected from retained raw counts
using sample-aware Seurat-v3 feature selection. Expression was normalized to
10,000 counts per cell and natural-log transformed. S-phase and G2/M scores
were computed using 43 and 53 matched reference genes, respectively. An
ordinary least-squares model containing an intercept and both scores was used
to obtain residual expression for the selected genes. Unregressed and residual
expression were separately SD-scaled and upper-clipped at 10, followed by
centered 50-component PCA, cosine 30-neighbor graph construction and UMAP
(min_dist0.3, spread1, seed20260903). Both embeddings contained the identical
retained cell population, and original count data were preserved.

**Do not yet write that the new HiCAT consensus, cluster number, MapMyCells
annotations or URD trajectories are finalized.** Their Methods/result prose
must be updated from completed output contracts and review records. The
HiCAT parameters above describe the submitted/running implementation;
MapMyCells details above describe method selection and reference validation.


## Operational concurrency history — 4 November 2026

The previous September production run allowed up to **10 simultaneous iteration
jobs**. It scheduled 100 iterations; the separately published 34-cluster review
used 98. That earlier 12-hour recovery snapshot included reused checkpoints and
is not a clean end-to-end runtime benchmark.

The current November run initially allowed **four simultaneous iteration jobs**
(array 65323908, tasks 1–99, after separately completed iteration 0). On
4 November the user requested increasing concurrency above ten urgently. The
array throttle was raised from **4 to 20** without restarting completed work or
changing the 100 iterations, seeds, sampling or scientific parameters. Each task
requests 8 CPUs and 150 GiB; 20 running tasks request 160 CPUs and 3,000 GiB
in total. The throttle is a maximum, not a promise of scheduler availability.

**Observed scheduler snapshot: 2026-11-04 13:39 UTC / 08:39 EST:**
20 iteration tasks actually RUNNING (032–051), 48 PENDING (052–099);
32 iterations completed in the immediately preceding checkpoint/accounting
check (000–031). Full initializer 65323907 completed; final consensus/DE
65323909 awaits the iteration array. No final cluster count is available.

Runtime change record: `RUN/hicat_min500_cycle/CONCURRENCY_UPDATE_20261104.json`.
No background polling or monitoring process was started. A subsequent user
message requested cancellation with an ambiguous target; clarification was
requested before cancelling any jobs. This snapshot does not claim cancellation.


### 2026-11-04 09:12 EST — temporary-storage failure and recovery

At 09:11 EST, five tasks failed (037, 038, 040, 042, 046). Logs show
Errno28/no space in inherited `/home/elcrespo/.codex/tmp`, contrary to the
intended node-local scratch behavior. Increasing concurrency exposed this
storage problem. The prior afternoon completion estimate is withdrawn.
Stopped original unfinished array 65323908 and dependent consensus 65323909;
all 32 completed iteration seals (000–031) preserved. Submitted remaining
68 iterations as recovery array **65352448**, maximum **20 concurrent**,
using a unique node-local `/tmp` directory per job with exit cleanup.
Replacement consensus/DE **65352449** depends on successful recovery.
Original frozen scientific source/config/seeds and wrapper are unchanged;
a separate recovery wrapper sets temporary-storage environment explicitly.
Runtime evidence: `hicat_min500_cycle/NODE_LOCAL_RECOVERY_20261104.json`;
`JOBS.json` now records replacement jobs and preserves original job IDs.
Recovery submitted; completion and actual concurrency require subsequent
scheduler/checkpoint verification. No final consensus count available.


## Restart correction — 4 November 2026, 10:05 AM EST

Scratch migration65352824 completed at 09:23:57 EST; all908 files checksum
verified. Array65352825 then failed on existing interrupted progress logs:
`Progress` refuses to overwrite `hicat_progress_events.jsonl`. No new fits
completed;32 sealed fits preserved. Failed/interrupted iteration directories
are moved intact to `hicat_min500_cycle/interrupted_attempts_before_retry_20261104/`
by preparation65353943 after old writers stop. Old array/consensus cancelled.
Qualification iteration032:65353944. Remaining033–099:65353945, max20,
dependent on successful qualification. Final consensus:65353946. This gate
prevents repeating the failure across the full array before one full fit passes.
No scientific settings changed. Scratch `SCRATCH_LOCATION.json` and HiCAT
`JOBS.json` record current IDs. Prior IDs above are historical. Completion
ETA withdrawn pending successful qualification and resumed throughput.


### Immediate parallel release — 4 November 2026

Removed single-fit dependency from array65353945 at user request.
Iteration032 (65353944) continues independently;033–099 may run up to20
at once (up to21 total while032 is active). Final65353946 explicitly awaits
both65353944 and65353945. No cancelled/restarted numerical work in this change.
Scratch outputs and sealed fits preserved. Current numerical progress already
passed original scratch path and progress-file errors. See runtime
`PARALLEL_RELEASE_20261104.json`.


### User concurrency correction — 4 November 2026

User reports prior testing supports only **10 jobs simultaneously**. This
supersedes requests/records for20. Array65353945 throttle reduced to**9**,
plus independent iteration032/65353944 gives at most10. No single-fit
dependency on array; jobs eligible immediately. Keep total active HiCAT fits
<=10. After032 completes, array may be raised to10 after verification.

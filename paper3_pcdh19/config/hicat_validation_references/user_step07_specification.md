# HiCAT Pilot Biological Validation and Stability Review

## Goal

We have completed an exploratory coarse-to-fine HiCAT pilot on 12,000 cells, sampled as 1,000 cells from each of 12 samples. The current fine clustering contains approximately 38–39 exploratory clusters. Before scaling HiCAT to the full 446,349-cell Step 02 dataset, perform a dedicated biological validity, nuisance-axis, reproducibility, and parameter-sensitivity review.

The goal is not to optimize parameters to obtain a preferred number of clusters and not to lock biological annotations.

The goal is to determine:

1. whether the coarse hierarchy captures recognizable biological compartments;
2. whether fine clusters represent reproducible biological states rather than stochastic graph partitions;
3. whether major boundaries are explained primarily by developmental identity, maturation, cell cycle, sample composition, stress/QC, or technical effects;
4. whether proposed provisional identities such as vRG-like, bRG/oRG-like, IPC, cycling neuroblast, postmitotic neuroblast, immature interneuron, MGE-like, LGE-like, etc. are supported by independent canonical marker programs;
5. which boundaries remain stable under an Allen-reference HiCAT parameterization;
6. whether the current clustering strategy is sufficiently credible to justify a full-data run.

Do not perform full-data HiCAT clustering in this step.

Do not use provisional biological labels to fit clusters or tune parameters. They are hypotheses to test after clustering.

## Inputs

Use the completed expanded coarse-to-fine HiCAT pilot as the primary clustering input.

Primary clustering run:

- 12,000 cells total
- exactly 1,000 cells per sample
- 12 samples
- 19,071 genes
- raw counts from approved Step 02
- existing coarse and fine cluster assignments
- existing parent-child relationships
- existing second-seed clustering
- existing top-gene statistics
- existing Step 06 UMAP coordinates for display only
- no batch correction
- no regression
- no integration

Also use the approved Step 02 / Step 06 full-data object containing 446,349 cells for display and marker-score validation only.

The full-data object must not be reclustered during this checkpoint.

Use the attached provisional annotation workbook only as a source of annotation hypotheses and candidate marker expectations. Do not treat its labels as ground truth.

## Overall output

Create one versioned review package:

results/hicat/02_hierarchy_validation/<RUN_ID>/

The package must contain:

outputs/
    figures/
    tables/
    marker_programs/
    stability/
    sample_composition/
    cell_cycle/
    full_data_projection/
    parameter_sensitivity/
    annotation_review/
    validation/

Produce:

1. one main multipanel PDF called:
hicat_hierarchy_validation_main.pdf

2. one detailed supplementary PDF called:
hicat_hierarchy_validation_detailed.pdf

3. individual PNG/PDF panels;

4. complete TSV/CSV source tables for every plotted value;

5. a machine-readable summary;

6. a human-readable:
REVIEW_README.md

The README should end with a concise table classifying each current coarse and fine cluster as:

SUPPORTED
PLAUSIBLE_BUT_UNSTABLE
CELL_CYCLE_DOMINATED
MATURATION_DOMINATED
SAMPLE_ASSOCIATED
TECHNICAL_OR_STRESS_ASSOCIATED
INSUFFICIENT_EVIDENCE

These are review categories, not biological annotations.

## Analysis 1 — canonical biological program scores

The first priority is to determine what major biological axis each cluster represents.

Compute per-cell scores for curated developmental and regional gene programs.

Do not derive these programs from the current cluster DEGs.

Use canonical literature-supported genes.

At minimum construct the following programs.

### A. Radial glia / neural progenitor program

Include markers such as:

Sox2
Sox1
Fabp7
Slc1a3
Vim
Nes
Hes1
Hes5
Notch1
Glul
Aldoc

Separate this from cell cycle.

### B. Apical / ventricular radial glia-like program

Candidate genes:

Sox2
Hes1
Hes5
Pard3
Cdh2
Prom1
Ezr
Scrib
Tjp1
Vim

Use only genes actually detected and appropriate for embryonic mouse tissue.

The purpose is not to claim true anatomical ventricular position from scRNA-seq alone.

Call this program:

apical_RG_like

rather than definitive "vRG" unless supported.

### C. Basal/oRG-like or delaminating progenitor program

Evaluate markers associated with basal radial glia / delaminating neurogenic progenitors.

Candidate genes may include:

Hopx
Ptprz1
Tnc
Fbln1
Fbln2
Igfbp7
Gpm6a
Ptn
Hmga2

However, this is embryonic mouse MGE and canonical cortical human oRG signatures may not translate directly.

Therefore:

- report whether an oRG/bRG-like program is supported;
- do not automatically call a cluster oRG based on one marker;
- explicitly report whether the evidence is weak, partial, or regionally inappropriate.

### D. IPC / neurogenic progenitor program

Candidate genes:

Ascl1
Insm1
Eomes
Neurog1
Neurog2
Ccnd2
Dlx1
Dlx2
Dlx5
Dlx6
Sp9

For MGE, ASCL1/INSM1/DLX programs may be more appropriate than cortical EOMES.

Do not require Eomes for an MGE IPC-like state.

### E. Cycling progenitor S-phase program

Use established S-phase genes, including:

Mcm2
Mcm3
Mcm4
Mcm5
Mcm6
Mcm7
Pcna
Tyms
Rrm1
Rrm2
Fen1

### F. Cycling progenitor G2/M program

Include:

Mki67
Top2a
Cenpf
Ccnb1
Ccnb2
Cdk1
Ube2c
Birc5
Nusap1
Tpx2
Aurkb

### G. Neuroblast / neuronal commitment program

Include genes such as:

Dcx
Tubb3
Stmn2
Stmn4
Sox11
Sox4
Tbr1
Elavl3
Elavl4
Gap43
Map1b

Adjust markers where appropriate for MGE lineage.

### H. Immature inhibitory neuron program

Include:

Gad1
Gad2
Slc32a1
Dlx1
Dlx2
Dlx5
Dlx6
Arx
Sox6
Erbb4

### I. MGE identity program

Include:

Nkx2-1
Lhx6
Lhx8
Sox6
Arx
Shh
Gli1
Olig2
Zic1
Zic4

Use developmental-stage judgment.

Do not require mature interneuron genes in progenitors.

### J. LGE / striatal lineage program

Include markers such as:

Gsx2
Meis2
Isl1
Ebf1
Foxp1
Foxp2
Sp8
Sp9
Penk

Interpret developmental-stage-specific expression carefully.

### K. CGE-like program

Include appropriate markers such as:

Nr2f1
Nr2f2
Prox1
Htr3a
Reln
Sp8
Sp9

Do not expect mature CGE markers to be strongly expressed in all embryonic progenitors.

### L. POA / preoptic-like program

Use literature-supported developmental POA markers where available.

Keep this program separate from MGE.

### M. neuronal maturation program

Capture transition from early neuroblast to increasingly differentiated neuron.

Candidate genes:

Dcx
Tubb3
Stmn2
Map1b
Gap43
Elavl3
Elavl4
Snap25
Syt1
Rbfox3
Map2

Distinguish early neuronal commitment from later maturation.

### N. stress / dissociation program

Include:

Fos
Jun
Junb
Egr1
Atf3
Ddit3
Hspa1a
Hspa1b
Hsp90aa1

### O. mitochondrial / low-quality expression program

Use mitochondrial fraction and appropriate mitochondrial genes as a diagnostic axis.

Do not interpret mitochondrial genes as biological identity.

### P. ribosomal / translation program

Calculate ribosomal fraction or an aggregate ribosomal score because some apparent clusters may be driven mainly by ribosomal abundance.

## Figure 1 — cluster × biological-program heatmap

Produce a heatmap with:

Rows:
all coarse clusters
all fine clusters

Columns:
RG/stemness
apical_RG_like
basal_RG_like
IPC/neurogenic
S-phase
G2M
neuroblast
immature inhibitory neuron
neuronal maturation
MGE
LGE
CGE
POA
stress
mitochondrial/QC
ribosomal

For each cluster report:

- mean score;
- median score;
- fraction of cells above an appropriate program threshold;
- cluster size.

Z-score program means across clusters only for visualization.

Retain absolute scores in source tables.

Separate coarse and fine heatmaps if necessary for readability.

Do not cluster the biological programs into an uninterpretable order.

Use a biologically meaningful column order.

## Figure 2 — canonical marker dotplot

Generate a carefully curated dotplot independent of top-20 DEGs.

Rows or columns should be clusters.

Markers should be grouped visibly into biological categories.

Use:

- color = mean normalized expression;
- point size = percentage of cells with detectable expression.

Include approximately 40–80 informative markers.

At minimum include representatives of:

Sox2
Fabp7
Slc1a3
Hes1
Hes5

Ascl1
Insm1
Ccnd2
Dlx1
Dlx2

Mki67
Top2a
Pcna
Mcm genes

Dcx
Tubb3
Sox11
Stmn2
Map1b

Gad1
Gad2
Slc32a1
Arx
Sox6

Nkx2-1
Lhx6
Lhx8
Olig2

Gsx2
Meis2
Isl1
Ebf1

Nr2f1
Nr2f2
Prox1

Snap25
Syt1
Rbfox3

Fos
Jun
Atf3
Hspa1a

Generate separate coarse and fine versions if a single figure becomes unreadable.

## Analysis 2 — determine whether cluster boundaries are primarily cell-cycle driven

Do not immediately regress cell cycle from the expression matrix.

First quantify its relationship to the clustering.

Assign or score every pilot cell for:

G1-like
S
G2/M

Use a standard mouse-compatible cell-cycle gene set.

Generate the following.

### Figure 3A — cell-cycle composition per cluster

For every coarse and fine cluster plot:

fraction G1
fraction S
fraction G2/M

as a 100% stacked bar.

### Figure 3B — S and G2/M score distributions

For every cluster show distributions of:

S score
G2/M score

using violin/box plots or equivalent.

### Figure 3C — identity versus cell-cycle scatter

At the cell level, plot:

x = progenitor → neuronal maturation score
y = cell-cycle score

or two versions:

y = S score
y = G2/M score

color cells by fine cluster.

The objective is to determine whether fine clusters occupy distinct developmental identities or simply different cell-cycle phases.

### Figure 3D — within-cell-cycle marker reassessment

For clusters suspected to differ because of cycling:

compare their canonical identity programs within matched cell-cycle states.

For example:

cluster A G1 vs cluster B G1
cluster A S vs cluster B S
cluster A G2/M vs cluster B G2/M

where sufficient cells exist.

Do not require formal DE if cell counts are too small.

The key output is whether regional/developmental identity differences remain after stratification.

For each fine-cluster pair within a parent, report:

identity_difference_persists_after_cell_cycle_stratification:
    yes
    partial
    no
    insufficient_cells

## Analysis 3 — sample representation and potential sample-driven clusters

Every cluster must be evaluated for sample composition.

Use all 12 samples.

Do not aggregate samples by genotype for the primary diagnostic.

### Figure 4A — sample composition of each cluster

For every fine cluster create a 100% stacked bar:

fraction from JZ-1
fraction from JZ-2
...
fraction from JZ-12

### Figure 4B — cluster composition within each sample

Plot the reciprocal view:

for each sample:

fraction of its 1,000 pilot cells assigned to each coarse/fine cluster

### Figure 4C — observed versus expected sample contribution

Because sampling was exactly balanced at 1,000 cells/sample, expected contribution under complete sample independence is approximately 1/12.

Calculate for every cluster and sample:

observed_fraction / expected_fraction

Create a heatmap.

Also calculate:

- maximum sample fraction;
- sample entropy;
- effective number of samples;
- chi-square statistic or another descriptive test of sample imbalance.

Do not use significance alone because large clusters can make tiny deviations significant.

Flag clusters descriptively if:

- one sample contributes disproportionately;
- only a few samples contribute most cells;
- the cluster is essentially absent from most samples.

Do not automatically remove these clusters.

## Analysis 4 — seed stability

We already have two fits using different random seeds on the same 12,000 cells.

Analyze them explicitly.

### Figure 5A — coarse seed overlap matrix

Rows:
baseline coarse clusters

Columns:
repeat-seed coarse clusters

Values:
number or fraction of shared cells

Show both raw counts and row-normalized fractions.

### Figure 5B — fine seed overlap matrix

Same analysis for fine clusters.

Order clusters to maximize correspondence while retaining the actual parent structure.

Do not rename repeat clusters as if labels were inherently matched.

### Quantitative stability metrics

For every baseline cluster calculate:

- best matching repeat cluster;
- Jaccard index;
- precision;
- recall;
- F1;
- fraction of cells retained together;
- number of repeat clusters receiving >10% of the baseline cluster;
- whether the baseline cluster splits;
- whether multiple baseline clusters merge in the repeat.

Also report global:

ARI
NMI

for coarse and fine partitions.

Assign descriptive categories:

HIGH_STABILITY
MODERATE_STABILITY
LOW_STABILITY

Use explicit thresholds in the config and report them.

Do not claim bootstrap stability because the same cells were used in both fits.

## Analysis 5 — Allen-reference parameter sensitivity

Run one additional controlled HiCAT comparison on the same exact 12,000 selected cells.

Keep all non-DE parameters fixed wherever possible.

Do not perform a broad parameter sweep.

Compare:

### Existing fine settings

q1 = 0.3
qdiff = 0.5
DE score = 150

against:

### Allen-reference fine settings

q1 = 0.4
qdiff = 0.7
DE score = 150

Do not change multiple unrelated parameters.

Use the same starting cells and fixed reproducible seed.

The purpose is not to determine an ideal K.

The purpose is to determine which fine boundaries survive a more canonical Allen-like separation requirement.

Generate:

### Figure 6A — current versus Allen-reference overlap matrix

Rows:
current fine clusters

Columns:
Allen-reference fine clusters

Values:
shared-cell fraction

### Figure 6B — boundary survival table

For every current fine cluster report:

current_cluster
parent
n_cells
best_Allen_match
Jaccard
fraction_retained
merges_with
splits_into
boundary_survives

Set:

boundary_survives = strong / partial / no

based on explicit documented rules.

### Figure 6C — parent-level summary

For every coarse parent report:

- number of current fine clusters;
- number of Allen-reference fine clusters;
- which boundaries disappear;
- which boundaries survive;
- which new boundaries appear.

This is particularly important for parents containing many fine clusters.

## Analysis 6 — developmental coherence

Use the provisional developmental-state hypothesis:

radial glia-like
    ↓
neurogenic progenitor / IPC-like
    ↓
cycling neuroblast
    ↓
postmitotic neuroblast
    ↓
immature inhibitory neuron
    ↓
more mature neuronal state

Do not force every cluster into this chain.

Some branches may represent:

- separate regional identity;
- parallel progenitor states;
- MGE versus LGE lineage;
- non-neural cells;
- stress or technical states.

### Figure 7 — developmental marker progression heatmap

Order clusters according to the best provisional developmental interpretation.

Plot biologically ordered markers/programs:

stemness
apical polarity
neurogenic commitment
cell cycle
neuronal commitment
migration/neuroblast
GABAergic identity
MGE identity
neuronal maturation

The purpose is to ask:

Does the proposed developmental ordering produce a coherent biological progression?

Expected patterns may include:

- decreasing Sox2/Fabp7/Hes programs;
- transient Ascl1/Insm1/Ccnd2;
- cycling programs peaking in proliferative intermediates;
- increasing Dcx/Tubb3/Sox11/Stmn2;
- increasing Gad1/Gad2/Arx/Sox6 in inhibitory neuronal populations;
- increasing Snap25/Syt1/Rbfox3 with maturation.

Do not interpret an imperfect monotonic trajectory as failure if the data show branching.

Explicitly identify apparent branches.

## Analysis 7 — coarse parent quality

Before worrying about all ~39 fine clusters, evaluate whether the coarse groups themselves make biological sense.

For every coarse cluster create a compact summary containing:

cluster ID
n cells
dominant developmental program
dominant regional program
cell-cycle composition
sample entropy
strongest canonical markers
strongest conflicting markers
provisional interpretation
confidence
major concern

Classify each coarse group as one of:

biologically coherent
mixed but interpretable
primarily maturation-defined
primarily cell-cycle-defined
sample-associated
poorly resolved

If the coarse structure itself is weak, clearly state that scaling the hierarchy is premature.

## Analysis 8 — fine cluster biological audit

For every fine cluster generate one compact review page.

Do not rely only on top-20 DEGs.

Each page should contain:

### Header

fine cluster
coarse parent
n cells
fraction of pilot

### Panel A

Canonical marker expression.

### Panel B

Biological program scores.

### Panel C

Cell-cycle composition.

### Panel D

Sample composition.

### Panel E

Seed stability.

### Panel F

Allen-parameter stability.

### Panel G

Top positive inspection genes.

### Text summary

Report:

provisional identity/state:
supporting genes:
supporting programs:
conflicting evidence:
cell-cycle influence:
sample influence:
seed stability:
Allen-parameter stability:
alternative interpretation:
confidence:
additional markers needed:

Use the annotation workbook's proposed names only as hypotheses.

If the data do not support the proposed identity, say so.

Do not force a replacement annotation.

## Analysis 9 — full 446,349-cell contextual validation

Use the full Step 06 unintegrated UMAP only for visualization.

Do not recluster the full data.

Do not transfer the 12k cluster labels to all 446k cells by nearest neighbors in this step.

Generate full-data UMAPs for major canonical genes and program scores.

At minimum show:

Sox2
Fabp7
Ascl1
Insm1
Mki67
Top2a
Dcx
Tubb3
Gad1
Gad2
Nkx2-1
Lhx6
Lhx8
Gsx2
Meis2
Sox11
Snap25

Also show:

RG score
neurogenic progenitor score
S-phase score
G2/M score
neuroblast score
MGE score
LGE score
neuronal maturation score
stress score

Use identical coordinates for all panels.

Use robust plotting limits so a few extreme cells do not dominate the color scale.

### Pilot-cell location overlay

On the same full-data UMAP, highlight the 12,000 cells used in the HiCAT pilot.

Then generate separate overlays for each coarse cluster.

For fine clusters, create paginated overlays rather than trying to show all ~39 simultaneously.

The question is:

Do the pilot-derived groups occupy coherent regions of the full-data manifold, or are they scattered broadly across it?

Do not use UMAP location alone as proof of identity.

## Analysis 10 — cluster relationship dendrogram versus actual hierarchy

We already have expression dendrograms.

Make the distinction explicit between:

HiCAT parent-child hierarchy

and:

expression-similarity dendrogram

Generate a comparison figure showing:

- actual coarse→fine parent tree;
- expression dendrogram of fine cluster means;
- major biological program annotations alongside both.

Flag cases where:

- two fine clusters from different parents are extremely similar;
- fine siblings are unexpectedly distant;
- cross-parent pairs fail final DE separation;
- expression relationships conflict strongly with the imposed coarse nesting.

These are particularly important before scaling.

## Main multipanel figure

Create a concise main review figure titled:

Biological and stability validation of the exploratory HiCAT hierarchy

The figure should contain approximately these panels:

A. Coarse → fine hierarchy with cell counts.
B. Cluster × biological-program heatmap.
C. Canonical developmental/regional marker dotplot.
D. Cell-cycle composition by fine cluster.
E. Sample composition / observed-versus-expected sample heatmap.
F. Seed-to-seed fine cluster overlap.
G. Current-versus-Allen-reference fine cluster overlap.
H. Full-data UMAP showing major developmental programs and pilot-cell locations.

The main figure should answer:

Are the major HiCAT boundaries biologically coherent, reproducible, sample-independent, and robust to a canonical Allen-like separation criterion?

## Summary decision table

Create:

cluster_validation_summary.tsv

One row per fine cluster.

Columns should include at least:

fine_cluster
coarse_parent
n_cells

provisional_identity
identity_confidence

dominant_developmental_program
dominant_regional_program

RG_score
IPC_score
neuroblast_score
MGE_score
LGE_score
maturation_score

S_score
G2M_score
cell_cycle_class
cell_cycle_concern

top_sample
top_sample_fraction
sample_entropy
sample_concern

seed_best_match
seed_jaccard
seed_stability

Allen_best_match
Allen_jaccard
Allen_boundary_survival

stress_score
technical_concern

canonical_supporting_markers
conflicting_markers
alternative_interpretation

overall_validation_status
recommended_action

Allowed recommended actions:

retain_candidate_cluster
retain_but_review
likely_merge
cell_cycle_boundary_review
sample_effect_review
insufficient_evidence

Do not actually merge clusters during this checkpoint.

## Parent-level decision table

Also create:

parent_validation_summary.tsv

For each coarse parent report:

parent
n_cells
n_fine_current
n_fine_repeat
n_fine_Allen_reference
dominant_identity
dominant_developmental_state
cell_cycle_structure
sample_structure
fine_boundary_stability
cross_parent_DE_issue
recommended_status

## Explicit questions the report must answer

The final README must directly answer the following.

1. Are the coarse groups biologically interpretable?
2. Which coarse groups are primarily distinguished by regional identity?
3. Which groups are primarily distinguished by progenitor-to-neuron maturation?
4. Which fine boundaries are strongly associated with cell cycle?
5. After cell-cycle stratification, do those clusters retain independent identity differences?
6. Are any fine clusters dominated by one or a few samples?
7. Which fine clusters reproduce under the second random seed?
8. Which fine boundaries survive the Allen-reference q1=0.4, qdiff=0.7, DEscore=150 condition?
9. Which current clusters disappear or merge under that condition?
10. Do the proposed vRG-like, bRG/oRG-like, IPC, cycling neuroblast, neuroblast and immature-neuron annotations form a coherent developmental organization?
11. Which provisional annotations are strongly supported, weakly supported, contradicted, or presently unresolvable?
12. Does the 12k pilot adequately sample the biological structures visible in the full 446k-cell manifold?
13. Are there major reasons not to proceed to full-data HiCAT?

## Important interpretation rules

Do not equate:

cluster = cell type

A cluster may represent:

cell identity
developmental state
maturation
cell cycle
stress
technical effect
sample effect
continuous transition

Explicitly separate these explanations.

Do not call a cycling cluster a distinct progenitor identity merely because it has strong DEGs.

Do not call a neuroblast subtype solely because of different cell-cycle status.

Do not call oRG/bRG solely from HOPX or a partial cortical oRG program.

Do not call MGE solely from one NKX2-1-positive gene.

Use combinations of regional and developmental markers.

Do not require mature markers in embryonic progenitors.

Do not interpret genotype enrichment as evidence that a cluster is biologically valid.

Do not tune clustering parameters based on whether they recover a desired biological label.

Do not remove clusters during this checkpoint.

Do not integrate or batch-correct.

Do not regress cell-cycle genes from the clustering during this checkpoint.

Cell-cycle stratification is diagnostic.

## Criteria for proceeding to the full-data run

At the end of the report, provide one of:

READY_FOR_FULL_DATA
READY_WITH_CAVEATS
NOT_READY

A recommendation of READY_FOR_FULL_DATA should require that:

1. the coarse hierarchy is biologically coherent;
2. major clusters are not obviously single-sample artifacts;
3. important fine boundaries are reasonably seed-stable;
4. biologically important boundaries are not exclusively cell-cycle driven;
5. a meaningful subset of fine boundaries survives the Allen-reference sensitivity condition;
6. full-data marker distributions support the broad developmental/regional structure seen in the pilot;
7. there is no major implementation or hierarchy inconsistency.

Do not require all fine clusters to be perfectly stable.

It is acceptable for the full-data run to resolve weak pilot boundaries differently.

The key question is whether the clustering framework itself is producing biologically interpretable and technically credible structure.

## Priority

If runtime or report complexity becomes limiting, prioritize in this order:

1. biological-program heatmap;
2. canonical marker dotplot;
3. cell-cycle dependence;
4. sample-composition diagnostics;
5. seed stability;
6. Allen-reference sensitivity;
7. full-data UMAP marker/program validation;
8. individual fine-cluster review pages.

Do not sacrifice the first six analyses in order to generate more decorative figures.

Stop after generating the complete review package.

Do not proceed automatically to the 446,349-cell HiCAT run.

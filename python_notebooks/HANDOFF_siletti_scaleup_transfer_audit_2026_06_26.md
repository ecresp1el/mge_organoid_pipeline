# Handoff: Siletti Scale-Up Transfer Audit

Date logged: 2026-06-26

This audit explains why the original Siletti MGE/Jia-style transfers looked
reasonable, but the scaled all-supercluster transfer failed or behaved oddly.

## Bottom Line

The all-supercluster v1 failure was a real pipeline bug. The all-supercluster
v2 run fixed that bug, but it is still not directly comparable to the original
MGE/Jia-style transfer because the biological question and label space changed.

The old run asked:

```text
Within a restricted adult MGE/LLC reference, which Jia-style adult MGE subtype
does each DIV90 cell resemble?
```

The scaled run asked:

```text
Across all 31 Siletti WHB superclusters, including Splatter and non-neuronal
classes, which broad source supercluster is nearest?
```

Those are not equivalent classification problems.

## Root Cause 1: Multi-H5AD Gene-Order Bug

The bridge exporter originally built one gene-to-column index from the first
H5AD in the scope and reused it for every other H5AD.

This was wrong because the Siletti/CELLxGENE supercluster H5ADs have the same
gene set but different `var["Gene"]` order.

Confirmed examples:

```text
First all-supercluster H5AD:
  Upper-layer intratelencephalic

GAD1 index in first H5AD:
  39638

GAD1 index in other H5ADs:
  MGE interneuron: 39967
  CGE interneuron: 48190
  LAMP5-LHX6 and Chandelier: 41673
  Splatter: 12845
```

Consequence:

```text
All-supercluster v1 had only the first H5AD, Upper-layer intratelencephalic,
properly gene-aligned. Other supercluster blocks were column-scrambled, so
fast-kNN collapsed all 16,206 DIV90 cells to Upper-layer intratelencephalic.
```

Fix:

```text
python_notebooks/scripts/export_siletti_div90_seurat_bridge.py
```

The exporter now builds per-H5AD gene indices and intersects unique genes across
all selected reference H5ADs plus DIV90 query.

## Why Older MGE Runs Looked Better

The original focused run:

```text
bridge: siletti_div90_seurat_bridge_v1/mge_llc
scope: MGE interneuron + LAMP5-LHX6 and Chandelier
label_column: candidate_jia_group
excluded label: Excluded / not assigned to Jia-style 9 groups
reference cells exported: 26,290
reference cells used after label filtering: 18,459
query cells: 16,206
```

Reference label distribution before filtering:

```text
Excluded / not assigned to Jia-style 9 groups: 7,831
Cortical PV+ basket neurons: 7,227
Cortical SST+ Mt neurons: 6,889
Cortical SST+ nMt neurons: 2,400
Cortical SST+ LRP neurons: 700
Cortical PV+ Chandelier neurons: 700
Subpallial SST+ neurons: 300
Subpallial SST+ LRP neurons: 143
Subpallial PV+ neurons: 100
```

Old MGE/LLC candidate-Jia assignments:

```text
Cortical PV+ basket neurons: 8,997
Subpallial SST+ neurons: 3,869
Cortical SST+ Mt neurons: 3,057
Cortical SST+ LRP neurons: 277
Cortical SST+ nMt neurons: 5
Subpallial PV+ neurons: 1
```

Important interpretation:

```text
The old run did not ask "which of all Siletti superclusters is this cell?"
It forced the query into a curated adult MGE/Jia-style label space after
removing the non-Jia bucket. That is why it gave MGE-like labels.
```

The old multi-H5AD focused scopes should still be treated cautiously because
they were generated before the per-H5AD gene-order fix. In `mge_llc`, the first
H5AD was MGE, so the MGE block was correctly indexed; the LLC block was likely
gene-order scrambled. This can still yield MGE-like calls because the usable
candidate-Jia reference was dominated by MGE-derived labels.

## Root Cause 2: Scaled v2 Uses a Different Label Space

Corrected all-supercluster v2:

```text
bridge: siletti_div90_all_supercluster_plot_bridge_v2/all_superclusters
transfer: siletti_div90_fast_knn_all_supercluster_source_supercluster_v2
label_column: source_supercluster
excluded labels: none
reference cells: 60,000
query cells: 16,206
```

The v2 bridge is gene-aligned correctly, but `source_supercluster` includes all
31 broad WHB labels, including `Splatter`, non-neuronal classes, and many cells
that were not part of the focused Jia-style label set.

All-supercluster v2 reference composition:

```text
Splatter: 20,181
MGE interneuron: 4,078
CGE interneuron: 3,392
Amygdala excitatory: 2,768
Upper-layer intratelencephalic: 2,740
Medium spiny neuron: 2,709
Deep-layer intratelencephalic: 2,558
Deep-layer corticothalamic and 6b: 2,122
Thalamic excitatory: 1,913
Miscellaneous: 1,807
```

The broad all-supercluster downsample is not balanced by supercluster. It is
downsampled by Siletti subcluster first, then capped at 60,000 total cells.
Because Splatter has many subclusters, it contributes about one third of the
all-supercluster reference.

All-supercluster v2 assignments:

```text
Splatter: 16,122
Miscellaneous: 31
MGE interneuron: 26
CGE interneuron: 24
Committed oligodendrocyte precursor: 2
Fibroblast: 1
```

Interpretation:

```text
This is not the same as the old MGE/Jia-style transfer. It says that when the
classifier is allowed to choose broad Siletti superclusters from an imbalanced
all-WHB reference, most DIV90 cells choose Splatter. That may reflect Splatter's
heterogeneous/ambiguous reference content and reference imbalance, not a useful
biological identity call.
```

## Root Cause 3: Feature Selection Changes With Scope

The fast-kNN workflow selects the top 3,000 variable genes across the stacked
reference-plus-query matrix.

In the old focused MGE/LLC run, those variable genes were selected within a
restricted inhibitory-neuron reference. This makes the feature space more
relevant to MGE/PV/SST-like distinctions.

In the all-supercluster run, feature selection spans excitatory neurons,
inhibitory neurons, glia, vascular cells, Splatter, and DIV90 query. The
selected features can be driven by broad adult brain class differences and
non-neuronal contrasts rather than the fetal/MGE distinctions that motivated
the original plot.

## Audit Test Submitted

To separate matrix correctness from label-space effects, an audit-only transfer
was submitted on the corrected all-supercluster v2 bridge using the old
candidate-Jia label/exclusion rules.

```text
job_id: 52413882
job_name: siletti-all-jia-audit-v2
final state: COMPLETED
elapsed: 00:01:01
max RSS: 11532112K
node: gl3358

bridge: siletti_div90_all_supercluster_plot_bridge_v2/all_superclusters
transfer run label: siletti_div90_fast_knn_all_supercluster_candidate_jia_v2_audit
label_column: candidate_jia_group
exclude labels:
  Excluded / not assigned to Jia-style 9 groups
  unassigned_jia_group
  unlabeled_or_na
```

Audit result:

```text
Reference after filtering to candidate-Jia labels: 3,928 cells
Query cells: 16,206
Labels available: 9

Predicted Jia group counts:
  Subpallial SST+ LRP neurons: 6,587
  Cortical SST+ Mt neurons: 3,763
  Cortical PV+ basket neurons: 3,621
  Subpallial Cholinergic neurons: 2,124
  Subpallial SST+ neurons: 110
  Cortical SST+ nMt neurons: 1

Score summary:
  mean max score: 0.399710
  median max score: 0.356882
  min max score: 0.199915
  max max score: 1.000000
```

Audit interpretation:

```text
The corrected all-supercluster v2 matrix can produce MGE/Jia-style assignments
when the old curated label column and exclusion rules are used. Therefore, the
continued Splatter dominance in the v2 `source_supercluster` run is not caused
by persistent matrix corruption. It is mainly a label-space/reference-composition
problem:

  - raw `source_supercluster` is too broad for the original MGE subtype question;
  - Splatter is included as a normal competing label;
  - the v2 all-supercluster reference is heavily enriched for Splatter;
  - the old run excluded non-Jia labels and forced classification among curated
    MGE/Jia-style adult groups.
```

## Practical Recommendation

Do not use `source_supercluster` across all 31 WHB superclusters as the primary
final biological label for DIV90 MGE organoids.

For final figure work, use one of these safer designs:

```text
1. Rerun the original focused MGE/LLC or MGE/CGE/LLC reference with the fixed
   per-H5AD gene indexing.

2. Transfer curated labels such as `candidate_jia_group` or
   `transferred_mtg_label`, not raw `source_supercluster`.

3. Exclude or isolate Splatter instead of letting it compete as a normal broad
   label.

4. If an all-supercluster context figure is needed, show it as a diagnostic
   reference screen, not as the main identity assignment.

5. Consider balancing by supercluster as well as by subcluster before
   source-supercluster voting.
```

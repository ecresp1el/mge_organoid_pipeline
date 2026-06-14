# Handoff: Siletti 2023 Adult Human Whole-Brain Reference Metadata

This handoff records the first metadata-only reconnaissance for:

```text
K. Siletti et al.
Transcriptomic diversity of cell types across the adult human brain
Science 382, eadd7046
```

This is a major future reference for the project. It is not staged as an
analysis object yet, and no expression matrices have been downloaded in this
repo/project run.

## Current Status

Date logged: 2026-06-14

Only lightweight public taxonomy metadata was inspected. This was done without
Slurm and without a compute node.

Do not treat this as a completed Siletti reference workflow. Treat it as the
initial metadata inventory that should guide a future Slurm-managed staging
workflow.

The generated metadata tables and report are saved under the project results
root, not inside the Git repo:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_metadata/siletti_2023_whb_metadata_inventory_v1
```

## Public Metadata Source

Allen Brain Cell Atlas public S3:

```text
s3://allen-brain-cell-atlas/metadata/WHB-taxonomy/20240330/
```

Files used for this metadata audit:

```text
cluster.csv
cluster_annotation_term.csv
cluster_to_cluster_annotation_membership.csv
```

The audit did not use the large WHB expression matrices.

## Results Metadata Log Files

Generated summary tables:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_metadata/siletti_2023_whb_metadata_inventory_v1/tables/siletti_whb_supercluster_summary.tsv
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_metadata/siletti_2023_whb_metadata_inventory_v1/tables/siletti_whb_nested_cluster_names_by_supercluster.tsv
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_metadata/siletti_2023_whb_metadata_inventory_v1/tables/siletti_whb_taxonomy_metadata_sources.tsv
```

Focused MGE/CGE/LAMP5-LHX6/chandelier tables:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_metadata/siletti_2023_whb_metadata_inventory_v1/tables/siletti_whb_focus_interneuron_supercluster_summary.tsv
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_metadata/siletti_2023_whb_metadata_inventory_v1/tables/siletti_whb_focus_interneuron_nested_cluster_names.tsv
```

Markdown report:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_metadata/siletti_2023_whb_metadata_inventory_v1/reports/siletti_whb_metadata_inventory_report.md
```

Interpretation of table fields:

```text
supercluster_name
  Top-level WHB taxonomy name.

supercluster_order
  WHB taxonomy order from cluster_annotation_term.csv.

supercluster_total_cells_from_term_table
  Total cells reported for that supercluster by the WHB taxonomy term table.

n_unique_nested_cluster_names
  Number of intermediate WHB cluster names nested inside the supercluster.

nested_cluster_name
  Intermediate WHB taxonomy cluster name, one level below supercluster.

n_subclusters
  Number of finest subcluster labels nested under that intermediate cluster.
```

## Supercluster Summary

There are 31 WHB superclusters in the taxonomy. Of these, 30 have more than one
intermediate nested cluster name. `Bergmann glia` is the only supercluster with
exactly one nested cluster name.

| Order | Supercluster | Cells | Nested cluster names |
|---:|---|---:|---:|
| 0 | Upper-layer intratelencephalic | 455,006 | 16 |
| 1 | Deep-layer intratelencephalic | 228,467 | 16 |
| 2 | Deep-layer near-projecting | 18,856 | 13 |
| 3 | Deep-layer corticothalamic and 6b | 78,396 | 17 |
| 4 | MGE interneuron | 222,434 | 26 |
| 5 | CGE interneuron | 227,671 | 21 |
| 6 | LAMP5-LHX6 and Chandelier | 45,118 | 12 |
| 7 | Miscellaneous | 25,071 | 22 |
| 8 | Hippocampal CA1-3 | 74,979 | 15 |
| 9 | Hippocampal CA4 | 10,654 | 9 |
| 10 | Hippocampal dentate gyrus | 67,533 | 7 |
| 11 | Amygdala excitatory | 109,452 | 20 |
| 12 | Medium spiny neuron | 152,189 | 18 |
| 13 | Eccentric medium spiny neuron | 40,144 | 14 |
| 14 | Splatter | 291,833 | 92 |
| 15 | Mammillary body | 16,602 | 11 |
| 16 | Thalamic excitatory | 85,546 | 17 |
| 17 | Midbrain-derived inhibitory | 126,782 | 12 |
| 18 | Upper rhombic lip | 137,162 | 6 |
| 19 | Cerebellar inhibitory | 14,411 | 10 |
| 20 | Lower rhombic lip | 52,650 | 8 |
| 21 | Oligodendrocyte | 490,246 | 8 |
| 22 | Committed oligodendrocyte precursor | 4,720 | 7 |
| 23 | Oligodendrocyte precursor | 105,734 | 5 |
| 24 | Astrocyte | 155,025 | 13 |
| 25 | Ependymal | 5,882 | 10 |
| 26 | Microglia | 91,838 | 9 |
| 27 | Vascular | 9,932 | 11 |
| 28 | Bergmann glia | 8,041 | 1 |
| 29 | Fibroblast | 9,156 | 8 |
| 30 | Choroid plexus | 7,689 | 7 |

## Inhibitory-Neuron-Relevant Top-Level Categories

The adult WHB taxonomy explicitly separates several inhibitory or
subpallial-relevant top-level categories:

| Supercluster | Cells | Nested cluster names | Initial interpretation |
|---|---:|---:|---|
| MGE interneuron | 222,434 | 26 | Forebrain interneuron class with MGE identity in the WHB taxonomy. |
| CGE interneuron | 227,671 | 21 | Forebrain interneuron class with CGE identity in the WHB taxonomy. |
| LAMP5-LHX6 and Chandelier | 45,118 | 12 | Adult interneuron class related to LAMP5-LHX6 and chandelier biology; keep separate from broad MGE/CGE until reviewed. |
| Medium spiny neuron | 152,189 | 18 | Subpallial basal nuclei/striatal neuron class; not the same thing as cortical MGE interneurons. |
| Eccentric medium spiny neuron | 40,144 | 14 | Subpallial/amygdala-basal forebrain-related spiny neuron class. |
| Midbrain-derived inhibitory | 126,782 | 12 | Inhibitory but not forebrain MGE/CGE; mostly midbrain/thalamic context in the WHB taxonomy. |
| Cerebellar inhibitory | 14,411 | 10 | Inhibitory but cerebellar/rhombic context; not a GE reference class. |

Focused isolated categories requested for the first reference decision:

| Supercluster | Cells | Nested cluster names | Actual nested cluster-name prefix |
|---|---:|---:|---|
| MGE interneuron | 222,434 | 26 | `MGE_*` |
| CGE interneuron | 227,671 | 21 | `CGE_*` |
| LAMP5-LHX6 and Chandelier | 45,118 | 12 | `LLC_*` |

The complete nested names, cell counts, subcluster counts, and neurotransmitter
labels for these three categories are in:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_metadata/siletti_2023_whb_metadata_inventory_v1/tables/siletti_whb_focus_interneuron_nested_cluster_names.tsv
```

## Operational Notes

This metadata audit is small enough for login-node inspection because it only
reads the public taxonomy CSVs and writes compact summaries. Any future Siletti
reference staging that downloads large metadata, expression matrices, AnnData
objects, or performs label transfer should follow the project Slurm pattern
used by the Schmitz and Shi reference workflows.

Planned next stage, once compute/data staging is available:

```text
Download/stage Siletti WHB data as a dedicated adult human whole-brain
reference and use it to integrate or label-transfer DIV90 cells against adult
MGE interneuron, CGE interneuron, and LAMP5-LHX6/chandelier reference classes.
```

Future Siletti staging should use its own workflow namespace, for example:

```text
results/siletti_2023_whb_reference_metadata/
results/siletti_2023_whb_reference_label_transfer/
```

Do not mix Siletti adult WHB taxonomy labels into the current DIV30/DIV90 URD
root/tip runs. Use it as an external adult human reference after a dedicated
staging and label-selection decision.

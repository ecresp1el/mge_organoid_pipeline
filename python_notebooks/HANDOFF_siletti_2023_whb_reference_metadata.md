# Handoff: Siletti 2023 Adult Human Whole-Brain Reference Metadata

This handoff records the first metadata-only reconnaissance for:

```text
K. Siletti et al.
Transcriptomic diversity of cell types across the adult human brain
Science 382, eadd7046
```

This is a major reference for the project. The initial metadata-only inventory
has now been extended into a Slurm-managed DIV90 label-transfer staging run using
the adult Siletti/CELLxGENE MGE, CGE, and LAMP5-LHX6/chandelier supercluster
H5AD files. Large generated data are saved under the project results root, not
inside this Git repo.

## Current Status

Date logged: 2026-06-14

The first pass only inspected lightweight public taxonomy metadata without Slurm
or a compute node. On 2026-06-14 this was extended into a Great Lakes
Slurm-managed staging and Seurat label-transfer workflow.

Do not treat this as a final biological call yet. Treat it as a reproducible
first Siletti adult-reference transfer sweep whose knobs and outputs are logged
below for review and follow-up plotting.

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

## Linnarsson Subcluster Annotation Workbook

Additional cross-reference metadata was downloaded from the Linnarsson lab
adult-human-brain repository:

```text
https://github.com/linnarsson-lab/adult-human-brain/blob/main/tables/subcluster_annotation.xlsx
```

Local source copy:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_metadata/siletti_2023_whb_metadata_inventory_v1/source/subcluster_annotation.xlsx
```

This workbook is important because it links WHB subcluster/cluster/supercluster
metadata to author cross-reference fields including:

```text
Subcluster
Cluster
Supercluster
Transferred MTG Label (Transferred from cluster level)
Class
Neurotransmitter
Neuropeptide
Top ROIGroupFine
Top ROI
Top enriched genes
Number of cells
Number of donors
```

Extracted workbook tables and reports:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_metadata/siletti_2023_whb_metadata_inventory_v1/tables/siletti_linnarsson_subcluster_annotation__sheet1.tsv
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_metadata/siletti_2023_whb_metadata_inventory_v1/tables/siletti_linnarsson_subcluster_annotation_focus_mge_cge_llc.tsv
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_metadata/siletti_2023_whb_metadata_inventory_v1/tables/siletti_linnarsson_focus_mge_cge_llc_mtg_label_summary.tsv
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_metadata/siletti_2023_whb_metadata_inventory_v1/reports/siletti_linnarsson_subcluster_annotation_inventory.md
```

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

Important subtype-label caveat from the ABC taxonomy metadata-only audit:

```text
The WHB taxonomy metadata inspected here does not explicitly label these
MGE/CGE/LLC nested cluster names as PV/PVALB, SST, or LRP classes. The available
metadata labels at this stage are broad supercluster names, intermediate names
such as MGE_261/CGE_288/LLC_264, finest subcluster IDs, and neurotransmitter
labels such as GABA or GABA VGLUT3.
```

Therefore, any future PV/SST/LRP call should be treated as a marker/expression
analysis or label-transfer result, not as something already supplied by this
metadata inventory.

However, the Linnarsson `subcluster_annotation.xlsx` workbook adds transferred
MTG cluster-level labels that can be used as an explicit cross-reference for
adult interneuron subtype interpretation:

| Supercluster | Transferred MTG label | Subclusters | Cells | Clusters |
|---|---|---:|---:|---:|
| MGE interneuron | Sst | 100 | 109,584 | 11 |
| MGE interneuron | Pvalb | 74 | 101,898 | 11 |
| MGE interneuron | unlabeled/NA | 37 | 10,952 | 4 |
| CGE interneuron | Vip | 97 | 122,280 | 10 |
| CGE interneuron | Lamp5 | 16 | 58,999 | 3 |
| CGE interneuron | Sncg | 14 | 20,701 | 2 |
| CGE interneuron | Pax6 | 7 | 7,177 | 1 |
| CGE interneuron | unlabeled/NA | 46 | 18,514 | 5 |
| LAMP5-LHX6 and Chandelier | Lamp5 Lhx6 | 27 | 26,400 | 7 |
| LAMP5-LHX6 and Chandelier | Chandelier | 8 | 10,920 | 2 |
| LAMP5-LHX6 and Chandelier | unlabeled/NA | 20 | 7,798 | 3 |

These transferred MTG labels should be used as cross-reference annotations, not
as evidence that the ABC WHB taxonomy files alone already encode PV/SST/LRP
labels.

## Jia-Style Adult Groups And Fetal Marker-Pair Validation

Date logged: 2026-06-14

The locally transferred `RUN_THIS_ANALYSIS.py` from the prior
`siletti_jia9_fetal_marker_pair_validation` work was adopted into the repo as a
reusable workflow script:

```text
python_notebooks/scripts/siletti_jia9_fetal_marker_pair_validation.py
```

This is still a lightweight workbook-only workflow. It does not download
Siletti WHB expression matrices. It was first smoke-tested directly with the
project Python env and then submitted through the normal Great Lakes Slurm
job-file pattern.

Repo Slurm template:

```text
slurm_templates/38_siletti_jia9_fetal_marker_pair_validation.sbatch.template
```

Copied job file:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/38_siletti_jia9_fetal_marker_pair_validation.sbatch
```

Slurm submission:

```text
job_id: 51778113
job_name: siletti-jia9-fetal
state: COMPLETED
exit_code: 0:0
elapsed: 00:00:09
max_rss: 133496K
```

Slurm logs:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/siletti-jia9-fetal-51778113.out
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/siletti-jia9-fetal-51778113.err
```

The project Python env used by the Slurm job:

```text
/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python
```

Underlying workflow command inside the Slurm job:

```bash
/home/elcrespo/miniconda3/envs/mge-organoid-python/bin/python \
  python_notebooks/scripts/siletti_jia9_fetal_marker_pair_validation.py \
  --overwrite \
  --previous-folder /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_metadata/local_siletti_jia9_fetal_marker_pair_validation_import/previous_jia9_gene_signature_v2_outputs
```

Output run:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_metadata/siletti_jia9_fetal_marker_pair_validation_v1
```

Output layout:

```text
tables/
plots/
reports/
source/
```

Key tables:

```text
tables/best_fetal_marker_pair_by_adult_candidate_group.csv
tables/adult_candidate_group_by_fetal_marker_pair_summary.csv
tables/row_level_fetal_marker_pair_scores_selected_only.csv
```

Key reports:

```text
reports/README_fetal_marker_pair_validation.md
reports/INTERPRETATION_fetal_marker_pair_validation.md
```

The imported Mac/Google Drive transfer is preserved separately as provenance,
not as the canonical workflow:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_metadata/local_siletti_jia9_fetal_marker_pair_validation_import
```

Run summary:

| Stage | Rows | Subclusters | Parent clusters | Cells |
|---|---:|---:|---:|---:|
| all workbook | 3,313 | 3,313 | 461 | 3,369,219 |
| candidate Jia-style 9 groups | 204 | 204 | 30 | 225,959 |
| excluded | 3,109 | 3,109 | 437 | 3,143,260 |

Best fetal marker-pair matches by adult candidate group:

| Adult candidate group | Best fetal marker pair | Cell-weighted score |
|---|---|---:|
| Subpallial Cholinergic neurons | LHX8/ISL1 fetal neurons | 0.745 |
| Subpallial PV+ neurons | LHX6/NFIA fetal neurons | 0.500 |
| Subpallial SST+ LRP neurons | LHX6/NFIA fetal neurons | 0.408 |
| Cortical PV+ basket neurons | LHX6/NFIA fetal neurons | 0.123 |
| Subpallial SST+ neurons | LHX6/NFIA fetal neurons | 0.055 |
| Cortical SST+ Mt neurons | LHX6/NFIA fetal neurons | 0.040 |
| Cortical SST+ nMt neurons | LHX6/NFIA fetal neurons | 0.021 |
| Cortical PV+ Chandelier neurons | EPHA5/MEF2C fetal neurons | 0.000 |
| Cortical SST+ LRP neurons | EPHA5/MEF2C fetal neurons | 0.000 |

Important interpretation notes:

- This is marker-pair overlap from the workbook `Top enriched genes` field, not
  label transfer.
- `Cortical PV+ basket neurons` is an operational adult group label assigned
  from `Pallial/cortical` anatomy, transferred MTG label `Pvalb`, and
  `MGE interneuron` supercluster. The word `basket` is not a literal
  Siletti/Linnarsson workbook label.
- A case-insensitive metadata search found no literal `basket` label in the
  workbook.
- The fetal marker-pair layer tests `EPHA5/MEF2C`, `LHX6/NFIA`,
  `CRABP1/ANGPT2`, `NR2F1/NR2F2`, and `LHX8/ISL1`.
- Absence of fetal marker-pair hits should not be interpreted as absence of a
  developmental relationship because adult endpoints may not retain fetal
  marker genes.

## DIV90 Jia-Like Siletti Label-Transfer Staging

Date logged: 2026-06-14

Purpose: run a Jia-like adult-reference label transfer test from Siletti adult
human brain interneuron reference classes onto the existing DIV90 neuron-only
query set. This is the first reproducible pass and is deliberately batched as a
knob sweep rather than a single opaque job.

Primary public source used for the large AnnData objects:

```text
CELLxGENE collection: 283d65eb-dd53-496d-adb7-7570c7caa443
K. Siletti et al., Transcriptomic diversity of cell types across the adult human brain
```

The Linnarsson README advertises large downloadable `Neurons.h5ad` and
`Nonneurons.h5ad` files, but the direct legacy Google Storage path tested on
Great Lakes returned `NoSuchBucket`. The workflow therefore stages the relevant
CELLxGENE per-supercluster H5AD assets instead.

Downloaded reference H5AD assets:

| Adult reference scope | CELLxGENE dataset id | Local file | Cells |
|---|---|---|---:|
| MGE interneuron | `e4710a02-8abc-48d5-a3e8-9ae7e9d79bdb` | `siletti_whb_mge_interneuron.h5ad` | 222,434 |
| CGE interneuron | `bdb26abd-f4ba-4ea3-8862-c2340e7a4f55` | `siletti_whb_cge_interneuron.h5ad` | 227,671 |
| LAMP5-LHX6 and Chandelier | `8f48c031-6692-4c82-b6c2-e4ad35cbc0aa` | `siletti_whb_lamp5_lhx6_and_chandelier.h5ad` | 45,118 |

Local staging root:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer
```

Repo workflow files added for this run:

```text
python_notebooks/scripts/fetch_siletti_cellxgene_supercluster_h5ads.py
python_notebooks/scripts/prepare_siletti_div90_transfer_inputs.py
python_notebooks/scripts/export_siletti_div90_seurat_bridge.py
scripts/31_run_siletti_div90_seurat_label_transfer.R
slurm_templates/39_fetch_siletti_cellxgene_supercluster_h5ads.sbatch.template
slurm_templates/40_prepare_siletti_div90_transfer_inputs.sbatch.template
slurm_templates/41_export_siletti_div90_seurat_bridge_array.sbatch.template
slurm_templates/42_siletti_div90_seurat_label_transfer_array.sbatch.template
```

DIV90 query source:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/varela_div90.h5ad
```

Neuron-only DIV90 clusters included for the query:

| Cluster | DIV90 label |
|---:|---|
| 0 | MGE Striatal/GP Fated |
| 1 | SST+, NPY+, Cortical Fated |
| 2 | CRABP1+/PV Precursors |
| 3 | PV precursors/Migrating cells/Cortical-fated |
| 5 | LHX8+ vMGE GABergic Striatal/GP fated 1 |
| 8 | LHX8+ vMGE GABergic Striatal/GP fated 2 |
| 11 | PV Precursors |

Excluded for this transfer pass: astrocyte/pre-astrocyte, OPC/pre-OPC, stressed,
and dividing clusters. Included query size: 16,206 cells.

Completed upstream jobs:

| Job id | Job name | State | Purpose |
|---:|---|---|---|
| 51779592 | `siletti-fetch-h5ad` | COMPLETED | Download CELLxGENE MGE, CGE, LLC H5ADs and write checksums/manifests. |
| 51779593 | `siletti-div90-prep` | COMPLETED | Build DIV90 neuron manifest, marker/DE summaries, and Figure A heatmaps. |
| 51779660 | `siletti-div90-bridge` | FAILED | First bridge attempt; path comparison bug collected no AnnData blocks. |
| 51779661 | `siletti-div90-xfer` | CANCELLED | Dependency-cancelled after failed bridge attempt. |
| 51779674 | `siletti-div90-bridge` | COMPLETED | Resubmitted bridge export after path-resolution fix. |
| 51779675 | `siletti-div90-xfer` | FAILED/CANCELLED | First transfer array; failed on Seurat 5.1.0 `FindTransferAnchors` API mismatch. |
| 51779847 | `siletti-div90-xfer` | TIMEOUT | Corrected v2 transfer array. All 8 tasks timed out at 8 hours with no final prediction or diagnostic tables. |
| 51793560 | `siletti-div90-xfer-long` | RUNNING as of 2026-06-16 13:50 EDT | Long-walltime v3 backup array released after v2 timeout. Still no final prediction/diagnostic tables at ~16:49 elapsed. |
| 51848144 | `siletti-div90-xfer-3d` | PENDING dependency as of 2026-06-16 13:50 EDT | 3-day v4 fallback array submitted with `afternotok:51793560_*` after Slurm refused live extension of 51793560. |

Bridge-export scopes:

| Scope | Adult reference included | Exported reference cells | Query cells | Shared genes |
|---|---|---:|---:|---:|
| `mge_llc` | MGE + LAMP5-LHX6/chandelier | 26,290 | 16,206 | 17,849 |
| `mge_cge_llc` | MGE + CGE + LAMP5-LHX6/chandelier | 43,916 | 16,206 | 17,849 |

Reference downsampling knobs used during bridge export:

```text
max reference cells per subcluster: 100
max reference cells total per scope: 60000
seed: 0
adult label column prepared for transfer: candidate_jia_group
```

Transfer array knobs:

```text
reference scopes: mge_llc, mge_cge_llc
reductions: rpca, cca
dims: 20, 50
nfeatures: 3000 variable shared features selected inside Seurat
npcs: 50
k.weight requested: 50
seed: 0
Seurat module: r-seurat/5.1.0-R-4.4.1-c3m7yfq
R script: scripts/31_run_siletti_div90_seurat_label_transfer.R
```

The first transfer array (`51779675`) failed because this installed Seurat
version accepts `reference.reduction` but not `query.reduction` in
`FindTransferAnchors`. The runner was patched to use `do.call()` with the
installed function's supported argument set and to explicitly select variable
shared genes instead of taking the first shared genes by file order. A corrected
v2 transfer array was submitted as job `51779847` with output run label:

```text
siletti_div90_seurat_label_transfer_sweep_v2
```

Current transfer-array status as of 2026-06-14 19:31 EDT:

```text
all 8 tasks in job 51779847 are RUNNING
elapsed: ~00:38:42 per task
completed transfer configs: 0/8
final prediction tables: not yet present
final transfer diagnostics: not yet present
files present so far: selected_transfer_features.tsv for each config
```

Live health check as of 2026-06-14 22:05 EDT:

```text
all 8 tasks in job 51779847 are still RUNNING
elapsed: ~03:12 per task
completed transfer configs: 0/8
final prediction tables: not yet present
final transfer diagnostics: not yet present
each R process is active at ~98-100% CPU, so the jobs are not sleeping
observed RSS: ~13.6 GB for mge_llc tasks and ~17.8-18.2 GB for mge_cge_llc tasks
important inefficiency: each task requested 8 CPUs but the Seurat/R process is
effectively using about one core during this stage
likely current stage: FindTransferAnchors/nearest-neighbor anchor search after
feature selection
```

v2 final status as of 2026-06-15:

```text
job_id: 51779847
state: TIMEOUT for all 8 array tasks
elapsed: ~08:00:26 per task
exit_code: 0:15
completed transfer configs: 0/8
final prediction tables: not present
final transfer diagnostics: not present
files present: selected_transfer_features.tsv for each config only
max RSS: ~13.9-14.2 GB for mge_llc tasks and ~18.1-19.8 GB for mge_cge_llc tasks
```

Backup long-walltime array prepared as of 2026-06-14 23:10 EDT, then released
on 2026-06-15 at 21:01 EDT:

```text
job_id: 51793560
job_name: siletti-div90-xfer-long
state: RUNNING as of 2026-06-15 21:01 EDT
release reason: v2 array 51779847 timed out
dependency after release: none
time_limit: 24:00:00
cpus_per_task: 8
memory: 180G
output run label: siletti_div90_seurat_label_transfer_sweep_v3_longtime
job file: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/42_siletti_div90_seurat_label_transfer_array_longtime_hold.sbatch
```

v3 live status as of 2026-06-16 13:50 EDT:

```text
job_id: 51793560
state: RUNNING for all 8 array tasks
elapsed: ~16:49 per task
time_limit: 24:00:00
time_remaining: ~07:11
completed transfer configs: 0/8
final prediction tables: not present
final transfer diagnostics: not present
files present: selected_transfer_features.tsv for each config only
live processes: active at ~99% CPU but effectively one R core per task
observed RSS: ~14.2-14.5 GB for mge_llc tasks and ~18.3-18.8 GB for mge_cge_llc tasks
attempted live extension to 3 days: failed with Slurm Access/permission denied
```

Automatic 3-day fallback submitted as of 2026-06-16 13:50 EDT:

```text
job_id: 51848144
job_name: siletti-div90-xfer-3d
state: PENDING
reason: Dependency
dependency: afternotok:51793560_*(unfulfilled)
time_limit: 3-00:00:00
cpus_per_task: 8
memory: 180G
output run label: siletti_div90_seurat_label_transfer_sweep_v4_3day
job file: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/42_siletti_div90_seurat_label_transfer_array_3day_fallback.sbatch
```

Key result locations:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/source_cellxgene_superclusters
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_neuron_prep_v1
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_seurat_bridge_v1
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_seurat_label_transfer_sweep_v2
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_seurat_label_transfer_sweep_v3_longtime
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_seurat_label_transfer_sweep_v4_3day
```

Figure/report status:

```text
Figure A heatmaps are generated in siletti_div90_neuron_prep_v1/plots.
Transfer prediction and diagnostic tables are expected from the Seurat array but
were not present at the 2026-06-14 19:31 EDT status check.
Integrated UMAP, river/alluvial plot, and stage/composition bar plots still need
the post-transfer plotting job after the selected transfer configuration is
chosen.
```

## Fast KNN Transfer Debug And Scale-Up

Date logged: 2026-06-16

Why this exists:

```text
The Seurat FindTransferAnchors-based workflow repeatedly reached feature
selection/PCA quickly and then spent many hours inside FindTransferAnchors
without writing final predictions. A capped Seurat pilot confirmed the bottleneck:
load/filter/normalize/scale/PCA were fast, but FindTransferAnchors remained the
wall even at 5,000 reference x 3,000 query cells.
```

Canonical fast-transfer code:

```text
python_notebooks/scripts/siletti_div90_fast_knn_label_transfer.py
```

Canonical fast-transfer Slurm template:

```text
slurm_templates/47_siletti_div90_fast_knn_label_transfer_pilot.sbatch.template
```

The script is intentionally not a Seurat anchor integration. It is a transparent
label-transfer classifier:

```text
bridge MatrixMarket counts
  -> log-normalize reference/query
  -> select high-variance genes
  -> sparse TruncatedSVD
  -> cosine kNN adult-reference vote
  -> per-cell DIV90 predicted adult label + score tables
```

Biology being transferred in the current fast-KNN jobs:

```text
reference scope: mge_llc
adult reference superclusters: MGE interneuron + LAMP5-LHX6/chandelier
query: DIV90 neuron-lineage cells
label column: candidate_jia_group
excluded label before transfer: Excluded / not assigned to Jia-style 9 groups
```

This does not change the biology relative to the original intent. It removes the
non-biological `Excluded / not assigned to Jia-style 9 groups` label before
training the transfer classifier. For the `candidate_jia_group` target,
`mge_cge_llc` is not the right first scale-up because the added CGE cells mostly
become this excluded bucket. If the target becomes broader adult interneuron
labels such as `transferred_mtg_label` (Vip/Lamp5/Sncg/Pax6 included), then
`mge_cge_llc` becomes biologically appropriate.

Cell-cap semantics:

```text
max_reference_cells and max_query_cells are compute/debug caps only.
They are applied after reference label filtering.
If a cap is >0, cells are uniformly sampled with the fixed seed.
If a cap is 0, all available cells are used.
The cap is not a biological subset rule.
```

Reference size for the current `candidate_jia_group` fast-transfer path:

| Stage | Cells |
|---|---:|
| Full staged Siletti MGE H5AD | 222,434 |
| Full staged Siletti LLC H5AD | 45,118 |
| Bridge-exported `mge_llc` reference | 26,290 |
| `mge_llc` after dropping excluded candidate_jia_group label | 18,459 |
| DIV90 neuron-lineage query | 16,206 |

Fast-transfer parameter ledger:

| Run label | Job id | Scope | Reference cap | Query cap | Features | SVD dims | k | Status |
|---|---:|---|---:|---:|---:|---:|---:|---|
| `siletti_div90_fast_knn_smoke_v1` | login smoke | `mge_llc` | 500 | 300 | 1,000 | 20 | 20 | COMPLETED in ~22 sec |
| `siletti_div90_fast_knn_transfer_pilot_v1` | 51852811 | `mge_llc` | 5,000 | 3,000 | 2,000 | 30 | 50 | COMPLETED in 27 sec |
| `siletti_div90_fast_knn_transfer_full_mge_llc_v1` | 51871420 | `mge_llc` | 0 = all usable 18,459 | 0 = all 16,206 | 3,000 | 50 | 50 | COMPLETED in 47 sec |

Fast-transfer output roots:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_fast_knn_smoke_v1
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_fast_knn_transfer_pilot_v1
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_fast_knn_transfer_full_mge_llc_v1
```

Expected fast-transfer outputs per run:

```text
fast_knn/run_progress.tsv
fast_knn/selected_transfer_features.tsv
fast_knn/*_predictions.tsv.gz
fast_knn/*_prediction_scores.tsv.gz
fast_knn/*_transfer_diagnostics.json
tables/*_query_obs_with_predictions.tsv.gz
tables/*_cluster_label_counts.tsv
```

Full `mge_llc` fast-KNN transfer completion:

```text
job_id: 51871420
state: COMPLETED
elapsed: 00:00:47
max_rss: 5,643,392K
n_reference_cells: 18,459
n_query_cells: 16,206
n_features: 3,000
n_components: 50
k: 50
diagnostics: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_fast_knn_transfer_full_mge_llc_v1/mge_llc/svd50_k50_ref0_query0/fast_knn/siletti_div90__candidate_jia_group__fast_knn__svd50__k50_transfer_diagnostics.json
cluster summary: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_fast_knn_transfer_full_mge_llc_v1/mge_llc/svd50_k50_ref0_query0/tables/siletti_div90__candidate_jia_group__fast_knn__svd50__k50_cluster_label_counts.tsv
```

## Provisional Update 2026-06-21: DIV90 Siletti Mapping Debug / Cholinergic Reference

This section is a working audit trail, not a final endorsed method. Several
scripts and result folders were added while debugging the Seurat/Jia-style
mapping path. They may be deleted or replaced later after the biology and
method choice are reviewed.

### Biological scope correction: adult cholinergic cells

The original `mge_llc` adult reference did **not** contain adult cholinergic
cells. It only included:

```text
MGE interneuron
LAMP5-LHX6 and Chandelier
```

The broader Siletti metadata has a Jia-style `Subpallial Cholinergic neurons`
candidate group, but those cells are in:

```text
Siletti supercluster: Splatter
Cluster: 400
Subclusters: 1634, 1635, 1636, 1637, 1638, 1640, 1641, 1642
Regions: mostly BasalForebrain / Amygdala
Marker rule in metadata script: Subpallial + NT-CHOL + CHAT or SLC5A7
```

To keep this explicit, a new adult bridge scope was added:

```text
mge_llc_cholinergic =
  MGE interneuron
  LAMP5-LHX6 and Chandelier
  Splatter restricted to cluster 400 cholinergic subclusters only
```

The full, non-subsampled bridge completed:

```text
job_id: 52082869_2
state: COMPLETED
run_label: siletti_div90_seurat_bridge_cholinergic_full_reference_v1
scope: mge_llc_cholinergic
max_ref_cells_per_subcluster: 0
max_ref_cells_total: 0
n_reference_cells_exported: 268,321
n_query_cells_exported: 16,206
n_shared_genes: 17,849
output: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_seurat_bridge_cholinergic_full_reference_v1/mge_llc_cholinergic
```

The corresponding full fast-KNN feature/transfer run completed:

```text
job_id: 52082870
state: COMPLETED
run_label: siletti_div90_fast_knn_transfer_full_mge_llc_cholinergic_full_reference_v1
scope: mge_llc_cholinergic
output: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_fast_knn_transfer_full_mge_llc_cholinergic_full_reference_v1/mge_llc_cholinergic/svd50_k50_ref0_query0
```

### Seurat RPCA/CCA issue

The initial Seurat RPCA mapping workflow stalled in `FindTransferAnchors`.

Cancelled job:

```text
job_id: 52081495
reason for cancel: ran ~26 hours but only accumulated ~26 minutes CPU
last progress line: find_transfer_anchors start
reference cells in that run: 18,459 usable old mge_llc cells
query cells: 16,206
features: 3,000
reduction: rpca
```

Microdebug jobs showed that this was not simply a full-data-size problem.
Four small tests all timed out at `FindTransferAnchors` after 20 minutes:

```text
job_id: 52083020_[0-3]
adult cells: ~975
DIV90 cells: 1,000
features: 500
dims: 1:10

A1 rpca, k.filter=100, nn.method=annoy -> timeout
A2 rpca, k.filter=NA,  nn.method=annoy -> timeout
A3 cca,  k.filter=100, nn.method=annoy -> timeout
A4 rpca, k.filter=100, nn.method=rann  -> timeout
```

Interpretation:

```text
The stall is not explained by cell count, feature count, k.filter, Annoy vs
RANN, or RPCA alone. It appears to be a Seurat 5.1.0 FindTransferAnchors
behavior for these RPCA/CCA-style calls in this environment/object setup.
```

### Object-integrity audit and pcaproject workaround

A hard object-integrity S0 audit was added before running any further anchor
calls. It validated the original and stripped objects for:

```text
Seurat/R version
Assays and DefaultAssay
RNA Assay5 layers
counts/data dimensions
non-empty RNA data layer with nonzero values
unique feature names
unique cell names
metadata rownames matching object cell names
shared feature length
selected anchor feature presence
reference PCA existence and alignment
requested dims within available PCs
finite PCA embeddings
query RNA data usability for selected features
```

S0 result:

```text
job_id: 52084521
state: COMPLETED
output: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/seurat_anchor_object_integrity_s0
Seurat: 5.1.0
SeuratObject: 5.0.2
R: 4.4.1
original objects: validation PASS
stripped RNA-only objects: validation PASS
stripped pcaproject FindTransferAnchors: COMPLETED
cells: 2,418 adult x 2,999 DIV90
features: 3,000
dims: 1:20
anchors found: 107
anchors retained after filtering: 39
elapsed FindTransferAnchors: 100.139 sec
```

Important interpretation:

```text
The Seurat objects are not obviously malformed by the hard audit.
Seurat FindTransferAnchors is not universally broken.
The validated path that worked was stripped RNA-only objects built from the RNA
data layer, using reduction='pcaproject' and reference.reduction='pca'.
This is a workaround, not yet a biologically finalized method.
```

### Current provisional Seurat pcaproject transfer

A production-style stripped RNA-only pcaproject transfer was launched against
the full cholinergic-aware adult reference. This is currently **in progress**
and should not be interpreted until it completes and outputs are reviewed.

```text
job_id: 52086244
state at handoff update: RUNNING on gl3090
run_label: seurat_pcaproject_transfer_stripped_cholinergic_full_reference
bridge: siletti_div90_seurat_bridge_cholinergic_full_reference_v1/mge_llc_cholinergic
transfer features: siletti_div90_fast_knn_transfer_full_mge_llc_cholinergic_full_reference_v1/mge_llc_cholinergic/svd50_k50_ref0_query0/fast_knn/selected_transfer_features.tsv
reference cells loaded: 268,321
reference cells after excluding non-Jia label: 223,436
query cells: 16,206
features: 3,000
dims: 1:20
method: stripped RNA-only pcaproject anchors + TransferData
output: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/seurat_pcaproject_transfer_stripped_cholinergic_full_reference
```

Current progress at time of handoff edit:

```text
2026-06-21 16:17:53 run start
2026-06-21 16:26:42 load_reference end, elapsed_sec=529.073
2026-06-21 16:26:57 load_query end, elapsed_sec=15.488
2026-06-21 16:26:57 object_sizes_initial reference=268321; query=16206; genes=17849
2026-06-21 16:27:25 filter_reference_label before=268321; after=223436
2026-06-21 16:27:25 prepare_original_reference start
```

Expected outputs if the run completes:

```text
seurat_pcaproject_per_cell_predictions.tsv
seurat_pcaproject_prediction_scores_by_class.tsv
seurat_pcaproject_anchor_summary.tsv
seurat_pcaproject_transfer_progress.tsv
tables/seurat_pcaproject_object_sizes.tsv
seurat/selected_transfer_features.tsv
```

Follow-up check:

```text
squeue -j 52086244
sacct -j 52086244 --format=JobID,JobName%35,State,ExitCode,Elapsed,MaxRSS,NodeList -P
tail -n 80 /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/seurat_pcaproject_transfer_stripped_cholinergic_full_reference/seurat_pcaproject_transfer_progress.tsv
cat /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/seurat_pcaproject_transfer_stripped_cholinergic_full_reference/seurat_pcaproject_anchor_summary.tsv
```

If this completes, review the number of anchors, TransferData score
distributions, and whether DIV90 LHX8/ISL1/cholinergic-like classes receive the
new adult `Subpallial Cholinergic neurons` label at plausible rates. Do not
promote this to the final Jia-style figure workflow until those checks are done.

### Pending RPCA rescue/debug array

An additional RPCA-specific debug array was launched to determine whether RPCA
itself works through Seurat's official integration path, or whether only
`FindTransferAnchors(reduction="rpca")` stalls.

```text
job_id: 52086386_[0-3]
run_label: seurat_rpca_rescue_debug
state at handoff update: PENDING/RUNNING depending on scheduler
output: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/seurat_rpca_rescue_debug
```

Array modes:

```text
0 integration_rpca
  FindIntegrationAnchors(object.list=list(reference, query),
                         reduction="rpca",
                         dims=1:20,
                         k.filter=NA)

1 transfer_rpca
  FindTransferAnchors(reduction="rpca",
                      reference.reduction="pca",
                      k.filter=NA)

2 transfer_rpca_l2_false
  Same transfer RPCA, plus l2.norm=FALSE

3 transfer_rpca_approx_false
  Same transfer RPCA, plus approx.pca=FALSE
```

All use stripped RNA-only objects, shared genes only, ~1,000 adult cells,
~1,000 DIV90 cells, 2,000 features, `npcs=20`, `dims=1:20`, `n.trees=10`, and a
20-minute timeout per task.

Follow-up check:

```text
squeue -j 52086386
sacct -j 52086386 --format=JobID,JobName%35,State,ExitCode,Elapsed,MaxRSS,NodeList -P
cat /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/seurat_rpca_rescue_debug/rpca_rescue_debug_summary.tsv
for d in /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/seurat_rpca_rescue_debug/*; do echo "### $d"; tail -n 40 "$d/progress.tsv" 2>/dev/null; cat "$d/summary.tsv" 2>/dev/null; done
```

Interpretation guide:

```text
FindIntegrationAnchors RPCA completes, transfer RPCA stalls:
  problem is specific to Seurat transfer RPCA.

Both integration RPCA and transfer RPCA stall:
  RPCA internals are likely broken/problematic in this environment or object
  setup.

transfer_rpca_approx_false completes:
  suspect truncated PCA / irlba / approximate PCA behavior.

transfer_rpca_l2_false completes:
  suspect L2 normalization / reduction handling.
```

### New/provisional scripts and templates added during this debug

These are useful for reproducing the audit, but they are not necessarily final
pipeline components:

```text
scripts/32_siletti_div90_seurat_rpca_jia_style_mapping.R
scripts/33_siletti_seurat_rpca_anchor_stall_debug.R
scripts/34_siletti_seurat_anchor_object_integrity_s0.R
scripts/35_siletti_seurat_pcaproject_transfer_stripped.R
scripts/36_siletti_seurat_rpca_rescue_debug.R

slurm_templates/53_siletti_div90_seurat_rpca_jia_style_mapping.sbatch.template
slurm_templates/54_siletti_seurat_rpca_anchor_stall_debug.sbatch.template
slurm_templates/55_siletti_seurat_anchor_parallel_microdebug.sbatch.template
slurm_templates/56_siletti_seurat_anchor_object_integrity_s0.sbatch.template
slurm_templates/57_siletti_seurat_pcaproject_transfer_stripped.sbatch.template
slurm_templates/58_siletti_seurat_rpca_rescue_debug.sbatch.template
```

Also modified during this debug:

```text
python_notebooks/scripts/fetch_siletti_cellxgene_supercluster_h5ads.py
python_notebooks/scripts/export_siletti_div90_seurat_bridge.py
slurm_templates/39_fetch_siletti_cellxgene_supercluster_h5ads.sbatch.template
slurm_templates/41_export_siletti_div90_seurat_bridge_array.sbatch.template
```

Core caution:

```text
Do not treat the stripped pcaproject output as the final Jia-style result until
we inspect anchor counts, prediction score distributions, cholinergic behavior,
and whether the mapping is biologically sensible for DIV90 classes.
```

## Full Staged-Supercluster Reference Asset Rerun

Date logged: 2026-06-26

Reason:

```text
The earlier Siletti DIV90 Jia-style figure run is incomplete for final-figure
use. It used the older `mge_llc` fast-kNN transfer output, which excluded CGE
and did not include the later cholinergic-aware Splatter reference correction.
Before moving to final figures, create a durable full-reference bridge asset so
future Siletti label-transfer/plot variants do not need to reread/re-export the
large CELLxGENE H5ADs unless the source H5ADs or DIV90 query set change.
```

Available staged H5AD superclusters for this rerun:

```text
MGE interneuron
CGE interneuron
LAMP5-LHX6 and Chandelier
Splatter
```

Code changes made for this durable rerun:

```text
python_notebooks/scripts/export_siletti_div90_seurat_bridge.py
  - added scope `mge_cge_llc_cholinergic`
  - added scope `mge_cge_llc_splatter`
  - generalized the Splatter cholinergic filter to any scope ending in
    `_cholinergic`

python_notebooks/scripts/siletti_div90_fast_knn_label_transfer.py
  - added `--exclude-labels`, a `||`-separated multi-label exclusion option
    for broad reference scopes that contain multiple non-biological/no-label
    buckets

slurm_templates/41_export_siletti_div90_seurat_bridge_array.sbatch.template
  - array now includes the two new scopes

slurm_templates/47_siletti_div90_fast_knn_label_transfer_pilot.sbatch.template
  - can pass either the legacy single `SILETTI_TRANSFER_EXCLUDE_LABEL` or the
    new multi-label `SILETTI_TRANSFER_EXCLUDE_LABELS`
```

Bridge asset submitted:

```text
job_id: 52393300_4
job_name: siletti-bridge
state at handoff update: RUNNING on gl3232
run_label: siletti_div90_full_staged_supercluster_bridge_v1
scope: mge_cge_llc_splatter
superclusters: MGE interneuron + CGE interneuron + LAMP5-LHX6 and Chandelier + Splatter
max_ref_cells_per_subcluster: 0
max_ref_cells_total: 0
query: existing DIV90 neuron-lineage query manifest
job file: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/41_export_siletti_div90_full_staged_supercluster_bridge_v1.sbatch
expected output: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_full_staged_supercluster_bridge_v1/mge_cge_llc_splatter
```

This bridge is the asset to preserve. Reuse its `seurat_bridge/` MatrixMarket
files and `tables/` metadata for future Siletti reruns. Do not rerun the H5AD
bridge export unless one of these changes:

```text
source CELLxGENE H5ADs
Siletti metadata/Jia-label mapping
DIV90 query object or query-cell manifest
shared-gene selection logic
desired source supercluster scope
```

Dependent fast-kNN runs submitted from the same bridge:

```text
job_id: 52393320
job_name: siletti-fullknn-jia
state at handoff update: PENDING dependency after 52393300
run_label: siletti_div90_fast_knn_full_staged_superclusters_candidate_jia_v1
scope: mge_cge_llc_splatter
label_column: candidate_jia_group
excluded labels:
  Excluded / not assigned to Jia-style 9 groups
  unassigned_jia_group
  unlabeled_or_na
max_reference_cells: 0
max_query_cells: 0
nfeatures: 3000
n_components: 50
k: 50

job_id: 52393321
job_name: siletti-fullknn-mtg
state at handoff update: PENDING dependency after 52393300
run_label: siletti_div90_fast_knn_full_staged_superclusters_mtg_label_v1
scope: mge_cge_llc_splatter
label_column: transferred_mtg_label
excluded labels:
  unlabeled_or_na
max_reference_cells: 0
max_query_cells: 0
nfeatures: 3000
n_components: 50
k: 50
```

Follow-up checks:

```text
squeue -j 52393300,52393320,52393321
sacct -j 52393300,52393320,52393321 --format=JobID,JobName%35,State,ExitCode,Elapsed,MaxRSS,NodeList -P

cat /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_full_staged_supercluster_bridge_v1/mge_cge_llc_splatter/tables/siletti_div90_seurat_bridge_config.json
cat /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_full_staged_supercluster_bridge_v1/mge_cge_llc_splatter/tables/siletti_reference_label_counts.tsv

cat /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_fast_knn_full_staged_superclusters_candidate_jia_v1/mge_cge_llc_splatter/svd50_k50_ref0_query0/fast_knn/siletti_div90__candidate_jia_group__fast_knn__svd50__k50_transfer_diagnostics.json
cat /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_fast_knn_full_staged_superclusters_mtg_label_v1/mge_cge_llc_splatter/svd50_k50_ref0_query0/fast_knn/siletti_div90__transferred_mtg_label__fast_knn__svd50__k50_transfer_diagnostics.json
```

Interpretation guardrail:

```text
The full `mge_cge_llc_splatter` bridge intentionally keeps all Splatter cells
because the goal is an archival broad staged-supercluster asset. For final
Jia-style biology, compare this broad run against the targeted
`mge_cge_llc_cholinergic` scope if all-Splatter neighbors dominate or obscure
MGE/CGE/LLC interpretation.
```

Update after submission:

```text
Full archival bridge completed:
  job_id: 52393300_4
  state: COMPLETED
  elapsed: 00:12:19
  max_rss: 85,869,300K
  output: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_full_staged_supercluster_bridge_v1/mge_cge_llc_splatter

Full fast-kNN jobs started:
  job_id: 52393320
  job_name: siletti-fullknn-jia
  run_label: siletti_div90_fast_knn_full_staged_superclusters_candidate_jia_v1

  job_id: 52393321
  job_name: siletti-fullknn-mtg
  run_label: siletti_div90_fast_knn_full_staged_superclusters_mtg_label_v1

Full Jia-style plot render submitted behind full Jia fast-kNN:
  job_id: 52395391
  job_name: siletti-full-jia-fig
  dependency: afterok:52393320
  run_label: siletti_div90_jia_style_figure_full_staged_superclusters_v1
  job file: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/48_siletti_full_staged_jia_style_figure_v1.sbatch
```

Downsampled plot-progress branch:

```text
Purpose:
  Keep a faster plotting branch available while full reference transfer/plotting
  runs. This uses the same broad `mge_cge_llc_splatter` scope, but exports at
  most 100 cells per Siletti subcluster and at most 60,000 total adult reference
  cells. Query cells remain uncapped.

Downsampled bridge:
  job_id: 52395346_4
  job_name: siletti-bridge-ds
  run_label: siletti_div90_plot_progress_downsample_bridge_v1
  scope: mge_cge_llc_splatter
  max_ref_cells_per_subcluster: 100
  max_ref_cells_total: 60000
  job file: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/41_export_siletti_div90_plot_progress_downsample_bridge_v1.sbatch
  expected output: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_plot_progress_downsample_bridge_v1/mge_cge_llc_splatter

Downsampled Jia fast-kNN:
  job_id: 52395347
  job_name: siletti-ds-knn-jia
  dependency: afterok:52395346
  run_label: siletti_div90_fast_knn_plot_progress_downsample_candidate_jia_v1

Downsampled Jia-style plot render:
  job_id: 52395348
  job_name: siletti-ds-jia-fig
  dependency: afterok:52395347
  run_label: siletti_div90_jia_style_figure_plot_progress_downsample_v1
  job file: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/48_siletti_plot_progress_downsample_jia_style_figure_v1.sbatch
```

## True All-Supercluster DIV90 Transfer Plots

Date logged: 2026-06-26

Correction:

```text
The `mge_cge_llc_splatter` runs above are broad Jia/Siletti staging runs, but
they are not all Siletti WHB superclusters. They include only the four H5ADs
that were staged at the time: MGE interneuron, CGE interneuron, LAMP5-LHX6 and
Chandelier, and Splatter.

The user requested analogous plots using all Siletti WHB superclusters. That
requires staging all 31 CELLxGENE supercluster H5ADs and transferring a label
that exists for every reference cell. The first all-supercluster pass therefore
uses `source_supercluster`, not `candidate_jia_group`.
```

Code changes for this corrected all-supercluster path:

```text
python_notebooks/scripts/export_siletti_div90_seurat_bridge.py
  - added `ALL_SUPERCLUSTERS` with the 31 WHB superclusters
  - added bridge scope `all_superclusters`
  - added generic H5AD filename resolution using the same safe-token convention
    as `fetch_siletti_cellxgene_supercluster_h5ads.py`

slurm_templates/41_export_siletti_div90_seurat_bridge_array.sbatch.template
  - added `all_superclusters` as array task 0

python_notebooks/scripts/plot_siletti_div90_all_supercluster_figure.py
  - new plotter for all-supercluster transfer results
  - renders adult-reference UMAP overlay, DIV90 predicted-supercluster UMAP,
    DIV90-class-to-Siletti-supercluster river plot, and sample composition plot
```

Submitted all-supercluster job chain:

```text
Fetch/stage all 31 supercluster H5ADs:
  job_id: 52396197
  job_name: siletti-fetch-all
  final state: COMPLETED
  elapsed: 00:20:36
  max RSS: 37247344K
  node: gl3343
  staged files: 31 / 31 complete H5ADs
  target superclusters: all 31 CELLxGENE datasets whose title starts with
    `Supercluster:`
  job file: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/39_fetch_siletti_all_superclusters_h5ads_v1.sbatch

All-supercluster downsampled bridge:
  job_id: 52396203_0
  job_name: siletti-all-bridge
  dependency: afterok:52396197
  run_label: siletti_div90_all_supercluster_plot_bridge_v1
  scope: all_superclusters
  max_ref_cells_per_subcluster: 100
  max_ref_cells_total: 60000
  final state: COMPLETED
  elapsed: 00:04:32
  max RSS: 17670324K
  node: gl3253
  expected output: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_all_supercluster_plot_bridge_v1/all_superclusters

All-supercluster fast-kNN transfer:
  job_id: 52396227
  job_name: siletti-all-knn
  dependency: afterok:52396203
  run_label: siletti_div90_fast_knn_all_supercluster_source_supercluster_v1
  scope: all_superclusters
  label_column: source_supercluster
  excluded labels: NONE
  max_reference_cells: 0
  max_query_cells: 0
  nfeatures: 3000
  n_components: 50
  k: 50
  final state: COMPLETED
  elapsed: 00:00:55
  max RSS: 7514800K
  node: gl3470

All-supercluster plots:
  job_id: 52396231
  job_name: siletti-all-plot
  dependency: afterok:52396227
  final state: COMPLETED
  elapsed: 00:03:09
  max RSS: 4406628K
  node: gl3470
  output: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_all_supercluster_source_supercluster_plots_v1
  completed plots:
    plots/figure_B_all_supercluster_reference_div90_overlay.png/pdf
    plots/figure_B_div90_predicted_siletti_superclusters.png/pdf
    plots/figure_C_div90_class_to_siletti_supercluster_river.png/pdf
    plots/figure_D_sample_predicted_siletti_supercluster_proportions.png/pdf
```

## All-Supercluster Assignment Workflow Details

Date logged: 2026-06-26

This section is the methods-level record for how DIV90 cells are being assigned
to Siletti labels in the true all-supercluster run. It is meant to be sufficient
for bioinformatics drafting without opening the scripts.

Workflow scripts:

```text
python_notebooks/scripts/export_siletti_div90_seurat_bridge.py
python_notebooks/scripts/siletti_div90_fast_knn_label_transfer.py
python_notebooks/scripts/plot_siletti_div90_all_supercluster_figure.py
```

Adult Siletti reference input:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/source_cellxgene_superclusters/h5ad
```

For the all-supercluster run, the bridge scope is `all_superclusters`, which is
the complete set of 31 WHB superclusters listed in the `Supercluster Summary`
section above. Each reference cell receives a `source_supercluster` value from
the H5AD it came from. This is the transferred label for the all-supercluster
plots.

DIV90 query input:

```text
DIV90 AnnData:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/python_anndata/varela_div90.h5ad

DIV90 neuron-lineage cell manifest:
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_neuron_prep_v1/tables/div90_siletti_query_neuron_cells.tsv
```

Bridge export:

```text
run_label: siletti_div90_all_supercluster_plot_bridge_v1
scope: all_superclusters
max_ref_cells_per_subcluster: 100
max_ref_cells_total: 60000
seed: 0
query cells: uncapped DIV90 neuron-lineage query, expected 16,206 cells
```

The bridge exporter reads all 31 Siletti supercluster H5ADs, keeps unique
reference genes from `var["Gene"]`, intersects them with unique DIV90 genes, and
writes matched sparse count matrices. MatrixMarket files are written as genes x
cells, with matching `*_genes.tsv`, `*_barcodes.tsv`, and `*_metadata.tsv.gz`.
The fast-kNN script reads these files and transposes them to cells x genes.

Bridge outputs:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_all_supercluster_plot_bridge_v1/all_superclusters/seurat_bridge/reference_counts.mtx
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_all_supercluster_plot_bridge_v1/all_superclusters/seurat_bridge/reference_metadata.tsv.gz
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_all_supercluster_plot_bridge_v1/all_superclusters/seurat_bridge/query_counts.mtx
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_all_supercluster_plot_bridge_v1/all_superclusters/seurat_bridge/query_metadata.tsv.gz
```

Assignment method:

```text
method: fast_knn_svd_cosine
label_column: source_supercluster
excluded labels: NONE
nfeatures: 3000
n_components: 50
k: 50
max_reference_cells: 0 after bridge export
max_query_cells: 0 after bridge export
seed: 0
```

This is not Seurat anchors and not Seurat `TransferData`. It is a transparent
SVD/cosine k-nearest-neighbor assignment:

```text
1. Read reference and query sparse count matrices from the bridge directory.
2. Confirm reference and query gene order match.
3. Keep all reference labels because `source_supercluster` exists for every
   adult reference cell.
4. Library-size normalize each cell to 10,000 counts and apply log1p.
5. Select the top 3,000 variable genes across the stacked normalized
   reference-plus-query matrix.
6. Fit TruncatedSVD on the adult reference matrix using 50 components.
7. Project DIV90 query cells into the same SVD space.
8. L2-normalize reference and query SVD coordinates.
9. For each DIV90 cell, find 50 adult reference nearest neighbors by cosine
   distance using brute-force sklearn NearestNeighbors.
10. Convert cosine distance to positive similarity weight with:
      weight = max(1 - distance, 0)
11. Normalize neighbor weights within each DIV90 cell.
12. Sum normalized weights by Siletti `source_supercluster`.
13. Assign `predicted.id` to the supercluster with the highest summed weight.
14. Record `prediction.score.max` as that winning summed weight.
```

Primary transfer output:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_fast_knn_all_supercluster_source_supercluster_v1/all_superclusters/svd50_k50_ref0_query0
```

Expected transfer tables:

```text
fast_knn/siletti_div90__source_supercluster__fast_knn__svd50__k50_predictions.tsv.gz
fast_knn/siletti_div90__source_supercluster__fast_knn__svd50__k50_prediction_scores.tsv.gz
fast_knn/siletti_div90__source_supercluster__fast_knn__svd50__k50_transfer_diagnostics.json
fast_knn/selected_transfer_features.tsv
fast_knn/run_progress.tsv
tables/siletti_div90__source_supercluster__fast_knn__svd50__k50_query_obs_with_predictions.tsv.gz
tables/siletti_div90__source_supercluster__fast_knn__svd50__k50_cluster_label_counts.tsv
```

Plotting/asset output:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_all_supercluster_source_supercluster_plots_v1/plots
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_all_supercluster_source_supercluster_plots_v1/tables
```

Expected saved plotting assets:

```text
plots/figure_B_all_supercluster_reference_div90_overlay.png/pdf
plots/figure_B_div90_predicted_siletti_superclusters.png/pdf
plots/figure_C_div90_class_to_siletti_supercluster_river.png/pdf
plots/figure_D_sample_predicted_siletti_supercluster_proportions.png/pdf
tables/figure_C_div90_class_to_siletti_supercluster_edges.tsv
tables/figure_D_sample_predicted_siletti_supercluster_proportions.tsv
tables/div90_query_cells_with_all_siletti_supercluster_assignments.tsv.gz
tables/siletti_reference_cells_with_all_supercluster_umap.tsv.gz
tables/plot_config.json
```

Completion validation:

```text
All four plot PNG/PDF pairs were written on 2026-06-26 at 12:53 EDT.
PNG dimensions:
  figure_B_all_supercluster_reference_div90_overlay.png: 3185 x 2096
  figure_B_div90_predicted_siletti_superclusters.png: 2894 x 2096
  figure_C_div90_class_to_siletti_supercluster_river.png: 3571 x 1817
  figure_D_sample_predicted_siletti_supercluster_proportions.png: 2635 x 1753
```

Post-completion issue found:

```text
Do not use the first all-supercluster plot outputs as final biological figures.

Issue 1: plotting merge bug.
  The all-supercluster plotter read the query-obs prediction table but merged
  predictions back onto the UMAP table using the barcode-like `cell_id` column
  instead of `seurat_cell_id`. This made the plotted/output assignment table say
  `Unassigned` even though the kNN prediction table contained labels.

  Fixed in:
    python_notebooks/scripts/plot_siletti_div90_all_supercluster_figure.py

  Plot-only rerun submitted:
    job_id: 52398318
    job_name: siletti-all-plot-fix
    state at issue note: PENDING

Issue 2: degenerate all-supercluster kNN assignment.
  The actual fast-kNN result is not unassigned, but it is biologically
  suspicious/degenerate: all 16,206 DIV90 query cells were assigned to
  `Upper-layer intratelencephalic`.

  Current assignment counts:
    Upper-layer intratelencephalic: 16,206 / 16,206 DIV90 cells

  By DIV90 broad class:
    MGE Striatal/GP Fated: 3,601 -> Upper-layer intratelencephalic
    SST+, NPY +, Cortical Fated: 3,548 -> Upper-layer intratelencephalic
    CRABP1+/PV Precursors: 3,503 -> Upper-layer intratelencephalic
    PV precursors/Migrating cells/Cortical-fated: 2,283 -> Upper-layer intratelencephalic
    LHX8+ vMGE GABergic Striatal/GP fated 1: 1,924 -> Upper-layer intratelencephalic
    LHX8+ vMGE GABergic Striatal/GP fated 2: 922 -> Upper-layer intratelencephalic
    PV Precursors: 425 -> Upper-layer intratelencephalic

  Score check:
    16,123 / 16,206 cells had Upper-layer score in [0.999, 1.001].
    Only 42 cells had any non-upper score > 0.1; Bergmann glia was the only
    runner-up label in those cases.

Interpretation:
  The original all-supercluster output is useful as a failed diagnostic and as
  proof that all 31 H5ADs/bridge assets can be staged, but not as a final figure.
  Next analysis should diagnose why all-supercluster SVD/cosine kNN collapses
  to the first/excitatory supercluster and should not package these plots as
  final biological evidence.
```

Root cause identified:

```text
The bridge exporter had a multi-H5AD gene-indexing bug.

Bug:
  `export_siletti_div90_seurat_bridge.py` built a gene-to-column index from
  only the first H5AD in the selected scope, then reused that same column index
  for every other H5AD. The 31 CELLxGENE supercluster H5ADs contain the same
  gene set but not the same `var["Gene"]` order.

Consequence:
  In the all-supercluster scope, the first H5AD is
  `Upper-layer intratelencephalic`, so only the upper-layer reference block was
  gene-aligned correctly. The other 30 reference blocks were column-scrambled.
  This explains why the kNN collapsed all DIV90 query cells to the first
  supercluster.

Why earlier focused runs could still look MGE-like:
  The earlier focused scopes started with the MGE H5AD and transferred
  MGE/Jia-style labels from a biologically restricted reference. That avoided
  the specific all-supercluster collapse to upper-layer, but the same
  cross-H5AD gene-order bug means any multi-H5AD bridge made before this fix
  should be treated cautiously and rerun before final use.

Fix:
  `export_siletti_div90_seurat_bridge.py` now builds a separate per-H5AD
  gene-to-column index and intersects unique genes across all selected reference
  H5ADs plus the DIV90 query. It also adds `source_supercluster` to the bridge
  label-count audit.
```

Corrected all-supercluster v2 rerun:

```text
Corrected bridge:
  job_id: 52398418_0
  job_name: siletti-all-bridge-v2
  run_label: siletti_div90_all_supercluster_plot_bridge_v2
  scope: all_superclusters
  output: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_all_supercluster_plot_bridge_v2/all_superclusters

Corrected fast-kNN:
  job_id: 52398421
  job_name: siletti-all-knn-v2
  dependency: afterok:52398418
  run_label: siletti_div90_fast_knn_all_supercluster_source_supercluster_v2
  label_column: source_supercluster
  excluded labels: NONE

Corrected plots:
  job_id: 52398423
  job_name: siletti-all-plot-v2
  dependency: afterok:52398421
  output: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_all_supercluster_source_supercluster_plots_v2

Canceled obsolete v1 plot-fix job:
  job_id: 52398318
```

Draft methods wording:

```text
Adult Siletti reference and DIV90 query count matrices were restricted to shared
unique genes, library-size normalized to 10,000 counts per cell, log1p
transformed, and embedded with TruncatedSVD using the top 3,000 variable shared
genes. DIV90 cells were assigned adult reference labels by cosine
k-nearest-neighbor voting in the SVD space using k = 50 neighbors. Neighbor
votes were weighted by normalized cosine similarity, and the label with the
highest summed weight was reported as the predicted Siletti supercluster.
```

Follow-up checks:

```text
squeue -j 52396197,52396203,52396227,52396231
sacct -j 52396197,52396203,52396227,52396231 --format=JobID,JobName%35,State,ExitCode,Elapsed,MaxRSS,NodeList -P
tail -n 80 /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/logs/siletti-fetch-all-52396197.out
find /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/source_cellxgene_superclusters/h5ad -maxdepth 1 -name 'siletti_whb_*.h5ad' -printf '%f\n' | sort | wc -l
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

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
| 51779847 | `siletti-div90-xfer` | RUNNING as of 2026-06-14 19:31 EDT | Corrected transfer array using installed Seurat 5.1.0 argument set. No array task had completed at the 19:31 status check. |

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

Backup long-walltime array prepared as of 2026-06-14 23:10 EDT:

```text
job_id: 51793560
job_name: siletti-div90-xfer-long
state: PENDING
reason: JobHeldUser
dependency: afternotok:51779847_*(unfulfilled)
time_limit: 24:00:00
cpus_per_task: 8
memory: 180G
output run label: siletti_div90_seurat_label_transfer_sweep_v3_longtime
job file: /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/42_siletti_div90_seurat_label_transfer_array_longtime_hold.sbatch
```

This backup array should remain held unless the current v2 array (`51779847`)
fails or times out. If it needs to run, release it with:

```bash
scontrol release 51793560
```

Key result locations:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/source_cellxgene_superclusters
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_neuron_prep_v1
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_seurat_bridge_v1
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siletti_2023_whb_reference_label_transfer/siletti_div90_seurat_label_transfer_sweep_v2
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

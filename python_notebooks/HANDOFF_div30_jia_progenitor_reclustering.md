# DIV30 Jia-Style Progenitor Reclustering Handoff

## Goal

Replicate Jia et al. Phase 1 progenitor reclustering on the current DIV30
progenitor compartment before attempting any lineage-committed progenitor
states.

This is intentionally a gate analysis. The only question for this phase is
whether the DIV30 progenitor-only compartment separates, without using Jia
scores for clustering, into:

```text
RGC1 / VZ-RGC
RGC2 / SVZ-RGC
IPC
```

Do not perform the second Jia progenitor reclustering until this Phase 1 result
has been reviewed.

## Position Relative To Current Work

This workflow is parallel to the active DIV30 and DIV90 URD work. It does not
modify existing URD objects or tree outputs.

Current nearby handoffs:

```text
python_notebooks/HANDOFF_div30_paper_clusters_urd.md
python_notebooks/HANDOFF_jia_program_div30_scoring.md
python_notebooks/HANDOFF_div90_jia_lineage_urd_plan.md
```

The recent branch history leading into this handoff is mostly DIV90
Jia-lineage URD smoke-run work, preceded by DIV30 URD/Jia validation commits.
This Phase 1 DIV30 progenitor reclustering is the missing upstream check before
using DIV30 progenitor substructure to steer later URD decisions.

## Inputs

Primary DIV30 Seurat object:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/varela_this_paper/varela_this_paper_seurat.rds
```

Post hoc Jia score sidecar:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/jia_program_div30_scoring/jia_program_div30_scoring_v1/tables/div30_jia_program_scores_obs.tsv
```

Input cells are selected only from the current DIV30 Seurat cluster column
(`seurat_clusters`, falling back to `RNA_snn_res.0.2`):

| DIV30 cluster | Current label |
|---:|---|
| 0 | Radial glia |
| 3 | Radial glia |
| 6 | Inhibitory progenitors |
| 7 | Radial glia |

The Jia score table is not used for selecting cells, clustering, PCA,
neighbors, UMAP, or marker detection. It is attached only after unsupervised
reclustering and `FindAllMarkers()` for validation overlays.

## Code Added

R runner:

```text
scripts/28_div30_jia_progenitor_reclustering.R
```

Slurm template:

```text
slurm_templates/35_div30_jia_progenitor_reclustering.sbatch.template
```

Handoff:

```text
python_notebooks/HANDOFF_div30_jia_progenitor_reclustering.md
```

## Reclustering Method

The runner creates a progenitor-only Seurat object from clusters `0,3,6,7`,
then applies the existing DIV30-style Seurat workflow:

```text
NormalizeData
FindVariableFeatures(selection.method = "vst", nfeatures = 3000)
ScaleData
RunPCA(npcs = 50)
FindNeighbors(dims = 1:30)
FindClusters(resolution = 0.8)
RunUMAP(dims = 1:30)
```

Default parameters can be overridden in Slurm with:

```text
NFEATURES=3000
NPCS=50
DIMS=30
RESOLUTION=0.8
SEED=7
MARKER_MIN_PCT=0.1
MARKER_LOGFC_THRESHOLD=0.25
TOP_N_MARKERS=50
```

## Marker Identification

The runner executes:

```text
FindAllMarkers()
```

on the new progenitor reclusters and exports both all markers and top markers.
The top-marker table keeps the key interpretation columns:

```text
cluster
gene
avg_log2FC
pct.1
pct.2
p_val_adj
```

## Biological Annotation

Tentative Phase 1 identities are assigned from candidate marker-set expression:

| Candidate identity | Markers |
|---|---|
| VZ-RGC / RGC1 | HES1, VIM, NES |
| SVZ-RGC / RGC2 | FBLN7, CACNA1E, DACH1 |
| IPC | DLX1, DLX2, ASCL1 |

This annotation is intentionally tentative and should be checked against both
the top markers and the post hoc Jia overlays.

## Validation

After reclustering and marker detection, the script overlays:

```text
jia_score_RGC1
jia_score_RGC2
jia_score_IPC
```

Validation outputs ask whether the unsupervised clusters correspond to:

```text
RGC1 / VZ-RGC
RGC2 / SVZ-RGC
IPC
```

## Expected Output Root

Default run label:

```text
div30_jia_progenitor_reclustering_phase1_v1
```

Default output root:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_jia_progenitor_reclustering/div30_jia_progenitor_reclustering_phase1_v1
```

Key tables:

```text
tables/div30_jia_progenitor_phase1_run_parameters.tsv
tables/div30_jia_progenitor_phase1_cluster_counts.tsv
tables/div30_jia_progenitor_phase1_all_markers.tsv.gz
tables/div30_jia_progenitor_phase1_top_markers.tsv
tables/div30_jia_progenitor_phase1_candidate_marker_panel.tsv
tables/div30_jia_progenitor_phase1_candidate_marker_expression_by_cluster.tsv
tables/div30_jia_progenitor_phase1_tentative_annotations.tsv
tables/div30_jia_progenitor_phase1_jia_score_summary_by_cluster.tsv
tables/div30_jia_progenitor_phase1_cells.tsv.gz
tables/div30_jia_progenitor_phase1_complete.tsv
```

Key plots:

```text
plots/div30_progenitor_umap_by_recluster.png
plots/div30_progenitor_umap_by_source_cluster.png
plots/div30_progenitor_candidate_marker_feature_grid.png
plots/div30_progenitor_candidate_marker_dotplot.png
plots/div30_progenitor_jia_score_umap_overlays.png
plots/div30_progenitor_jia_score_violins_by_recluster.png
```

Saved object:

```text
div30_jia_progenitor_phase1_seurat.rds
```

## Submission

Copy the template into the project jobs folder and submit:

```bash
cp slurm_templates/35_div30_jia_progenitor_reclustering.sbatch.template /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/35_div30_jia_progenitor_reclustering.sbatch
sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/35_div30_jia_progenitor_reclustering.sbatch
```

The template requests:

```text
time: 08:00:00
cpus: 8
mem: 160G
partition: standard
account: parent0
```

## Stop Rule

Stop after reviewing:

```text
div30_jia_progenitor_phase1_tentative_annotations.tsv
div30_jia_progenitor_phase1_top_markers.tsv
div30_jia_progenitor_phase1_jia_score_summary_by_cluster.tsv
plots/div30_progenitor_jia_score_umap_overlays.png
```

Only proceed to lineage-committed progenitor reclustering if the Phase 1
clusters convincingly reproduce the expected RGC1/VZ-RGC, RGC2/SVZ-RGC, and
IPC structure.

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

Default assay:

```text
RNA
```

The runner opens this object directly with `readRDS()`. It does not read the
AnnData/Matrix Market export and it does not create a new integrated object.

Post hoc Jia score sidecar:

```text
/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/jia_program_div30_scoring/jia_program_div30_scoring_v1/tables/div30_jia_program_scores_obs.tsv
```

Expected Jia score columns:

```text
cell_id
jia_score_RGC1
jia_score_RGC2
jia_score_IPC
```

The join key for Jia validation is:

```text
Seurat cell name / colnames(progenitor_object) == jia_scores$cell_id
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

## Exact Programmatic Logic

The R runner is:

```text
scripts/28_div30_jia_progenitor_reclustering.R
```

The exact default command invoked by the Slurm template is equivalent to:

```bash
Rscript scripts/28_div30_jia_progenitor_reclustering.R \
  --project-root /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder \
  --seurat-rds /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/varela_this_paper/varela_this_paper_seurat.rds \
  --jia-scores /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/jia_program_div30_scoring/jia_program_div30_scoring_v1/tables/div30_jia_program_scores_obs.tsv \
  --outdir /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div30_jia_progenitor_reclustering/div30_jia_progenitor_reclustering_phase1_v1 \
  --run-label div30_jia_progenitor_reclustering_phase1_v1 \
  --assay RNA \
  --source-cluster-col "" \
  --progenitor-clusters 0,3,6,7 \
  --nfeatures 3000 \
  --npcs 50 \
  --dims 30 \
  --resolution 0.8 \
  --seed 7 \
  --min-pct 0.1 \
  --logfc-threshold 0.25 \
  --top-n-markers 50
```

Programmatic sequence:

1. Load the DIV30 Seurat RDS.
2. Check that the requested assay exists; default is `RNA`.
3. Choose the source cluster column:
   `seurat_clusters` if present, otherwise `RNA_snn_res.0.2`.
4. Convert source cluster values to character and keep only cells whose source
   cluster is exactly one of `0,3,6,7`.
5. Create the progenitor-only Seurat object with `subset(obj, cells = keep_cells)`.
6. Add metadata:
   `div30_parent_cluster = source_cluster_col`.
7. Add metadata:
   `div30_parent_progenitor_label = "Inhibitory progenitors"` for parent
   cluster `6`, otherwise `"Radial glia"`.
8. Set `DefaultAssay(prog) = "RNA"`.
9. Run `NormalizeData(prog, assay = "RNA")`.
10. Run `FindVariableFeatures(selection.method = "vst", nfeatures = 3000)`.
11. Run `ScaleData()` with no regression variables.
12. Run `RunPCA(npcs = 50, seed.use = 7)`.
13. Run `FindNeighbors(reduction = "pca", dims = 1:30,
    graph.name = c("progenitor_nn", "progenitor_snn"))`.
14. Run `FindClusters(graph.name = "progenitor_snn", resolution = 0.8,
    random.seed = 7, cluster.name = "div30_progenitor_cluster")`.
15. Set identities to `div30_progenitor_cluster`.
16. Run `RunUMAP(reduction = "pca", dims = 1:30,
    reduction.name = "progenitor_umap", reduction.key = "PROGUMAP_",
    seed.use = 7)`.
17. Run `FindAllMarkers()` on the new `div30_progenitor_cluster` identities
    with `only.pos = TRUE`, `min.pct = 0.1`, and `logfc.threshold = 0.25`.
18. Export all markers and the top 50 positive markers per new cluster.
19. Summarize candidate marker expression from the normalized `RNA` data layer.
20. Tentatively annotate each new recluster by the candidate marker set with
    the highest mean expression score.
21. Only after steps 1-20, read the Jia score table and attach
    `jia_score_RGC1`, `jia_score_RGC2`, and `jia_score_IPC` by `cell_id`.
22. Export Jia score summaries and Jia score UMAP/violin validation plots.
23. Save the progenitor-only reclustered Seurat object.

What is not done:

```text
No Jia score is used to select cells.
No Jia score is used in NormalizeData, variable-gene selection, ScaleData, PCA, neighbors, clustering, or UMAP.
No Jia score is used in FindAllMarkers.
No new integration workflow is run.
No URD object is changed.
No lineage-committed progenitor reclustering is run.
```

## Programmatic Annotation Rule

Candidate marker genes are searched case-insensitively against the progenitor
object gene names:

| Candidate identity | Marker-set name in output | Markers |
|---|---|---|
| VZ-RGC / RGC1 | `VZ_RGC_RGC1` | HES1, VIM, NES |
| SVZ-RGC / RGC2 | `SVZ_RGC_RGC2` | FBLN7, CACNA1E, DACH1 |
| IPC | `IPC` | DLX1, DLX2, ASCL1 |

For each found marker gene and each new recluster, the script computes:

```text
mean_expression
median_expression
pct_expressed
```

using the normalized Seurat `RNA` data layer after `NormalizeData()`.

For each marker set within each new recluster, the script computes:

```text
marker_set_score = mean(mean_expression across found genes in that marker set)
```

The tentative identity is the marker set with the highest
`marker_set_score`. This is a first-pass programmatic label only. The final
biological call should be made by reviewing the top markers, candidate marker
table, and post hoc Jia score overlays together.

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

## Next-Step Note: Batch Correction

The Phase 1 reclustering intentionally used the existing DIV30 workflow without
creating a new integration workflow. However, the post hoc metadata UMAPs should
be reviewed for sample/cell-line structure before interpreting biology.

If reclusters track strongly with `orig.ident` / sample ID or with derived cell
line labels (`H9`, `79B`, `2E`), the next reclustering pass should explicitly
test sample-aware correction/integration on the same progenitor-only input
cells before accepting the three Jia-style progenitor groups.

Required comparison before next biological decision:

```text
Uncorrected progenitor UMAP/reclusters
vs.
Sample/batch-corrected progenitor UMAP/reclusters
```

The batch-corrected pass should still preserve the same biological guardrails:

```text
Input cells: DIV30 clusters 0,3,6,7 only
Do not use Jia scores for cell selection
Do not use Jia scores for clustering
Use Jia scores only after clustering for validation
Stop before lineage-committed progenitor reclustering
```

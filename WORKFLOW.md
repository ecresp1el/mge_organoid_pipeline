# Walsh & Bershteyn Single-cell Repro Workflow (Great Lakes)

## Roots
- Repo (code/config only): `/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline`
- Runtime workspace (data/jobs/results): `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder` (PROJECT_ROOT)

## Key Data Locations
- Raw GEO downloads:
  - Walsh: `PROJECT_ROOT/data/raw/walsh_2025_geo_files/`
  - Bershteyn: `PROJECT_ROOT/data/raw/bershteyn_2025_geo_files/`
- Raw GEO downloads (new):
  - Xiang_2018 (GSE97882): `PROJECT_ROOT/data/raw/xiang_2018_geo_files/`
  - Samarasinghe_2021 (GSE165577): `PROJECT_ROOT/data/raw/samarasinghe_2021_geo_files/`
- NeMO landing-page downloads (Siebert 2026, `nemo:dat-htzat9t`):
  - `PROJECT_ROOT/data/raw/siebert_2026_nemo/` (metadata, BDBags, subset fetch + downloads)
  - Status: **TBD** until NeMO publishes manifests/files.
- Processed / results:
  - Walsh: `PROJECT_ROOT/results/walsh_day75/`
  - Bershteyn: `PROJECT_ROOT/results/bershteyn/`
  - Samarasinghe_2021: `PROJECT_ROOT/results/samarasinghe_2021/` (Seurat + UMAP, TBD if rerun)
  - Samarasinghe_2021 LIGER (paper-matched): `PROJECT_ROOT/results/samarasinghe_2021_liger/`
- Logs: `PROJECT_ROOT/logs/`
- Job scripts: `PROJECT_ROOT/jobs/`

## Main Seurat Objects
- Walsh final (post-stress, annotated): `results/walsh_day75/walsh_day75_final_annotated.rds`
- Walsh final (post-stress, unannotated): `results/walsh_day75/walsh_day75_final.rds`
- Bershteyn Seurat object (converted from GEO-supplied serialized RDS): `results/bershteyn/bershteyn_2025_seurat.rds`

## Checkpoints (Walsh)
`results/walsh_day75/checkpoints/`
- `walsh_merged_raw.rds` — merged GSM7979671 dFB + GSM7979672 vFB
- `walsh_postqc.rds` — QC-filtered (1000–5000 genes, <15% MT)
- `walsh_norm_hvg.rds` — LogNormalize + HVG=5000
- `walsh_cellcycle.rds` — after CellCycleScoring
- `walsh_scaled.rds` — ScaleData (regress S.Score, G2M.Score)
- `walsh_pca.rds` — PCA (pre-stress)
- `walsh_pca_poststress.rds` — PCA after stress-cluster removal

## Plots (Walsh)
- Main UMAPs (selected PCs & 20-PC comparison): `results/walsh_day75/plots/`
- Elbows/JackStraw (fallback to Elbow): `results/walsh_day75/plots/`
- k sweeps: `results/walsh_day75/k_sweep_plots/`
- dims sweeps: `results/walsh_day75/dims_sweep_plots/`
- resolution sweeps: `results/walsh_day75/resolution_sweep_plots/`
- k+resolution sweep (e.g., 14/15 cluster views): `results/walsh_day75/kres_sweep_plots/`
- Walsh-group annotated UMAP: `results/walsh_day75/annotation_plots/umap_by_walsh_group.{png,pdf}`

## Tables (Walsh)
- Main run report: `results/walsh_day75/walsh_day75_report.txt`
- Stress scores: `results/walsh_day75/walsh_day75_stress_scores_by_cluster.csv`
- k sweep summary: `results/walsh_day75/walsh_day75_k_sweep_summary.tsv`
- dims sweep summary: `results/walsh_day75/walsh_day75_dims_sweep_summary.tsv`
- resolution sweep summary: `results/walsh_day75/walsh_day75_resolution_sweep_summary.tsv`
- k+resolution sweep summary: `results/walsh_day75/walsh_day75_kres_sweep_summary.tsv`
- Cluster→Walsh-group annotation: `results/walsh_day75/walsh_cluster_annotation.tsv`
- Domain composition by Walsh group: `results/walsh_day75/walsh_domain_composition_by_group.tsv`

## Scripts (repo)
- NeMO Siebert landing page/BDBag: `scripts/01b_download_siebert_nemo.sh`
- Extra GEO downloads (Xiang_2018, Samarasinghe_2021): `scripts/01c_download_extra_geo.sh`
- Samarasinghe 2021 Seurat/UMAP: `scripts/05_samarasinghe_2021_seurat.R`
- Samarasinghe 2021 LIGER (paper settings): `scripts/05b_samarasinghe_2021_liger.R`
- Main pipeline: `scripts/02_walsh_day75_seurat.R` (methods-locked QC/normalize/HVG/CC/Scale/PCA/stress removal/UMAP/clustering; checkpoints; Elbow fallback)
- Sweeps: `scripts/02b_walsh_k_sweep.R`, `02c_walsh_dims_sweep.R`, `02d_walsh_resolution_sweep.R`, `02e_walsh_kres_sweep.R`
- Annotation: `scripts/02f_walsh_annotation.R` (adds walsh_group labels, tables, annotated UMAP)
- Bershteyn: Seurat object already provided; no processing script run here.

## Canonical plot paths to share
- 15 clusters (labels 0–14): `results/walsh_day75/kres_sweep_plots/umap_by_cluster_res0.8_k30.png`
- 14 clusters (labels 0–13): `results/walsh_day75/kres_sweep_plots/umap_by_cluster_res0.8_k50.png`
- Walsh-group annotated UMAP: `results/walsh_day75/annotation_plots/umap_by_walsh_group.png`

## Notes
- All runs use R 4.1.1 + Seurat 4.1.1 (module r-seurat/4.1.1-R-4.1.1-qyci4bo); Python UMAP is not used.
- Methods-locked parameters retained: QC 1000–5000 genes, <15% MT (^MT-); LogNormalize (scale.factor=1e4); HVG=5000 (vst); regress S.Score & G2M.Score; dims=1:20 for neighbors/UMAP; clustering resolution default 2.0 (sweeps explore alternatives); stress removal clusters mean score >0.5.
- Checkpoints allow reruns without redoing QC/normalization/scaling; sweeps start from post-stress PCA checkpoint.

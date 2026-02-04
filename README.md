# mge_organoid_pipeline

Status snapshot (Feb 2026)
- Walsh: processed; Seurat objects and plots in `results/walsh_day75/`.
- Bershteyn 2025 (GSE283775): supplied Seurat object in `results/bershteyn/`.
- Xiang 2018 (GSE98201 scRNA 10x): Seurat + UMAP generated — see `results/xiang_2018/` (plots/umap_by_cluster.{png,pdf}).
- Bershteyn 2023 (GSE208672): supplied Seurat RDS plotted — see `results/bershteyn_2023/plots/umap_by_cluster.{png,pdf}`.
- Samarasinghe 2021 (GSE165577): counts downloaded; Seurat script ready; LIGER run pending package installs.
- Siebert 2026 (nemo:dat-htzat9t): no files published yet (metadata only).

Key scripts
- Downloads: `scripts/01c_download_extra_geo.sh` (Xiang_2018, Samarasinghe_2021, Bershteyn_2023)
- Xiang UMAP: `scripts/05c_xiang_2018_seurat.R`
- Bershteyn 2023 UMAP from provided RDS: `scripts/05d_bershteyn_2023_seurat_plot.R`

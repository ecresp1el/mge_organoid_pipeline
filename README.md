# mge_organoid_pipeline

Status snapshot (Feb 2026)
- Walsh (GSE250482): processed; Seurat objects and plots in `results/walsh_day75/`.
- Bershteyn 2025 (GSE283775): provided Seurat object downloaded to `data/raw/bershteyn_2025_geo_files/suppl/`; UMAP plotted to `results/bershteyn/plots/`.
- Xiang 2018 (GSE98201 scRNA 10x): Seurat + UMAP generated — see `results/xiang_2018/` (`plots/umap_by_cluster.{png,pdf}`).
- Bershteyn 2023 (GSE208672): provided Seurat object plotted + copied — see `results/bershteyn_2023/` (`plots/umap_by_cluster.{png,pdf}`).
- Samarasinghe 2021 (GSE165577): counts downloaded (`data/raw/samarasinghe_2021_geo_files/suppl/`); Seurat script/template ready; results not generated yet; LIGER run pending.
- Siebert 2026 (NeMO `nemo:dat-htzat9t`): no files published yet (metadata only).

Key scripts
- Downloads: `scripts/01c_download_extra_geo.sh` (Xiang_2018, Samarasinghe_2021, Bershteyn_2023)
- Xiang UMAP: `scripts/05c_xiang_2018_seurat.R`
- Bershteyn 2023 UMAP from provided RDS: `scripts/05d_bershteyn_2023_seurat_plot.R`
- Status audit: `scripts/00_audit_studies.sh` (prints a Markdown table)

Docs
- Workflow + directory conventions: `WORKFLOW.md`

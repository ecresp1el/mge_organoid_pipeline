# mge_organoid_pipeline

Status snapshot (Feb 2026)
- Walsh (GSE250482): processed; Seurat objects and plots in `results/walsh_day75/`.
- Bershteyn 2025 (GSE283775): canonical Seurat object + plots now in `results/bershteyn_2025/` (`bershteyn_2025_seurat.rds`, `plots/umap_by_cluster.{png,pdf}`); raw GEO object remains in `data/raw/bershteyn_2025_geo_files/suppl/`.
- Xiang 2018 (GSE98201 scRNA 10x): Seurat + UMAP generated — see `results/xiang_2018/` (`plots/umap_by_cluster.{png,pdf}`).
- Bershteyn 2023 (GSE208672): provided Seurat object plotted + copied — see `results/bershteyn_2023/` (`plots/umap_by_cluster.{png,pdf}`).
- Samarasinghe 2021 (GSE165577): counts downloaded (`data/raw/samarasinghe_2021_geo_files/suppl/`); Seurat script/template ready; results not generated yet; LIGER run pending.
- Siebert 2026 (NeMO `nemo:dat-htzat9t`): no files published yet (metadata only).

Path migration note (2026-02-08)
- Bershteyn 2025 canonical result path moved from legacy `results/bershteyn/` logic to `results/bershteyn_2025/`.
- Panel/config scripts now target `results/bershteyn_2025/bershteyn_2025_seurat.rds`.
- Legacy `results/bershteyn/plots/` files were left in place for backward compatibility because external files/jobs may still reference that older path.

Key scripts
- Downloads: `scripts/01c_download_extra_geo.sh` (Xiang_2018, Samarasinghe_2021, Bershteyn_2023)
- Xiang UMAP: `scripts/05c_xiang_2018_seurat.R`
- Bershteyn 2023 UMAP from provided RDS: `scripts/05d_bershteyn_2023_seurat_plot.R`
- Cross-study Panel B markers (figure assembly only): `scripts/06_cross_study_panelB_markers.R`
- Slurm template for Panel B (Seurat v4 runtime): `slurm_templates/06_cross_study_panelB_markers.sbatch.template`
- Slurm template for Panel B (Seurat v5 runtime / Assay5-aware): `slurm_templates/06_cross_study_panelB_markers_seurat5.sbatch.template`
- Status audit: `scripts/00_audit_studies.sh` (prints a Markdown table)

Run Panel B interactively (no Slurm required)
- Edit config template: `config/panel_b_cross_study_config.example.R`
- Optional config field: `studies$feature_map_path` for symbol->feature remapping (used for Xiang Ensembl features).
- Run:
  - `module load Bioinformatics`
  - `module load r-seurat/4.1.1-R-4.2.0-5z5hgo7`
  - `Rscript scripts/06_cross_study_panelB_markers.R --config config/panel_b_cross_study_config.example.R --project-root /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder --run-label panel_b_cross_study_v1 --detailed-log true --write-prepared true`
- Note: for Assay5 studies (for example Varela), use a Seurat v5 runtime so the script can use `LayerData` directly.
- Note: detailed logs now include assay slot inventory, assay layer inventory, feature namespace detection, cell overlap checks, and metadata structure summaries.
- Example config now includes both Bershteyn studies as separate columns:
  - `Bershteyn 2025`
  - `Bershteyn 2023`
- Plot layout is enforced programmatically before plotting:
  - Varela is always the left-most study column when present.
  - ON-target markers are assembled in the top block.
  - OFF-target markers are assembled in the bottom block.
  - Study column labels include plotted cell counts (`n=<cells>`) for quick auditing.
  - Figure subtitle includes `Plotted cells: <study>=<n> | ...` in left-to-right order.
- The script now writes a reusable input bundle for downstream scripts:
  - `panel_b_prepared_inputs.rds` (per-study UMAP coords + marker matrix + status/metadata fields)
  - This is Seurat-runtime independent and is written by default (`--write-prepared true`).
- Marker set used by Panel B (fixed order):
  - ON-target: `DCX,GAD2,DLX5,LHX6,MAF,SST,LHX8,SP8`
  - OFF-target: `PAX6,NEUROD2,ISL1,ACHE`
- Outputs:
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_cross_study_markers.png`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_cross_study_markers.pdf`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_cross_study_markers.svg`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_study_status.tsv`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_marker_presence.tsv`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_assay_slot_summary.tsv`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_reduction_summary.tsv`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_metadata_summary.tsv`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_metadata_columns.tsv`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_ident_counts.tsv`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_feature_space.tsv`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_prepared_inputs.rds`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_row_summary.tsv`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_issues.tsv`

How to interpret marker availability
- `panel_b_marker_presence.tsv` is the per-study, per-gene truth table.
- A row is `present=TRUE` when the marker exists in the study feature namespace and expression could be pulled for plotting.
- A row is `present=FALSE` with reason:
  - `Missing object`: study file is missing.
  - `Cell mismatch`: expression matrix cells do not match UMAP cells.
  - `Gene not found`: marker symbol absent from that study feature namespace.
  - `Gene not found (feature IDs are Ensembl-like)`: study uses Ensembl IDs, so symbol markers need mapping.
- `panel_b_study_status.tsv` now includes:
  - `feature_map_path`, `feature_map_resolved`, `feature_map_source`, `n_feature_map_pairs`
  - `n_marker_genes_mapped`, `marker_gene_feature_map`

Xiang Ensembl ID note
- Xiang (`results/xiang_2018/xiang_2018_seurat.rds`) stores features as Ensembl-like IDs (for example `ENSG...`), not symbols.
- The Panel B script now supports dynamic symbol remapping:
  - markers come from `GENE_ORDER` (so add/remove markers freely in code)
  - per-study optional `feature_map_path` in config is used to map symbols to feature IDs
  - mappings are applied automatically at runtime (no hardcoded per-marker logic in the extraction path)
- Because Panel B markers are symbols (`DCX`, `LHX6`, etc.), Xiang needs symbol-to-Ensembl mapping to avoid false "Gene not found" rows.
- See `panel_b_feature_space.tsv` for `feature_id_type` and `feature_id_examples`.
- Current symbol-to-Ensembl mappings for the Panel B markers (verified against Xiang object):
  - `DCX` -> `ENSG00000077279`
  - `GAD2` -> `ENSG00000136750`
  - `DLX5` -> `ENSG00000105880`
  - `LHX6` -> `ENSG00000106852`
  - `MAF` -> `ENSG00000178573`
  - `SST` -> `ENSG00000157005`
  - `LHX8` -> `ENSG00000162624`
  - `SP8` -> `ENSG00000164651`
  - `PAX6` -> `ENSG00000007372`
  - `NEUROD2` -> `ENSG00000171532`
  - `ISL1` -> `ENSG00000016082`
  - `ACHE` -> `ENSG00000087085`

Run Panel B on Slurm
- Seurat v4 runtime template:
  - `cp slurm_templates/06_cross_study_panelB_markers.sbatch.template /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/06_cross_study_panelB_markers.sbatch`
  - `sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/06_cross_study_panelB_markers.sbatch`
- Seurat v5 runtime template (recommended when config includes Assay5 objects):
  - `cp slurm_templates/06_cross_study_panelB_markers_seurat5.sbatch.template /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/06_cross_study_panelB_markers_seurat5.sbatch`
  - `export SEURAT5_MODULE=<your-seurat-v5-module>`
  - `sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/06_cross_study_panelB_markers_seurat5.sbatch`

Docs
- Workflow + directory conventions: `WORKFLOW.md`

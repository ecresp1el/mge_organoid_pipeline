# mge_organoid_pipeline

Status snapshot (Feb 2026)
- Walsh (GSE250482): processed; Seurat objects and plots in `results/walsh_day75/`.
- Bershteyn 2025 (GSE283775): canonical Seurat object + plots now in `results/bershteyn_2025/` (`bershteyn_2025_seurat.rds`, `plots/umap_by_cluster.{png,pdf}`); raw GEO object remains in `data/raw/bershteyn_2025_geo_files/suppl/`.
- Xiang 2018 (GSE98201 scRNA 10x): Seurat + UMAP generated — see `results/xiang_2018/` (`plots/umap_by_cluster.{png,pdf}`).
- Bershteyn 2023 (GSE208672): provided Seurat object plotted + copied — see `results/bershteyn_2023/` (`plots/umap_by_cluster.{png,pdf}`).
- Samarasinghe 2021 (GSE165577): counts downloaded (`data/raw/samarasinghe_2021_geo_files/suppl/`); Seurat script/template ready; results not generated yet; LIGER run pending.
- Siebert 2026 (NeMO `nemo:dat-htzat9t`): canonical Seurat landing path is `results/siebert_2026/siebert_2026_seurat.rds` (with NeMO metadata under `data/raw/siebert_2026_nemo/`).
- He et al HNOCA full V2 (Zenodo `14160929`): raw source path `data/raw/he_et_al_zenodo/suppl/hnoca_allmeta.h5ad`; SCN8A-scoped Seurat output path `results/he_et_al/he_et_al_scn8a_seurat.rds`.

Register a cleaned Siebert `.rds` (from local laptop to Great Lakes)
- Create canonical destination directory on Great Lakes:
  - `ssh elcrespo@gl-login1.arc-ts.umich.edu "mkdir -p /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siebert_2026"`
- Copy local `.rds` into canonical study path (run from local terminal):
  - `rsync -avh --progress "/Users/elcrespo/Downloads/<your_file>.rds" "elcrespo@gl-login1.arc-ts.umich.edu:/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siebert_2026/siebert_2026_seurat.rds"`
- Optional verify on Great Lakes:
  - `ls -lh /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siebert_2026/siebert_2026_seurat.rds`

Path migration note (2026-02-08)
- Bershteyn 2025 canonical result path moved from legacy `results/bershteyn/` logic to `results/bershteyn_2025/`.
- Panel/config scripts now target `results/bershteyn_2025/bershteyn_2025_seurat.rds`.
- Legacy `results/bershteyn/plots/` files were left in place for backward compatibility because external files/jobs may still reference that older path.

Key scripts
- Downloads: `scripts/01c_download_extra_geo.sh` (Xiang_2018, Samarasinghe_2021, Bershteyn_2023)
- He et al full V2 download: `scripts/01d_download_he_et_al_zenodo.sh`
- He et al SCN8A slice extraction (Python): `scripts/05e_extract_he_et_al_scn8a_slice.py`
- He et al SCN8A Seurat build: `scripts/05f_he_et_al_scn8a_seurat.R`
- Cross-study multi-gene panel (log1p): `scripts/06_cross_study_gene_panel_log.R`
- He-vs-Varela gene-panel config (LHX6, NKX2.1): `config/gene_panel_he_vs_varela_config.example.R`
- He LHX6/NKX2.1 extraction Slurm template: `slurm_templates/05e_extract_he_et_al_lhx6_nkx21_slice.sbatch.template`
- He-vs-Varela gene-panel Slurm template: `slurm_templates/06_cross_study_gene_panel_he_vs_varela_seurat5.sbatch.template`
- He-vs-Varela completion checker: `scripts/00b_check_he_vs_varela_gene_panel.sh`
- Xiang UMAP: `scripts/05c_xiang_2018_seurat.R`
- Bershteyn 2023 UMAP from provided RDS: `scripts/05d_bershteyn_2023_seurat_plot.R`
- Cross-study Panel B markers (figure assembly only): `scripts/06_cross_study_panelB_markers.R`
- SCN8A He-vs-Varela config: `config/scn8a_he_vs_varela_config.example.R`
- SCN8A He-vs-Varela Slurm template: `slurm_templates/06_cross_study_scn8a_log_he_vs_varela_seurat5.sbatch.template`
- Slurm template for Panel B (Seurat v4 runtime): `slurm_templates/06_cross_study_panelB_markers.sbatch.template`
- Slurm template for Panel B (Seurat v5 runtime / Assay5-aware): `slurm_templates/06_cross_study_panelB_markers_seurat5.sbatch.template`
- Human GE comparison placeholder (next stage draft): `scripts/07_compare_human_developing_ge_tbd.R`
- Slurm template for Human GE placeholder: `slurm_templates/07_compare_human_developing_ge_tbd.sbatch.template`
- Status audit: `scripts/00_audit_studies.sh` (prints a Markdown table)

He et al full V2 SCN8A workflow (additive; existing stages unchanged)
- Stage 1 (download full h5ad, 48.9 GB):
  - `cp slurm_templates/01d_download_he_et_al_zenodo.sbatch.template /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/01d_download_he_et_al_zenodo.sbatch`
  - `sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/01d_download_he_et_al_zenodo.sbatch`
- Stage 2 (extract SCN8A + UMAP slice from full h5ad):
  - `cp slurm_templates/05e_extract_he_et_al_scn8a_slice.sbatch.template /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/05e_extract_he_et_al_scn8a_slice.sbatch`
  - `sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/05e_extract_he_et_al_scn8a_slice.sbatch`
- Stage 3 (build SCN8A-ready Seurat object):
  - `cp slurm_templates/05f_he_et_al_scn8a_seurat.sbatch.template /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/05f_he_et_al_scn8a_seurat.sbatch`
  - `sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/05f_he_et_al_scn8a_seurat.sbatch`
- Stage 4 (run SCN8A panel: He et al vs Varela):
  - `cp slurm_templates/06_cross_study_scn8a_log_he_vs_varela_seurat5.sbatch.template /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/06_cross_study_scn8a_log_he_vs_varela_seurat5.sbatch`
  - `sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/06_cross_study_scn8a_log_he_vs_varela_seurat5.sbatch`
- Output:
  - `PROJECT_ROOT/results/panel_b_scn8a_he_vs_varela_v1/plots/scn8a_cross_study_log.png`

He-vs-Varela marker panel for `LHX6` and `NKX2.1` (grey->purple)
- Stage A (extract required genes from full h5ad):
  - `cp slurm_templates/05e_extract_he_et_al_lhx6_nkx21_slice.sbatch.template /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/05e_extract_he_et_al_lhx6_nkx21_slice.sbatch`
  - `sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/05e_extract_he_et_al_lhx6_nkx21_slice.sbatch`
- Stage B (build He Seurat from extracted slice):
  - `cp slurm_templates/05f_he_et_al_scn8a_seurat.sbatch.template /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/05f_he_et_al_scn8a_seurat.sbatch`
  - `sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/05f_he_et_al_scn8a_seurat.sbatch`
- Stage C (render cross-study marker panel):
  - `cp slurm_templates/06_cross_study_gene_panel_he_vs_varela_seurat5.sbatch.template /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/06_cross_study_gene_panel_he_vs_varela_seurat5.sbatch`
  - `sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/06_cross_study_gene_panel_he_vs_varela_seurat5.sbatch`
- Completion check:
  - `scripts/00b_check_he_vs_varela_gene_panel.sh --project-root /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder`
- Output image:
  - `PROJECT_ROOT/results/gene_panel_he_vs_varela_v1/plots/gene_panel_cross_study_log.png`

Run Panel B interactively (no Slurm required)
- Edit config template: `config/panel_b_cross_study_config.example.R`
- Optional config field: `studies$feature_map_path` for symbol->feature remapping (used for Xiang Ensembl features).
- Run:
  - `module load Bioinformatics`
  - `module load r-seurat/5.1.0-R-4.4.1-c3m7yfq`
  - `Rscript scripts/06_cross_study_panelB_markers.R --config config/panel_b_cross_study_config.example.R --project-root /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder --run-label panel_b_cross_study_v1_seurat5 --detailed-log true --write-prepared true --write-study-objects true --prepared-objects-root results/panel_b_prepared_objects`
- Note: Assay5 studies (for example both Varela objects) require a Seurat v5 runtime so the script can use `LayerData` directly.
- Note: detailed logs now include assay slot inventory, assay layer inventory, feature namespace detection, cell overlap checks, and metadata structure summaries.
- Current config (`config/panel_b_cross_study_config.example.R`) includes both required Varela study object paths:
  - `varela_this_paper`: `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/varela_this_paper/varela_this_paper_seurat.rds`
  - `varela_div90`: `/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/umap_props_output/clustered_day90_with_cluster_names_2.rds`
- Plot layout is enforced programmatically before plotting:
  - Studies are rendered as rows (top->bottom), genes are rendered as columns (left->right).
  - The first two study rows are required and ordered: `varela_this_paper` (DIV30), then `varela_div90` (DIV90).
  - ON-target markers are assembled in the top block and OFF-target markers in the bottom block.
  - Study row labels include plotted cell counts (`n=<cells>`) for quick auditing.
  - Figure subtitle includes `Plotted cells: <study>=<n> | ...` in top-to-bottom study order.
- The script now writes a reusable input bundle for downstream scripts:
  - `panel_b_prepared_inputs.rds` (per-study UMAP coords + marker matrix + status/metadata fields)
  - This is Seurat-runtime independent and is written by default (`--write-prepared true`).
- Panel B PNG exports are now written at high print quality (`dpi=600`).
- In addition to the combined ON+OFF figure, the script now saves separate ON-only and OFF-only figures with matching formatting.
- The script now also publishes validated per-study Seurat objects to a stable path for downstream pipeline stages:
  - `results/panel_b_prepared_objects/studies/<study_id>_panelb_ready_seurat.rds`
  - `results/panel_b_prepared_objects/panel_b_prepared_object_paths.tsv` (canonical manifest)
  - `<run_label>/plots/panel_b_prepared_object_paths.tsv` (run-scoped manifest snapshot)
- Marker set used by Panel B (fixed order):
  - ON-target: `DCX,GAD2,DLX5,LHX6,MAF,SST,ERBB4,MEF2C,MAFB,LHX8,NKX2-1`
  - OFF-target: `SP8,PAX6,NEUROD2,ISL1,ACHE,NKX6-2,MKI67`
- Outputs:
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_cross_study_markers.png`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_cross_study_markers.pdf`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_cross_study_markers.svg`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_cross_study_markers_on_target.png`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_cross_study_markers_on_target.pdf`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_cross_study_markers_on_target.svg`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_cross_study_markers_off_target.png`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_cross_study_markers_off_target.pdf`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_cross_study_markers_off_target.svg`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_study_status.tsv`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_marker_presence.tsv`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_assay_slot_summary.tsv`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_reduction_summary.tsv`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_metadata_summary.tsv`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_metadata_columns.tsv`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_ident_counts.tsv`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_feature_space.tsv`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_prepared_inputs.rds`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_prepared_object_paths.tsv`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_row_summary.tsv`
  - `PROJECT_ROOT/results/<run_label>/plots/panel_b_issues.tsv`
  - `PROJECT_ROOT/results/panel_b_prepared_objects/studies/<study_id>_panelb_ready_seurat.rds`
  - `PROJECT_ROOT/results/panel_b_prepared_objects/panel_b_prepared_object_paths.tsv`

Next stage placeholder: Human developing GE comparison
- Edit config template: `config/human_ge_comparison_tbd_config.example.R`
- Run:
  - `module load Bioinformatics`
  - `module load r-seurat/5.1.0-R-4.4.1-c3m7yfq`
  - `Rscript scripts/07_compare_human_developing_ge_tbd.R --config config/human_ge_comparison_tbd_config.example.R --project-root /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder --run-label human_ge_comparison_tbd_v1`
- Outputs:
  - `PROJECT_ROOT/results/<run_label>/human_ge_comparison_tbd/human_ge_input_manifest.tsv`
  - `PROJECT_ROOT/results/<run_label>/human_ge_comparison_tbd/human_ge_missing_required_studies.tsv`
  - `PROJECT_ROOT/results/<run_label>/human_ge_comparison_tbd/human_ge_comparison_tbd_notes.txt`

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
- Current symbol-to-Ensembl mappings for a subset of Panel B markers (verified against Xiang object):
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
- Additional Panel B markers (`ERBB4`, `MEF2C`, `MAFB`, `NKX2-1`, `NKX6-2`, `MKI67`) are remapped dynamically at runtime when present in the configured feature map; audit run-specific remaps in `panel_b_study_status.tsv` via `marker_gene_feature_map`.

Run Panel B on Slurm
- Seurat v4 runtime template:
  - `cp slurm_templates/06_cross_study_panelB_markers.sbatch.template /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/06_cross_study_panelB_markers.sbatch`
  - `sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/06_cross_study_panelB_markers.sbatch`
- Seurat v5 runtime template (recommended when config includes Assay5 objects):
  - `cp slurm_templates/06_cross_study_panelB_markers_seurat5.sbatch.template /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/06_cross_study_panelB_markers_seurat5.sbatch`
  - `export SEURAT5_MODULE=<your-seurat-v5-module>`
  - `sbatch /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/jobs/06_cross_study_panelB_markers_seurat5.sbatch`

Copy Panel B outputs to local machine (run from local terminal, not Great Lakes)
- Create local destination:
  - `mkdir -p "/Users/ecrespo/Downloads/panel_b_cross_study_v1_seurat5_plots"`
- Copy PNG/PDF/SVG outputs for the current Seurat v5 run label:
  - `rsync -avh --progress --prune-empty-dirs --include='*/' --include='*.png' --include='*.pdf' --include='*.svg' --exclude='*' "elcrespo@gl-login1.arc-ts.umich.edu:/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/panel_b_cross_study_v1_seurat5/plots/" "/Users/ecrespo/Downloads/panel_b_cross_study_v1_seurat5_plots/"`
- Tip:
  - replace `panel_b_cross_study_v1_seurat5` with any other run label you want to pull.

Docs
- Workflow + directory conventions: `WORKFLOW.md`

# Walsh & Bershteyn Single-cell Repro Workflow (Great Lakes)

## Roots
- Repo (code/config only): `/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline`
- Runtime workspace (data/jobs/results): `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder` (PROJECT_ROOT)

## Directory Conventions (PROJECT_ROOT)
Canonical layout for each study:

```
PROJECT_ROOT/
  data/
    raw/
      <study>_geo_files/          # GEO downloads (series matrix, SOFT, MINiML, supp files)
        matrix/
        miniml/
        soft/
        suppl/
        download_manifest.tsv
      <study>_nemo/               # NeMO landing-page/BDBag downloads
        metadata.json
        manifest_download.tsv
    processed/
      <study_or_run_label>/       # Optional: intermediate extracts (e.g., raw tar extraction)
  results/
    <study_or_run_label>/
      <study_or_run_label>_seurat.rds   # If generated/copied
      plots/
        umap_by_cluster.png
        umap_by_cluster.pdf
      checkpoints/                # Optional (Walsh pipeline)
  logs/                           # Slurm logs
  jobs/                           # Rendered sbatch files actually submitted
```

## Key Data Locations
- Raw GEO downloads:
  - Walsh: `PROJECT_ROOT/data/raw/walsh_2025_geo_files/`
  - Bershteyn: `PROJECT_ROOT/data/raw/bershteyn_2025_geo_files/`
- Raw GEO downloads (additional):
  - Xiang_2018 (GSE97882; scRNA SubSeries GSE98201 10x): `PROJECT_ROOT/data/raw/xiang_2018_geo_files/`
    - 10x matrix/genes/barcodes: `.../suppl/GSE98201_{matrix,genes,barcodes}.tsv.gz`
  - Bershteyn_2023 (GSE208672; supplied Seurat RDS): `PROJECT_ROOT/data/raw/bershteyn_2023_geo_files/`
  - Samarasinghe_2021 (GSE165577; counts CSVs): `PROJECT_ROOT/data/raw/samarasinghe_2021_geo_files/`
  - Shi_2019 (GSE135827; GEO raw-count matrix + week suffix barcodes): `PROJECT_ROOT/data/raw/shi_2019_geo_files/`
  - Liu_2025_hnbMO (GSE286235; Cell Ranger raw H5s): `PROJECT_ROOT/data/raw/liu_2025_hnbmo_geo_files/`
    - Healthy samples only by default: GSM8721440, GSM8721441, GSM8721442; DS sample GSM8721443 retained for optional runs.
  - He et al HNOCA full V2 (Zenodo 14160929): `PROJECT_ROOT/data/raw/he_et_al_zenodo/suppl/hnoca_allmeta.h5ad`
- NeMO landing-page downloads (Siebert 2026, `nemo:dat-htzat9t`):
  - `PROJECT_ROOT/data/raw/siebert_2026_nemo/` (metadata, BDBags, subset fetch + downloads)
  - Canonical cleaned Seurat landing path: `PROJECT_ROOT/results/siebert_2026/siebert_2026_seurat.rds`
- Processed / results:
  - Walsh: `PROJECT_ROOT/results/walsh_day75/`
  - Bershteyn_2025: `PROJECT_ROOT/results/bershteyn_2025/` (canonical Seurat copy + UMAP plots)
  - Xiang_2018: `PROJECT_ROOT/results/xiang_2018/` (Seurat + UMAP generated from GSE98201 10x)
  - Bershteyn_2023: `PROJECT_ROOT/results/bershteyn_2023/` (provided Seurat object copied + UMAP plotted)
  - Samarasinghe_2021: **not generated yet** (expected: `PROJECT_ROOT/results/samarasinghe_2021/`)
  - Samarasinghe_2021 LIGER (paper-matched): **not generated yet** (expected: `PROJECT_ROOT/results/samarasinghe_2021_liger/`)
  - Shi_2019: `PROJECT_ROOT/results/shi_2019/` (standalone Seurat + UMAP built from GSE135827 GEO count table)
  - Liu_2025_hnbMO: `PROJECT_ROOT/results/liu_2025_hnbmo/` (healthy-only Seurat + UMAP built from GSE286235 Cell Ranger raw H5s)
  - Siebert_2026: `PROJECT_ROOT/results/siebert_2026/` (canonical cleaned Seurat copy path)
  - Varela (this paper): `PROJECT_ROOT/results/varela_this_paper/varela_this_paper_seurat.rds` (true copy of `/nfs/turbo/umms-parent/mgeo_scRNAseq/day30_old/Day30.rds`)
  - He et al SCN8A slice (intermediate): `PROJECT_ROOT/data/processed/he_et_al_scn8a_slice/`
  - He et al SCN8A-ready Seurat: `PROJECT_ROOT/results/he_et_al/he_et_al_scn8a_seurat.rds`
- Logs: `PROJECT_ROOT/logs/`
- Job scripts: `PROJECT_ROOT/jobs/`

## Main Seurat Objects / Plots
- Walsh final (post-stress, annotated): `results/walsh_day75/walsh_day75_final_annotated.rds`
- Walsh final (post-stress, unannotated): `results/walsh_day75/walsh_day75_final.rds`
- Bershteyn 2025 (GSE283775) canonical Seurat: `results/bershteyn_2025/bershteyn_2025_seurat.rds`
  - Raw source (downloaded): `data/raw/bershteyn_2025_geo_files/suppl/GSE283775_Seurat_scRNA_seq.rds.gz`
  - UMAP: `results/bershteyn_2025/plots/umap_by_cluster.{png,pdf}`
- Bershteyn 2023 (GSE208672) Seurat: `results/bershteyn_2023/bershteyn_2023_seurat.rds`; UMAP `results/bershteyn_2023/plots/umap_by_cluster.{png,pdf}`
- Xiang 2018 Seurat (GSE98201 10x): `results/xiang_2018/xiang_2018_seurat.rds`; UMAP `results/xiang_2018/plots/umap_by_cluster.{png,pdf}`
- Samarasinghe 2021 (GSE165577) expected outputs (not generated yet):
  - Seurat: `results/samarasinghe_2021/samarasinghe_2021_seurat.rds`
  - UMAP: `results/samarasinghe_2021/plots/umap_by_cluster.{png,pdf}`
- Shi 2019 (GSE135827): `results/shi_2019/shi_2019_seurat.rds`; UMAP `results/shi_2019/plots/umap_by_cluster.{png,pdf}` and week UMAP `results/shi_2019/plots/umap_by_week.{png,pdf}`
- Liu 2025 hnbMO (GSE286235) expected outputs:
  - Healthy-only Seurat: `results/liu_2025_hnbmo/liu_2025_hnbmo_healthy_seurat.rds`
  - UMAP: `results/liu_2025_hnbmo/plots/umap_by_cluster.{png,pdf}`
- Siebert 2026 canonical Seurat path: `results/siebert_2026/siebert_2026_seurat.rds` (UMAP plot may be generated later under `results/siebert_2026/plots/`)
- He et al full V2 SCN8A-ready Seurat: `results/he_et_al/he_et_al_scn8a_seurat.rds`; sanity plot `results/he_et_al/plots/scn8a_umap.{png,pdf}`

## Study Status Audit (2026-02-04)

To regenerate this table from the current state of `PROJECT_ROOT`, run:
`scripts/00_audit_studies.sh --project-root /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder`

| Study | Raw inputs available | Seurat object available | UMAP PNG available | Notes |
| --- | --- | --- | --- | --- |
| Walsh (GSE250482) | Yes (GEO RAW tar) | Yes (`results/walsh_day75/walsh_day75_final*.rds`) | Yes (`results/walsh_day75/plots/umap_by_cluster.png`) | Run label is `walsh_day75` (not `walsh_2025`). |
| Bershteyn 2025 (GSE283775) | Yes (provided Seurat `.rds.gz`) | Yes (`results/bershteyn_2025/bershteyn_2025_seurat.rds`) | Yes (`results/bershteyn_2025/plots/umap_by_cluster.png`) | Canonical location is now `results/bershteyn_2025`; legacy `results/bershteyn/plots` retained for compatibility. |
| Xiang 2018 (GSE98201) | Yes (10x matrix/genes/barcodes) | Yes (`results/xiang_2018/xiang_2018_seurat.rds`) | Yes (`results/xiang_2018/plots/umap_by_cluster.png`) | Built from GSE98201 trio under `.../suppl/`. |
| Bershteyn 2023 (GSE208672) | Yes (provided Seurat `.rds.gz`) | Yes (`results/bershteyn_2023/bershteyn_2023_seurat.rds`) | Yes (`results/bershteyn_2023/plots/umap_by_cluster.png`) | Raw `suppl/` contains extra variants (`.rds`, `.gz2`) from download/read debugging. |
| Samarasinghe 2021 (GSE165577) | Yes (filtered/normalized counts CSVs) | Not yet (script ready) | Not yet | Needs Seurat run (`scripts/05_samarasinghe_2021_seurat.R`). LIGER step pending. |
| Shi 2019 (GSE135827) | Expected at `data/raw/shi_2019_geo_files/suppl/GSE135827_GE_mat_raw_count_with_week_info.txt.gz` | Expected at `results/shi_2019/shi_2019_seurat.rds` | Expected at `results/shi_2019/plots/umap_by_cluster.png` | Standalone GEO workflow with per-cell week metadata from barcode suffix (`-GWxx`). |
| Liu 2025 hnbMO (GSE286235) | Expected at `data/raw/liu_2025_hnbmo_geo_files/suppl/GSM872144*_raw_feature_bc_matrix.h5` | Expected at `results/liu_2025_hnbmo/liu_2025_hnbmo_healthy_seurat.rds` | Expected at `results/liu_2025_hnbmo/plots/umap_by_cluster.png` | Uses healthy samples only by default; raw Cell Ranger H5s require explicit cell filtering. |
| Varela (this paper, Day30) | Yes (`/nfs/turbo/umms-parent/mgeo_scRNAseq/day30_old/Day30.rds`) | Yes (`results/varela_this_paper/varela_this_paper_seurat.rds`) | Included in Panel B assembly | Canonical Panel B path is under `PROJECT_ROOT/results/...` and now uses a true copied file. |
| Siebert 2026 (NeMO `nemo:dat-htzat9t`) | Metadata available (`data/raw/siebert_2026_nemo/`) | Canonical path: `results/siebert_2026/siebert_2026_seurat.rds` | Optional (when plotted) | Cleaned Seurat object should be copied to canonical results path for pipeline use. |
| He et al HNOCA full V2 (Zenodo 14160929) | Yes (`data/raw/he_et_al_zenodo/suppl/hnoca_allmeta.h5ad`) | Yes (`results/he_et_al/he_et_al_scn8a_seurat.rds`) | Yes (`results/he_et_al/plots/scn8a_umap.png`) | SCN8A-only extracted slice for cross-study SCN8A panel. |

## Known Deviations / Cleanup Items
- If you ever see `PROJECT_ROOT/data/raw/{matrix,miniml,soft,suppl}/` at the top level (they should not exist), they are empty leftovers from an earlier download script bug; safe to delete.
- Bershteyn 2025 canonical output path is `results/bershteyn_2025/`; legacy `results/bershteyn/plots/` is intentionally left in place for backward compatibility with older references.
- Bershteyn 2023 raw `suppl/` contains multiple variants (`.rds`, `.rds.gz`, `.rds.gz2`) from download/read debugging; the canonical copy is `results/bershteyn_2023/bershteyn_2023_seurat.rds`.

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
- Audit current study status: `scripts/00_audit_studies.sh` (prints a Markdown table)
- NeMO Siebert landing page/BDBag: `scripts/01b_download_siebert_nemo.sh`
- Extra GEO downloads (Xiang_2018, Samarasinghe_2021, Bershteyn_2023): `scripts/01c_download_extra_geo.sh`
- He et al full V2 download (Zenodo 14160929): `scripts/01d_download_he_et_al_zenodo.sh`
- Shi et al GEO download (GSE135827): `scripts/01e_download_shi_geo.sh`
- Liu 2025 hnbMO GEO download (GSE286235): `scripts/01f_download_gse286235_geo.sh`
- Samarasinghe 2021 Seurat/UMAP: `scripts/05_samarasinghe_2021_seurat.R`
- Samarasinghe 2021 LIGER (paper settings): `scripts/05b_samarasinghe_2021_liger.R`
- Xiang 2018 scRNA 10x (GSE98201) Seurat/UMAP: `scripts/05c_xiang_2018_seurat.R`
- Bershteyn 2023 Seurat plotting: `scripts/05d_bershteyn_2023_seurat_plot.R`
- He et al SCN8A extraction from full h5ad: `scripts/05e_extract_he_et_al_scn8a_slice.py`
- He et al SCN8A Seurat build: `scripts/05f_he_et_al_scn8a_seurat.R`
- Shi et al standalone Seurat + UMAP: `scripts/05g_shi_2019_seurat.R`
- Shi Table S3 xlsx->tsv conversion (no extra deps): `scripts/05h_shi_table_s3_xlsx_to_tsv.py`
- Shi Table S3-based cluster annotation mapping: `scripts/05i_shi_2019_annotate_from_table_s3.R`
- Liu 2025 hnbMO healthy-only Seurat + UMAP: `scripts/05j_liu_2025_hnbmo_seurat.R`
- Cross-study marker Panel B figure assembly (no analysis/recompute): `scripts/06_cross_study_panelB_markers.R`
- Cross-study multi-gene UMAP panel (log1p): `scripts/06_cross_study_gene_panel_log.R`
- Cross-study SCN8A He-vs-Varela config: `config/scn8a_he_vs_varela_config.example.R`
- Cross-study LHX6/NKX2.1 He-vs-Varela config: `config/gene_panel_he_vs_varela_config.example.R`
- Slurm template for SCN8A He-vs-Varela (Seurat v5): `slurm_templates/06_cross_study_scn8a_log_he_vs_varela_seurat5.sbatch.template`
- Slurm template for He LHX6/NKX2.1 extraction: `slurm_templates/05e_extract_he_et_al_lhx6_nkx21_slice.sbatch.template`
- Slurm template for LHX6/NKX2.1 He-vs-Varela panel (Seurat v5): `slurm_templates/06_cross_study_gene_panel_he_vs_varela_seurat5.sbatch.template`
- Completion checker for He-vs-Varela gene panel: `scripts/00b_check_he_vs_varela_gene_panel.sh`
- Slurm template for cross-study Panel B (Seurat v4 runtime): `slurm_templates/06_cross_study_panelB_markers.sbatch.template`
- Slurm template for cross-study Panel B (Seurat v5 runtime / Assay5-aware): `slurm_templates/06_cross_study_panelB_markers_seurat5.sbatch.template`
- Slurm template for Shi Table S3 xlsx->tsv conversion: `slurm_templates/05h_shi_table_s3_xlsx_to_tsv.sbatch.template`
- Slurm template for Shi Table S3 annotation mapping: `slurm_templates/05i_shi_2019_annotate_from_table_s3.sbatch.template`
- Human developing GE comparison placeholder (next stage draft): `scripts/07_compare_human_developing_ge_tbd.R`
- Slurm template for human GE placeholder: `slurm_templates/07_compare_human_developing_ge_tbd.sbatch.template`

- Main pipeline: `scripts/02_walsh_day75_seurat.R` (methods-locked QC/normalize/HVG/CC/Scale/PCA/stress removal/UMAP/clustering; checkpoints; Elbow fallback)
- Sweeps: `scripts/02b_walsh_k_sweep.R`, `02c_walsh_dims_sweep.R`, `02d_walsh_resolution_sweep.R`, `02e_walsh_kres_sweep.R`
- Annotation: `scripts/02f_walsh_annotation.R` (adds walsh_group labels, tables, annotated UMAP)

## Register Siebert RDS
Run from your local terminal (Mac/Linux):
- `ssh elcrespo@gl-login1.arc-ts.umich.edu "mkdir -p /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siebert_2026"`
- `rsync -avh --progress "/Users/elcrespo/Downloads/<your_file>.rds" "elcrespo@gl-login1.arc-ts.umich.edu:/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siebert_2026/siebert_2026_seurat.rds"`
- Verify on Great Lakes:
  - `ls -lh /nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siebert_2026/siebert_2026_seurat.rds`

## Centralized Table: Final Seurat Object Paths

| Study                  | Seurat Object Path                                                                                      | Script/Step Producing It                        | UMAP Included? | Notes                                      |
|------------------------|--------------------------------------------------------------------------------------------------------|-------------------------------------------------|---------------|---------------------------------------------|
| Walsh (GSE250482)      | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/walsh_day75/walsh_day75_final.rds`  | scripts/02_walsh_day75_seurat.R                 | Yes           | Also: `walsh_day75_final_annotated.rds`     |
| Bershteyn 2025         | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/bershteyn_2025/bershteyn_2025_seurat.rds` | Downloaded, then processed/checked              | Yes           | Provided as .rds.gz, copied to results      |
| Bershteyn 2023         | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/bershteyn_2023/bershteyn_2023_seurat.rds` | Downloaded, then processed/checked              | Yes           | Provided as .rds.gz, copied to results      |
| Xiang 2018             | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/xiang_2018/xiang_2018_seurat.rds`   | scripts/05c_xiang_2018_seurat.R                 | Yes           | Built from 10x matrix/genes/barcodes        |
| Siebert 2026           | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/siebert_2026/siebert_2026_seurat.rds` | Registered/copy step                            | Yes           | Canonical cleaned Seurat object             |
| Varela (this paper, DIV30) | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/varela_this_paper/varela_this_paper_seurat.rds` | Copied from `/nfs/turbo/umms-parent/mgeo_scRNAseq/day30_old/Day30.rds` | Yes           | Used in Panel B assembly                    |
| Varela (DIV90)         | `/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/umap_props_output/clustered_day90_with_cluster_names_2.rds` | External analysis output | Yes | DIV90 timepoint, ventral SOSRS             |
| He et al SCN8A         | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/he_et_al/he_et_al_scn8a_seurat.rds`  | scripts/05f_he_et_al_scn8a_seurat.R             | Yes           | Built from extracted slice                  |
| Shi 2019 (GSE135827)   | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_2019/shi_2019_seurat.rds` | scripts/05g_shi_2019_seurat.R | Yes | Standalone GEO matrix workflow; metadata exports under `results/shi_2019/` |
| Shi 2019 (Table S3-annotated) | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/shi_2019_paper_qc/table_s3_annotation/shi_2019_seurat_annotated_table_s3_res0_11.rds` | scripts/05i_shi_2019_annotate_from_table_s3.R | Yes | Cluster labels mapped to paper major cell types via Table S3 DEGs |
| Liu 2025 hnbMO (GSE286235) | `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/liu_2025_hnbmo/liu_2025_hnbmo_healthy_seurat.rds` | scripts/05j_liu_2025_hnbmo_seurat.R | Yes | Healthy-only by default; DS sample optional with `--include-ds true` |
| Samarasinghe 2021      | *(Not present yet)*                                                                                    | scripts/05_samarasinghe_2021_seurat.R           | No            | Script ready, output not found              |

---

## Workflow Summary & R Handover Points

- Each study’s pipeline produces a final `.rds` Seurat object in the results directory above.
- UMAPs are generated and embedded in these objects by the R scripts listed.
- Once the final `.rds` is present, you can load it in R for downstream analysis and plotting—no need to rerun earlier steps.
- Comparative plots and cross-study analyses use these `.rds` files as input.

**Example R usage:**
```R
seu <- readRDS("/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/walsh_day75/walsh_day75_final.rds")
```

---

## Shared vs. Unique Features Across Studies

**Shared features:**
- All Seurat objects contain UMAP embeddings, cluster assignments, and basic metadata (cell/sample IDs, domains, etc.).

**Unique features:**
- Some studies include additional metadata (e.g., stress scores, domain labels, custom annotations).
- See the scripts and config files for details per study.

---

## Pipeline Workflow: Inputs, Outputs, and Scripts

Below is a high-level workflow for each major study, with scripts and key input/output files:

### Walsh (GSE250482)
- **Script:** `scripts/02_walsh_day75_seurat.R`
- **Inputs:** Raw GEO files, gene lists (hypoxia/glycolysis)
- **Outputs:** `walsh_day75_final.rds`, UMAP plots, annotation tables
- **R handover:** After `walsh_day75_final.rds` is created, all downstream analysis/plotting can be done in R

### Bershteyn 2025/2023
- **Script:** Download/copy, then processed/checked
- **Inputs:** Provided Seurat `.rds.gz` files
- **Outputs:** `bershteyn_2025_seurat.rds`, `bershteyn_2023_seurat.rds`, UMAP plots
- **R handover:** After `.rds` is copied to results, ready for R analysis

### Xiang 2018
- **Script:** `scripts/05c_xiang_2018_seurat.R`
- **Inputs:** 10x matrix/genes/barcodes
- **Outputs:** `xiang_2018_seurat.rds`, UMAP plots
- **R handover:** After `.rds` is created, ready for R analysis

### Siebert 2026
- **Script:** Registered/copy step
- **Inputs:** Cleaned Seurat object
- **Outputs:** `siebert_2026_seurat.rds`
- **R handover:** After `.rds` is copied, ready for R analysis

### Varela (this paper)
- **Script:** Copied from legacy path
- **Inputs:** `/nfs/turbo/umms-parent/mgeo_scRNAseq/day30_old/Day30.rds`
- **Outputs:** `varela_this_paper_seurat.rds`
- **R handover:** After `.rds` is copied, ready for R analysis

### He et al SCN8A
- **Script:** `scripts/05f_he_et_al_scn8a_seurat.R`
- **Inputs:** Extracted slice files
- **Outputs:** `he_et_al_scn8a_seurat.rds`
- **R handover:** After `.rds` is created, ready for R analysis

### Samarasinghe 2021
- **Plain Seurat script:** `scripts/05_samarasinghe_2021_seurat.R`
- **Paper-style LIGER script:** `scripts/05b_samarasinghe_2021_liger.R`
- **Inputs:** `GSE165577_Filtered_counts_all_samples.csv.gz`; normalized GEO counts are downloaded but not used by these scripts
- **Current plain Seurat output:** `results/samarasinghe_2021/samarasinghe_2021_seurat.rds`
- **Pending LIGER outputs:** `results/samarasinghe_2021_liger/samarasinghe_2021_liger.rds` and UMAP plots
- **LIGER module note:** `05b` uses module/prebuilt-library packages only; it intentionally does not install packages in the job. The runtime must provide `SeuratWrappers` and `rliger` in addition to Seurat.
- **R handover:** Use the plain Seurat object for non-paper UMAPs; use the LIGER object, once generated, for paper-style all-sample integration and WT-control-only UMAP views on the integrated manifold.

### Liu 2025 hnbMO (GSE286235)
- **Scripts:** `scripts/01f_download_gse286235_geo.sh`, `scripts/05j_liu_2025_hnbmo_seurat.R`
- **Inputs:** Cell Ranger `raw_feature_bc_matrix.h5` files; healthy samples GSM8721440, GSM8721441, GSM8721442 are used by default
- **Outputs:** `liu_2025_hnbmo_healthy_seurat.rds`, UMAP/QC plots, sample and QC tables
- **R handover:** After `.rds` is created, ready for cross-study marker plotting or label transfer

---

## Comparative Plotting: On/Off Target Markers

- Comparative plots use the final `.rds` objects above as input.
- Panel B and cross-study marker scripts (e.g., `scripts/06_cross_study_panelB_markers.R`) expect these paths.
- All paths have been confirmed to exist (except Samarasinghe 2021, pending generation).

---
- 15 clusters (labels 0–14): `results/walsh_day75/kres_sweep_plots/umap_by_cluster_res0.8_k30.png`
- 14 clusters (labels 0–13): `results/walsh_day75/kres_sweep_plots/umap_by_cluster_res0.8_k50.png`
- Walsh-group annotated UMAP: `results/walsh_day75/annotation_plots/umap_by_walsh_group.png`

## Notes
- Great Lakes modules control R/Seurat (no in-job installs). Common patterns:
  - Walsh pipeline: `module load Bioinformatics` + `module load r-seurat/4.1.1-R-4.1.1-qyci4bo`
  - Plotting / Xiang / Bershteyn_2023: `module load Bioinformatics` + `module load r-seurat/4.1.1-R-4.2.0-5z5hgo7`
  - Samarasinghe LIGER: `slurm_templates/05b_samarasinghe_2021_liger.sbatch.template` loads `SEURAT_MODULE` plus optional `LIGER_EXTRA_MODULES`; the currently checked Seurat modules do not include `SeuratWrappers`/`rliger`, so this needs a site/user module stack that provides them before submission will pass preflight.
  - Panel B with Assay5 studies (e.g., Varela): submit with the Seurat v5 template and set `SEURAT5_MODULE` to an available v5 module.
  - Current Panel B Seurat v5 template uses:
    - `CONFIG_PATH=/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline/config/panel_b_cross_study_config.example.R`
    - `RUN_LABEL=panel_b_cross_study_v1_seurat5`
    - outputs under `/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/panel_b_cross_study_v1_seurat5/plots/`
  - Panel B marker list (fixed order; defined in `scripts/06_cross_study_panelB_markers.R`) is:
    - ON-target: `DCX,GAD2,DLX5,LHX6,MAF,SST,ERBB4,MEF2C,MAFB,LHX8,NKX2-1`
    - OFF-target: `SP8,PAX6,NEUROD2,ISL1,ACHE,NKX6-2,MKI67` (`ACHE`, not `ACHF`)
  - Xiang Panel B note: Xiang features are Ensembl-like IDs; symbol markers require symbol-to-Ensembl mapping for direct lookup.
  - Panel B symbol remap is dynamic: set optional `feature_map_path` per study in `config/panel_b_cross_study_config.example.R` and the script remaps current `GENE_ORDER` markers automatically.
  - Panel B combined figure is emitted as PNG/PDF/SVG with fixed layout logic: study rows (top->bottom), gene columns (left->right), ON-target block on top, OFF-target block on bottom.
  - Panel B now requires both Varela studies and fixed first-two row order:
    - `varela_this_paper` (DIV30), then `varela_div90` (DIV90).
  - Panel B also emits separate ON-only and OFF-only PNG/PDF/SVG files; PNG exports are high print quality (`dpi=600`).
  - Panel B figure labels now include plotted cell counts per study (`n=<cells>`), and the config example includes both `Bershteyn 2025` and `Bershteyn 2023`.
  - Panel B now also writes `results/<run_label>/plots/panel_b_prepared_inputs.rds` for downstream scripts (matched coords + marker matrix per study).
  - Panel B publishes validated per-study Seurat objects to a stable canonical path for downstream steps:
    `results/panel_b_prepared_objects/studies/<study_id>_panelb_ready_seurat.rds`
    with manifest `results/panel_b_prepared_objects/panel_b_prepared_object_paths.tsv`.
  - Stage 07 placeholder consumes that manifest and writes run-scoped planning artifacts under:
    `results/<run_label>/human_ge_comparison_tbd/`.
  - To download Panel B outputs to local Downloads (run from local terminal, not Great Lakes login node), for example:
    `rsync -avh --progress --prune-empty-dirs --include='*/' --include='*.png' --include='*.pdf' --include='*.svg' --exclude='*' "elcrespo@gl-login1.arc-ts.umich.edu:/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/panel_b_cross_study_v1_seurat5/plots/" "/Users/ecrespo/Downloads/panel_b_cross_study_v1_seurat5_plots/"`
- Slurm does **not** expand shell variables in `#SBATCH --output/--error`; templates use absolute log paths to keep runtime artifacts out of the repo.
- Methods-locked parameters retained: QC 1000–5000 genes, <15% MT (^MT-); LogNormalize (scale.factor=1e4); HVG=5000 (vst); regress S.Score & G2M.Score; dims=1:20 for neighbors/UMAP; clustering resolution default 2.0 (sweeps explore alternatives); stress removal clusters mean score >0.5.
- Checkpoints allow reruns without redoing QC/normalization/scaling; sweeps start from post-stress PCA checkpoint.

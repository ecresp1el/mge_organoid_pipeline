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
- NeMO landing-page downloads (Siebert 2026, `nemo:dat-htzat9t`):
  - `PROJECT_ROOT/data/raw/siebert_2026_nemo/` (metadata, BDBags, subset fetch + downloads)
  - Status: **TBD** until NeMO publishes manifests/files.
- Processed / results:
  - Walsh: `PROJECT_ROOT/results/walsh_day75/`
  - Bershteyn_2025: `PROJECT_ROOT/results/bershteyn_2025/` (canonical Seurat copy + UMAP plots)
  - Xiang_2018: `PROJECT_ROOT/results/xiang_2018/` (Seurat + UMAP generated from GSE98201 10x)
  - Bershteyn_2023: `PROJECT_ROOT/results/bershteyn_2023/` (provided Seurat object copied + UMAP plotted)
  - Samarasinghe_2021: **not generated yet** (expected: `PROJECT_ROOT/results/samarasinghe_2021/`)
  - Samarasinghe_2021 LIGER (paper-matched): **not generated yet** (expected: `PROJECT_ROOT/results/samarasinghe_2021_liger/`)
  - Varela (this paper): `PROJECT_ROOT/results/varela_this_paper/varela_this_paper_seurat.rds` (true copy of `/nfs/turbo/umms-parent/mgeo_scRNAseq/day30_old/Day30.rds`)
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
| Varela (this paper, Day30) | Yes (`/nfs/turbo/umms-parent/mgeo_scRNAseq/day30_old/Day30.rds`) | Yes (`results/varela_this_paper/varela_this_paper_seurat.rds`) | Included in Panel B assembly | Canonical Panel B path is under `PROJECT_ROOT/results/...` and now uses a true copied file. |
| Siebert 2026 (NeMO `nemo:dat-htzat9t`) | Not published (metadata only) | No | No | Waiting on NeMO manifests/files; placeholder only. |

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
- Samarasinghe 2021 Seurat/UMAP: `scripts/05_samarasinghe_2021_seurat.R`
- Samarasinghe 2021 LIGER (paper settings): `scripts/05b_samarasinghe_2021_liger.R`
- Xiang 2018 scRNA 10x (GSE98201) Seurat/UMAP: `scripts/05c_xiang_2018_seurat.R`
- Bershteyn 2023 Seurat plotting: `scripts/05d_bershteyn_2023_seurat_plot.R`
- Cross-study marker Panel B figure assembly (no analysis/recompute): `scripts/06_cross_study_panelB_markers.R`
- Slurm template for cross-study Panel B (Seurat v4 runtime): `slurm_templates/06_cross_study_panelB_markers.sbatch.template`
- Slurm template for cross-study Panel B (Seurat v5 runtime / Assay5-aware): `slurm_templates/06_cross_study_panelB_markers_seurat5.sbatch.template`
- Main pipeline: `scripts/02_walsh_day75_seurat.R` (methods-locked QC/normalize/HVG/CC/Scale/PCA/stress removal/UMAP/clustering; checkpoints; Elbow fallback)
- Sweeps: `scripts/02b_walsh_k_sweep.R`, `02c_walsh_dims_sweep.R`, `02d_walsh_resolution_sweep.R`, `02e_walsh_kres_sweep.R`
- Annotation: `scripts/02f_walsh_annotation.R` (adds walsh_group labels, tables, annotated UMAP)

## Canonical plot paths to share
- 15 clusters (labels 0–14): `results/walsh_day75/kres_sweep_plots/umap_by_cluster_res0.8_k30.png`
- 14 clusters (labels 0–13): `results/walsh_day75/kres_sweep_plots/umap_by_cluster_res0.8_k50.png`
- Walsh-group annotated UMAP: `results/walsh_day75/annotation_plots/umap_by_walsh_group.png`

## Notes
- Great Lakes modules control R/Seurat (no in-job installs). Common patterns:
  - Walsh pipeline: `module load Bioinformatics` + `module load r-seurat/4.1.1-R-4.1.1-qyci4bo`
  - Plotting / Xiang / Bershteyn_2023: `module load Bioinformatics` + `module load r-seurat/4.1.1-R-4.2.0-5z5hgo7`
  - Panel B with Assay5 studies (e.g., Varela): submit with the Seurat v5 template and set `SEURAT5_MODULE` to an available v5 module.
  - Panel B marker list includes `ACHE` (not `ACHF`).
  - Xiang Panel B note: Xiang features are Ensembl-like IDs; symbol markers require symbol-to-Ensembl mapping for direct lookup.
  - Panel B symbol remap is dynamic: set optional `feature_map_path` per study in `config/panel_b_cross_study_config.example.R` and the script remaps current `GENE_ORDER` markers automatically.
  - Panel B final figure is emitted as PNG/PDF/SVG with fixed layout logic: Varela left-most column when present, ON-target block on top, OFF-target block on bottom.
  - Panel B figure labels now include plotted cell counts per study (`n=<cells>`), and the config example includes both `Bershteyn 2025` and `Bershteyn 2023`.
  - Panel B now also writes `results/<run_label>/plots/panel_b_prepared_inputs.rds` for downstream scripts (matched coords + marker matrix per study).
- Slurm does **not** expand shell variables in `#SBATCH --output/--error`; templates use absolute log paths to keep runtime artifacts out of the repo.
- Methods-locked parameters retained: QC 1000–5000 genes, <15% MT (^MT-); LogNormalize (scale.factor=1e4); HVG=5000 (vst); regress S.Score & G2M.Score; dims=1:20 for neighbors/UMAP; clustering resolution default 2.0 (sweeps explore alternatives); stress removal clusters mean score >0.5.
- Checkpoints allow reruns without redoing QC/normalization/scaling; sweeps start from post-stress PCA checkpoint.

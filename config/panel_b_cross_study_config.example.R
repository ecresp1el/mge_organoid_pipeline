list(
  # Optional if PROJECT_ROOT env or --project-root is provided
  project_root = "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder",
  run_label = "panel_b_cross_study_v1",
  studies = data.frame(
    study_id = c(
      "walsh",
      "bershteyn_2025",
      "xiang_2018",
      "samarasinghe_2021",
      "varela_this_paper"
    ),
    study_label = c(
      "Walsh",
      "Bershteyn",
      "Xiang",
      "Samarasinghe",
      "Varela (This paper)"
    ),
    # Paths are explicit and resolved relative to project_root if not absolute.
    # No path guessing and no globbing are used.
    object_path = c(
      "results/walsh_day75/walsh_day75_final_annotated.rds",
      "data/raw/bershteyn_2025_geo_files/suppl/GSE283775_Seurat_scRNA_seq.rds.gz",
      "results/xiang_2018/xiang_2018_seurat.rds",
      "results/samarasinghe_2021/samarasinghe_2021_seurat.rds",
      "results/varela_this_paper/varela_this_paper_seurat.rds"
    ),
    # Note: Walsh final object uses "umap_sel"/"umap20" rather than "umap".
    reduction = c("umap_sel", "umap", "umap", "umap", "umap"),
    assay = c("RNA", "RNA", "RNA", "RNA", "RNA"),
    # Optional. Script defaults to "data" then falls back to "counts" if needed.
    expression_slot = c("data", "data", "data", "data", "data"),
    stringsAsFactors = FALSE
  )
)

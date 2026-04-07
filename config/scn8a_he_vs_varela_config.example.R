list(
  # Optional if PROJECT_ROOT env or --project-root is provided.
  project_root = "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder",
  run_label = "panel_b_scn8a_he_vs_varela_v1",
  studies = data.frame(
    study_id = c(
      "he_et_al",
      "varela_this_paper"
    ),
    study_label = c(
      "He et al (HNOCA Full V2)",
      "Varela (This paper)"
    ),
    object_path = c(
      "results/he_et_al/he_et_al_scn8a_seurat.rds",
      "results/varela_this_paper/varela_this_paper_seurat.rds"
    ),
    reduction = c("umap", "umap"),
    assay = c("RNA", "RNA"),
    # He object is generated from SCN8A-only slice; data slot is populated via NormalizeData.
    expression_slot = c("data", "data"),
    stringsAsFactors = FALSE
  )
)

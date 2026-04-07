list(
  project_root = "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder",
  run_label = "gene_panel_he_vs_varela_v1",
  genes = c("LHX6", "NKX2.1"),
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
    expression_slot = c("data", "data"),
    stringsAsFactors = FALSE
  )
)

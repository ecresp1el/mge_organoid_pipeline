# Interactive driver for Panel B plotting in RStudio.
# Open this file in RStudio and run line-by-line.
# For Assay5 objects (for example Varela), run in a Seurat v5-capable R session.

# 1) Point RStudio at the repo.
setwd("/home/elcrespo/Desktop/githubprojects/mge_organoid_pipeline")

# 2) Load functions (does not auto-run the pipeline).
source("scripts/06_cross_study_panelB_markers.R")

# 3) Edit these values if needed.
project_root <- "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"
config_path <- "config/panel_b_cross_study_config.example.R"
run_label <- paste0("panel_b_rstudio_", format(Sys.time(), "%Y%m%d_%H%M%S"))

# 4) Run locally and preview each row plot in RStudio Plots pane.
res <- run_panel_b_local(
  config_path = config_path,
  project_root = project_root,
  run_label = run_label,
  retain_seurat = TRUE,
  export_global = TRUE,
  show_progress_plots = FALSE
)

# 5) Inspect outputs in RStudio.
View(panel_b_studies)
View(panel_b_metadata_summary)
View(panel_b_metadata_columns)
View(panel_b_assay_slot_summary)
View(panel_b_reduction_summary)
View(panel_b_ident_counts)
View(panel_b_feature_space)
View(panel_b_rows)
View(panel_b_issues)

panel_b_row_plots[["DCX"]]
panel_b_final_plot

res$output_paths$pdf
res$output_paths$svg
file.exists(res$output_paths$pdf)
file.exists(res$output_paths$svg)




####### trouble shooting ######
# Grab the object 
cfg <-dget("config/panel_b_cross_study_config.example.R")
cfg$studies

# pick the row for walsh
p <- file.path(cfg$project_root, cfg$studies$object_path[cfg$studies$study_id=="walsh"])
obj <- read_rds_any(p)

names(obj@assays)
DefaultAssay(obj)
names(obj@reductions)

#

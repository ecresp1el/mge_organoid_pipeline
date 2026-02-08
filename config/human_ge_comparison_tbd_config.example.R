list(
  # Optional if PROJECT_ROOT env or --project-root is provided.
  project_root = "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder",
  run_label = "human_ge_comparison_tbd_v1",

  # Canonical output from scripts/06_cross_study_panelB_markers.R
  prepared_object_manifest_path = "results/panel_b_prepared_objects/panel_b_prepared_object_paths.tsv",
  prepared_object_studies_dir = "results/panel_b_prepared_objects/studies",

  # Studies expected for the GE-comparison stage.
  required_studies = c(
    "walsh",
    "bershteyn_2025",
    "bershteyn_2023",
    "xiang_2018",
    "varela_this_paper"
  ),

  # Placeholder path for forthcoming human GE reference data integration.
  human_ge_reference_path = "metadata/human_developing_ge/TBD_reference_path.tsv"
)

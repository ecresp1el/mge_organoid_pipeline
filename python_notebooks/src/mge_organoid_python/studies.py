"""Study definitions used by the Python notebook workflow."""

from dataclasses import dataclass
from pathlib import Path

from .paths import default_anndata_dir


@dataclass(frozen=True)
class StudySpec:
    """Description of one Seurat object to convert to AnnData."""

    study_id: str
    label: str
    seurat_path: str
    assay: str = "RNA"
    reduction: str = "umap"
    expression_layer: str = "data"
    output_filename: str = ""

    @property
    def source_path(self):
        return Path(self.seurat_path).expanduser()

    def h5ad_path(self, project_root=None, output_dir=None):
        filename = self.output_filename or "{}.h5ad".format(self.study_id)
        base = Path(output_dir).expanduser() if output_dir else default_anndata_dir(project_root)
        return base / filename


def default_studies():
    """Return the canonical studies for the notebook entry point."""
    return [
        StudySpec(
            study_id="shi_2019_paper_qc",
            label="Shi 2019 paper QC",
            seurat_path=(
                "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/"
                "results/shi_2019_paper_qc/shi_2019_seurat.rds"
            ),
            assay="RNA",
            reduction="umap",
            expression_layer="data",
            output_filename="shi_2019_paper_qc.h5ad",
        ),
        StudySpec(
            study_id="varela_div30",
            label="Varela DIV30",
            seurat_path=(
                "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/"
                "results/varela_this_paper/varela_this_paper_seurat.rds"
            ),
            assay="RNA",
            reduction="umap",
            expression_layer="data",
            output_filename="varela_div30.h5ad",
        ),
        StudySpec(
            study_id="varela_div90",
            label="Varela DIV90",
            seurat_path=(
                "/nfs/turbo/umms-parent/Manny_test/ventral_sosrs_output/"
                "umap_props_output/clustered_day90_with_cluster_names_2.rds"
            ),
            assay="RNA",
            reduction="umap",
            expression_layer="data",
            output_filename="varela_div90.h5ad",
        ),
    ]

"""Validation and reporting helpers for AnnData conversion."""

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class AnnDataReport:
    study_id: str
    label: str
    h5ad_path: str
    n_obs: int
    n_vars: int
    has_umap: bool
    n_obs_columns: int
    n_var_columns: int

    def as_dict(self):
        return {
            "study_id": self.study_id,
            "label": self.label,
            "h5ad_path": self.h5ad_path,
            "n_obs": self.n_obs,
            "n_vars": self.n_vars,
            "has_umap": self.has_umap,
            "n_obs_columns": self.n_obs_columns,
            "n_var_columns": self.n_var_columns,
        }


def validate_source_paths(studies):
    """Return a list of missing Seurat source paths."""
    missing = []
    for study in studies:
        path = Path(study.seurat_path).expanduser()
        if not path.exists():
            missing.append((study.study_id, str(path)))
    return missing


def validate_anndata(study, adata, h5ad_path):
    """Validate minimal AnnData structure needed for downstream notebooks."""
    errors = []
    if adata.n_obs <= 0:
        errors.append("AnnData has zero cells")
    if adata.n_vars <= 0:
        errors.append("AnnData has zero genes/features")
    if "X_umap" not in adata.obsm:
        errors.append('AnnData missing obsm["X_umap"]')
    else:
        umap = adata.obsm["X_umap"]
        if getattr(umap, "shape", (0, 0))[0] != adata.n_obs:
            errors.append('obsm["X_umap"] row count does not match n_obs')
        if getattr(umap, "shape", (0, 0))[1] < 2:
            errors.append('obsm["X_umap"] has fewer than two columns')
    if adata.obs is None or adata.obs.shape[0] != adata.n_obs:
        errors.append("obs metadata is missing or has the wrong row count")
    if adata.var is None or adata.var.shape[0] != adata.n_vars:
        errors.append("var metadata is missing or has the wrong row count")
    if errors:
        raise ValueError("{} validation failed: {}".format(study.study_id, "; ".join(errors)))

    return AnnDataReport(
        study_id=study.study_id,
        label=study.label,
        h5ad_path=str(Path(h5ad_path)),
        n_obs=int(adata.n_obs),
        n_vars=int(adata.n_vars),
        has_umap=True,
        n_obs_columns=int(adata.obs.shape[1]),
        n_var_columns=int(adata.var.shape[1]),
    )

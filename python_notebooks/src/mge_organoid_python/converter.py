"""Seurat RDS to AnnData conversion through a Python-driven R bridge."""

from pathlib import Path
from datetime import datetime
import os
import platform
import subprocess
import sys
import tempfile

from .paths import default_anndata_dir, ensure_under_path, resolve_project_root
from .validation import validate_anndata


def _timestamp():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


class SeuratToAnnDataConverter:
    """Convert Seurat `.rds` objects to cached AnnData `.h5ad` files."""

    def __init__(self, project_root=None, output_dir=None, overwrite=False, verbose=True, repo_root=None):
        self.project_root = resolve_project_root(project_root)
        self.output_dir = Path(output_dir).expanduser() if output_dir else default_anndata_dir(self.project_root)
        self.output_dir = ensure_under_path(self.output_dir, self.project_root)
        self.overwrite = bool(overwrite)
        self.verbose = bool(verbose)
        self.repo_root = Path(repo_root).expanduser().resolve() if repo_root else self._infer_repo_root()
        self._r_convert = None

    def _infer_repo_root(self):
        here = Path(__file__).resolve()
        # converter.py -> mge_organoid_python -> src -> python_notebooks -> repo
        return here.parents[3]

    def log(self, message):
        if self.verbose:
            print("[{}] {}".format(_timestamp(), message), flush=True)

    def output_path(self, study):
        path = study.h5ad_path(project_root=self.project_root, output_dir=self.output_dir)
        return ensure_under_path(path, self.project_root)

    def needs_conversion(self, study):
        source = Path(study.seurat_path).expanduser()
        target = self.output_path(study)
        if self.overwrite or not target.exists():
            return True
        if source.exists() and source.stat().st_mtime > target.stat().st_mtime:
            return True
        return False

    def convert(self, study, overwrite=None):
        """Convert one study if needed and return `(adata, report)`."""
        target = self.convert_file(study, overwrite=overwrite)

        self.log("Study {}: loading H5AD into AnnData".format(study.study_id))
        adata = self._read_h5ad(target)
        self.log("Study {}: validating AnnData".format(study.study_id))
        report = validate_anndata(study, adata, target)
        self.log(
            "Study {study_id}: ready n_obs={n_obs:,} n_vars={n_vars:,} has_umap={has_umap}".format(
                study_id=study.study_id,
                n_obs=report.n_obs,
                n_vars=report.n_vars,
                has_umap=report.has_umap,
            )
        )
        return adata, report

    def convert_file(self, study, overwrite=None):
        """Convert one study if needed and return the cached `.h5ad` path.

        This does not load the H5AD into memory. Use this for large studies.
        """
        self._refuse_login_node_conversion()
        source = Path(study.seurat_path).expanduser()
        if not source.exists():
            raise FileNotFoundError("Missing Seurat source for {}: {}".format(study.study_id, source))

        target = self.output_path(study)
        target.parent.mkdir(parents=True, exist_ok=True)

        do_overwrite = self.overwrite if overwrite is None else bool(overwrite)
        needs_conversion = do_overwrite or self.needs_conversion(study)
        self.log(
            "Study {study_id}: source={source} target={target} needs_conversion={needs}".format(
                study_id=study.study_id,
                source=source,
                target=target,
                needs=needs_conversion,
            )
        )
        if needs_conversion:
            self.log("Study {}: starting RDS -> H5AD conversion".format(study.study_id))
            self._run_rscript_converter(source, target, study, do_overwrite)
            self.log("Study {}: finished conversion write".format(study.study_id))
        else:
            self.log("Study {}: using existing cached H5AD".format(study.study_id))

        return target

    def _refuse_login_node_conversion(self):
        host = platform.node()
        allow = os.environ.get("MGE_ALLOW_LOGIN_NODE_CONVERSION", "").lower() in {"1", "true", "yes"}
        if host.startswith("gl-login") and not allow:
            raise RuntimeError(
                "Refusing Seurat -> AnnData conversion on login node '{}'. "
                "Connect VS Code/Jupyter to an allocated compute node first. "
                "This guard can be bypassed only by setting MGE_ALLOW_LOGIN_NODE_CONVERSION=true, "
                "which is not recommended.".format(host)
            )

    def convert_many(self, studies, overwrite=None):
        """Convert studies and return `(adatas, reports)` dictionaries/lists."""
        studies = list(studies)
        self.log("Starting conversion/load for {} studies".format(len(studies)))
        adatas = {}
        reports = []
        for idx, study in enumerate(studies, start=1):
            self.log("({}/{}) {}".format(idx, len(studies), study.study_id))
            adata, report = self.convert(study, overwrite=overwrite)
            adatas[study.study_id] = adata
            reports.append(report)
        self.log("All studies complete")
        return adatas, reports

    def convert_many_files(self, studies, overwrite=None):
        """Convert studies and return a dict of study id to cached `.h5ad` path."""
        studies = list(studies)
        self.log("Starting file conversion for {} studies".format(len(studies)))
        paths = {}
        for idx, study in enumerate(studies, start=1):
            self.log("({}/{}) {}".format(idx, len(studies), study.study_id))
            paths[study.study_id] = self.convert_file(study, overwrite=overwrite)
        self.log("All file conversions complete")
        return paths

    def _run_rscript_converter(self, source, target, study, overwrite):
        script = self.repo_root / "python_notebooks" / "scripts" / "seurat_export_for_anndata.R"
        if not script.exists():
            raise FileNotFoundError("Missing R conversion helper: {}".format(script))

        export_dir = Path(
            tempfile.mkdtemp(
                prefix="{}_".format(study.study_id),
                dir=str(target.parent),
            )
        )
        cmd = [
            "Rscript",
            str(script),
            "--seurat",
            str(source),
            "--outdir",
            str(export_dir),
            "--assay",
            study.assay,
            "--reduction",
            study.reduction,
            "--expression_layer",
            study.expression_layer,
        ]
        self.log("Running Rscript subprocess: {}".format(" ".join(cmd)))

        env = os.environ.copy()
        env.setdefault("PROJECT_ROOT", str(self.project_root))
        proc = subprocess.Popen(
            cmd,
            cwd=str(self.repo_root),
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            print(line, end="", flush=True)
        returncode = proc.wait()
        if returncode != 0:
            raise RuntimeError(
                "Rscript Seurat -> H5AD conversion failed for {} with exit code {}".format(
                    study.study_id,
                    returncode,
                )
            )
        self._write_h5ad_from_export(export_dir, target, study)

    def _write_h5ad_from_export(self, export_dir, target, study):
        self.log("Study {}: reading R export into Python".format(study.study_id))
        try:
            import anndata as ad
            import pandas as pd
            from scipy import io as scipy_io
        except ImportError as exc:
            raise RuntimeError(
                "Missing Python packages for H5AD writing. Need anndata, pandas, and scipy."
            ) from exc

        matrix_path = export_dir / "matrix.mtx"
        features_path = export_dir / "features.tsv"
        barcodes_path = export_dir / "barcodes.tsv"
        obs_path = export_dir / "obs.tsv"
        umap_path = export_dir / "umap.tsv"
        manifest_path = export_dir / "manifest.tsv"

        x = scipy_io.mmread(matrix_path).tocsr().transpose().tocsr()
        features = pd.read_csv(features_path, sep="\t")
        barcodes = pd.read_csv(barcodes_path, sep="\t")
        obs = pd.read_csv(obs_path, sep="\t", dtype=str)
        umap = pd.read_csv(umap_path, sep="\t")
        manifest = pd.read_csv(manifest_path, sep="\t")

        obs = obs.set_index("cell_id", drop=False)
        var = features.set_index("feature_id", drop=False)
        if list(obs.index) != list(barcodes["cell_id"]):
            obs = obs.loc[list(barcodes["cell_id"])]
        umap = umap.set_index("cell_id").loc[obs.index]

        adata = ad.AnnData(X=x, obs=obs, var=var)
        adata.obsm["X_umap"] = umap[["UMAP_1", "UMAP_2"]].to_numpy()
        adata.uns["source_seurat_path"] = str(study.seurat_path)
        adata.uns["seurat_assay"] = study.assay
        adata.uns["seurat_reduction"] = study.reduction
        adata.uns["conversion_manifest"] = dict(zip(manifest["key"], manifest["value"]))

        self.log(
            "Study {}: writing H5AD with n_obs={:,} n_vars={:,}".format(
                study.study_id,
                adata.n_obs,
                adata.n_vars,
            )
        )
        adata.write_h5ad(target)
        self.log("Study {}: finished Python H5AD write: {}".format(study.study_id, target))

    def _read_h5ad(self, path):
        try:
            import anndata as ad
        except ImportError as exc:
            raise RuntimeError(
                "Missing Python package 'anndata'. Activate the mge-organoid-python conda env."
            ) from exc
        return ad.read_h5ad(path)

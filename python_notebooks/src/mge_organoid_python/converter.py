"""Seurat RDS to AnnData conversion through a Python-driven R bridge."""

from pathlib import Path

from .paths import default_anndata_dir, ensure_under_path, resolve_project_root
from .validation import validate_anndata


class SeuratToAnnDataConverter:
    """Convert Seurat `.rds` objects to cached AnnData `.h5ad` files."""

    def __init__(self, project_root=None, output_dir=None, overwrite=False):
        self.project_root = resolve_project_root(project_root)
        self.output_dir = Path(output_dir).expanduser() if output_dir else default_anndata_dir(self.project_root)
        self.output_dir = ensure_under_path(self.output_dir, self.project_root)
        self.overwrite = bool(overwrite)
        self._r_convert = None

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
        source = Path(study.seurat_path).expanduser()
        if not source.exists():
            raise FileNotFoundError("Missing Seurat source for {}: {}".format(study.study_id, source))

        target = self.output_path(study)
        target.parent.mkdir(parents=True, exist_ok=True)

        do_overwrite = self.overwrite if overwrite is None else bool(overwrite)
        if do_overwrite or self.needs_conversion(study):
            r_convert = self._get_r_converter()
            r_convert(
                str(source),
                str(target),
                study.assay,
                study.reduction,
                study.expression_layer,
                bool(do_overwrite),
            )

        adata = self._read_h5ad(target)
        report = validate_anndata(study, adata, target)
        return adata, report

    def convert_many(self, studies, overwrite=None):
        """Convert studies and return `(adatas, reports)` dictionaries/lists."""
        adatas = {}
        reports = []
        for study in studies:
            adata, report = self.convert(study, overwrite=overwrite)
            adatas[study.study_id] = adata
            reports.append(report)
        return adatas, reports

    def _read_h5ad(self, path):
        try:
            import anndata as ad
        except ImportError as exc:
            raise RuntimeError(
                "Missing Python package 'anndata'. Activate the mge-organoid-python conda env."
            ) from exc
        return ad.read_h5ad(path)

    def _get_r_converter(self):
        if self._r_convert is not None:
            return self._r_convert

        try:
            from rpy2 import robjects
        except ImportError as exc:
            raise RuntimeError(
                "Missing Python package 'rpy2'. Activate the mge-organoid-python conda env."
            ) from exc

        robjects.r(
            r'''
            .mge_convert_seurat_to_h5ad <- function(
                seurat_path,
                h5ad_path,
                assay,
                reduction,
                expression_layer,
                overwrite
            ) {
                suppressPackageStartupMessages({
                    library(Seurat)
                    library(SeuratObject)
                    library(SingleCellExperiment)
                    library(SummarizedExperiment)
                    library(zellkonverter)
                })

                if (!file.exists(seurat_path)) {
                    stop("Seurat source does not exist: ", seurat_path, call. = FALSE)
                }
                if (file.exists(h5ad_path) && !isTRUE(overwrite)) {
                    return(normalizePath(h5ad_path, mustWork = TRUE))
                }

                obj <- readRDS(seurat_path)
                if (!inherits(obj, "Seurat")) {
                    stop("Object is not a Seurat object: ", seurat_path, call. = FALSE)
                }
                if (!(assay %in% Seurat::Assays(obj))) {
                    stop(
                        "Missing assay '", assay, "'. Available assays: ",
                        paste(Seurat::Assays(obj), collapse = ", "),
                        call. = FALSE
                    )
                }
                if (!(reduction %in% Seurat::Reductions(obj))) {
                    stop(
                        "Missing reduction '", reduction, "'. Available reductions: ",
                        paste(Seurat::Reductions(obj), collapse = ", "),
                        call. = FALSE
                    )
                }

                Seurat::DefaultAssay(obj) <- assay
                if (
                    inherits(obj[[assay]], "Assay5") &&
                    exists("JoinLayers", where = asNamespace("SeuratObject"), inherits = FALSE)
                ) {
                    obj <- SeuratObject::JoinLayers(obj, assay = assay)
                }
                sce <- Seurat::as.SingleCellExperiment(obj, assay = assay)

                emb <- Seurat::Embeddings(obj, reduction = reduction)
                if (is.null(dim(emb)) || nrow(emb) == 0 || ncol(emb) < 2) {
                    stop("Reduction '", reduction, "' does not contain a two-dimensional embedding.", call. = FALSE)
                }
                missing_emb <- setdiff(colnames(sce), rownames(emb))
                if (length(missing_emb) > 0) {
                    stop(
                        "Reduction '", reduction, "' is missing embeddings for ",
                        length(missing_emb), " cells.",
                        call. = FALSE
                    )
                }
                emb <- emb[colnames(sce), , drop = FALSE]
                emb <- as.matrix(emb[, seq_len(2), drop = FALSE])
                colnames(emb) <- c("UMAP_1", "UMAP_2")
                SingleCellExperiment::reducedDim(sce, "X_umap") <- emb

                assay_names <- names(SummarizedExperiment::assays(sce))
                x_name <- NULL
                if (identical(expression_layer, "data") && "logcounts" %in% assay_names) {
                    x_name <- "logcounts"
                } else if (expression_layer %in% assay_names) {
                    x_name <- expression_layer
                } else if ("counts" %in% assay_names) {
                    x_name <- "counts"
                } else if (length(assay_names) > 0) {
                    x_name <- assay_names[[1]]
                }

                dir.create(dirname(h5ad_path), showWarnings = FALSE, recursive = TRUE)
                zellkonverter::writeH5AD(sce, file = h5ad_path, X_name = x_name)
                normalizePath(h5ad_path, mustWork = TRUE)
            }
            '''
        )
        self._r_convert = robjects.globalenv[".mge_convert_seurat_to_h5ad"]
        return self._r_convert

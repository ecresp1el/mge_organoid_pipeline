# Safely inspect a serialized RDS object without triggering computation
install.packages("Seurat")
library(Seurat)
library(ggplot2)

#version of seurat loaded
packageVersion("Seurat")


## input path to the SeuratObject ##
input_rds_path <- "/nfs/turbo/umms-parent/mgeo_scRNAseq/day30_old/Day30.rds"

# load the Seurat object 
mge_organoid_div30 <- readRDS(input_rds_path)

# version of the object 
mge_organoid_div30@version

# basic identity and structure 
mge_organoid_div30

# extract meta data of the Seurat object, and print class to console
md = mge_organoid_div30@meta.data
head(md)
ncol(md) # should be 11 for this case

# make the umap if stored in the Seurat Object 
DimPlot(
  mge_organoid_div30, 
  reduction= "umap",
  label = TRUE, 
  repel = TRUE, 
  group.by = "seurat_clusters"
  )

#### function to make #### 
plot_gene_with_umaps_and_clusters <- function(
  obj,
  gene,
  slot = "data",
  pt.size = 0.4
) {
  stopifnot(inherits(obj, "Seurat"))

  # Safety checks
  if (!(gene %in% rownames(obj))) {
    stop(paste("Gene", gene, "not found in DefaultAssay"))
  }

  if (!("seurat_clusters" %in% colnames(obj@meta.data))) {
    stop("seurat_clusters not found in metadata")
  }

  if (!("umap" %in% Reductions(obj))) {
    stop("UMAP reduction not found in object")
  }

  ## 1️⃣ UMAP with clusters
  p_clusters <- DimPlot(
    obj,
    reduction = "umap",
    group.by = "seurat_clusters",
    label = TRUE,
    repel = TRUE,
    pt.size = pt.size
  ) +
    ggplot2::ggtitle("UMAP: Seurat clusters")

  ## 2️⃣ UMAP with gene expression
  p_gene <- FeaturePlot(
    obj,
    features = gene,
    reduction = "umap",
    slot = slot,
    pt.size = pt.size
  ) +
    ggplot2::ggtitle(paste(gene, "expression (UMAP)"))

  ## 3️⃣ Gene expression per cluster
  p_cluster_expr <- VlnPlot(
    obj,
    features = gene,
    group.by = "seurat_clusters",
    slot = slot,
    pt.size = 0
  ) +
    ggplot2::ggtitle(paste(gene, "expression per cluster")) +
    ggplot2::theme_classic()

  # Combine all three
  p_clusters + p_gene + p_cluster_expr
}

plot_gene_with_umaps_and_clusters <- function(
    obj,
    gene,
    slot = "data",
    pt.size = 0.4
) {
  stopifnot(inherits(obj, "Seurat"))
  
  # Safety checks
  if (!(gene %in% rownames(obj))) {
    stop(paste("Gene", gene, "not found in DefaultAssay"))
  }
  
  if (!("seurat_clusters" %in% colnames(obj@meta.data))) {
    stop("seurat_clusters not found in metadata")
  }
  
  if (!("umap" %in% Reductions(obj))) {
    stop("UMAP reduction not found in object")
  }
  
  ## 1️⃣ UMAP with clusters
  p_clusters <- DimPlot(
    obj,
    reduction = "umap",
    group.by = "seurat_clusters",
    label = TRUE,
    repel = TRUE,
    pt.size = pt.size
  ) +
    ggplot2::ggtitle("UMAP: Seurat clusters")
  
  ## 2️⃣ UMAP with gene expression
  p_gene <- FeaturePlot(
    obj,
    features = gene,
    reduction = "umap",
    slot = slot,
    pt.size = pt.size
  ) +
    ggplot2::ggtitle(paste(gene, "expression (UMAP)"))
  
  ## 3️⃣ Gene expression per cluster
  p_cluster_expr <- VlnPlot(
    obj,
    features = gene,
    group.by = "seurat_clusters",
    slot = slot,
    pt.size = 0
  ) +
    ggplot2::ggtitle(paste(gene, "expression per cluster")) +
    ggplot2::theme_classic()
  
  # Combine all three
  p_clusters + p_gene + p_cluster_expr
}




### sanity check 
stopifnot(
  inherits(mge_organoid_div30, "Seurat"),
  "umap" %in% Reductions(mge_organoid_div30),
  "seurat_clusters" %in% colnames(mge_organoid_div30@meta.data)
)

### 
plot_gene_umap_and_clusters(
  obj  = mge_organoid_div30,
  gene = "LHX6"
)

### 
plot_gene_with_umaps_and_clusters(
  mge_organoid_div30,
  gene = "LHX6"
)




############### TEST 
############################################################
## R EXPLORATION HELPERS FOR SEURAT v5
##
## PURPOSE:
##   - Visualize UMAPs + gene expression reproducibly
##   - Save high-quality figures automatically
##   - Avoid Seurat v5 API pitfalls
##   - Be explicit, readable, and debuggable
##
## ASSUMPTIONS:
##   - Seurat >= 5.0
##   - UMAP already computed and stored in object
##   - seurat_clusters exists in meta.data
##   - DefaultAssay() is already set correctly by user
##
## AUTHORSHIP MODEL:
##   - Treat Seurat as a data container
##   - Avoid mutating object state
##   - Fail loudly if expectations are not met
############################################################


##############################
## REQUIRED LIBRARIES
##############################

# Seurat for scRNA-seq object + plotting
library(Seurat)

# ggplot2 for plot composition and saving
library(ggplot2)

# patchwork is automatically used by Seurat for "+"
# (loaded implicitly, no need to library() it)


############################################################
## HELPER 1: OUTPUT DIRECTORY CREATION
############################################################
## PURPOSE:
##   - Ensure a single, consistent output directory exists
##   - Avoid scattered plot outputs
##
## RETURNS:
##   - Absolute path to output directory
############################################################

get_output_dir <- function() {
  
  # Hard-coded, explicit path (per your request)
  out_dir <- "/home/elcrespo/Desktop/r_exploration_output"
  
  # Create directory if it does not exist
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  return(out_dir)
}


############################################################
## HELPER 2: SAVE PLOTS AT HIGH QUALITY
############################################################
## PURPOSE:
##   - Save plots in BOTH PNG and SVG formats
##   - PNG: high-resolution raster for quick viewing
##   - SVG: vector graphics for Illustrator / Inkscape
##
## ARGUMENTS:
##   p              : ggplot / patchwork object
##   filename_base  : base filename WITHOUT extension
##   width, height  : inches (controls aspect + readability)
##   dpi            : raster resolution (PNG only)
############################################################

save_plot <- function(p,
                      filename_base,
                      width = 16,
                      height = 6,
                      dpi = 300) {
  
  # Resolve output directory
  out_dir <- get_output_dir()
  
  # ---- PNG (high-resolution raster) ----
  ggsave(
    filename = file.path(out_dir, paste0(filename_base, ".png")),
    plot     = p,
    width    = width,
    height   = height,
    dpi      = dpi,
    limitsize = FALSE
  )
  
  # ---- SVG (vector graphics) ----
  ggsave(
    filename = file.path(out_dir, paste0(filename_base, ".svg")),
    plot     = p,
    width    = width,
    height   = height,
    limitsize = FALSE
  )
}


############################################################
## FUNCTION 1: GENE UMAP + CLUSTER EXPRESSION (2 PANELS)
############################################################
## PURPOSE:
##   - Visualize expression of ONE gene
##   - Panel 1: UMAP colored by gene expression
##   - Panel 2: Gene expression per Seurat cluster
##
## THIS FUNCTION DOES NOT:
##   - Modify Idents()
##   - Modify DefaultAssay()
##   - Assume cluster names exist
##
## ARGUMENTS:
##   obj      : Seurat object
##   gene     : gene symbol (character)
##   slot     : "data" (normalized), "counts", or "scale.data"
##   pt.size  : point size for UMAP
##   save     : TRUE/FALSE — whether to save to disk
##
## RETURNS:
##   - patchwork plot object (still printed to screen)
############################################################

plot_gene_umap_and_clusters <- function(
    obj,
    gene,
    slot = "data",
    pt.size = 0.4,
    save = TRUE
) {
  
  ## ---- Safety checks ----
  
  # Confirm object type
  stopifnot(inherits(obj, "Seurat"))
  
  # Confirm gene exists in ACTIVE assay
  if (!(gene %in% rownames(obj))) {
    stop(
      paste(
        "Gene", gene,
        "not found in DefaultAssay:",
        DefaultAssay(obj)
      )
    )
  }
  
  # Confirm clusters exist
  if (!("seurat_clusters" %in% colnames(obj@meta.data))) {
    stop("seurat_clusters not found in meta.data")
  }
  
  # Confirm UMAP exists
  if (!("umap" %in% Reductions(obj))) {
    stop("UMAP reduction not found in object")
  }
  
  ## ---- Panel 1: UMAP colored by gene expression ----
  p_gene_umap <- FeaturePlot(
    obj,
    features  = gene,
    reduction = "umap",
    slot      = slot,
    pt.size   = pt.size
  ) +
    ggtitle(paste(gene, "expression (UMAP)"))
  
  ## ---- Panel 2: Gene expression per cluster ----
  p_cluster_expr <- VlnPlot(
    obj,
    features = gene,
    group.by = "seurat_clusters",
    slot     = slot,
    pt.size  = 0
  ) +
    ggtitle(paste(gene, "expression per cluster")) +
    theme_classic()
  
  ## ---- Combine panels ----
  p <- p_gene_umap + p_cluster_expr
  
  ## ---- Save to disk if requested ----
  if (isTRUE(save)) {
    save_plot(
      p,
      filename_base = paste0("gene_", gene, "_umap_and_clusters"),
      width  = 12,
      height = 5
    )
  }
  
  return(p)
}


############################################################
## FUNCTION 2: CLUSTER UMAP + GENE UMAP + CLUSTER EXPRESSION
##             (3 PANELS)
############################################################
## PURPOSE:
##   - Provide full spatial + expression context
##   - Panel 1: UMAP colored by Seurat clusters
##   - Panel 2: UMAP colored by gene expression
##   - Panel 3: Gene expression per cluster
##
## THIS FUNCTION DOES NOT:
##   - Assign biological names
##   - Modify object state
##
## ARGUMENTS:
##   obj      : Seurat object
##   gene     : gene symbol (character)
##   slot     : "data", "counts", or "scale.data"
##   pt.size  : point size for UMAP
##   save     : TRUE/FALSE — whether to save to disk
##
## RETURNS:
##   - patchwork plot object
############################################################

plot_gene_with_umaps_and_clusters <- function(
    obj,
    gene,
    slot = "data",
    pt.size = 0.4,
    save = TRUE
) {
  
  ## ---- Safety checks ----
  stopifnot(inherits(obj, "Seurat"))
  
  if (!(gene %in% rownames(obj))) {
    stop(
      paste(
        "Gene", gene,
        "not found in DefaultAssay:",
        DefaultAssay(obj)
      )
    )
  }
  
  if (!("seurat_clusters" %in% colnames(obj@meta.data))) {
    stop("seurat_clusters not found in meta.data")
  }
  
  if (!("umap" %in% Reductions(obj))) {
    stop("UMAP reduction not found in object")
  }
  
  ## ---- Panel 1: UMAP colored by clusters ----
  p_clusters <- DimPlot(
    obj,
    reduction = "umap",
    group.by  = "seurat_clusters",
    label     = TRUE,
    repel     = TRUE,
    pt.size   = pt.size
  ) +
    ggtitle("UMAP: Seurat clusters")
  
  ## ---- Panel 2: UMAP colored by gene expression ----
  p_gene_umap <- FeaturePlot(
    obj,
    features  = gene,
    reduction = "umap",
    slot      = slot,
    pt.size   = pt.size
  ) +
    ggtitle(paste(gene, "expression (UMAP)"))
  
  ## ---- Panel 3: Gene expression per cluster ----
  p_cluster_expr <- VlnPlot(
    obj,
    features = gene,
    group.by = "seurat_clusters",
    slot     = slot,
    pt.size  = 0
  ) +
    ggtitle(paste(gene, "expression per cluster")) +
    theme_classic()
  
  ## ---- Combine all panels ----
  p <- p_clusters + p_gene_umap + p_cluster_expr
  
  ## ---- Save to disk if requested ----
  if (isTRUE(save)) {
    save_plot(
      p,
      filename_base = paste0("gene_", gene, "_clusters_umap_expression"),
      width  = 18,
      height = 6
    )
  }
  
  return(p)
}

# Make sure DefaultAssay is correct ONCE
DefaultAssay(mge_organoid_div30) <- "RNA"

# 2-panel version
plot_gene_umap_and_clusters(
  mge_organoid_div30,
  gene = "LHX6"
)

# 3-panel version
plot_gene_with_umaps_and_clusters(
  mge_organoid_div30,
  gene = "LHX6"
)




# Safely inspect a serialized RDS object without triggering computation
install.packages("Seurat")
library(Seurat)


## input path to the SeuratObject ##
input_rds_path <- "/nfs/turbo/umms-parent/mgeo_scRNAseq/day30_old/Day30.rds"

# load the Seurat object 
mge_organoid_div30 <- readRDS(input_rds_path)

# basic identity and structure 
mge_organoid_div30mge

# extract meta data of the Seurat object, and print class to console
md = mge_organoid_div30@meta.data
head(md)
ncol(md) # should be 11 for this case

# make the umap if stored in the Seurat Object 
DimPlot(
  mge_organoid_div30, 
  reduction= "umap", 
  group.by = "seurat_clusters"
)

#!/usr/bin/env Rscript
# Bounded feature export from the existing official processed object.
# Does not alter source objects, perform clustering, or discover gates.
suppressPackageStartupMessages({library(SeuratObject);library(Matrix);library(jsonlite)})
root <- '/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results'
source_path <- file.path(root,'samarasinghe_2021_zenodo_processed_object/samarasinghe_2021_zenodo_seurat.rds')
base <- file.path(root,'div30_antecedent_identity_v1')
out <- file.path(base,'external_extra/samarasinghe')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
requested <- jsonlite::fromJSON(file.path(base,'external/cache/genes.json'))
before <- file.info(source_path)[,c('size','mtime')]
message('Loading official processed Samarasinghe RDS read-only: ',source_path)
obj <- readRDS(source_path)
stopifnot('RNA' %in% names(obj@assays))
assay <- obj@assays[['RNA']]
counts <- if ('counts' %in% slotNames(assay)) methods::slot(assay,'counts') else LayerData(obj,assay='RNA',layer='counts')
stopifnot(inherits(counts,'dgCMatrix'))
message('Loaded RNA counts: ',nrow(counts),' genes x ',ncol(counts),' cells; ',length(counts@x),' stored values')
stopifnot(!anyDuplicated(rownames(counts)),!anyDuplicated(colnames(counts)))
stopifnot(all(is.finite(counts@x)),all(counts@x>=0))
fraction_noninteger <- mean(abs(counts@x-round(counts@x))>1e-6)
if(fraction_noninteger>0)stop('RNA counts contain noninteger nonzero values; refusing unverified renormalization')
total <- Matrix::colSums(counts)
ngenes <- if(all(counts@x>0))diff(counts@p) else Matrix::colSums(counts>0)
stopifnot(all(total>0))
meta <- obj@meta.data
stopifnot(identical(rownames(meta),colnames(counts)))
meta$source_cell_id <- colnames(counts)
meta$cell_id <- colnames(counts)
meta$raw_total_counts <- total
meta$raw_n_genes <- ngenes
if('orig.ident' %in% colnames(meta))meta$sample <- as.character(meta$orig.ident)
metadata <- gzfile(file.path(out,'cells.tsv.gz'),'wt',compression=1)
write.table(meta,metadata,sep='\t',quote=TRUE,row.names=FALSE,na='')
close(metadata)
index <- match(requested,rownames(counts));present <- !is.na(index)
availability <- data.frame(gene=requested,assayed=present,source_row_1based=index,stringsAsFactors=FALSE)
write.table(availability,file.path(out,'gene_availability.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='')
selected <- counts[index[present],,drop=FALSE]
selected@x <- log1p(selected@x * rep.int(1e4/total,diff(selected@p)))
stopifnot(all(is.finite(selected@x)),all(selected@x>=0))
existing_data_check <- NA_real_
if('data' %in% slotNames(assay)){
  old <- methods::slot(assay,'data')
  if(identical(dim(old),dim(counts))){
    gi <- which(present)[seq_len(min(20,sum(present)))];ci <- seq_len(min(200,ncol(counts)))
    expected <- as.matrix(counts[index[gi],ci,drop=FALSE]);expected <- log1p(sweep(expected,2,total[ci],'/')*1e4)
    existing_data_check <- max(abs(expected-as.matrix(old[index[gi],ci,drop=FALSE])))
  }
}
message('Exporting ',nrow(selected),'/',length(requested),' requested features, normalized using all ',nrow(counts),' source genes')
writeBin(selected@x,file.path(out,'expression_values.f32'),size=4,endian='little')
writeBin(selected@i,file.path(out,'expression_indices.i32'),size=4,endian='little')
writeBin(selected@p,file.path(out,'expression_indptr.i32'),size=4,endian='little')
write_json(rownames(selected),file.path(out,'present_genes.json'),auto_unbox=FALSE,pretty=TRUE)
for(nm in intersect(c('orig.ident','Time','Genotype','sample','nCount_RNA','nFeature_RNA'),colnames(meta))){
  if(length(unique(meta[[nm]]))<=100){
    z<-as.data.frame(table(meta[[nm]],useNA='ifany'));colnames(z)<-c(nm,'n_cells');write.table(z,file.path(out,paste0('metadata_counts_',nm,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
  }
}
after <- file.info(source_path)[,c('size','mtime')]
stopifnot(identical(before,after))
provenance <- list(status='exported_sparse_arrays',source=source_path,source_size=unname(before$size),source_mtime=as.character(before$mtime),source_unchanged_size_mtime=TRUE,
  assay='RNA',layer='counts',count_matrix_class=class(counts),counts_shape=dim(counts),selected_shape=dim(selected),selected_nnz=length(selected@x),fraction_noninteger_counts=fraction_noninteger,
  normalization='log1p(raw RNA gene counts / total raw RNA counts over ALL source count features * 10000); performed once; requested absent genes remain unassayed/NaN in final NPY',
  original_data_max_abs_difference_small_check=existing_data_check,requested_genes=length(requested),assayed_requested_genes=sum(present),unassayed_genes=requested[!present],
  raw_total_counts_range=range(total),raw_n_genes_range=range(ngenes),metadata_columns=colnames(meta),cell_order='Original RNA counts columns and official object metadata rownames, verified identical',
  cell_filter='None during export; all official processed cells preserved. Downstream control-only primary comparison must use original Genotype metadata.',
  no_large_pipeline=TRUE,no_gate_search=TRUE,source_objects_modified=FALSE,R_version=R.version.string,SeuratObject_version=as.character(packageVersion('SeuratObject')),Matrix_version=as.character(packageVersion('Matrix')))
write_json(provenance,file.path(out,'provenance.json'),auto_unbox=TRUE,pretty=TRUE,na='null')
message('Bounded official-object sparse feature export complete')

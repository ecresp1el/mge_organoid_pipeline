#!/usr/bin/env Rscript
# Fit actual R Hypergate models on exported expression; never fit gates per sample.
suppressPackageStartupMessages(library(hypergate))
suppressPackageStartupMessages(library(jsonlite))
args <- commandArgs(trailingOnly=TRUE)
input <- args[1]; jobs_path <- args[2]; output <- args[3]
runtime_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset="/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder"), mustWork=TRUE)
output <- normalizePath(output, mustWork=FALSE)
if(!startsWith(output, paste0(runtime_root, "/"))) {
 stop("Hypergate fit outputs must be under the external PROJECT_ROOT runtime workspace.")
}
d <- read.delim(gzfile(input), check.names=FALSE, stringsAsFactors=FALSE)
jobs <- read.delim(jobs_path, check.names=FALSE, stringsAsFactors=FALSE)
dir.create(dirname(output), recursive=TRUE, showWarnings=FALSE)
if(file.exists(output)) stop("Refusing to overwrite existing fit checkpoint: ",output)
for(i in seq_len(nrow(jobs))) {
 j <- jobs[i,]; genes <- strsplit(j$features,";",fixed=TRUE)[[1]]
 keep <- !is.na(d[[j$label_column]]) & d[[j$label_column]] != "ambiguous"
 xp <- as.matrix(d[keep,genes,drop=FALSE]); storage.mode(xp)<-"double"
 y <- d[[j$label_column]][keep] == j$target
 start <- proc.time()[3]
 result <- tryCatch({
   hg <- hypergate(xp, as.integer(y), level=1, beta=j$beta, delta_add=0, verbose=FALSE)
   capture <- subset_matrix_hg(hg,xp)
   # Read exact package parameters; independent Python evaluation checks these.
   pars <- tail(hg$pars.history,1)[1,]
   rules <- lapply(hg$active_channels,function(ch) list(
      gene=substr(ch,1,nchar(ch)-4),
      op=if(substr(ch,nchar(ch)-3,nchar(ch))=="_min") ">=" else "<=",
      threshold=unname(pars[ch])))
   list(rules=toJSON(rules,auto_unbox=TRUE,digits=16),tp=sum(capture & y),
        fp=sum(capture & !y),fn=sum(!capture & y),tn=sum(!capture & !y),error="")
 }, error=function(e) list(rules="[]",tp=NA,fp=NA,fn=NA,tn=NA,error=conditionMessage(e)))
 row <- cbind(j,as.data.frame(result,stringsAsFactors=FALSE),elapsed_seconds=proc.time()[3]-start)
 write.table(row,output,sep="\t",row.names=FALSE,col.names=!file.exists(output),
             append=file.exists(output),quote=TRUE,na="")
 if(i%%10==0 || i==nrow(jobs)) {cat(i,"/",nrow(jobs),"completed\n");flush.console()}
}

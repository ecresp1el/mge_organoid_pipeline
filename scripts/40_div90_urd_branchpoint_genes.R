#!/usr/bin/env Rscript
# Downstream only: immutable URD topology, identities, pseudotime and source RDS.
# Native markersAUCPR ranks specificity; markersBinom tests detection frequency.
# Runtime artifacts must be outside the checkout and outside the source run.

parse_args <- function(a) {
  if (identical(a, "--self-test")) return(list(self_test = TRUE))
  if (length(a) %% 2L) stop("Expected --name value pairs")
  out <- as.list(a[seq(2, length(a), 2)])
  names(out) <- sub("^--", "", a[seq(1, length(a), 2)])
  out
}
log_msg <- function(...) { cat(format(Sys.time(), "%FT%T%z"), "pid",Sys.getpid(), ..., "\n"); flush.console() }
needed <- c("URD", "Matrix", "data.table", "ggplot2", "cowplot", "svglite", "digest")
for (p in needed) if (!requireNamespace(p, quietly = TRUE)) stop("Missing installed package: ", p)
suppressPackageStartupMessages({library(URD); library(Matrix); library(data.table); library(ggplot2)})
options(stringsAsFactors = FALSE)

read_settings <- function(path) {
  z <- fread(path); out <- as.list(setNames(z$value, z$parameter))
  for (n in setdiff(names(out), c("pseudotime", "window_fractions", "anchor_fractions"))) out[[n]] <- as.numeric(out[[n]])
  for (n in c("window_fractions", "anchor_fractions")) out[[n]] <- as.numeric(strsplit(out[[n]], ",", fixed = TRUE)[[1]])
  out
}
safe_name <- function(s) gsub("[^A-Za-z0-9_.-]", "_", s)
bind <- function(x) rbindlist(x, fill = TRUE, use.names = TRUE)
write_tsv <- function(x, p) fwrite(x, p, sep = "\t", na = "NA")

# No topology changes: locally select only original exclusive segment memberships.
select_window <- function(pt, a, b, lo, hi, cfg, balanced = FALSE) {
  overlap <- intersect(a, b); a <- setdiff(a, overlap); b <- setdiff(b, overlap)
  # Shift -w followed by +w can leave a 1-ULP numerical sliver at the split.
  # It is a zero-width biological interval, not a testable set of bins.
  tolerance <- 64 * .Machine$double.eps * max(1,abs(lo),abs(hi))
  if (!is.finite(lo) || !is.finite(hi) || hi-lo <= tolerance) return(list(status = "no_common_interval", a = character(), b = character(), overlap = length(overlap)))
  a <- a[is.finite(pt[a]) & pt[a] >= lo & pt[a] <= hi]
  b <- b[is.finite(pt[b]) & pt[b] >= lo & pt[b] <= hi]
  edges <- seq(lo, hi, length.out = cfg$matching_bins + 1L)
  ia <- cut(pt[a], edges, include.lowest = TRUE, labels = FALSE)
  ib <- cut(pt[b], edges, include.lowest = TRUE, labels = FALSE)
  good <- which(tabulate(ia, cfg$matching_bins) >= cfg$minimum_bin_cells & tabulate(ib, cfg$matching_bins) >= cfg$minimum_bin_cells)
  aa <- a[ia %in% good]; bb <- b[ib %in% good]
  if (balanced && length(good)) {
    set.seed(cfg$seed)
    aa <- bb <- character()
    for (k in good) {
      ca <- a[ia == k]; cb <- b[ib == k]; n <- min(length(ca), length(cb))
      aa <- c(aa, ca[sample.int(length(ca), n)])
      bb <- c(bb, cb[sample.int(length(cb), n)])
    }
  }
  status <- if (length(good) < cfg$minimum_joint_bins || min(length(aa), length(bb)) < cfg$minimum_cells) "insufficient_cells_or_bins" else "supported"
  list(status = status, a = aa, b = bb, overlap = length(overlap), joint_bins = length(good), edges = edges)
}

# Recover the enriched-direction native binomial p-value without outcome filters.
# Equal detection has no enriched direction: raw p is NA; BH input is 1.
# Direction selection is conservatively corrected by multiplying native p by 2.
# Native mean/count summaries can drop dimensions on one returned gene; pad only
# that case with a duplicate row in a temporary object, then discard the pad.
native_detection <- function(object, a, b, genes, fa, fb, pn) {
  different <- genes[fa != fb]
  answer <- data.table(gene = genes, native_p_enriched = NA_real_, native_fdr = NA_real_, native_padding = FALSE)
  if (!length(different)) return(answer)
  work <- object
  selected <- unique(c(a, b))
  work@logupx.data <- object@logupx.data[different, selected, drop = FALSE]
  work@count.data <- object@count.data[different, selected, drop = FALSE]
  work@group.ids <- data.frame(bp_native_group=ifelse(selected %in% a,"A","B"),row.names=selected)
  padded <- length(different) == 1L
  if (padded) {
    pad <- "__URD_NATIVE_DIMENSION_PADDING__"
    stopifnot(!pad %in% genes)
    work@logupx.data <- rbind(work@logupx.data, work@logupx.data)
    work@count.data <- rbind(work@count.data, work@count.data)
    rownames(work@logupx.data) <- rownames(work@count.data) <- c(different, pad)
  }
  tested <- URD::markersBinom(work, pseudotime = pn, clust.1 = "A", clust.2 = "B", clustering="bp_native_group", cells.1 = a, cells.2 = b,
    genes.use = rownames(work@logupx.data), effect.size = .Machine$double.eps, p.thresh = 1, frac.must.express = 0)
  result <- data.table(gene = rownames(tested), native_p_enriched = tested$pval, native_fdr = tested$pval.fdr)
  stopifnot(all(different %in% result$gene), !anyDuplicated(result$gene))
  j <- match(answer$gene, result$gene)
  answer[, `:=`(native_p_enriched = result$native_p_enriched[j], native_fdr = result$native_fdr[j], native_padding = padded)]
  answer
}

native_statistics <- function(object, a, b, cfg) {
  started <- Sys.time()
  log_msg("COMPARISON_START", "cells_A",length(a),"cells_B",length(b))
  expr <- object@logupx.data
  fa <- Matrix::rowMeans(expr[, a, drop = FALSE] > 0)
  fb <- Matrix::rowMeans(expr[, b, drop = FALSE] > 0)
  genes <- rownames(expr)[pmax(fa, fb) >= cfg$minimum_detection]
  if (!length(genes)) return(data.table())
  log_msg("GENE_SCREEN_DONE","eligible_genes",length(genes))
  fa <- fa[genes]; fb <- fb[genes]
  # URD logupx is log2(x+1); match native mean.of.logs (linear mean, logged).
  mean_linear <- function(cells) {
    m <- expr[genes, cells, drop = FALSE]; m@x <- 2^m@x - 1
    Matrix::rowMeans(m)
  }
  ma <- mean_linear(a); mb <- mean_linear(b)
  z <- data.table(gene = genes, mean_normalized_A = ma, mean_normalized_B = mb,
    mean_log2upx_A = log2(ma + 1), mean_log2upx_B = log2(mb + 1),
    log2fc_A_vs_B = log2((ma + 1e-6)/(mb + 1e-6)), detection_A = fa, detection_B = fb,
    detection_difference_A_vs_B = fa - fb, n_A = length(a), n_B = length(b))
  log_msg("NATIVE_BINOM_START",length(genes),"genes")
  de <- native_detection(object, a, b, genes, fa, fb, cfg$pseudotime)
  log_msg("NATIVE_BINOM_DONE","elapsed_seconds",round(as.numeric(difftime(Sys.time(),started,units="secs")),2))
  z <- merge(z, de, by = "gene", sort = FALSE)
  z[, detection_test_status := fifelse(is.na(native_p_enriched), "equal_detection_no_directional_test", "native_enriched_direction_times_two")]
  z[, detection_p_direction_corrected := fifelse(is.na(native_p_enriched), 1, pmin(1, 2*native_p_enriched))]
  z[, local_fdr := p.adjust(detection_p_direction_corrected, "BH")]
  for (direction in c("A", "B")) {
    aa <- if (direction == "A") a else b; bb <- if (direction == "A") b else a
    native_fc <- log2((2^round(log2(ma+1),3)-1)/(2^round(log2(mb+1),3)-1))
    eligible <- genes[round(fa,3)!=round(fb,3) & if(direction=="A") native_fc>0 else native_fc<0]
    if(length(eligible)) {
      batch_size <- cfg$gene_batch_size %||% 128L
      batches <- split(eligible,ceiling(seq_along(eligible)/batch_size))
      chunks <- vector("list",length(batches))
      selected <- unique(c(aa,bb))
      for(i in seq_along(batches)) {
        block <- batches[[i]]
        log_msg("NATIVE_AUCPR_BATCH_START","direction",direction,"batch",i,"of",length(batches),"genes",length(block))
        # Native marker calculations are gene-independent. Restricting the
        # transient matrix avoids repeatedly indexing all 20,049 cells and
        # 15,239 sparse rows for each gene; the native method is unchanged.
        small <- methods::new("URD")
        small@logupx.data <- expr[block,selected,drop=FALSE]
        chunks[[i]] <- URD::markersAUCPR(small, cells.1 = aa, cells.2 = bb, genes.use = block,
          effect.size = 0, frac.must.express = 0, frac.min.diff = 0, auc.factor = 0, max.auc.threshold = 1)
        log_msg("NATIVE_AUCPR_BATCH_DONE","direction",direction,"batch",i,"of",length(batches),
          "elapsed_seconds",round(as.numeric(difftime(Sys.time(),started,units="secs")),2))
      }
      au <- do.call(rbind,chunks)
    } else au <- data.frame(AUCPR=numeric(),AUCPR.ratio=numeric())
    idx <- match(z$gene, rownames(au)); baseline <- length(aa)/(length(aa) + length(bb))
    z[, (paste0("auc_", direction)) := au$AUCPR[idx]]
    z[, (paste0("auc_ratio_", direction)) := au$AUCPR.ratio[idx]]
    z[, (paste0("auc_baseline_", direction)) := baseline]
  }
  z
}

self_test <- function() {
  cfg <- list(matching_bins = 5, minimum_bin_cells = 2, minimum_joint_bins = 3, minimum_cells = 10, seed = 7,
    pseudotime = "pt", minimum_detection = .1)
  ca <- paste0("a", 1:50); cb <- paste0("b", 1:50); pt <- setNames(rep(seq(.01,.99,length.out=50),2),c(ca,cb))
  s <- select_window(pt, ca, cb, 0, 1, cfg, TRUE)
  stopifnot(s$status == "supported", length(s$a) == length(s$b), !length(intersect(s$a,s$b)))
  stopifnot(select_window(pt, ca, cb, 1, 0, cfg)$status == "no_common_interval")
  stopifnot(select_window(pt, ca, cb, .272363130443725, .272363130443725+1e-16, cfg)$status == "no_common_interval")
  m <- matrix(0, 4, 100, dimnames = list(c("earlyA","earlyB","equal","absent"),c(ca,cb)))
  m[1,1:40] <- 6; m[1,51:55] <- 1; m[2,1:5] <- 1; m[2,51:90] <- 6; m[3,c(1:25,51:75)] <- 2
  obj <- methods::new("URD"); obj@count.data <- as(m,"dgCMatrix"); obj@logupx.data <- as(log2(m+1),"dgCMatrix")
  out <- native_statistics(obj, ca, cb, cfg)
  stopifnot(out[gene=="earlyA", log2fc_A_vs_B] > 0, out[gene=="earlyB", log2fc_A_vs_B] < 0,
    out[gene=="earlyA", detection_p_direction_corrected] < .05,
    out[gene=="equal", detection_p_direction_corrected] == 1, !"absent" %in% out$gene)
  cfg$gene_batch_size <- 1
  batched <- native_statistics(obj,ca,cb,cfg)
  comparison <- all.equal(as.data.frame(out),as.data.frame(batched),check.attributes=FALSE)
  if(!isTRUE(comparison))stop("Batch equivalence check: ",paste(comparison,collapse="; "))
  direct <- URD::markersAUCPR(obj,cells.1=ca,cells.2=cb,genes.use="earlyA",effect.size=0,frac.must.express=0,frac.min.diff=0,auc.factor=0)
  stopifnot(isTRUE(all.equal(out[gene=="earlyA",auc_A],direct["earlyA","AUCPR"])))
  fa <- setNames(.8,"earlyA"); fb <- setNames(.1,"earlyA")
  one <- native_detection(obj, ca, cb, "earlyA", fa, fb, "pt")
  stopifnot(is.finite(one$native_p_enriched), one$native_padding)
  expected <- pbinom(40,50,.1,lower.tail=FALSE)+dbinom(40,50,.1)
  stopifnot(isTRUE(all.equal(one$native_p_enriched,expected)))
  obj@logupx.data[,cb] <- obj@logupx.data[,ca]
  obj@count.data[,cb] <- obj@count.data[,ca]
  equal <- native_statistics(obj,ca,cb,cfg)
  stopifnot(all(equal$detection_p_direction_corrected==1),all(is.na(equal$auc_A)),all(is.na(equal$auc_B)))
  cat("SELF_TEST PASS: window matching, directionality, equality, sparse single-gene native edge case, native p-value recovery\n")
}

main <- function(opt) {
  if (isTRUE(opt$self_test)) {self_test(); return(invisible(NULL))}
  required <- c("tree-rds","outdir","config","tf-reference")
  stopifnot(all(required %in% names(opt)))
  source <- normalizePath(opt[["tree-rds"]], mustWork=TRUE); outdir <- normalizePath(opt$outdir,mustWork=TRUE)
  # Runtime must be a separate Turbo package, never the protected run or repo.
  runtime_parent <- "/nfs/turbo/umms-parent/mgeo_neuron_scrnaseq_projectfolder/results/div90_urd_branchpoint_genes/"
  stopifnot(startsWith(paste0(outdir,"/"),runtime_parent), !startsWith(source,paste0(outdir,"/")))
  cfg <- read_settings(opt$config); cores <- as.integer(opt$cores %||% "1")
  if (file.exists(file.path(outdir,"COMPLETED.txt"))) stop("Refusing completed output overwrite")
  for (d in c("tables","plots","checkpoints","provenance")) dir.create(file.path(outdir,d),showWarnings=FALSE)
  set.seed(cfg$seed)
  log_msg("Reading protected saved tree",source)
  object <- readRDS(source)
  stopifnot(inherits(object,"URD"),identical(dimnames(object@count.data),dimnames(object@logupx.data)))
  genes <- rownames(object@logupx.data); cells <- colnames(object@logupx.data)
  stopifnot(!is.null(genes),!is.null(cells), !anyDuplicated(genes),!anyDuplicated(cells))
  pt <- setNames(object@pseudotime[cells,cfg$pseudotime],cells)
  stopifnot(all(is.finite(pt)),max(abs(object@tree$pseudotime[cells]-pt)) < 1e-10)
  joins <- as.data.table(object@tree$segment.joins); joins[, `:=`(parent=as.character(parent),child=as.character(child),pseudotime=as.numeric(pseudotime))]
  memberships <- lapply(object@tree$cells.in.segment, function(x) intersect(as.character(x),cells))
  get_children <- function(s) as.character(URD::segChildren(object,s))
  descendants <- function(s) unique(c(s,as.character(URD::segChildrenAll(object,s))))
  terminal <- setdiff(joins$child,joins$parent)
  # Terminal names come from the immutable source run's mapping table.
  map_path <- file.path(dirname(source),"tables","tree_tip_mapping.tsv")
  tipmap <- fread(map_path,colClasses="character")
  labels <- setNames(tipmap$paper_cluster_annotation,tipmap$tip_id)
  leaves <- function(s) intersect(descendants(s),terminal)
  label <- function(s) paste(labels[leaves(s)],collapse=" + ")
  limits <- object@tree$segment.pseudotime.limits
  limit <- function(s,k) as.numeric(limits[s,k])
  parents <- unique(joins[order(pseudotime),parent]); R <- diff(range(pt))
  focused <- isTRUE(cfg$single_comparison==1)
  analysis_parents <- if(focused) as.character(cfg$focus_parent) else parents
  branchpoints <- bind(lapply(parents,function(p) {
    kids <- get_children(p); stopifnot(setequal(kids,joins[parent==p,child]))
    data.table(branchpoint=paste0("BP",p),parent=p,daughter=kids,split_time=unique(joins[parent==p,pseudotime]),
      daughter_lineage=vapply(kids,label,""),n_native_cells=vapply(kids,function(k) length(memberships[[k]]),integer(1)),
      parent_cells=length(memberships[[p]]))
  }))
  write_tsv(branchpoints,file.path(outdir,"tables/branchpoints.tsv"))
  write_tsv(data.table(cell=cells,sample=object@meta[cells,"orig.ident"],assigned_to_tree=cells %in% unique(unlist(memberships))),file.path(outdir,"tables/cell_assignment_audit.tsv"))
  ref <- fread(opt[["tf-reference"]],check.names=FALSE)
  stopifnot(all(c("HGNC symbol","Is TF?","DBD") %in% names(ref)))
  tf <- unique(ref[,.(gene=`HGNC symbol`,is_tf=tolower(`Is TF?`) %in% c("yes","true"),tf_family=DBD)],by="gene")
  tf[,tf_reference := "Lambert_HumanTFs_v1.01"]
  write_tsv(merge(data.table(gene=genes),tf,by="gene",all.x=TRUE),file.path(outdir,"provenance/tf_gene_mapping.tsv"))
  requests <- list(); selections <- list(); support <- list(); cell_rows <- list()
  add_request <- function(p,aa,bb,leaf_a,leaf_b,stage,wname,anchor,lo,hi,requested_lo,requested_hi,ca,cb) {
    for (variant in if(focused) "native_all_cells" else c("native_all_cells","balanced_bins")) {
      selected <- select_window(pt,ca,cb,lo,hi,cfg,variant=="balanced_bins")
      id <- paste0("C",sprintf("%04d",length(requests)+1L))
      key <- if(selected$status=="supported") digest::digest(list(sort(selected$a),sort(selected$b)),algo="sha256") else NA_character_
      row <- data.table(comparison_id=id,selection_id=key,branchpoint=paste0("BP",p),parent=p,
        daughter_A=aa,daughter_B=bb,daughter_A_label=label(aa),daughter_B_label=label(bb),
        leaf_A=leaf_a,leaf_B=leaf_b,stage=stage,window=wname,anchor_fraction=anchor,variant=variant,
        split_time=unique(joins[parent==p,pseudotime]),requested_low=requested_lo,requested_high=requested_hi,
        time_low=lo,time_high=hi,status=selected$status,n_selected_A=length(selected$a),n_selected_B=length(selected$b),
        joint_bins=selected$joint_bins %||% 0L,excluded_overlap=selected$overlap,
        observed_min_A=if(length(selected$a))min(pt[selected$a]) else NA_real_,
        observed_max_A=if(length(selected$a))max(pt[selected$a]) else NA_real_,
        observed_min_B=if(length(selected$b))min(pt[selected$b]) else NA_real_,
        observed_max_B=if(length(selected$b))max(pt[selected$b]) else NA_real_)
      requests[[length(requests)+1L]] <<- row
      if (selected$status=="supported" && is.null(selections[[key]])) selections[[key]] <<- selected[c("a","b")]
      for (side in c("a","b")) {
        cc <- selected[[side]]
        if(length(cc)) cell_rows[[length(cell_rows)+1L]] <<- data.table(comparison_id=id,side=toupper(side),cell=cc,pseudotime=pt[cc],sample=object@meta[cc,"orig.ident"])
      }
    }
  }
  path_segments <- function(start,leaf) {
    path <- leaf; cursor <- leaf
    while(cursor != start) {
      cursor <- joins[child==cursor,parent]
      if(length(cursor)!=1L) stop("Invalid path")
      path <- c(cursor,path)
    }; path
  }
  log_msg("Enumerating original daughter branches and window support")
  for(p in analysis_parents) {
    b <- unique(joins[parent==p,pseudotime]); kids <- get_children(p)
    for(pair in combn(kids,2,simplify=FALSE)) {
      aa <- pair[1]; bb <- pair[2]
      if(focused && !setequal(pair,as.character(c(cfg$focus_daughter_A,cfg$focus_daughter_B))))next
      local_end <- min(limit(aa,2),limit(bb,2))
      for(anchor in cfg$anchor_fractions) for(w in cfg$window_fractions) {
        pivot <- b+anchor*R; lo <- max(b,pivot,limit(aa,1),limit(bb,1)); hi <- min(pivot+w*R,local_end)
        add_request(p,aa,bb,NA_character_,NA_character_,"local",paste0("w",w),anchor,lo,hi,pivot,pivot+w*R,memberships[[aa]],memberships[[bb]])
      }
      if(focused)next
      # Late tests remain on separate daughter-to-tip paths; never pool fates.
      late_start <- min(b+max(cfg$window_fractions)*R,local_end)
      for(la in leaves(aa)) for(lb in leaves(bb)) {
        ca <- unique(unlist(memberships[path_segments(aa,la)])); cb <- unique(unlist(memberships[path_segments(bb,lb)]))
        overlap <- intersect(ca,cb);ca <- setdiff(ca,overlap);cb <- setdiff(cb,overlap)
        common_end <- min(max(pt[ca]),max(pt[cb]))
        if(common_end<=late_start) next
        edges <- seq(late_start,common_end,length.out=cfg$late_bins+1L)
        for(k in seq_len(cfg$late_bins)) add_request(p,aa,bb,la,lb,"late",paste0("late",k),0,
          edges[k],edges[k+1L],edges[k],edges[k+1L],ca,cb)
      }
    }
  }
  requests <- bind(requests)
  requests[,duplicate_selection:=!is.na(selection_id) & duplicated(selection_id)]
  write_tsv(requests,file.path(outdir,"tables/window_support.tsv"))
  memberships_export <- bind(cell_rows);fwrite(memberships_export,file.path(outdir,"tables/comparison_cells.tsv.gz"),sep="\t")
  write_tsv(memberships_export[,.(n_cells=.N),by=.(comparison_id,side,sample)],file.path(outdir,"tables/comparison_sample_counts.tsv"))
  rm(cell_rows,memberships_export);gc()
  log_msg("Supported requests",sum(requests$status=="supported"),"unique native comparisons",length(selections))
  work <- function(key) {
    cache <- file.path(outdir,"checkpoints",paste0(key,".rds"))
    if(file.exists(cache)) return(readRDS(cache))
    s <- selections[[key]]; started <- Sys.time()
    z <- native_statistics(object,s$a,s$b,cfg)
    z[,selection_id:=key]
    saveRDS(z,cache)
    log_msg("Completed native comparison",substr(key,1,10),"genes",nrow(z),"seconds",round(as.numeric(difftime(Sys.time(),started,units="secs"))))
    z
  }
  pieces <- parallel::mclapply(names(selections),work,mc.cores=cores,mc.preschedule=FALSE)
  if(any(vapply(pieces,inherits,logical(1),"try-error"))) stop("Native comparison failed; checkpoints retained")
  stats <- bind(pieces); rm(pieces);gc()
  stopifnot(nrow(stats)>0)
  stats[,global_fdr:=p.adjust(detection_p_direction_corrected,"BH")]
  saveRDS(stats,file.path(outdir,"checkpoints/unique_gene_statistics.rds"))
  # Two direction rows share one direction-corrected significance test per gene.
  detail <- merge(requests[status=="supported"],stats,by="selection_id",allow.cartesian=TRUE)
  results <- bind(lapply(c("A","B"),function(side) {
    d <- copy(detail);sgn <- if(side=="A") 1 else -1
    d[,`:=`(enriched_side=side,reported_daughter=get(paste0("daughter_",side)),
      log2fc=sgn*log2fc_A_vs_B,detection_difference=sgn*detection_difference_A_vs_B,
      AUCPR=get(paste0("auc_",side)),AUCPR_ratio=get(paste0("auc_ratio_",side)),AUCPR_baseline=get(paste0("auc_baseline_",side)))]
    d
  }));rm(detail);gc()
  results <- merge(results,tf,by="gene",all.x=TRUE,sort=FALSE)
  results[is.na(is_tf),`:=`(is_tf=FALSE,tf_family="unmapped_or_nonTF",tf_reference="Lambert_HumanTFs_v1.01")]
  results[,`:=`(effect_supported=log2fc>=cfg$minimum_log2fc & detection_difference>=cfg$minimum_detection_difference & !is.na(AUCPR_ratio) & AUCPR_ratio>=cfg$minimum_auc_ratio,
    detection_significant=global_fdr<=cfg$maximum_fdr)]
  results[,candidate:=effect_supported & detection_significant]
  # Summaries use distinct original-anchor widths, avoiding duplicated clipped windows.
  local <- unique(results[stage=="local" & anchor_fraction==0 & variant=="native_all_cells"],by=c("branchpoint","daughter_A","daughter_B","reported_daughter","gene","selection_id"))
  group <- c("branchpoint","daughter_A","daughter_B","reported_daughter","gene")
  timing <- local[,.(supported_widths=.N,early_effect_widths=sum(effect_supported),early_candidate_widths=sum(candidate),
    smallest_width_effect=effect_supported[which.min(as.numeric(sub("w","",window)))],
    max_early_log2fc=max(log2fc),best_local_auc_ratio=suppressWarnings(max(AUCPR_ratio,na.rm=TRUE))),by=group]
  timing[!is.finite(best_local_auc_ratio),best_local_auc_ratio:=NA_real_]
  late <- results[stage=="late" & variant=="native_all_cells"]
  late_summary <- late[,.(later_effect_bins=sum(effect_supported),later_candidate_bins=sum(candidate)),
    by=c(group,"leaf_A","leaf_B")][,.(later_effect_bins=max(later_effect_bins),later_candidate_bins=max(later_candidate_bins)),by=group]
  timing <- merge(timing,late_summary,by=group,all=TRUE)
  for(n in c("supported_widths","early_effect_widths","early_candidate_widths","later_effect_bins","later_candidate_bins")) set(timing,which(is.na(timing[[n]])),n,0)
  timing[,timing_class:="window_sensitive_or_unresolved"]
  timing[early_effect_widths>=2 & smallest_width_effect==TRUE,timing_class:="early_associated"]
  timing[timing_class=="early_associated" & later_effect_bins>=2,timing_class:="early_and_sustained"]
  timing[supported_widths>=2 & max_early_log2fc<cfg$minimum_log2fc & later_effect_bins>=2,timing_class:="late_emerging"]
  if(focused) timing[early_effect_widths>=1,timing_class:="near_split_difference_single_window"]
  sensitivity <- results[stage=="local",.(anchor_effect_fraction=mean(effect_supported),balanced_effect_fraction=mean(effect_supported[variant=="balanced_bins"])),by=group]
  timing <- merge(timing,sensitivity,by=group,all.x=TRUE)
  onset <- results[effect_supported==TRUE & variant=="native_all_cells" & anchor_fraction==0,
    .(first_supported_low=min(time_low),first_supported_high=time_high[which.min(time_low)]),by=group]
  timing <- merge(timing,onset,by=group,all.x=TRUE)
  write_tsv(timing,file.path(outdir,"tables/timing_window_stability.tsv"))
  results <- merge(results,timing,by=group,all.x=TRUE,sort=FALSE)
  # All-sibling labels require consistent evidence against every actual sibling.
  required_widths <- if(focused)1L else 2L
  sibling <- unique(results[stage=="local" & anchor_fraction==0 & variant=="native_all_cells" & candidate==TRUE & early_candidate_widths>=required_widths,
    .(branchpoint,reported_daughter,gene,daughter_A,daughter_B)])
  sibling <- sibling[,.(siblings_beaten=.N),by=.(branchpoint,reported_daughter,gene)]
  sibling[,siblings_required:=vapply(branchpoint,function(bp)nrow(branchpoints[branchpoint==bp])-1L,integer(1))]
  sibling[,all_siblings:=siblings_beaten>=siblings_required]
  results <- merge(results,sibling,by=c("branchpoint","reported_daughter","gene"),all.x=TRUE,sort=FALSE)
  results[is.na(all_siblings),all_siblings:=FALSE]
  ranking <- unique(results[anchor_fraction==0 & variant=="native_all_cells" & log2fc>0],by=c(group,"selection_id"))
  setorderv(ranking,c(group[1:4],"early_candidate_widths","early_effect_widths","AUCPR_ratio","log2fc"),c(rep(1,4),-1,-1,-1,-1),na.last=TRUE)
  ranking <- unique(ranking,by=group)
  rank_group <- group[1:4]
  ranking[,all_gene_rank:=seq_len(.N),by=rank_group]
  ranking[is_tf==TRUE,tf_rank:=seq_len(.N),by=rank_group]
  ranks <- ranking[,c(group,"all_gene_rank","tf_rank"),with=FALSE]
  results <- merge(results,ranks,by=group,all.x=TRUE,sort=FALSE)
  write_tsv(ranking,file.path(outdir,"tables/ranked_divergence_genes.tsv"))
  write_tsv(ranking[is_tf==TRUE],file.path(outdir,"tables/ranked_candidate_tfs.tsv"))
  write_tsv(results[stage=="late" & candidate==TRUE],file.path(outdir,"tables/ranked_late_genes.tsv"))
  # Final timestamps and links are shared by CSV and detailed tables.
  stamp <- format(Sys.time(),"%FT%T%z");run_id <- basename(outdir)
  results[,`:=`(run_id=run_id,generated_at=stamp,URD_version=as.character(packageVersion("URD")),
    tree_plot="plots/urd_tree_branchpoints.pdf",trajectory_plot=paste0("plots/",branchpoint,"/gene_trajectories.pdf"),
    heatmap_plot=paste0("plots/",branchpoint,"/divergence_heatmap.pdf"))]
  insufficient <- requests[status!="supported"]
  insufficient[,`:=`(gene=NA_character_,run_id=run_id,generated_at=stamp)]
  results <- bind(list(results,insufficient))
  fwrite(results,file.path(outdir,"DIV90_URD_gene_results.csv"),na="NA")
  # Spreadsheet-sized entry point: one best representative per ranked
  # gene/branchpoint/daughter contrast, retaining timing/window robustness.
  summary_csv <- copy(ranking)
  summary_csv[,`:=`(run_id=run_id,generated_at=stamp,full_results="DIV90_URD_gene_results.csv",
    trajectory_plot=paste0("plots/",branchpoint,"/gene_trajectories.pdf"),
    heatmap_plot=paste0("plots/",branchpoint,"/divergence_heatmap.pdf"))]
  fwrite(summary_csv,file.path(outdir,"DIV90_URD_gene_summary.csv"),na="NA")
  fwrite(results,file.path(outdir,"tables/branchpoint_gene_statistics.tsv.gz"),sep="\t",na="NA")
  saveRDS(list(branchpoints=branchpoints,ranking=ranking,requests=requests,timing=timing),file.path(outdir,"checkpoints/summary_tables.rds"))
  log_msg("Master CSV written",nrow(results),"rows; preparing selected gene trajectories")
  # Observed expression along EACH original segment; do not pool future descendants.
  selected_by_bp <- lapply(unique(branchpoints$branchpoint),function(bp) {
    rr <- ranking[branchpoint==bp]
    unique(c(head(unique(rr[order(all_gene_rank),gene]),cfg$top_plot_genes),
      head(unique(rr[is_tf==TRUE][order(tf_rank),gene]),cfg$top_plot_tfs)))
  });names(selected_by_bp)<-unique(branchpoints$branchpoint)
  trajectories <- list()
  for(bp in names(selected_by_bp)) {
    parent <- sub("BP","",bp);gs <- selected_by_bp[[bp]]
    if(!length(gs)) next
    segments <- unique(c(parent,unlist(lapply(get_children(parent),descendants))))
    for(segment in segments) {
      cc <- memberships[[segment]];if(length(cc)<2L)next
      if(diff(range(pt[cc]))<=0)next
      edges <- seq(min(pt[cc]),max(pt[cc]),length.out=cfg$trajectory_bins+1L)
      bins <- cut(pt[cc],edges,include.lowest=TRUE,labels=FALSE)
      for(k in sort(unique(bins))) {
        ids<-cc[bins==k]; if(length(ids)<cfg$minimum_bin_cells)next
        for(sample in c("ALL",sort(unique(object@meta[ids,"orig.ident"])))) {
          use<-if(sample=="ALL")ids else ids[object@meta[ids,"orig.ident"]==sample]
          if(length(use)<2L)next
          mm<-object@logupx.data[gs,use,drop=FALSE]
          trajectories[[length(trajectories)+1L]]<-data.table(branchpoint=bp,segment=segment,sample=sample,gene=gs,
            pseudotime=mean(pt[use]),time_low=min(pt[use]),time_high=max(pt[use]),n_cells=length(use),
            mean_log2upx=Matrix::rowMeans(mm),detection=Matrix::rowMeans(mm>0))
        }
      }
    }
  }
  trajectories<-bind(trajectories)
  fwrite(trajectories,file.path(outdir,"tables/expression_trajectories.tsv.gz"),sep="\t")
  plot_assets <- list(tree=object@tree$tree.layout,segments=object@tree$segment.layout,
    cell_layout=object@tree$cell.layout,branchpoints=branchpoints,joins=joins,limits=limits,
    labels=labels,memberships=memberships,pt=pt,ranking=ranking,trajectories=trajectories,requests=requests,
    cfg=cfg,run_id=run_id,stamp=stamp,source=source,analysis_branchpoints=paste0("BP",analysis_parents))
  saveRDS(plot_assets,file.path(outdir,"checkpoints/plot_assets.rds"))
  render_plots(plot_assets,outdir)
  summary <- c(paste("Completed",stamp),paste("Run",run_id),paste("Source",source),
    paste("Analyzed branchpoints",length(analysis_parents),"supported requests",sum(requests$status=="supported"),"unique comparisons",length(selections)),
    paste("Master CSV rows",nrow(results),"bytes",file.info(file.path(outdir,"DIV90_URD_gene_results.csv"))$size),
    if(nrow(results)>1048575)"WARNING: master CSV exceeds Excel's worksheet row limit; filter/import with a data tool." else "Master CSV fits Excel's worksheet row limit.",
    "Native binomial p-values test detection; enriched-direction p multiplied by two before local/global BH.",
    "Global BH family counts each unique cell selection/gene once, including balanced and shifted windows; mirrored direction rows share that test.",
    "AUCPR is marker specificity, not a p-value. Equal detection: raw native p NA, conservative BH input 1.",
    "Statistical support is exploratory cell-level evidence; not causal or donor-level validation.",
    "Timing rules use descriptive effect support; inspect global FDR, sample curves and window support separately.",
    "No fitted impulse curves are required: trajectories show observed binned expression; later lineages stay separate.")
  for(bp in paste0("BP",analysis_parents)) {
    tf_hits<-unique(ranking[branchpoint==bp & is_tf==TRUE & candidate==TRUE,gene])
    summary<-c(summary,paste(bp,"TFs passing candidate thresholds:",if(length(tf_hits))paste(head(tf_hits,12),collapse=", ")else"none; see ranked TF table for descriptive results"))
  }
  if(focused)summary<-c(summary,"ONE comparison only: native BP12 OPC (segment 7) versus neuronal branch (segment 11), 5% pseudotime window. No window/anchor/balancing/late DE sweep. Later curves are descriptive; single-window findings do not establish robust early/late timing.")
  writeLines(summary,file.path(outdir,"analysis_summary.txt"))
  writeLines(capture.output(sessionInfo()),file.path(outdir,"provenance/sessionInfo.txt"))
  writeLines(c(stamp,"Analysis and plots completed; source hash check pending in job wrapper."),file.path(outdir,"ANALYSIS_COMPLETE.txt"))
  log_msg("Completed",outdir)
}

`%||%` <- function(a,b) if(is.null(a)) b else a

render_plots <- function(a,outdir) {
  cfg<-a$cfg
  subtitle<-paste(a$run_id,a$stamp,sep=" | ")
  save_set<-function(p,path,w=12,h=8) {
    dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
    ggsave(paste0(path,".png"),p,width=w,height=h,dpi=cfg$png_dpi,bg="white",limitsize=FALSE)
    ggsave(paste0(path,".pdf"),p,width=w,height=h,device=grDevices::cairo_pdf,limitsize=FALSE)
    ggsave(paste0(path,".svg"),p,width=w,height=h,device=svglite::svglite,limitsize=FALSE)
  }
  layout<-as.data.table(a$cell_layout)
  assigned<-bind(lapply(names(a$memberships),function(s)data.table(cell=a$memberships[[s]],segment=s)))
  layout<-merge(layout,assigned,by="cell",allow.cartesian=TRUE)
  segs<-sort(unique(assigned$segment));palette<-setNames(grDevices::hcl.colors(length(segs),"Dark 3"),segs)
  legend_labels<-setNames(paste0("s",segs,": ",ifelse(segs %in% names(a$labels),gsub("^tip_","",a$labels[segs]),"shared segment")),segs)
  # Cell coordinates come directly from the saved object; joins are drawn only
  # between original segment endpoints; no embedding/layout is recomputed.
  centers<-layout[,.(x=median(x),y=median(y)),by=segment]
  tree<-ggplot()+geom_segment(data=as.data.table(a$tree),aes(x=x1,y=y1,xend=x2,yend=y2),color="grey35",linewidth=.45)+
    geom_point(data=layout,aes(x,y,color=segment),size=.3,alpha=.4)+
    scale_color_manual(values=palette)+theme_void()+theme(legend.position="none")+
    labs(title=if(isTRUE(cfg$single_comparison==1))"Saved seven-tip URD context: only BP12 OPC versus neuronal branch tested" else "DIV90 saved URD tree: four inferred branchpoints",subtitle=subtitle)
  label_rows<-bind(lapply(unique(a$branchpoints$parent),function(s) {
    kids<-a$branchpoints[parent==s,daughter]
    # Find the existing parent/daughter junction from layout boundary cells.
    cc<-layout[segment %in% kids]
    cc[,pt:=a$pt[cell]]
    q<-cc[pt<=quantile(pt,.1)]
    # Use an actual saved edge endpoint at the inferred split. The cell
    # boundary is only a tiebreaker when several endpoints share that time.
    edges<-as.data.table(a$tree)
    endpoints<-unique(rbind(data.table(x=edges$x1,y=edges$y1),data.table(x=edges$x2,y=edges$y2)))
    split<-unique(a$branchpoints[parent==s,split_time])
    axis<-if(abs(cor(layout$x,a$pt[layout$cell],use="complete.obs"))>abs(cor(layout$y,a$pt[layout$cell],use="complete.obs")))"x" else "y"
    candidates<-endpoints[abs(get(axis)-split)<1e-7]
    if(!nrow(candidates)) {
      # Saved edges use node-mean times and need not contain an endpoint at
      # exactly the inferred split. Place a callout at its native segment x
      # and inferred time without moving any edge/cell.
      seg<-as.data.table(a$segments)
      candidates<-data.table(x=as.numeric(seg[as.character(segment)==s,x]),y=split)
    }
    nearest<-which.min((candidates$x-median(q$x))^2+(candidates$y-median(q$y))^2)
    data.table(x=candidates$x[nearest],y=candidates$y[nearest],label=paste0("BP",s,"\n",paste(kids,collapse=", ")))
  }))
  tree<-tree+geom_label(data=label_rows,aes(x,y,label=label),inherit.aes=FALSE,size=3,color="black",fill="white",alpha=.85)
  key<-a$branchpoints[,paste0(branchpoint," at ",round(split_time,4),": segment ",daughter," = ",daughter_lineage)]
  tree<-tree+labs(caption=paste(c(key,"Labels mark inferred splits, not biological commitment times. Original cell coordinates preserved."),collapse="\n"))+
    theme(plot.caption=element_text(hjust=0,size=7))
  save_set(tree,file.path(outdir,"plots/urd_tree_branchpoints"),14,11)
  grDevices::cairo_pdf(file.path(outdir,"DIV90_URD_gene_overview.pdf"),width=14,height=11,onefile=TRUE)
  on.exit(grDevices::dev.off(),add=TRUE);print(tree)
  for(bp in a$analysis_branchpoints %||% unique(a$branchpoints$branchpoint)) {
    dir<-file.path(outdir,"plots",bp);rr<-a$ranking[branchpoint==bp]
    split<-unique(a$branchpoints[branchpoint==bp,split_time])
    selected<-head(unique(rr[order(all_gene_rank),gene]),cfg$top_plot_genes)
    selected_tf<-head(unique(rr[is_tf==TRUE][order(tf_rank),gene]),cfg$top_plot_tfs)
    if(!length(selected)) {
      p<-ggplot()+annotate("text",x=0,y=0,label=paste(bp,"No supported local genes; inspect window_support.tsv"))+theme_void()+labs(title=bp,subtitle=subtitle)
      print(p);save_set(p,file.path(dir,"insufficient_support"));next
    }
    curve<-function(gs,title) {
      d<-a$trajectories[branchpoint==bp & sample=="ALL" & gene %in% gs]
      ggplot(d,aes(pseudotime,mean_log2upx,color=segment,group=segment))+
        annotate("rect",xmin=rep(split,length(cfg$window_fractions)),xmax=split+cfg$window_fractions*diff(range(a$pt)),
          ymin=-Inf,ymax=Inf,alpha=.035,fill="grey30")+
        geom_line(linewidth=.5)+geom_point(aes(size=n_cells),alpha=.7)+
        geom_vline(xintercept=split,linetype=2,color="black")+
        scale_color_manual(values=palette,labels=legend_labels)+scale_size_continuous(range=c(.5,2))+
        facet_wrap(~gene,scales="free_y",ncol=4)+theme_bw(base_size=9)+
        labs(title=title,subtitle=subtitle,x="Original pseudotime",y="Mean log2(UPX + 1)",
          caption="Observed bin means; colors denote original segments. Parent cells are shared context, not assigned future fates.")
    }
    p<-curve(selected,paste(bp,"divergence-associated genes (ranked; inspect FDR)"));print(p);save_set(p,file.path(dir,"gene_trajectories"),14,11)
    if(length(selected_tf)){p<-curve(selected_tf,paste(bp,"candidate DNA-binding transcription factors"));print(p);save_set(p,file.path(dir,"tf_trajectories"),14,11)}
    d<-copy(a$trajectories[branchpoint==bp & sample=="ALL" & gene %in% selected])
    d[,z:=if(sd(mean_log2upx)>0)as.numeric(scale(mean_log2upx)) else 0,by=gene]
    d[,z:=pmax(-2,pmin(2,z))]
    p<-ggplot(d,aes(pseudotime,gene,fill=z))+geom_tile()+facet_grid(~segment,scales="free_x",space="free_x")+
      scale_fill_gradient2(low="#2166ac",mid="white",high="#b2182b",limits=c(-2,2))+
      theme_bw(base_size=9)+labs(title=paste(bp,"expression across original daughter segments"),subtitle=subtitle,
        x="Original pseudotime",y="Gene",fill="Gene z score",caption="One shared gene-wise scale across segments; clipped at ±2. Color is descriptive, not significance.")
    print(p);save_set(p,file.path(dir,"divergence_heatmap"),14,9)
    # Sample curves expose composition effects rather than imply donor validation.
    sd<-a$trajectories[branchpoint==bp & sample!="ALL" & gene %in% head(selected,8)]
    p<-ggplot(sd,aes(pseudotime,mean_log2upx,color=segment,group=interaction(segment,sample)))+
      geom_line(alpha=.5)+facet_wrap(~gene,scales="free_y",ncol=4)+scale_color_manual(values=palette,labels=legend_labels)+theme_bw(base_size=9)+
      labs(title=paste(bp,"per-sample trajectories"),subtitle=subtitle,y="Mean log2(UPX + 1)",x="Original pseudotime")
    print(p);save_set(p,file.path(dir,"sample_trajectories"),14,9)
    rw<-rr[gene %in% head(selected,16),.(early_effect_widths=max(early_effect_widths)),by=.(gene,reported_daughter)]
    p<-ggplot(rw,aes(reorder(gene,early_effect_widths),early_effect_widths,fill=reported_daughter))+
      geom_col(position="dodge")+coord_flip()+theme_bw(base_size=9)+
      labs(title=paste(bp,"distinct local windows supporting expression differences"),subtitle=subtitle,
        x="Gene",y="Distinct original-anchor widths",caption="Descriptive robustness; significance is reported separately in the master CSV.")
    print(p);save_set(p,file.path(dir,"window_robustness"),12,8)
  }
}

if (sys.nframe() == 0L) main(parse_args(commandArgs(trailingOnly=TRUE)))

# Input counts are already summed within independent biological samples/libraries.
# Nuclei/cells are not treated as independent replicates here.
pseudobulk_edgeR <- function(counts, metadata, control) {
  stopifnot(!anyNA(counts), all(counts>=0), all(counts==floor(counts)),
            !anyDuplicated(rownames(counts)), !anyDuplicated(colnames(counts)),
            !anyDuplicated(metadata$sample), setequal(colnames(counts),metadata$sample))
  metadata <- metadata[match(colnames(counts),metadata$sample),,drop=FALSE]
  levels <- c(control,setdiff(unique(as.character(metadata$group)),control))
  metadata$group <- factor(metadata$group,levels=levels)
  stopifnot(length(levels)>1, all(table(metadata$group)>=2))
  design <- model.matrix(~0+group,metadata)
  colnames(design) <- levels
  y <- edgeR::DGEList(counts)
  keep <- edgeR::filterByExpr(y,design=design)
  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y,design)
  fit <- edgeR::glmQLFit(y,design,robust=TRUE)
  results <- lapply(seq.int(2,length(levels)),function(i) {
    contrast <- rep(0,length(levels));contrast[1]<- -1;contrast[i]<-1
    z <- edgeR::topTags(edgeR::glmQLFTest(fit,contrast=contrast),n=Inf)$table
    z$gene<-rownames(z);z$comparison<-paste(levels[i],'vs',control);z
  })
  list(results=do.call(rbind,results),logCPM=edgeR::cpm(y,log=TRUE,prior.count=2))
}


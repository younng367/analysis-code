# Supply ranked numeric gene scores with unique names in the same identifier system
# as TERM2GENE (columns term and gene). GO/KEGG mappings must be obtained separately.
ranked_enrichment <- function(ranks, term2gene, min_size, max_size, seed) {
  stopifnot(is.numeric(ranks), !anyNA(ranks), all(is.finite(ranks)),
            !is.null(names(ranks)), !anyDuplicated(names(ranks)),
            all(c('term','gene') %in% names(term2gene)))
  set.seed(seed)
  z <- clusterProfiler::GSEA(sort(ranks,decreasing=TRUE),
      TERM2GENE=unique(term2gene[,c('term','gene')]),
      minGSSize=min_size,maxGSSize=max_size,pvalueCutoff=1,
      pAdjustMethod='BH',verbose=FALSE,by='fgsea')
  as.data.frame(z) # Retain all returned terms and the sign of NES.
}

overrepresentation <- function(selected, tested_universe, term2gene, min_size, max_size) {
  stopifnot(all(selected %in% tested_universe))
  as.data.frame(clusterProfiler::enricher(unique(selected),universe=unique(tested_universe),
      TERM2GENE=unique(term2gene[,c('term','gene')]),minGSSize=min_size,
      maxGSSize=max_size,pvalueCutoff=1,qvalueCutoff=1,pAdjustMethod='BH'))
}

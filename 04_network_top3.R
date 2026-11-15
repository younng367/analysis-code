# Supply complete Cytoscape/cytoHubba exports: algorithm, gene, score, export_order.
# Exactly three per algorithm. Ties keep recorded export order, not rank <= 3.
network_top3 <- function(rankings) {
  algorithms <- c('DMNC','MCC','MNC')
  stopifnot(all(c('algorithm','gene','score','export_order') %in% names(rankings)),
            setequal(unique(rankings$algorithm), algorithms),
            !anyNA(rankings), all(is.finite(rankings$score)))
  sets <- setNames(lapply(algorithms, function(a) {
    x <- rankings[rankings$algorithm==a, , drop=FALSE]
    stopifnot(nrow(x)>=3, !anyDuplicated(x$gene), !anyDuplicated(x$export_order))
    x <- x[order(-x$score, x$export_order), , drop=FALSE]
    head(x, 3)
  }), algorithms)
  list(top3=do.call(rbind, sets), intersection=Reduce(intersect,lapply(sets, `[[`, 'gene')))
}
# No upstream differential-expression recalculation is performed by this function.
# It does not establish robustness to alternative network thresholds.


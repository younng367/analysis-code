# Input: normalized log-scale gene-by-sample expression and named gene-set list.
# Supply licensed gene sets separately. No marker collection is distributed here.
ssgsea_scores <- function(expression, signatures) {
  stopifnot(is.matrix(expression), !anyNA(expression), all(is.finite(expression)),
            !anyDuplicated(rownames(expression)), !is.null(names(signatures)))
  signatures <- lapply(signatures, function(x) intersect(unique(x),rownames(expression)))
  stopifnot(all(lengths(signatures)>0))
  parameter <- GSVA::ssgseaParam(expression, signatures, normalize=TRUE)
  GSVA::gsva(parameter,verbose=FALSE)
}


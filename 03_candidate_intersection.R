# All input vectors must use the same species and identifier convention.
# A marker list is distinct from a highly-variable-gene list.
candidate_intersection <- function(bulk_genes, marker_genes, mechanics_sets) {
  stopifnot(is.list(mechanics_sets), length(mechanics_sets)>0)
  vectors <- c(list(bulk_genes, marker_genes), mechanics_sets)
  stopifnot(all(vapply(vectors, function(z) is.character(z) && !anyNA(z), logical(1))))
  mechanics <- unique(unlist(mechanics_sets, use.names=FALSE))
  common <- Reduce(intersect, list(unique(bulk_genes), unique(marker_genes), mechanics))
  list(candidates=sort(common), by_category=lapply(mechanics_sets, function(z)
    Reduce(intersect, list(unique(bulk_genes), unique(marker_genes), unique(z)))))
}


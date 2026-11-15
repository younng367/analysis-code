# Single-section templates: descriptive spatial comparisons, no animal-level P value.
spatial_summary <- function(spots) {
  stopifnot(all(c('spot','section','region','target','value') %in% names(spots)),
            !anyNA(spots),!anyDuplicated(spots[,c('spot','section','target')]))
  aggregate(value~section+region+target,spots,
            function(x)c(n_spots=length(x),mean=mean(x),SD=sd(x)))
}

spatial_root_refits <- function(coordinates, clusters, roots) {
  stopifnot(is.matrix(coordinates),!anyNA(coordinates),
            nrow(coordinates)==length(clusters),!is.null(rownames(coordinates)),
            all(roots %in% clusters),!is.null(names(roots)))
  # Explicitly a refit on supplied coordinates, not a recovered historical object.
  fits <- lapply(roots,function(root) slingshot::slingshot(
      coordinates,clusterLabels=clusters,start.clus=root))
  lapply(fits,function(fit) list(pseudotime=slingshot::slingPseudotime(fit),
                               curves=slingshot::slingCurves(fit)))
}
# Spots from one tissue section are not independent animals.


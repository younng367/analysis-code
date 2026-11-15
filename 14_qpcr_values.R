# Input is one previously normalized expression value per animal and target.
# This template does not reconstruct Cq, delta-Ct or delta-delta-Ct.
qpcr_values <- function(values) {
  stopifnot(all(c('animal','group','target','value') %in% names(values)),
            !anyNA(values[,c('animal','group','target','value')]),
            !anyDuplicated(values[,c('animal','target')]),
            all(is.finite(values$value)),all(values$value>=0))
  values$assay <- 'qPCR'
  values[,c('assay','target','animal','group','value')]
}
# Document normalization and technical-repeat aggregation in the input provenance.
# Source 13_animal_statistics.R, then test at the animal level.


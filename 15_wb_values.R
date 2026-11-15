# Generic normalized ratios from paired target/reference band metrics on the same membrane.
# band_metric definition is an external input; this template does not process images.
wb_values <- function(bands) {
  required <- c('animal','group','target','membrane','target_metric','reference_metric')
  stopifnot(all(required %in% names(bands)),!anyNA(bands[,required]),
            !anyDuplicated(bands[,c('animal','target')]),
            all(is.finite(bands$target_metric)),all(bands$target_metric>=0),
            all(is.finite(bands$reference_metric)),all(bands$reference_metric>0))
  # Each row must pair the target with its OWN same-membrane loading control.
  normalized <- data.frame(assay='WB',target=bands$target,animal=bands$animal,
                          group=bands$group,value=bands$target_metric/bands$reference_metric)
  reference <- data.frame(assay=paste0('reference_',bands$membrane),target=bands$target,
                         animal=bands$animal,group=bands$group,value=bands$reference_metric)
  list(normalized=normalized,reference=reference)
}
# Separate reference comparisons do not prove control equivalence or linear detection.


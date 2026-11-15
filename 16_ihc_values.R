# Input: recorded field IOD and measured area, not images or newly selected fields.
# Keep the area denominator definition and OD calibration consistent across all rows.
ihc_values <- function(fields, fields_per_animal, display_multiplier=1) {
  required <- c('animal','group','target','field','IOD','area')
  stopifnot(all(required %in% names(fields)),!anyNA(fields[,required]),
            all(is.finite(fields$IOD)),all(fields$IOD>=0),
            all(is.finite(fields$area)),all(fields$area>0),
            !anyDuplicated(fields[,c('animal','target','field')]))
  size <- aggregate(field~animal+group+target,fields,length)
  stopifnot(all(size$field==fields_per_animal))
  fields$value <- fields$IOD/fields$area * display_multiplier
  values <- aggregate(value~animal+group+target,fields,mean)
  stopifnot(!anyDuplicated(values[,c('animal','target')]))
  values$assay <- 'IHC'
  list(fields=fields,animals=values[,c('assay','target','animal','group','value')])
}
# Mean of field MODs is not the pooled sum(IOD)/sum(area) when field areas differ.
# No H-score, positivity threshold or cell segmentation is fabricated.


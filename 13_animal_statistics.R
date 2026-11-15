# One row per animal, assay and target; no technical-repeat or field pseudoreplication.
animal_statistics <- function(data, control, treatment) {
  required <- c('assay','target','animal','group','value')
  stopifnot(all(required %in% names(data)),!anyNA(data[,required]),
            all(is.finite(data$value)),
            setequal(unique(data$group),c(control,treatment)),
            !anyDuplicated(data[,c('assay','target','animal')]))
  keys <- unique(data[,c('assay','target')])
  out <- lapply(seq_len(nrow(keys)),function(i) {
    z <- data[data$assay==keys$assay[i] & data$target==keys$target[i],]
    x <- z$value[z$group==control];y<-z$value[z$group==treatment]
    stopifnot(length(x)>=2,length(y)>=2)
    w <- t.test(y,x,var.equal=FALSE,alternative='two.sided',conf.level=.95)
    exact <- !anyDuplicated(c(x,y)) && length(x)<50 && length(y)<50
    m <- wilcox.test(y,x,alternative='two.sided',exact=exact,correct=!exact)
    u <- unname(m$statistic)
    data.frame(assay=keys$assay[i],target=keys$target[i],
      control=control,treatment=treatment,n_control=length(x),n_treatment=length(y),
      mean_control=mean(x),SD_control=sd(x),mean_treatment=mean(y),SD_treatment=sd(y),
      difference=mean(y)-mean(x),CI_low=w$conf.int[1],CI_high=w$conf.int[2],
      Welch_t=unname(w$statistic),Welch_df=unname(w$parameter),Welch_P=w$p.value,
      Mann_Whitney_U=u,Mann_Whitney_P=m$p.value,Mann_Whitney_exact=exact,
      rank_biserial=2*u/(length(x)*length(y))-1)
  })
  out <- do.call(rbind,out)
  out$BH_Welch_P <- p.adjust(out$Welch_P,'BH')
  out
}
# BH applies across every row supplied in one call; define the test family in advance.
# Direction of difference and rank-biserial effect is treatment minus control.


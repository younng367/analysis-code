# Plot animal-level values with mean +/- SD; colours supplied explicitly by the caller.
animal_barplots <- function(data, statistics, group_order, colours) {
  stopifnot(setequal(unique(data$group),group_order),
            all(group_order %in% names(colours)))
  keys <- unique(data[,c('assay','target')])
  plots <- lapply(seq_len(nrow(keys)),function(i) {
    z <- data[data$assay==keys$assay[i]&data$target==keys$target[i],]
    z$group <- factor(z$group,levels=group_order)
    s <- data.frame(group=factor(group_order,levels=group_order),
      mean=vapply(group_order,function(g)mean(z$value[z$group==g]),numeric(1)),
      SD=vapply(group_order,function(g)sd(z$value[z$group==g]),numeric(1)))
    stat <- statistics[statistics$assay==keys$assay[i]&statistics$target==keys$target[i],]
    stopifnot(nrow(stat)==1)
    ggplot2::ggplot(s,ggplot2::aes(group,mean,fill=group)) +
      ggplot2::geom_col(width=.6,colour='black',linewidth=.3) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin=mean-SD,ymax=mean+SD),width=.15) +
      ggplot2::geom_point(data=z,ggplot2::aes(group,value),inherit.aes=FALSE,
        position=ggplot2::position_jitter(width=.07,height=0,seed=1),size=2) +
      ggplot2::scale_fill_manual(values=colours) + ggplot2::theme_classic() +
      ggplot2::theme(legend.position='none') +
      ggplot2::labs(title=paste(keys$assay[i],keys$target[i]),x=NULL,y='Value (input units)',
        subtitle=sprintf("Welch P = %.4g; mean +/- SD",stat$Welch_P))
  })
  names(plots) <- paste(keys$assay,keys$target,sep='_')
  plots
}
# Save a returned plot with ggplot2::ggsave; keep the assay's actual units on the y-axis.

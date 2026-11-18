# Recompute the six original assay-gene comparisons from supplied animal values.
# Usage: Rscript reproduce_six_from_tables.R /path/to/repository_additions /path/to/output_directory
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2) stop("Supply the supplementary folder and an output directory")
input <- normalizePath(args[1]); dest <- args[2]
dir.create(dest, recursive=TRUE, showWarnings=FALSE)
a <- read.csv(file.path(input, "data/animals/animal_values.csv"))
s<-list();k<-0
for(assay in c('qPCR','WB','IHC')) for(gene in c('CORO1A','LCP1')) {
 k<-k+1;z<-a[a$assay==assay&a$gene==gene,];x<-z$value[z$group=='Sham'];y<-z$value[z$group=='TBI'];stopifnot(length(x)==5,length(y)==5,length(unique(c(x,y)))==10)
 w<-t.test(y,x,var.equal=FALSE); m<-wilcox.test(y,x,exact=TRUE,correct=FALSE,alternative='two.sided');u<-unname(m$statistic)
 s[[k]]<-data.frame(assay=assay,gene=gene,n_Sham=5,n_TBI=5,mean_Sham=mean(x),SD_Sham=sd(x),mean_TBI=mean(y),SD_TBI=sd(y),difference_TBI_minus_Sham=mean(y)-mean(x),CI95_low=w$conf.int[1],CI95_high=w$conf.int[2],Welch_t=unname(w$statistic),Welch_df=unname(w$parameter),Welch_P=w$p.value,Mann_Whitney_U_TBI=u,Mann_Whitney_U_min=min(u,25-u),Mann_Whitney_exact_P=m$p.value,rank_biserial_TBI_minus_Sham=2*u/25-1)
}
s<-do.call(rbind,s);s$BH_adjusted_Welch_P_6<-p.adjust(s$Welch_P,'BH')

expected <- read.csv(file.path(input, "results/animals/animal_statistics.csv"))
stopifnot(identical(as.character(s$assay), as.character(expected$assay)),
          identical(as.character(s$gene), as.character(expected$gene)))
for(n in names(s)[sapply(s,is.numeric)]) {
  if(max(abs(s[[n]]-expected[[n]])) > 1e-10) stop(paste("Mismatch:",n))
}
write.csv(s,file.path(dest,"six_endpoint_statistics.csv"),row.names=FALSE)
capture.output(sessionInfo(), file=file.path(dest,"R_sessionInfo.txt"))
writeLines("All numeric columns reproduced within 1e-10.",file.path(dest,"verification.txt"))
print(s[,c("assay","gene","Welch_P","Mann_Whitney_exact_P","BH_adjusted_Welch_P_6")],row.names=FALSE)

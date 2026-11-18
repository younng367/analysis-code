args<-commandArgs(FALSE);root<-normalizePath(file.path(dirname(sub('^--file=','',args[grepl('^--file=',args)])),'..'))
.libPaths(c(file.path(root,'work/Rlib'),.libPaths()))
suppressPackageStartupMessages({library(clusterProfiler);library(org.Mm.eg.db);library(AnnotationDbi)})
out<-file.path(root,'results');dir.create(file.path(out,'enrichment'),showWarnings=FALSE)
d<-read.csv(file.path(out,'bulk_normalized_counts.csv'),row.names=1,check.names=FALSE);d<-log2(as.matrix(d)+1)
mapping<-AnnotationDbi::select(org.Mm.eg.db,keys=rownames(d),keytype='SYMBOL',columns='ENTREZID');mapping<-mapping[!is.na(mapping$ENTREZID),]
write.csv(mapping,file.path(out,'enrichment/symbol_to_entrez_mapping.csv'),row.names=FALSE)
for(g in c('Lcp1','Coro1a')){
  co<-apply(d,1,function(v)if(sd(v)>0)cor(v,d[g,],method='pearson') else NA_real_)
  co<-co[names(co)!=g&!is.na(co)]
  z<-merge(data.frame(SYMBOL=names(co),correlation=co),mapping,by='SYMBOL')
  z<-z[order(-abs(z$correlation),z$SYMBOL),];z<-z[!duplicated(z$ENTREZID),]
  ranks<-setNames(z$correlation,z$ENTREZID);ranks<-sort(ranks,decreasing=TRUE)
  write.csv(z,file.path(out,paste0('enrichment/',g,'_correlation_ranking.csv')),row.names=FALSE)
  set.seed(20261114)
  go<-gseGO(ranks,OrgDb=org.Mm.eg.db,ont='ALL',keyType='ENTREZID',minGSSize=10,maxGSSize=500,pvalueCutoff=1,pAdjustMethod='BH',eps=0,verbose=FALSE)
  write.csv(as.data.frame(go),file.path(out,paste0('enrichment/',g,'_GSEA_GO.csv')),row.names=FALSE)
  set.seed(20261114)
  ke<-tryCatch(gseKEGG(ranks,organism='mmu',minGSSize=10,maxGSSize=500,pvalueCutoff=1,pAdjustMethod='BH',eps=0,verbose=FALSE),error=function(e){writeLines(conditionMessage(e),file.path(out,paste0('enrichment/',g,'_KEGG_ERROR.txt')));NULL})
  if(!is.null(ke))write.csv(as.data.frame(ke),file.path(out,paste0('enrichment/',g,'_GSEA_KEGG.csv')),row.names=FALSE)
}
sets<-read.csv(file.path(out,'mechanics_gene_sets_as_used.csv'))
de<-read.csv(file.path(out,'bulk_DESeq2_corrected.csv'));alltested<-de$gene[!is.na(de$pvalue)]
for(cat in c('MS','MM','CF')){
 genes<-de$gene[de$selected&toupper(de$gene)%in%sets$gene[sets$category==cat]]
 go<-enrichGO(gene=genes,universe=alltested,OrgDb=org.Mm.eg.db,keyType='SYMBOL',ont='BP',pAdjustMethod='BH',pvalueCutoff=1,qvalueCutoff=1,minGSSize=10,maxGSSize=500)
 write.csv(as.data.frame(go),file.path(out,paste0('enrichment/bulk_',cat,'_GO.csv')),row.names=FALSE)
}
writeLines(capture.output(sessionInfo()),file.path(out,'enrichment/sessionInfo.txt'))

args <- commandArgs(trailingOnly=FALSE)
script <- sub('^--file=', '', args[grepl('^--file=',args)])
root <- normalizePath(file.path(dirname(script),'..'))
.libPaths(c(file.path(root,'work/Rlib'),.libPaths()))
suppressPackageStartupMessages({library(DESeq2);library(edgeR);library(GSVA);library(readxl);library(jsonlite)})
bio <- normalizePath(file.path(root,'../../TBI_es_补'))
out <- file.path(root,'results')
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))

# Keep symbols and counts in one table during sorting and duplicate handling.
raw <- read.csv(file.path(bio,'BULK/bulk/GSE79441_gene_exp_filtered.csv'),check.names=FALSE)
samplecols <- grep('_COUNT$',names(raw),value=TRUE)
stopifnot(length(samplecols)==6,all(as.matrix(raw[,samplecols])%%1==0))
raw$original_row <- seq_len(nrow(raw))
raw$mean_count <- rowMeans(raw[,samplecols])
raw <- raw[order(raw$gene,-raw$mean_count,raw$original_row),]
dup <- raw[duplicated(raw$gene)|duplicated(raw$gene,fromLast=TRUE),]
write.csv(dup,file.path(out,'bulk_duplicate_resolution.csv'),row.names=FALSE)
raw <- raw[!duplicated(raw$gene),]
cts <- as.matrix(raw[,samplecols]);rownames(cts)<-raw$gene
meta <- data.frame(row.names=samplecols,group=factor(ifelse(grepl('TBI',samplecols),'TBI','Control'),levels=c('Control','TBI')))
dds <- DESeqDataSetFromMatrix(cts,meta,~group)
dds <- dds[rowSums(counts(dds))>=9,]
dds <- DESeq(dds,quiet=TRUE)
res <- as.data.frame(results(dds,contrast=c('group','TBI','Control'),alpha=.05))
res$gene<-rownames(res);res$selected<-!is.na(res$padj)&res$padj<.05&abs(res$log2FoldChange)>1
write.csv(res,file.path(out,'bulk_DESeq2_corrected.csv'),row.names=FALSE)
norm <- counts(dds,normalized=TRUE)
write.csv(data.frame(gene=rownames(norm),norm,check.names=FALSE),file.path(out,'bulk_normalized_counts.csv'),row.names=FALSE)
saveRDS(dds,file.path(root,'work/bulk_corrected_dds.rds'))

# Original marker list is inspected separately; it is not silently relabelled HVGs.
e<-new.env();load(file.path(bio,'Microglia_DEG.Rdata'),e)
writeLines(capture.output(str(as.list(e),max.level=2)),file.path(out,'microglia_marker_object_structure.txt'))
markers<-e$deg
write.csv(markers,file.path(out,'microglia_original_marker_results.csv'),row.names=FALSE)
microset<-unique(toupper(markers$gene[markers$avg_log2FC>1&markers$p_val_adj<.05]))
sets<-list()
for (pair in list(c('MS','Mechanosensation & Signaling'),c('MM','Mechanical Microenvironment'),c('CF','Cytoskeleton & Force Generation'))){
  z<-unlist(read_excel(file.path(bio,'BULK/bulk/MECHANIC GENE.xlsx'),sheet=pair[2]),use.names=FALSE)
  sets[[pair[1]]]<-unique(toupper(trimws(z[!is.na(z)&trimws(z)!=''])))
}
write.csv(do.call(rbind,lapply(names(sets),function(k)data.frame(category=k,gene=sets[[k]]))),file.path(out,'mechanics_gene_sets_as_used.csv'),row.names=FALSE)
old<-read.csv(file.path(bio,'BULK/common_genes.csv'))[[1]]
allmech<-unique(unlist(sets));sen<-list();candidates<-list()
for(fc in c(0,.5,1,1.5))for(q in c(.01,.05,.1)){
  sel<-!is.na(res$padj)&res$padj<q&abs(res$log2FoldChange)>fc
  bulkset<-toupper(res$gene[sel]);common<-intersect(intersect(bulkset,allmech),microset)
  sen[[length(sen)+1]]<-data.frame(abs_log2FC=fc,FDR=q,bulk_n=sum(sel),MS_n=length(intersect(bulkset,sets$MS)),MM_n=length(intersect(bulkset,sets$MM)),CF_n=length(intersect(bulkset,sets$CF)),micro_marker_n=length(microset),shared_n=length(common),LCP1='LCP1'%in%common,CORO1A='CORO1A'%in%common,shared_genes=paste(sort(common),collapse=';'))
  if(fc==1&&q==.05)candidates<-common
}
write.csv(do.call(rbind,sen),file.path(out,'bulk_screen_threshold_sensitivity.csv'),row.names=FALSE)
audit<-res[toupper(res$gene)%in%union(old,candidates),]
audit$original_16<-toupper(audit$gene)%in%old;audit$corrected_selected_candidate<-toupper(audit$gene)%in%candidates
write.csv(audit,file.path(out,'candidate_selection_audit.csv'),row.names=FALSE)

# Recompute immune scores on the correct six-sample count matrix using current
# GSVA ssGSEA with the supplied marker signatures. This is a revision rerun.
cm<-read.csv(file.path(bio,'BULK/GSVA/cellMarker.csv'),check.names=FALSE)
names(cm)[2]<-'celltype'; sig<-lapply(split(cm$Metagene,cm$celltype),unique)
x<-log2(norm+1);rownames(x)<-toupper(rownames(x))
x<-x[!duplicated(rownames(x)),,drop=FALSE]
scores<-gsva(ssgseaParam(x,sig,normalize=TRUE),verbose=FALSE)
write.csv(data.frame(signature=rownames(scores),scores,check.names=FALSE),file.path(out,'bulk_immune_ssGSEA_scores.csv'),row.names=FALSE)
write.csv(cm,file.path(out,'immune_marker_panel.csv'),row.names=FALSE)

# edgeR quasi-likelihood pseudobulk analysis: library is the experimental unit.
pb<-read.csv(gzfile(file.path(out,'microglia_pseudobulk_counts.csv.gz')),row.names=1,check.names=FALSE)
md<-read.csv(file.path(out,'microglia_library_metadata.csv'))
md<-md[match(colnames(pb),md$sample),];stopifnot(all(md$sample==colnames(pb)))
md$group<-factor(md$group,levels=c('Control','6H_DPI','2D_DPI','4D_DPI'))
design<-model.matrix(~0+group,md);colnames(design)<-levels(md$group)
y<-DGEList(pb);keep<-filterByExpr(y,design=design);y<-y[keep,,keep.lib.sizes=FALSE]
y<-calcNormFactors(y);y<-estimateDisp(y,design);fit<-glmQLFit(y,design,robust=TRUE)
allres<-list()
for(k in 2:4){cvec<-rep(0,4);cvec[1]<- -1;cvec[k]<-1;t<-topTags(glmQLFTest(fit,contrast=cvec),n=Inf)$table;t$gene<-rownames(t);t$comparison<-paste0(colnames(design)[k],'_vs_Control');allres[[k-1]]<-t}
pr<-do.call(rbind,allres);write.csv(pr,file.path(out,'microglia_pseudobulk_all_tests.csv'),row.names=FALSE)
tr<-pr[tolower(pr$gene)%in%c('lcp1','coro1a'),];tr$BH_six_target_tests<-p.adjust(tr$PValue,'BH')
write.csv(tr,file.path(out,'microglia_pseudobulk_targets.csv'),row.names=FALSE)
lcpm<-cpm(y,log=TRUE,prior.count=2)
write.csv(data.frame(gene=rownames(lcpm),lcpm,check.names=FALSE),file.path(out,'microglia_pseudobulk_logCPM.csv'),row.names=FALSE)
print(res[tolower(res$gene)%in%c('lcp1','coro1a'),]);print(do.call(rbind,sen));print(tr)

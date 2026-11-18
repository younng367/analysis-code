args<-commandArgs(FALSE);root<-normalizePath(file.path(dirname(sub('^--file=','',args[grepl('^--file=',args)])),'..'))
.libPaths(c(file.path(root,'work/Rlib'),.libPaths()))
suppressPackageStartupMessages({library(qs);library(Biobase);library(igraph)})
x<-qread(file.path(root,'../../TBI_es_补/monocle/micro_monocle.qs'));a<-attributes(x)
pd<-attributes(a$phenoData)$data;aux<-get('DDRTree',a$auxOrderingData);g<-upgrade_graph(aux$pr_graph_cell_proj_tree)
original<-as.numeric(distances(g,v=aux$root_cell,to=rownames(pd)))
stopifnot(max(abs(pd$Pseudotime-original))<1e-8)
roots<-c(original_stored_root=aux$root_cell)
ctrl<-which(pd$group=='Control'&pd$cell_type=='Homeostatic 1')
roots<-c(roots,Control_Homeostatic1=rownames(pd)[ctrl[which.min(pd$Pseudotime[ctrl])]])
for(state in c('Homeostatic 2','DAM','Pro-inflammatory')){
 idx<-which(pd$cell_type==state);roots<-c(roots,setNames(rownames(pd)[idx[which.max(pd$Pseudotime[idx])]],paste0('distal_',gsub(' ','_',state))))
}
o<-data.frame(nucleus=rownames(pd),sample=pd$orig.ident,timepoint=pd$group,state=pd$cell_type)
summary<-list()
for(n in names(roots)){
 d<-as.numeric(distances(g,v=roots[n],to=o$nucleus));o[[n]]<-d
 summary[[n]]<-data.frame(orientation=n,root_nucleus=roots[n],root_state=pd[roots[n],'cell_type'],root_timepoint=pd[roots[n],'group'],spearman_with_original=cor(original,d,method='spearman'),maximum_distance=max(d),n_nuclei=length(d),n_libraries=length(unique(pd$orig.ident)))
}
write.csv(o,file.path(root,'results/microglia_alternative_orderings.csv'),row.names=FALSE)
write.csv(do.call(rbind,summary),file.path(root,'results/microglia_root_sensitivity.csv'),row.names=FALSE)
write.csv(as_data_frame(g,what='edges'),file.path(root,'results/microglia_projected_tree_edges.csv'),row.names=FALSE)
write.csv(data.frame(nucleus=rownames(pd),t(a$reducedDimS)),file.path(root,'results/microglia_trajectory_coordinates.csv'),row.names=FALSE)
print(do.call(rbind,summary)[,c('orientation','root_state','root_timepoint','spearman_with_original')])

# Focused revision analysis of the archived network. No STRING re-query or DEG re-analysis.
args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[1] else 'revision_20261114/integration_20261114')
src <- normalizePath(if(length(args)>=2)args[2] else file.path(root, '../../../TBI_es_补/BULK'))
out <- if(length(args)>=3)args[3] else file.path(root,'work/network_sensitivity/results');dir.create(out,recursive=TRUE,showWarnings=FALSE)
edges <- read.delim(file.path(src,'String/string612.tsv'),check.names=FALSE)
names(edges)[1] <- 'node1'
candidates <- read.csv(file.path(src,'common_genes.csv'))[[1]]
known <- unique(c(edges$node1,edges$node2)); mapped <- known[match(toupper(candidates),toupper(known))]
nodes <- ifelse(is.na(mapped),candidates,mapped)
stopifnot(length(nodes)==16,length(unique(nodes))==16,all(known %in% nodes))
metrics <- c('DMNC','MCC','MNC')
exports <- setNames(lapply(metrics,function(m){
 d<-read.csv(file.path(src,'String',paste0(tolower(m),'.csv')),skip=1)
 d<-d[!duplicated(d$Name),];d<-d[order(-d$Score,seq_len(nrow(d))),]
 d$export_position<-seq_len(nrow(d));d
}),metrics)
joins<-list();details<-list()
for(k in 2:5){
 sets<-lapply(exports,function(d)head(d$Name,k));common<-Reduce(intersect,sets)
 joins[[length(joins)+1]]<-data.frame(top_n=k,rule='Exactly N entries; score descending; archived order for ties',DMNC=paste(sets$DMNC,collapse=';'),MCC=paste(sets$MCC,collapse=';'),MNC=paste(sets$MNC,collapse=';'),intersection=paste(sort(common),collapse=';'),n_intersection=length(common),Lcp1='Lcp1'%in%common,Coro1a='Coro1a'%in%common)
 for(m in metrics){d<-head(exports[[m]],k);d$algorithm<-m;d$top_n<-k;details[[length(details)+1]]<-d}
}
write.csv(do.call(rbind,joins),file.path(out,'TopN_archived_intersections.csv'),row.names=FALSE)
write.csv(do.call(rbind,details),file.path(out,'TopN_selected_entries.csv'),row.names=FALSE)
components<-function(a){
 left<-seq_len(nrow(a));out<-list()
 while(length(left)){
  reached<-left[1];repeat {new<-unique(c(reached,which(colSums(a[reached,,drop=FALSE])>0)));if(length(new)==length(reached))break;reached<-new}
  out[[length(out)+1]]<-reached;left<-setdiff(left,reached)
 };out
}
score_graph<-function(a){
 # Enumeration is feasible for this small graph; MCC sums (clique size - 1)! over maximal cliques.
 active<-which(rowSums(a)>0);cliques<-list()
 if(length(active)>=2)for(k in 2:min(length(active),max(rowSums(a))+1)){
  cmb<-combn(active,k,simplify=FALSE)
  for(s in cmb)if(all(a[s,s][upper.tri(a[s,s])])){
   outside<-setdiff(seq_len(nrow(a)),s)
   extendable<-length(outside)&&any(colSums(a[s,outside,drop=FALSE])==length(s))
   if(!extendable)cliques[[length(cliques)+1]]<-s
  }
 }
 do.call(rbind,lapply(seq_len(nrow(a)),function(i){
  nb<-which(a[i,]);n<-0;e<-0
  if(length(nb)) {sub<-a[nb,nb,drop=FALSE];cc<-components(sub);largest<-cc[[which.max(lengths(cc))]];n<-length(largest);e<-sum(sub[largest,largest,drop=FALSE])/2}
  data.frame(gene=rownames(a)[i],degree=sum(a[i,]),DMNC=if(n)e/n^1.7 else 0,MCC=sum(vapply(cliques,function(s)if(i%in%s)factorial(length(s)-1)else 0,numeric(1))),MNC=n)
 }))
}
all_scores<-list();summaries<-list();validations<-list()
for(threshold in c(.4,.7)){
 ee<-edges[edges$combined_score>=threshold,];a<-matrix(FALSE,16,16,dimnames=list(nodes,nodes))
 for(i in seq_len(nrow(ee))){a[ee$node1[i],ee$node2[i]]<-TRUE;a[ee$node2[i],ee$node1[i]]<-TRUE}
 d<-score_graph(a);d$confidence_threshold<-threshold
 for(m in metrics){d[[paste0(m,'rank_min')]]<-rank(-d[[m]],ties.method='min');d[[paste0(m,'tie_size')]]<-vapply(d[[m]],function(v)sum(abs(d[[m]]-v)<1e-12),integer(1))}
 all_scores[[length(all_scores)+1]]<-d
 summaries[[length(summaries)+1]]<-data.frame(confidence_threshold=threshold,candidate_nodes=16,edge_supported_nodes=sum(rowSums(a)>0),isolated_candidates=sum(rowSums(a)==0),edges=nrow(ee),observed_additional_nodes=length(setdiff(known,nodes)),DMNC_distinct_scores=length(unique(d$DMNC)),Lcp1_degree=d$degree[d$gene=='Lcp1'],Coro1a_degree=d$degree[d$gene=='Coro1a'],Top3_intersection=if(threshold==.4)'Coro1a;Lcp1 (archived fixed-three rule)'else 'Not uniquely defined: all DMNC scores are zero')
 write.csv(ee,file.path(out,paste0('Edges_score_',threshold,'.csv')),row.names=FALSE)
 if(threshold==.4)for(m in metrics){
  q<-exports[[m]];actual<-d[[m]][match(q$Name,d$gene)];stopifnot(all(abs(actual-q$Score)<1e-12))
  validations[[length(validations)+1]]<-data.frame(algorithm=m,gene=q$Name,archived_score=q$Score,recomputed_score=actual,absolute_difference=abs(actual-q$Score))
 }
}
write.csv(do.call(rbind,all_scores),file.path(out,'Confidence_threshold_all_scores.csv'),row.names=FALSE)
write.csv(do.call(rbind,summaries),file.path(out,'Confidence_threshold_summary.csv'),row.names=FALSE)
write.csv(do.call(rbind,validations),file.path(out,'Baseline_score_validation.csv'),row.names=FALSE)
write.csv(data.frame(input_position=seq_along(nodes),original_symbol=candidates,network_label=nodes,appears_in_edge_export=nodes%in%known),file.path(out,'Candidate_node_inventory.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
print(do.call(rbind,joins)[,c('top_n','intersection','n_intersection')])
print(do.call(rbind,summaries))

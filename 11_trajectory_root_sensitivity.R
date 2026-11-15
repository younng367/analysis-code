# Input: an already fitted, undirected projected tree with named cell vertices.
# This changes distance-to-root orientation; it does not refit Monocle or infer time.
trajectory_root_sensitivity <- function(tree, roots, cell_order, reference) {
  stopifnot(!igraph::is_directed(tree),igraph::is_tree(tree),
            !is.null(names(roots)), reference %in% names(roots),
            all(roots %in% igraph::V(tree)$name),
            all(cell_order %in% igraph::V(tree)$name))
  distance <- sapply(roots,function(root) as.numeric(igraph::distances(
      tree,v=root,to=cell_order,mode='all')))
  stopifnot(all(is.finite(distance)))
  if(is.null(dim(distance))) distance<-matrix(distance,ncol=1,
                                            dimnames=list(NULL,names(roots)))
  list(orderings=data.frame(cell=cell_order,distance,check.names=FALSE),
       comparison=data.frame(root=names(roots),spearman_to_reference=vapply(
         names(roots),function(k) cor(distance[,reference],distance[,k],method='spearman'),numeric(1))))
}


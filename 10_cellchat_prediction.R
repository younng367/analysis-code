# Generic demonstration for one prepared subset. This file is not an executed rerun.
# Input: normalized gene-by-cell expression, metadata with labels, CellChat database.
cellchat_prediction <- function(expression, metadata, database, min_cells, seed) {
  stopifnot(identical(colnames(expression),rownames(metadata)),
            'labels' %in% names(metadata),!anyNA(metadata$labels))
  cc <- CellChat::createCellChat(object=expression,meta=metadata,group.by='labels')
  cc@DB <- database
  cc <- CellChat::subsetData(cc)
  cc <- CellChat::identifyOverExpressedGenes(cc)
  cc <- CellChat::identifyOverExpressedInteractions(cc)
  # raw.use=TRUE means unprojected signaling expression, NOT raw sequencing counts.
  cc <- CellChat::computeCommunProb(cc,type='triMean',raw.use=TRUE,
                                   population.size=TRUE,seed.use=seed)
  cc <- CellChat::filterCommunication(cc,min.cells=min_cells)
  cc <- CellChat::computeCommunProbPathway(cc)
  CellChat::aggregateNet(cc)
}

pathway_information_flow <- function(cellchat) {
  p <- cellchat@netP$prob
  values <- apply(p,3,sum)
  total <- sum(values)
  data.frame(pathway=names(values),information_flow=as.numeric(values),
             relative_flow=if(total>0) as.numeric(values/total) else NA_real_)
}
# Pooled predictions are not evidence of physical signaling or sample-level validation.


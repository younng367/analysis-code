# Template for a prepared, annotated Seurat object; QC/integration are input requirements.
# This marker contrast is NOT a replicate-aware treatment-effect test.
cell_type_markers <- function(object, annotation_column, target_cell_type,
                             min_pct, logfc_cutoff, adjusted_p_cutoff) {
  stopifnot(annotation_column %in% colnames(object[[]]))
  SeuratObject::Idents(object) <- annotation_column
  result <- Seurat::FindMarkers(object, ident.1=target_cell_type,
      ident.2=NULL, test.use='wilcox', only.pos=TRUE,
      min.pct=min_pct, logfc.threshold=0, verbose=FALSE)
  stopifnot(all(c('avg_log2FC','p_val_adj') %in% names(result)))
  result$gene <- rownames(result)
  result$selected <- !is.na(result$p_val_adj) & result$p_val_adj < adjusted_p_cutoff &
                     result$avg_log2FC > logfc_cutoff
  result
}


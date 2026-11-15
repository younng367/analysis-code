# General template: supply an already validated gene-by-sample integer count matrix.
# No study-specific count extraction, identifier repair or spreadsheet parsing.
bulk_differential_expression <- function(counts, metadata, control, treatment,
                                        min_total, fdr_cutoff, abs_log2fc_cutoff) {
  stopifnot(is.matrix(counts), !anyNA(counts), all(counts >= 0),
            all(counts == floor(counts)), !is.null(rownames(counts)),
            !anyDuplicated(rownames(counts)), !anyDuplicated(colnames(counts)),
            !anyDuplicated(metadata$sample),
            setequal(colnames(counts), metadata$sample))
  metadata <- metadata[match(colnames(counts), metadata$sample), , drop=FALSE]
  stopifnot(setequal(unique(metadata$group), c(control, treatment)))
  metadata$group <- factor(metadata$group, levels=c(control, treatment))
  rownames(metadata) <- metadata$sample
  keep <- rowSums(counts) >= min_total
  dds <- DESeq2::DESeqDataSetFromMatrix(counts[keep,,drop=FALSE], metadata, ~group)
  dds <- DESeq2::DESeq(dds, quiet=TRUE)
  result <- as.data.frame(DESeq2::results(dds, contrast=c('group',treatment,control),
                                        alpha=fdr_cutoff))
  result$gene <- rownames(result)
  result$selected <- !is.na(result$padj) & result$padj < fdr_cutoff &
                     abs(result$log2FoldChange) > abs_log2fc_cutoff
  list(results=result, normalized_counts=DESeq2::counts(dds, normalized=TRUE))
}


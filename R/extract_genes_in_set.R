#' Extract genes in set from \code{data.frame} containing all gene p-values.
#'
#' @param genes_in_set \code{vector} containing names of genes in set of
#' interest.
#' @param gene_pvalues \code{data.frame} containing gene names p-values and
#' rank_statistic.
#' @param randomize_gene_labels \code{boolean}. If \code{TRUE}, randomize gene
#'  labels among all genes.
#' @return \code{data.frame} containing subset of \code{gene_pvalues} found in
#' \code{genes_in_set}.
#' @export
extract_genes_in_set <- function(genes_in_set, gene_pvalues,
                                 randomize_gene_labels = FALSE) {
  if (randomize_gene_labels) {
    gene_pvalues$Gene <- gene_pvalues[sample(1:nrow(gene_pvalues),
                                             nrow(gene_pvalues)), "Gene"]
  }

  pvalues_for_genes_in_set <- gene_pvalues[(gene_pvalues$Gene %in% genes_in_set),
    ]

  # Remove duplicate genes
  pvalues_for_genes_in_set <- pvalues_for_genes_in_set[
    order(pvalues_for_genes_in_set$Pvalue),
  ]
  pvalues_for_genes_in_set <- pvalues_for_genes_in_set[
    !duplicated(pvalues_for_genes_in_set$Gene),
  ]

  gene_set_mapped <- pvalues_for_genes_in_set$Gene
  gene_set_mapped <- gene_set_mapped[order(gene_set_mapped)]
  n_genes_mapped <- length(gene_set_mapped)

  gene_set_not_mapped <- genes_in_set[!genes_in_set %in% gene_set_mapped]
  gene_set_not_mapped <- gene_set_not_mapped[order(gene_set_not_mapped)]
  n_genes_not_mapped <- length(gene_set_not_mapped)

  total_genes <- length(genes_in_set)

  return(list(pvalues = pvalues_for_genes_in_set,
              gene_set_mapped = paste(gene_set_mapped,
                                      collapse = ":"),
              gene_set_not_mapped = paste(gene_set_not_mapped, collapse = ":"),
              n_mapped = n_genes_mapped, n_not_mapped = n_genes_not_mapped,
              n_total = total_genes))
}


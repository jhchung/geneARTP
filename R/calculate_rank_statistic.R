#' Calculate rank statistic
#'
#' \deqn{r_i = rank of gene_i / K}
#'
#' Where K = total number of genes
#'
#' @param gene_pvalues \code{data.frame} containing gene name and p-value
#' @param pvalue_col \code{character} or \code{integer} defining the column containing p-values
#' @return \code{data.frame} contining original \code{gene_pvalues} with an additional column of the rank statistic
#' @export
calculate_rank_statistic <- function(gene_pvalues, pvalue_col = "Pvalue"){
  gene_pvalues <- gene_pvalues[order(gene_pvalues[, pvalue_col]),]
  gene_pvalues$rank_statistic <- seq(1:nrow(gene_pvalues))/nrow(gene_pvalues)
  return(gene_pvalues)
}

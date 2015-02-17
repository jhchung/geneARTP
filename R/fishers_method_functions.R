#' Calculate Fishers method test statistic.
#'
#' @param gene_rank_statistic \code{vector} containing rank statistic of genes
#'  of interest.
#' @export
calculate_fishers_method_test_statistic <- function(gene_rank_statistic) {
  # Calculate test statistic
  fishers_statistic <- -2 * (sum(log(gene_rank_statistic)))
  return(fishers_statistic)
}

#' Calculate empirical p-value for Fisher's method.
#'
#' \enumerate{
#'  \item{First calculate rank statistic for all genes being tested.}
#'  \item{Extract genes in set.}
#'  \item{Calculate fishers statistic for the observed data.}
#'  \item{Calculate fishers statistics for permuted data.}
#'  \item{Calculate empirical p-values using observed and permuted data.}
#' }
#'
#' @param genes_in_set \code{vector} containing genes in the set of interest.
#' @param gene_pvalues \code{data.frame} containing gene names and rank
#' statistic calculated by \code{\link{calculate_rank_statistic}}.
#' @param n_perm \code{integer} defining the number of permutations to perform
#' @export
fishers_method_permutation <- function(gene_set, gene_pvalues, n_perm = 1000) {
  require(plyr)
  gene_pvalues <- calculate_rank_statistic(gene_pvalues)
  gene_set_pvalues <- extract_genes_in_set(gene_set, gene_pvalues)
  observed_fm <- calculate_fishers_method_test_statistic(
    gene_set_pvalues$pvalues$rank_statistic
  )

  all_perm_fm <- c()
  for (i in 1:n_perm) {
    perm_gene_set_pvalues <- extract_genes_in_set(gene_set, gene_pvalues,
                                                  randomize_gene_labels = TRUE)
    perm_fm <- calculate_fishers_method_test_statistic(
      perm_gene_set_pvalues$pvalues$rank_statistic
    )
    all_perm_fm <- c(all_perm_fm, perm_fm)
  }

  # Get number of permutations greater than test statistic
  count_perm_greater_than_observed <- count_a_less_than_b(observed_fm,
                                                          all_perm_fm)

  # Calculate empirical p-value
  empirical_fm_pvalue <- (count_perm_greater_than_observed + 1)/(n_perm + 1)

  gene_set_pvalues$empirical_fisher_pvalue <- empirical_fm_pvalue
  return(gene_set_pvalues)
}

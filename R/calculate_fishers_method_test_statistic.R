#' Calculate Fishers method test statistic.
#'
#' @param gene_rank_statistic \code{vector} containing rank statistic of genes of interest.
#' @export
calculate_fishers_method_test_statistic <- function(gene_rank_statistic){
  # Calculate test statistic
  fishers_statistic <- -2*(sum(log(gene_rank_statistic)))
  return(fishers_statistic)
}

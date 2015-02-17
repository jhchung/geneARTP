#' Calculate empirical p-values from permutation data.
#'
#' @param perm_data \code{list} containing permutation data where length is the number of simulations.
#' @return \code{numeric} of the empirical p-value.
#' @export
calculate_simulated_pvalues <- function(perm_data){
  count_sim_less_than_all <- c()
  for (sim_number in 1:length(perm_data)){
    count_sim_less_than_all[[sim_number]] <- count_a_less_than_b(perm_data, perm_data[sim_number])
  }
  sum_of_counts <- llply(count_sim_less_than_all, function(x) sum(as.integer(unlist(x))))
  sum_of_counts <- unlist(sum_of_counts)
  estimated_pvalues <- sum_of_counts / (length(sum_of_counts))
  return(estimated_pvalues)
}

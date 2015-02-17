#' Calculate rank truncation product
#'
#' @param gene_rank_statistics \code{vector} containing gene rank statistics calculated by \code{\link{calculate_rank_statistic}}
#' @param num_genes \code{integer} defining the truncation point to test
#' @return \code{numeric} giving the RTP for the given truncation point
#' @export
calculate_rtp <- function(num_genes, gene_rank_statistics){
  gene_rank_statistics <- gene_rank_statistics[order(gene_rank_statistics)]
  rtp <- sum(log(gene_rank_statistics[1:num_genes]))
  return(rtp)
}

#' Calculate rank truncated product for all truncation points.
#'
#' @param rank_statistics \code{vector} containing rank_statistics_to_test.
#' @param gene_set \code{vector} contianing genes in set to test.
#' @param gene_pvalues \code{data.frame} containing gene names and p-values.
#' @param pvalue_col \code{character} giving name of column containing gene p-values.
#' @param randomize_gene_labels \code{boolean}. If \code{TRUE}, randomize gene labels in \code{gene_pvalues} to perform permutations.
#' @return \code{vector} of rtp for all truncation points in gene set.
#' @export
calculate_all_truncation_rtp <- function(rank_statistics){
  require(plyr)

  # calculate RTP for all truncation points
  all_rtp <- llply(.data = seq(1:length(rank_statistics)),
                   .fun = calculate_rtp,
                   gene_rank_statistics = rank_statistics)
  all_rtp <- unlist(all_rtp)
  # Calcualate RTP for each truncation point
  #   for (i in 1:nrow(gene_set_pvalues$pvalues)){
  #     test_rtp <- calculate_rtp(
  #       gene_rank_statistics = gene_set_pvalues$pvalues$rank_statistic,
  #       num_genes = i)
  #     all_rtp[[i]] <- test_rtp
  #     names(all_rtp[[i]]) <- paste("truncation_point_", i, sep = "")
  #   }
  #   all_rtp <- unlist(all_rtp)
  #
  #   gene_set_pvalues$all_rtp <- all_rtp
  return(all_rtp)
}

#' Calculate empirical ARTP p-values using permutations.
#'
#' @param gene_set \code{vector} contianing genes in set to test.
#' @param gene_pvalues \code{data.frame} containing gene names and p-values.
#' @param n_perm \code{integer} number of permutations to perform.
#' @return \code{data.frame} containing empirical ARTP p-values.
#' @export
calculate_artp_perm <- function(gene_set, gene_pvalues, n_perm = 1000){
  require(plyr)

  gene_pvalues <- calculate_rank_statistic(gene_pvalues)


  # Calcualte RTP statistic for observed dataset
  observed_gene_set <- extract_genes_in_set(genes_in_set = gene_set,
                                            gene_pvalues = gene_pvalues,
                                            randomize_gene_labels = FALSE)
  observed_rank_statistic <- observed_gene_set$pvalues$rank_statistic
  observed_rtp_results <- calculate_all_truncation_rtp(observed_rank_statistic)

  simulated_rtp_results <- list()
  # Calculate RTP statistic for each simulated dataset
  for(i in 1:n_perm){
    simulated_gene_set <- extract_genes_in_set(genes_in_set = gene_set,
                                               gene_pvalues = gene_pvalues,
                                               randomize_gene_labels = TRUE)
    simulated_rank_statistic <- simulated_gene_set$pvalues$rank_statistic
    simulated_rtp_results[[i]] <- calculate_all_truncation_rtp(simulated_rank_statistic)
  }

  # Combine observed and simulated RTP
  all_rtp <- c(list(observed_rtp_results), simulated_rtp_results)
  all_rtp <- as.data.frame(do.call(rbind, all_rtp))
  rownames(all_rtp) <- paste("b", seq(0, (nrow(all_rtp) - 1)), sep = "")
  names(all_rtp) <- paste("J", seq(1:ncol(all_rtp)), sep = "")

  # estimate p-values for each simulated dataset
  simulated_pvalues <- apply(all_rtp, 2, calculate_simulated_pvalues)
  rownames(simulated_pvalues) <- rownames(all_rtp)

  # Select minimum truncation point p-value for each simulation
  min_truncation_pvalues <- unlist(alply(simulated_pvalues, 1, min))
  count_pvalues_less_than_observed <- sum(as.integer(min_truncation_pvalues <=
                                                       min_truncation_pvalues[1]))
  artp_pvalue <- count_pvalues_less_than_observed / (length(min_truncation_pvalues))

  return(artp_pvalue)
}

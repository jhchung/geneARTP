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

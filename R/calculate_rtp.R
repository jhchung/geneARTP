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

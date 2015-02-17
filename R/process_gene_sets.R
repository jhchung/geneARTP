#' Filter gene set on size.
#'
#' @param all_gene_sets List containing the gene sets.
#' @param ... Parameters to pass to \code{check_gene_set_size}
#' @export
filter_gene_set_on_size <- function(all_gene_sets, ...){
  require(plyr)
  filtered_sets <- llply(all_gene_sets, check_gene_set_size, ...)
  filtered_sets <- filtered_sets[!is.na(filtered_sets)]
  return(filtered_sets)
}

#' Check gene set size.
#'
#' Check gene set size to see if it is within the size boundaries.
#'
#' @param gene_set Gene set to test.
#' @param nMin Integer. Minimum number of genes.
#' @param nMax Integer. Maximum number of genes.
#' @export
check_gene_set_size <- function(gene_set, nMin, nMax){
  if (length(gene_set) >= nMin & length(gene_set) <= nMax){
    return(gene_set)
  } else {
    return(NA)
  }
}

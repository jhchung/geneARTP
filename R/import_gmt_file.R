#' Import GMT formatted file.
#'
#' Import GMT (Gene Matrix Transposed) file. A tab delimited file where each row
#' is a gene set. The first column is the gene set name, second column contains
#' a brief description, and the remaining columns contain genes in the gene set.
#'
#' @param gmt_file Path to input \code{.gmt} file.
#' @return list containing gene sets.
#' @export
read_gmt_file <- function(gmt_file) {
  message("Reading .gmt file")
  gmt_data <- scan(gmt_file, what = "character", sep = "\n", fill = TRUE,
                   quote = "")
  gmt_data <- strsplit(gmt_data, "\t")
  gmt_data <- format_gmt(gmt_data)
  return(gmt_data)
}

#' Format \code{.gmt} data.
#'
#' @param gmt_data \code{.gmt} data as read using \code{scan}.
format_gmt <- function(gmt_data) {
  require(plyr)
  formatted_gmt_data <- llply(gmt_data, function(x) as.array(x[3:length(x)]))
  set_names <- laply(gmt_data, function(x) x[[1]])

  names(formatted_gmt_data) <- set_names
  return(formatted_gmt_data)
}

#' Count the number of elements in \code{a} are less than \code{b}
#'
#' @param a thing that less
#' @param b thing that is more
#' @return \code{integer} count of times \code{a} is less than \code{b}
#' @export
count_a_less_than_b <- function(a, b){
  count_of_a_less_than_b <- a <= b
  count_of_a_less_than_b <- sum(as.integer(count_of_a_less_than_b))
  return(count_of_a_less_than_b)
}

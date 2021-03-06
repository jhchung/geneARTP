% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/import_gmt_file.R
\name{read_gmt_file}
\alias{read_gmt_file}
\title{Import GMT formatted file.}
\usage{
read_gmt_file(gmt_file)
}
\arguments{
\item{gmt_file}{Path to input \code{.gmt} file.}
}
\value{
list containing gene sets.
}
\description{
Import GMT (Gene Matrix Transposed) file. A tab delimited file where each row
is a gene set. The first column is the gene set name, second column contains
a brief description, and the remaining columns contain genes in the gene set.
}


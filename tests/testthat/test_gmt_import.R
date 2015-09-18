library(geneARTP)

test_that("Test if GMT file is imported correctly",{
  test_gmt_line <- "gene_set1\tGene set description\tGene1\tGene2\tGene3"

  expect_match(names(format_gmt(test_gmt_line)), "gene_set1")
  expect_equal(length(format_gmt(test_gmt_line)), 1)
  expect_equal(length(format_gmt(test_gmt_line)[[1]]), 3)
})

pvals <- c(0.05, 0.05, 0.05)
test_genes <- data.frame(Gene = c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5"),
                         Pvalue = c(0.001, 0.01, 0.05, 0.00001, NA))
gene_set <- c("Gene1", "Gene3", "Gene4", "Gene5")

test_that("calculate_rank_statistic gives proper rank", {
  test_genes_rank <- calculate_rank_statistic(test_genes)
  expect_equal(nrow(test_genes_rank), 4)
  expect_equivalent(test_genes_rank$rank_statistic,
                    c(0.25, 0.50, 0.75, 1.00))
})

test_that("calculate_fishers_method_test_statistic gives correct statistic",{
  expect_equal(round(calculate_fishers_method_test_statistic(pvals), 5),
               17.97439)
})

test_that("fishers_method_permutation works", {
  set.seed(1)
  fisher_perm_results <- fishers_method_permutation(gene_set, test_genes,
                                                    n_perm = 100)
  expect_equal(length(fisher_perm_results), 7)
  expect_equal(fisher_perm_results$n_mapped, 3)
  expect_equal(fisher_perm_results$n_not_mapped, 1)
  expect_equal(round(fisher_perm_results$empirical_fisher_pvalue, 7),
               0.5445545)
})

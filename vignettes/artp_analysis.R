library("optparse")
library("plyr")
library("dplyr")
library("geneARTP")
options(stringsAsFactors = FALSE)

# Set up command line options ----
option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              default = "C:/Lab/projects/WES/results/wes_wgs/metaskat/wes_wgs.deleterious.maf_0.01.skat_burden.txt",
              help = "Path to input gene p-value file"),
  make_option(c("-g", "--gmt"), type = "character",
              default = "C:/Lab/projects/WES/data/gene_lists/tbx1_and_noncanonical_wnt.gmt",
              help = "Path to GMT file containing gene set definitions"),
  make_option(c("-o", "--output_dir"), type = "character",
              default = "C:/Lab/projects/WES/results/wes_wgs/artp",
              help = "Output file prefix"),
  make_option(c("-p", "--permutations"), type = "integer", default = 1000,
              help = "Number of permutations to run. [default]"),
  make_option("--set_min", type = "integer", default = 20,
              help = "Minimum number of genes in set. [default]"),
  make_option("--set_max", type = "integer", default = 250,
              help = "Maximum number of genes in set [default]"),
  make_option("--adjust_method", type = "character", default = "fdr",
              help = "Method to adjust p-values c('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY','fdr', 'none')"),
  make_option("--pvalue_col", type = "character",
              default = "MetaSKAT",
              help = "Name of the p-value column"),
  make_option("--gene_col", type = "character", default = "Gene",
              help = "Name of the gene column"),
  make_option("--gene_pval_cutoff",
              type = "numeric",
              default = 0.5,
              help = "Maximum p-value for genes")
)

#-------------------------------------------------------------------------------
# Begin script
#-------------------------------------------------------------------------------
args <- parse_args(object = OptionParser(option_list = option_list),
                   args = commandArgs(trailingOnly = TRUE))

# Import data ----
pvalues <- read.table(args$input, header = TRUE) %>%
  dplyr::select_(args$gene_col, args$pvalue_col) %>%
  dplyr::rename_(Pvalue = args$pvalue_col)

gene_sets <- read_gmt_file(args$gmt)

# Filter genes ----
pvalues <- pvalues %>%
  dplyr::filter(Pvalue < args$gene_pval_cutoff)

# Filter gene sets ----
gene_sets <- filter_gene_set_on_size(gene_sets,
                                     nMin = args$set_min,
                                     nMax = args$set_max)

message("Testing ", length(gene_sets), " gene sets with: ", args$set_min,
        " and ", args$set_max, " genes")

# Calculate Fisher's p-values ----
message("Calculate fisher p-values")
fisher_pvalues <- plyr::llply(.data        = gene_sets,
                              .fun         = fishers_method_permutation,
                              gene_pvalues = pvalues,
                              n_perm       = args$permutations,
                              .progress    = "text")

# Convert list to data.frame
fisher_pvalues <- ldply(fisher_pvalues, function(x){
  as.data.frame(x[c("empirical_fisher_pvalue",
                    "gene_set_mapped",
                    "gene_set_not_mapped",
                    "n_mapped",
                    "n_not_mapped",
                    "n_total")])
})

# Adjust p-values for multiple testing
fisher_pvalues <- fisher_pvalues %>%
  dplyr::mutate(fisher_pvalue_adjust = p.adjust(
    fisher_pvalues$empirical_fisher_pvalue, method = args$adjust_method)
  ) %>%
  dplyr::rename(set           = .id,
                fisher_pvalue = empirical_fisher_pvalue)

# Calculate ARTP p-values ----
message("Calculate ARTP p-values")
artp_pvalues <- llply(.data        = gene_sets,
                      .fun         = calculate_artp_perm,
                      gene_pvalues = pvalues,
                      n_perm       = args$permutations,
                      .progress    = "text")

# Convert to data.frame
artp_pvalues <- ldply(artp_pvalues) %>%
  dplyr::rename(set         = .id,
                artp_pvalue = V1) %>%
  dplyr::mutate(artp_pvalue_adjust = p.adjust(artp_pvalue,
                                              method = args$adjust_method))

# Combine Fisher's method and ARTP ----
merged_pvalues <- fisher_pvalues %>%
  dplyr::inner_join(artp_pvalues, by = "set") %>%
  dplyr::arrange(artp_pvalue) %>%
  dplyr::select(one_of(c("set", "fisher_pvalue", "artp_pvalue",
                         "fisher_pvalue_adjust", "artp_pvalue_adjust",
                         "n_mapped", "n_not_mapped", "n_total",
                         "gene_set_mapped", "gene_set_not_mapped")))

# Export results ----
merged_output_file <-gsub(".txt$", "", basename(args$input)) %>%
  paste("merged_pvalues.txt", sep = ".")
merged_output_path <- file.path(args$output_dir, merged_output_file)

dir.create(args$output_dir, recursive=TRUE)
write.table(merged_pvalues,
            file = merged_output_path,
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

message("Finished")

message("Arguments")
print(args)
cat("\n\n")
print(sessionInfo())

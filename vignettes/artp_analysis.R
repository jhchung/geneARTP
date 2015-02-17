library("GWASAnalysis")
library("optparse")
library("plyr")
library("GenomicFeatures")
library(org.Hs.eg.db)
options(stringsAsFactors = FALSE)

#' Load or create gene ranges for each isoform using GenomicFeatures and UCSC
#' browser track.
#'
#' Genes with multiple isoforms are collapsed into the largest continuous
#' interval
#'
#' @param genome_build \code{character} describing the genome build to use.
#' @param coordinate_file \code{character} file path to premade coordinate file
#'  or location to save the coordinate file.
#' @return \code{data.frame} containing gene coordinates. Contains the following
#'  columns:
#'    entrez: entrez ID.
#'    chr: Chromsome.
#'    start: Start position of the largest interval in bp.
#'    end: End position of the largest interval in bp.
#'    width: Width of the largest interval in bp.
#'    strand: Strand the gene is located on.
#'  @export
load_genomic_ranges <- function(genome_build = "hg18", coordinate_file = NA){
  require(GenomicFeatures)

  if (file.exists(coordinate_file)){
    gene_coordinates <- read.table(coordinate_file, header = TRUE)
  } else {
    message("Download UCSC refGene track")
    gene_transcripts <- makeTranscriptDbFromUCSC(genome = genome_build,
                                                 tablename = "refGene")

    message("Extracting gene transcripts")
    gene_transcripts <- transcriptsBy(gene_transcripts, by = "gene")

    message("Collapsing gene isoforms")
    gene_bounds <- seqapply(gene_transcripts, range)

    message("Extracting chromosome positions")
    gene_coordinates <- as.data.frame(gene_bounds)
    gene_coordinates$seqnames <- gsub("chr", "", as.character(gene_coordinates$seqnames))
    names(gene_coordinates) <- c("entrez", "chr", "start", "end", "width", "strand")

    message("Saving gene coordinates file to:\n\t", cooridnate_file)
    write.table(gene_coordinates, file = coordinate_file, col.names = TRUE,
                row.names = FALSE, quote = FALSE, sep = "\t")
  }
  return(gene_coordinates)
}

#' Attempt to resolve gene ID's by mapping gene alias' to entrez ID and then to
#' current gene symbols.
#'
#' Order of matching gene names:
#'  \enumerate{
#'    \item If there is a single entrez entry in \code{gene_data}, return that
#'      entry. Otherwise...
#'    \item Get gene coordinates for the entrez entries.
#'    \item Filter gene coordinates based on chromosome.
#'    \item If the entrez entries are not found in \code{gene_coordinates}, try
#'      and match gene names using \code{org.Hs.eg.db}. Else...
#'    \item Check if gene start is within \code{flank} distance from the entrez
#'      start.
#'    \item Check for an exact match with \code{gene_data} gene name and the
#'      current symbol.
#'    \item If unable to resolve gene symbol, set as \code{NA}.
#'    \item If no match in \code{entrez_coordinates} set as \code{NA}.
#'  }
#'
#' @param gene_data row in a \code{data.frame} that contains the gene of
#'  interest.
#' @param flank distance to match start positions.
#' @export
resolve_gene_id <- function(gene_data, flank, gene_coordinates){
  entrez_id <- gene_data$entrez[[1]]

  matched_id <- NULL

  # If there is already only a single match, set matched_id to the current entrez
  if (length(gene_data$entrez[[1]]) == 1){
    matched_id <- gene_data$entrez
  }

  if ((length(entrez_id) > 1) & is.null(matched_id)){
    # Get chromosome locations for entrez genes
    entrez_coordinates <- gene_coordinates[as.character(gene_coordinates$entrez) %in% entrez_id, ]

    # Check chromosome
    entrez_coordinates <- subset(entrez_coordinates,
                                 subset = (chr == gene_data$Chr))

    if (nrow(entrez_coordinates) == 0){
      # Double check if gene name matches current Symbol
      gene_name <- select(org.Hs.eg.db,
                          keys = as.character(entrez_id),
                          columns = "SYMBOL",
                          keytype = "ENTREZID")[["SYMBOL"]]
      gene_name <- grep(gene_data$Gene, gene_name, fixed = TRUE, value = TRUE)

      if (length(gene_name) == 1){
        matched_id <- get(gene_name, org.Hs.egSYMBOL2EG)
      }
    }

    # Filter further based on start positions
    if (nrow(entrez_coordinates) > 1){
      matched_entrez <- subset(entrez_coordinates,
                               subset = (abs(start - gene_data$Start) <= flank))
      if(nrow(matched_entrez) == 1){
        entrez_coordinates <- matched_entrez
      }
    }

    # Last resort to match based on gene grep of gene name. This may not produce
    # the correct mapping because the same symbol may have been used for a
    # different gene in the past.
    if(nrow(entrez_coordinates) > 1){
      gene_name <- select(org.Hs.eg.db,
                          keys = as.character(entrez_coordinates$entrez),
                          columns = "SYMBOL",
                          keytype = "ENTREZID")[["SYMBOL"]]
      gene_name_match <- grep(paste("^", gene_data$Gene, "$", sep = ""),
                              gene_name,
                              perl = TRUE)
      if(length(gene_name_match) == 1){
        entrez_coordinates <- entrez_coordinates[gene_name_match,]
      }
    }

    # Return the first row of entrez_coordinates if nothing else worked
    if (nrow(entrez_coordinates) > 1){
      warning(paste("Could not resolve entrez IDs for:", gene_data$Gene,
                    "\n\t\tTrying to resolve IDs:",
                    paste(entrez_coordinates$entrez, collapse = ", "),
                    "\n\tSetting to NA\n"))
      entrez_coordinates <- entrez_coordinates[1,]
      entrez_coordinates$entrez <- NA
    }

    # If in the end, nothing works, return NA
    if(nrow(entrez_coordinates) == 0){
      warning(paste("Could not resolve entrez IDs for:", gene_data$Gene,
                    "\n\t\tTrying to resolve IDs:",
                    paste(entrez_coordinates$entrez, collapse = ", "),
                    "\n\tSetting to NA\n"))
      matched_id <- NA
    } else {
      matched_id <- entrez_coordinates$entrez
    }
  }
  gene_data$entrez <- unlist(matched_id)
  return(gene_data)
}

# Set up command line options
option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              default = "C:/Lab/projects/gwas/results/TOF6/vegas/filtered.BMGWAS4BCXBD.TOF6.nominal_eigencorr.assoc.logistic.top10.out",
              help = "Path to input gene p-value file"),
  make_option(c("-g", "--gmt"), type = "character",
              default = "C:/Lab/projects/gwas/data/gmt_files/c2.cp.v4.0.symbols.gmt",
              help = "Path to GMT file containing gene set definitions"),
  make_option(c("-o", "--output_prefix"), type = "character",
              default = "C:/Lab/projects/gwas/results/TOF6_maf05/artp/filtered.BMGWAS4BCXBD.TOF6.nominal_eigencorr.top10",
              help = "Output file prefix"),
  make_option(c("-p", "--permutations"), type = "integer", default = 9,
              help = "Number of permutations to run. [default]"),
  make_option("--min", type = "integer", default = 20,
              help = "Minimum number of genes in set. [default]"),
  make_option("--max", type = "integer", default = 250,
              help = "Maximum number of genes in set [default]"),
  make_option("--adjust_method", type = "character", default = "bonferroni",
              help = "Method to adjust p-values c('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY','fdr', 'none')"),
  make_option(c("-u", "--update_gene_symbol"), type = "logical",
              default = TRUE,
              help = "Update gene symbols by converting alias' to entrez ID and then to Symbols"),
  make_option("--genome_build", type = "character",
              default = "hg18",
              help = "Genome build used to match gene coordinates"),
  make_option(c("-c", "--convert_to_entrez"), type = "logical",
              default = FALSE,
              help = "Convert input gene names to entrez id"),
  make_option("--pvalue_col", type = "character",
              default = "Top-0.1-pvalue",
              help = "Name of the p-value column"),
  make_option("--gene_col", type = "character",
              default = "Gene",
              help = "Name of the gene column"),
  make_option("--gene_flank", type = "integer",
              default = 100,
              help = "Flanking distance when trying to map genes to entrez ID"),
  make_option("--coordinate_file", type = "character",
              default = "C:/Lab/projects/gwas/data/hg18_gene_coordinates.txt",
              help = "data.frame containing the gene coordinates used for resolving entrez ID")
)

args <- parse_args(object = OptionParser(option_list = option_list),
                   args = commandArgs(trailingOnly = TRUE))

pvalues <- read.table(args$input, header = TRUE)

if (args$update_gene_symbol){
  # Get list of all gene alias'
  alias2eg <- org.Hs.egALIAS2EG
  mapped_genes <- mappedkeys(alias2eg)
  alias2eg <- as.list(alias2eg[mapped_genes])

  # Check if the genes in pvalues data.frame are among known alias'
  message("Removing genes with unknown alias")
  pvalues <- pvalues[(pvalues[, args$gene_col] %in% names(alias2eg)), ]

  # Add entrez ID
  message("Converting to ENTREZ ID")
  pvalues$entrez <- mget(pvalues[, args$gene_col], org.Hs.egALIAS2EG)

  # Resolve genes that have multiple entrez ID's
  message("\tResolving genes with multiple ENTREZ ID")
  gene_coordinates <- load_genomic_ranges(genome_build = args$genome_build,
                                          coordinate_file = args$coordinate_file)

  pb <- txtProgressBar(min = 0, max = nrow(pvalues), style = 3)
  for (i in 1:nrow(pvalues)){
    setTxtProgressBar(pb, i)
    pvalues[i, ] <- resolve_gene_id(pvalues[i, ],
                                    args$gene_flank,
                                    gene_coordinates)
  }
  close(pb)

  pvalues$entrez  <- unlist(pvalues$entrez)

  # Remove NAs
  pvalues <- pvalues[complete.cases(pvalues$entrez), ]

  # Add current gene symbol
  pvalues$symbol <- select(org.Hs.eg.db,
                           keys = as.character(pvalues$entrez),
                           columns = "SYMBOL",
                           keytype = "ENTREZID")[["SYMBOL"]]

  if (args$convert_to_entrez){
    pvalues <- pvalues[, c("entrez", "Pvalue")]
  } else {
    pvalues <- pvalues[, c("symbol", "Pvalue")]
  }

  names(pvalues) <- c("Gene", "Pvalue")
} else {
  pvalues <- pvalues[, c(args$gene_col, args$pvalue_col)]
  names(pvalues) <- c("Gene", "Pvalue")
}

gene_sets <- read_gmt_file(args$gmt)
gene_sets <- filter_gene_set_on_size(gene_sets, nMin = args$min, nMax = args$max)
message("Testing ", length(gene_sets), " with: ", args$min, " and ", args$max, " genes")

message("Calculate fisher p-values")
fisher_pvalues <- llply(gene_sets,
                        fishers_method_permutation,
                        gene_pvalues = pvalues,
                        n_perm = args$permutations,
                        .progress = "text")

fisher_pvalues <- llply(fisher_pvalues, function(x) x[c("empirical_fisher_pvalue",
                                                        "gene_set_mapped",
                                                        "gene_set_not_mapped",
                                                        "n_mapped",
                                                        "n_not_mapped",
                                                        "n_total")])

fisher_pvalues <- ldply(fisher_pvalues, function(x) as.data.frame(x))

fisher_pvalues <- fisher_pvalues[order(fisher_pvalues$empirical_fisher_pvalue), ]
fisher_pvalues$fisher_pvalue_fdr <- p.adjust(fisher_pvalues$empirical_fisher_pvalue,
                                             method = args$adjust_method)
names(fisher_pvalues)[1:2] <- c("set", "fisher_pvalue")
# fisher_output_file <- paste(args$output_prefix, "fisher_pvalues.txt", sep = ".")

#######
# Calculate ARTP p-values
message("Calculate ARTP p-values")
artp_pvalues <- llply(.data = gene_sets,
                      .fun = calculate_artp_perm,
                      .progress = "text",
                      gene_pvalues = pvalues,
                      n_perm = args$permutations)

artp_pvalues <- ldply(artp_pvalues, function(x) c(names(x), x))
names(artp_pvalues) <- c("set", "artp_pvalue")
artp_pvalues <- artp_pvalues[order(artp_pvalues$artp_pvalue), ]
artp_pvalues$artp_pvalue_fdr <- p.adjust(artp_pvalues$artp_pvalue,
                                         method = args$adjust_method)
artp_output_file <- paste(args$output_prefix, "artp_pvalues.txt", sep = ".")

merged_pvalues <- merge(fisher_pvalues, artp_pvalues, by.x = "set", by.y = "set")
merged_pvalues <- merged_pvalues[order(merged_pvalues$artp_pvalue), ]
merged_pvalues <- merged_pvalues[, c("set",
                                     "fisher_pvalue",
                                     "artp_pvalue",
                                     "fisher_pvalue_fdr",
                                     "artp_pvalue_fdr",
                                     "n_mapped",
                                     "n_not_mapped",
                                     "n_total",
                                     "gene_set_mapped",
                                     "gene_set_not_mapped")]

merged_output_file <- paste(args$output_prefix, "merged_pvalues.txt", sep = ".")
dir.create(dirname(merged_output_file), recursive=TRUE)
write.table(merged_pvalues,
            file = merged_output_file,
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = "\t")

message("Finished")

message("Arguments")
print(args)
cat("\n\n")
print(sessionInfo())

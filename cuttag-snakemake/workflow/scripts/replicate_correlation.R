#!/usr/bin/env Rscript
#' CUT&Tag Replicate Correlation Analysis
#' 
#' Based on: https://yezhengstat.github.io/CUTTag_tutorial/
#' Generates correlation heatmap for replicate reproducibility assessment.

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(corrplot)
})

# Get parameters from command line or Snakemake
args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 3) {
  # Command line usage: Rscript replicate_correlation.R bin_file1,bin_file2,... sample1,sample2,... output.pdf
  bin_files <- strsplit(args[1], ",")[[1]]
  samples <- strsplit(args[2], ",")[[1]]
  output_pdf <- args[3]
} else if (exists("snakemake")) {
  # Snakemake usage
  bin_files <- snakemake@input[["bin_counts"]]
  samples <- snakemake@params[["samples"]]
  output_pdf <- snakemake@output[["pdf"]]
} else {
  stop("Please provide input files, sample names, and output file path.")
}

# Read and merge fragment count data
fragCount <- NULL

for (i in seq_along(bin_files)) {
  sample <- samples[i]
  
  # Read bin count file
  tmp <- read.table(bin_files[i], header = FALSE, 
                    col.names = c("chrom", "bin", sample))
  
  if (is.null(fragCount)) {
    fragCount <- tmp
  } else {
    fragCount <- full_join(fragCount, tmp, by = c("chrom", "bin"))
  }
}

# Calculate correlation matrix (log2 transformed)
countMatrix <- fragCount %>%
  select(-c("chrom", "bin")) %>%
  as.matrix()

# Replace NA with 0 and add pseudocount for log transformation
countMatrix[is.na(countMatrix)] <- 0
countMatrix <- countMatrix + 1  # Add pseudocount

# Log2 transform
logMatrix <- log2(countMatrix)

# Calculate Pearson correlation
M <- cor(logMatrix, use = "complete.obs")

# Generate correlation plot
pdf(output_pdf, width = 10, height = 10)

corrplot(M, 
         method = "color", 
         outline = TRUE, 
         addgrid.col = "darkgray", 
         order = "hclust", 
         addrect = 3, 
         rect.col = "black", 
         rect.lwd = 3,
         cl.pos = "b", 
         tl.col = "indianred4", 
         tl.cex = 1, 
         cl.cex = 1, 
         addCoef.col = "black", 
         number.digits = 2, 
         number.cex = 0.8, 
         col = colorRampPalette(c("midnightblue", "white", "darkred"))(100),
         title = "Replicate Correlation (log2 Fragment Counts)",
         mar = c(0, 0, 2, 0))

dev.off()

message("Correlation plot saved to: ", output_pdf)

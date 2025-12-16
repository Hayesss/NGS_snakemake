#!/usr/bin/env Rscript
#' CUT&Tag Differential Peak Analysis
#' 
#' Based on: https://yezhengstat.github.io/CUTTag_tutorial/
#' Performs differential analysis using DESeq2 between conditions.

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(GenomicRanges)
  library(chromVAR)
  library(DESeq2)
})

# Get parameters from command line
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  cat("Usage: Rscript differential_analysis.R peak_dir bam_dir conditions samples output_prefix\n")
  cat("  peak_dir: Directory containing SEACR peak files\n")
  cat("  bam_dir: Directory containing BAM files\n")
  cat("  conditions: Comma-separated list of conditions (e.g., K27me3,K27me3,K4me3,K4me3)\n")
  cat("  samples: Comma-separated list of sample names\n")
  cat("  output_prefix: Prefix for output files\n")
  quit(status = 1)
}

peak_dir <- args[1]
bam_dir <- args[2]
conditions <- strsplit(args[3], ",")[[1]]
samples <- strsplit(args[4], ",")[[1]]
output_prefix <- args[5]

# Filter out IgG samples
non_igg <- conditions != "IgG"
conditions <- conditions[non_igg]
samples <- samples[non_igg]

if (length(unique(conditions)) < 2) {
  cat("Need at least 2 different conditions for differential analysis.\n")
  quit(status = 0)
}

# Create master peak list from all samples
cat("Creating master peak list...\n")
mPeak <- GRanges()

for (sample in samples) {
  peak_file <- file.path(peak_dir, paste0(sample, "_seacr_control.peaks.stringent.bed"))
  if (file.exists(peak_file)) {
    peakRes <- read.table(peak_file, header = FALSE, fill = TRUE)
    peakGR <- GRanges(seqnames = peakRes$V1, 
                      IRanges(start = peakRes$V2, end = peakRes$V3), 
                      strand = "*")
    mPeak <- append(mPeak, peakGR)
  }
}

# Merge overlapping peaks
masterPeak <- reduce(mPeak)
cat(paste("Master peak list contains", length(masterPeak), "peaks\n"))

# Get fragment counts for each peak
cat("Counting fragments in peaks...\n")
countMat <- matrix(NA, length(masterPeak), length(samples))

for (i in seq_along(samples)) {
  sample <- samples[i]
  bam_file <- file.path(bam_dir, paste0(sample, "_bowtie2.mapped.bam"))
  
  if (file.exists(bam_file)) {
    fragment_counts <- getCounts(bam_file, masterPeak, paired = TRUE, 
                                  by_rg = FALSE, format = "bam")
    countMat[, i] <- counts(fragment_counts)[, 1]
  }
}

colnames(countMat) <- samples

# Filter low count peaks
selectR <- which(rowSums(countMat, na.rm = TRUE) > 5)
dataS <- countMat[selectR, , drop = FALSE]
cat(paste("Retained", nrow(dataS), "peaks after filtering\n"))

# Run DESeq2
cat("Running DESeq2 differential analysis...\n")
condition_factor <- factor(conditions)

dds <- DESeqDataSetFromMatrix(countData = dataS,
                               colData = DataFrame(condition = condition_factor),
                               design = ~ condition)

DDS <- DESeq(dds)

# Get normalized counts
normDDS <- counts(DDS, normalized = TRUE)
colnames(normDDS) <- paste0(colnames(normDDS), "_norm")

# Get results
res <- results(DDS, independentFiltering = FALSE, altHypothesis = "greaterAbs")

# Combine results
countMatDiff <- cbind(as.data.frame(dataS), 
                      as.data.frame(normDDS), 
                      as.data.frame(res))

# Add peak coordinates
peak_coords <- as.data.frame(masterPeak[selectR])
countMatDiff <- cbind(peak_coords[, c("seqnames", "start", "end")], countMatDiff)

# Write results
output_file <- paste0(output_prefix, "_differential_peaks.tsv")
write.table(countMatDiff, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(paste("Results saved to:", output_file, "\n"))

# Summary statistics
sig_peaks <- sum(res$padj < 0.05, na.rm = TRUE)
cat(paste("\nSignificant peaks (padj < 0.05):", sig_peaks, "\n"))

# Save DESeq2 object
rds_file <- paste0(output_prefix, "_deseq2.rds")
saveRDS(DDS, rds_file)
cat(paste("DESeq2 object saved to:", rds_file, "\n"))

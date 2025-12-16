#!/usr/bin/env Rscript
#' CUT&Tag Alignment QC Visualization
#' 
#' Based on: https://yezhengstat.github.io/CUTTag_tutorial/
#' Generates boxplots for sequencing depth, alignment rate, duplication rate, etc.

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(viridis)
  library(ggpubr)
})

# Get parameters from command line or Snakemake
args <- commandArgs(trailingOnly = TRUE)

if (length(args) >= 2) {
  # Command line usage
  summary_file <- args[1]
  output_pdf <- args[2]
} else if (exists("snakemake")) {
  # Snakemake usage
  summary_file <- snakemake@input[["summary"]]
  output_pdf <- snakemake@output[["pdf"]]
} else {
  stop("Please provide input and output file paths.")
}

# Read summary data
alignSummary <- read.table(summary_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Extract numeric values from percentage strings
alignSummary <- alignSummary %>%
  mutate(
    AlignmentRate_hg38_num = as.numeric(gsub("%", "", AlignmentRate_hg38)),
    DuplicationRate_num = as.numeric(gsub("%", "", DuplicationRate))
  )

# Get unique histones
histList <- unique(alignSummary$Histone)
alignSummary$Histone <- factor(alignSummary$Histone, levels = histList)

# Extract replicate info from sample name
alignSummary <- alignSummary %>%
  mutate(Replicate = gsub(".*_", "", Sample))

# Figure A: Sequencing Depth
fig_A <- alignSummary %>%
  ggplot(aes(x = Histone, y = SequencingDepth / 1000000, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15), size = 2) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 14) +
  ylab("Sequencing Depth (Million)") +
  xlab("") +
  ggtitle("A. Sequencing Depth")

# Figure B: Mapped Fragments
fig_B <- alignSummary %>%
  ggplot(aes(x = Histone, y = MappedFragments_hg38 / 1000000, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15), size = 2) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 14) +
  ylab("Mapped Fragments (Million)") +
  xlab("") +
  ggtitle("B. Alignable Fragments (hg38)")

# Figure C: Alignment Rate
fig_C <- alignSummary %>%
  ggplot(aes(x = Histone, y = AlignmentRate_hg38_num, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15), size = 2) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 14) +
  ylab("Alignment Rate (%)") +
  xlab("") +
  ggtitle("C. Alignment Rate (hg38)")

# Figure D: Duplication Rate
fig_D <- alignSummary %>%
  ggplot(aes(x = Histone, y = DuplicationRate_num, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15), size = 2) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 14) +
  ylab("Duplication Rate (%)") +
  xlab("") +
  ggtitle("D. Duplication Rate")

# Figure E: Estimated Library Size
fig_E <- alignSummary %>%
  ggplot(aes(x = Histone, y = EstimatedLibrarySize, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15), size = 2) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 14) +
  ylab("Estimated Library Size") +
  xlab("") +
  ggtitle("E. Estimated Library Size")

# Figure F: Unique Fragments
fig_F <- alignSummary %>%
  ggplot(aes(x = Histone, y = UniqueFragments / 1000000, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15), size = 2) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 14) +
  ylab("Unique Fragments (Million)") +
  xlab("") +
  ggtitle("F. Unique Fragments")

# Combine plots and save
combined_plot <- ggarrange(fig_A, fig_B, fig_C, fig_D, fig_E, fig_F,
                           ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")

# Save to PDF
pdf(output_pdf, width = 15, height = 10)
print(combined_plot)
dev.off()

message("Alignment QC plots saved to: ", output_pdf)

#!/usr/bin/env Rscript
#' CUT&Tag Fragment Length Distribution Plot
#' 
#' Based on: https://yezhengstat.github.io/CUTTag_tutorial/
#' Generates violin plot and line plot of fragment length distributions.

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(viridis)
  library(ggpubr)
})

# Get parameters from Snakemake
fraglen_files <- snakemake@input[["fraglen"]]
samples <- snakemake@params[["samples"]]
histones <- snakemake@params[["histones"]]
output_pdf <- snakemake@output[["pdf"]]

# Collect fragment length data from all samples
fragLen <- data.frame()

for (i in seq_along(samples)) {
  sample <- samples[i]
  histone <- histones[i]
  
  # Read fragment length file
  frag_data <- read.table(fraglen_files[i], header = FALSE, 
                          col.names = c("fragLen", "fragCount"))
  
  # Add sample info
  frag_data <- frag_data %>%
    mutate(
      fragLen = as.numeric(fragLen),
      fragCount = as.numeric(fragCount),
      Weight = fragCount / sum(fragCount),
      Histone = histone,
      Sample = sample
    )
  
  fragLen <- rbind(fragLen, frag_data)
}

# Get unique histones for factor levels
histList <- unique(histones)
fragLen$Histone <- factor(fragLen$Histone, levels = histList)
fragLen$Sample <- factor(fragLen$Sample, levels = samples)

# Generate violin plot (fragment size distribution)
fig_violin <- fragLen %>%
  ggplot(aes(x = Sample, y = fragLen, weight = Weight, fill = Histone)) +
  geom_violin(bw = 5) +
  scale_y_continuous(breaks = seq(0, 800, 50)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Fragment Length (bp)") +
  xlab("") +
  ggtitle("A. Fragment Length Distribution (Violin)")

# Generate line plot (fragment size distribution)
fig_line <- fragLen %>%
  ggplot(aes(x = fragLen, y = fragCount, color = Histone, group = Sample)) +
  geom_line(linewidth = 0.8, alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 14) +
  xlab("Fragment Length (bp)") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500)) +
  ggtitle("B. Fragment Length Distribution (Line)")

# Combine plots and save
combined_plot <- ggarrange(fig_violin, fig_line, ncol = 2, common.legend = TRUE, legend = "bottom")

# Save to PDF
pdf(output_pdf, width = 14, height = 6)
print(combined_plot)
dev.off()

message("Fragment length plot saved to: ", output_pdf)

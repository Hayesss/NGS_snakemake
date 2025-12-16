#!/usr/bin/env Rscript
#' CUT&Tag Peak Calling QC and FRiP Calculation
#' 
#' Based on: https://yezhengstat.github.io/CUTTag_tutorial/
#' Calculates peak statistics, reproducibility, and FRiP scores.

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(viridis)
  library(ggpubr)
  library(GenomicRanges)
  library(chromVAR)
})

# Get parameters from command line
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  cat("Usage: Rscript peak_qc.R peak_dir bam_dir samples histones output_prefix\n")
  quit(status = 1)
}

peak_dir <- args[1]
bam_dir <- args[2]
samples <- strsplit(args[3], ",")[[1]]
histones <- strsplit(args[4], ",")[[1]]
output_prefix <- args[5]

# Filter out IgG samples
non_igg <- histones != "IgG"
target_samples <- samples[non_igg]
target_histones <- histones[non_igg]

histList <- unique(target_histones)

# 1. Count peaks and calculate width
cat("Analyzing peak statistics...\n")
peakN <- data.frame()
peakWidth <- data.frame()
peakType <- c("control", "top0.01")

for (i in seq_along(target_samples)) {
  sample <- target_samples[i]
  hist <- target_histones[i]
  rep <- gsub(".*_", "", sample)
  
  for (type in peakType) {
    peak_file <- file.path(peak_dir, paste0(sample, "_seacr_", type, ".peaks.stringent.bed"))
    
    if (file.exists(peak_file)) {
      peakInfo <- read.table(peak_file, header = FALSE, fill = TRUE)
      peakInfo$width <- abs(peakInfo$V3 - peakInfo$V2)
      
      peakN <- rbind(peakN, data.frame(
        peakN = nrow(peakInfo),
        peakType = type,
        Histone = hist,
        Replicate = rep,
        Sample = sample
      ))
      
      peakWidth <- rbind(peakWidth, data.frame(
        width = peakInfo$width,
        peakType = type,
        Histone = hist,
        Replicate = rep
      ))
    }
  }
}

peakN$Histone <- factor(peakN$Histone, levels = histList)
peakWidth$Histone <- factor(peakWidth$Histone, levels = histList)

# 2. Calculate peak reproducibility
cat("Calculating peak reproducibility...\n")
repL <- unique(gsub(".*_", "", target_samples))
peakOverlap <- data.frame()

for (type in peakType) {
  for (hist in histList) {
    overlap.gr <- GRanges()
    hist_samples <- target_samples[target_histones == hist]
    
    for (sample in hist_samples) {
      peak_file <- file.path(peak_dir, paste0(sample, "_seacr_", type, ".peaks.stringent.bed"))
      
      if (file.exists(peak_file)) {
        peakInfo <- read.table(peak_file, header = FALSE, fill = TRUE)
        peakInfo.gr <- GRanges(peakInfo$V1, 
                               IRanges(start = peakInfo$V2, end = peakInfo$V3), 
                               strand = "*")
        
        if (length(overlap.gr) > 0) {
          overlap.gr <- overlap.gr[findOverlaps(overlap.gr, peakInfo.gr)@from]
        } else {
          overlap.gr <- peakInfo.gr
        }
      }
    }
    
    peakOverlap <- rbind(peakOverlap, data.frame(
      peakReprod = length(overlap.gr),
      Histone = hist,
      peakType = type
    ))
  }
}

# Join with peak numbers
peakReprod <- left_join(peakN, peakOverlap, by = c("Histone", "peakType")) %>%
  mutate(peakReprodRate = peakReprod / peakN * 100)

# 3. Calculate FRiP (Fragment proportion in Peaks)
cat("Calculating FRiP scores...\n")
inPeakData <- data.frame()

for (i in seq_along(target_samples)) {
  sample <- target_samples[i]
  hist <- target_histones[i]
  rep <- gsub(".*_", "", sample)
  
  peak_file <- file.path(peak_dir, paste0(sample, "_seacr_control.peaks.stringent.bed"))
  bam_file <- file.path(bam_dir, paste0(sample, "_bowtie2.mapped.bam"))
  
  if (file.exists(peak_file) && file.exists(bam_file)) {
    peakRes <- read.table(peak_file, header = FALSE, fill = TRUE)
    peak.gr <- GRanges(seqnames = peakRes$V1, 
                       IRanges(start = peakRes$V2, end = peakRes$V3), 
                       strand = "*")
    
    fragment_counts <- getCounts(bam_file, peak.gr, paired = TRUE, 
                                  by_rg = FALSE, format = "bam")
    inPeakN <- sum(counts(fragment_counts)[, 1])
    
    inPeakData <- rbind(inPeakData, data.frame(
      inPeakN = inPeakN,
      Histone = hist,
      Replicate = rep,
      Sample = sample
    ))
  }
}

# Save peak summary
peak_summary <- peakN %>%
  select(Sample, Histone, Replicate, peakType, peakN)
write.table(peak_summary, paste0(output_prefix, "_peak_summary.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Save reproducibility summary
reprod_summary <- peakReprod %>%
  select(Sample, Histone, Replicate, peakType, peakN, peakReprod, peakReprodRate)
write.table(reprod_summary, paste0(output_prefix, "_reproducibility.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Save FRiP data
write.table(inPeakData, paste0(output_prefix, "_frip.tsv"), 
            sep = "\t", row.names = FALSE, quote = FALSE)

# Generate plots
cat("Generating plots...\n")

# Peak number plot
fig_A <- peakN %>%
  ggplot(aes(x = Histone, y = peakN, fill = Histone)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  facet_grid(~peakType) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 14) +
  ylab("Number of Peaks") +
  xlab("") +
  ggtitle("A. Number of Peaks")

# Peak width plot
fig_B <- peakWidth %>%
  ggplot(aes(x = Histone, y = width, fill = Histone)) +
  geom_violin() +
  facet_grid(Replicate ~ peakType) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_y_continuous(trans = "log10") +
  theme_bw(base_size = 14) +
  ylab("Width of Peaks (bp)") +
  xlab("") +
  ggtitle("B. Peak Width Distribution")

# Reproducibility plot
fig_C <- peakReprod %>%
  ggplot(aes(x = Histone, y = peakReprodRate, fill = Histone)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = round(peakReprodRate, 1)), vjust = -0.5, size = 3) +
  facet_grid(Replicate ~ peakType) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  theme_bw(base_size = 14) +
  ylab("% of Peaks Reproduced") +
  xlab("") +
  ggtitle("C. Peak Reproducibility")

# FRiP plot (if data available)
if (nrow(inPeakData) > 0) {
  inPeakData$Histone <- factor(inPeakData$Histone, levels = histList)
  
  fig_D <- inPeakData %>%
    ggplot(aes(x = Histone, y = inPeakN / 1000000, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 14) +
    ylab("Fragments in Peaks (Million)") +
    xlab("") +
    ggtitle("D. Fragments in Peaks")
  
  combined <- ggarrange(fig_A, fig_B, fig_C, fig_D, ncol = 2, nrow = 2, 
                        common.legend = TRUE, legend = "bottom")
} else {
  combined <- ggarrange(fig_A, fig_B, fig_C, ncol = 2, nrow = 2, 
                        common.legend = TRUE, legend = "bottom")
}

# Save plot
pdf(paste0(output_prefix, "_peak_qc.pdf"), width = 14, height = 12)
print(combined)
dev.off()

cat(paste("Peak QC analysis complete. Output files saved with prefix:", output_prefix, "\n"))

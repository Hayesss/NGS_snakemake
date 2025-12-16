# CUT&Tag Data Processing and Analysis Snakemake Pipeline

A Snakemake implementation of the [CUT&Tag Data Processing and Analysis Tutorial](https://yezhengstat.github.io/CUTTag_tutorial/) by Ye Zheng, Kami Ahmad, and Steven Henikoff.

## Overview

This pipeline processes CUT&Tag data from raw FASTQ files through:

1. **Quality Control** - FastQC analysis
2. **Alignment** - Bowtie2 mapping to reference genome (hg38)
3. **Spike-in Calibration** - E. coli alignment for normalization
4. **Duplicate Analysis** - Picard MarkDuplicates
5. **Fragment Analysis** - Size distribution assessment
6. **Peak Calling** - SEACR peak calling
7. **Visualization** - BigWig generation and heatmaps (deepTools)
8. **Reproducibility** - Replicate correlation analysis
9. **QC Summaries** - Alignment/peak QC plots and FRiP

## Requirements

### Software Dependencies

- **Python** >= 3.8
- **Snakemake** >= 6.0
- **FastQC** >= 0.11.9
- **Bowtie2** >= 2.3.4.3
- **samtools** >= 1.10
- **bedtools** >= 2.29.1
- **Picard** >= 2.18.29
- **SEACR** >= 1.3
- **deepTools** >= 2.0

### R Packages

- dplyr
- stringr
- ggplot2
- viridis
- GenomicRanges
- chromVAR
- DESeq2
- ggpubr
- corrplot

### Installation

1. **Create conda environment:**

```bash
conda env create -f workflow/envs/cuttag.yaml
conda activate cuttag
```

2. **Prepare reference genomes:**

```bash
bash workflow/scripts/prepare_references.sh /path/to/references 8
```

This will download and prepare:
- hg38 reference genome and Bowtie2 index
- E. coli (U00096.3) spike-in genome and Bowtie2 index
- Gene annotations
- SEACR peak caller

3. **Download tutorial data (optional):**

```bash
bash workflow/scripts/download_tutorial_data.sh .
```

## Configuration

### 1. Edit `config/config.yaml`

Update the following paths:

```yaml
# Sample information
samples: "config/samples.tsv"

# Input FASTQ directory
fastq_dir: "data/fastq"

# Reference genome paths
bowtie2_index: "/path/to/bowtie2Index/hg38"
chrom_sizes: "/path/to/hg38.chrom.sizes"
gene_bed: "/path/to/hg38_genes.bed"

# Spike-in reference
spikein:
  enabled: true
  bowtie2_index: "/path/to/bowtie2Index/Ecoli"
  multiplier: 10000

# SEACR path
seacr_path: "/path/to/SEACR/SEACR_1.3.sh"
```

### 2. Edit `config/samples.tsv`

Sample information file (tab-separated):

| sample | histone | replicate | control |
|--------|---------|-----------|---------|
| K27me3_rep1 | K27me3 | rep1 | IgG_rep1 |
| K27me3_rep2 | K27me3 | rep2 | IgG_rep2 |
| K4me3_rep1 | K4me3 | rep1 | IgG_rep1 |
| K4me3_rep2 | K4me3 | rep2 | IgG_rep2 |
| IgG_rep1 | IgG | rep1 | IgG_rep1 |
| IgG_rep2 | IgG | rep2 | IgG_rep2 |

### 3. FASTQ File Naming

Place FASTQ files in `data/fastq/` with naming convention:
- `{sample}_R1.fastq.gz`
- `{sample}_R2.fastq.gz`

Example:
- `K27me3_rep1_R1.fastq.gz`
- `K27me3_rep1_R2.fastq.gz`

## Usage

### Run the complete pipeline

```bash
# Dry run to check workflow
snakemake -n

# Run with 8 cores
snakemake --cores 8

# Run on cluster (SLURM example)
snakemake --cluster "sbatch -p queue -n {threads}" --jobs 10
```

### Run specific steps

```bash
# Only FastQC
snakemake --cores 4 results/fastqc/{sample}_R1_fastqc.html

# Only alignment
snakemake --cores 8 results/alignment/bam/{sample}_bowtie2.mapped.bam

# Only peak calling
snakemake --cores 4 results/peakCalling/SEACR/{sample}_seacr_control.peaks.stringent.bed
```

### Generate DAG visualization

```bash
snakemake --dag | dot -Tpdf > dag.pdf
```

## Output Structure

```
results/
├── fastqc/                          # FastQC reports
│   ├── {sample}_R1_fastqc.html
│   └── {sample}_R2_fastqc.html
├── alignment/
│   ├── sam/
│   │   ├── bowtie2_summary/         # Alignment statistics
│   │   │   ├── {sample}_bowtie2.txt
│   │   │   └── {sample}_bowtie2_spikeIn.txt
│   │   └── fragmentLen/             # Fragment length distributions
│   │       └── {sample}_fragmentLen.txt
│   ├── bam/                         # BAM files
│   │   ├── {sample}_bowtie2.mapped.bam
│   │   └── {sample}.sorted.bam
│   ├── bed/                         # BED files
│   │   ├── {sample}_bowtie2.fragments.bed
│   │   └── {sample}_bowtie2.fragmentsCount.bin500.bed
│   ├── bedgraph/                    # Normalized coverage
│   │   └── {sample}_bowtie2.fragments.normalized.bedgraph
│   ├── bigwig/                      # BigWig for visualization
│   │   └── {sample}_raw.bw
│   └── removeDuplicate/
│       └── picard_summary/          # Duplication metrics
│           └── {sample}_picard.dupMark.txt
├── peakCalling/
│   └── SEACR/                       # SEACR peaks
│       ├── {sample}_seacr_control.peaks.stringent.bed
│       └── {sample}_seacr_top0.01.peaks.stringent.bed
├── visualization/                   # Heatmaps
│   └── Histone_gene_heatmap.png
└── qc/                              # Summary reports
    ├── alignment_summary.tsv
    ├── alignment_qc_plots.pdf
    ├── fragment_length_summary.pdf
    ├── replicate_correlation.pdf
    └── peak/
        ├── peak_qc_peak_summary.tsv
        ├── peak_qc_reproducibility.tsv
        ├── peak_qc_frip.tsv
        └── peak_qc.pdf
```

## Pipeline Workflow

```
                    ┌──────────────────────────────────┐
                    │          Raw FASTQ Files         │
                    └──────────────┬───────────────────┘
                                   │
                    ┌──────────────▼───────────────────┐
                    │            FastQC                │
                    └──────────────┬───────────────────┘
                                   │
          ┌────────────────────────┴────────────────────────┐
          │                                                  │
┌─────────▼─────────┐                            ┌──────────▼──────────┐
│  Bowtie2 (hg38)   │                            │ Bowtie2 (E. coli)   │
└─────────┬─────────┘                            └──────────┬──────────┘
          │                                                  │
┌─────────▼─────────┐                            ┌──────────▼──────────┐
│  Picard MarkDup   │                            │   Spike-in Depth    │
└─────────┬─────────┘                            └──────────┬──────────┘
          │                                                  │
┌─────────▼─────────┐                                        │
│  SAM → BAM → BED  │                                        │
└─────────┬─────────┘                                        │
          │                                                  │
┌─────────▼─────────┐                            ┌──────────▼──────────┐
│  Fragment BED     │◄───────────────────────────│  Normalization      │
└─────────┬─────────┘                            └──────────┬──────────┘
          │                                                  │
          │                               ┌──────────────────┘
          │                               │
┌─────────▼─────────┐            ┌────────▼────────┐
│ Fragment Length   │            │ Normalized      │
│ Distribution      │            │ BedGraph        │
└─────────┬─────────┘            └────────┬────────┘
          │                               │
          │                      ┌────────▼────────┐
          │                      │  SEACR Peak     │
          │                      │  Calling        │
          │                      └────────┬────────┘
          │                               │
┌─────────▼─────────────────────┬────────▼────────┐
│      QC Summary Reports       │    BigWig &     │
│                               │    Heatmaps     │
└───────────────────────────────┴─────────────────┘
```

## Key Parameters (from tutorial)

### Bowtie2 Alignment

- `--end-to-end --very-sensitive` - Sensitive end-to-end alignment
- `--no-mixed --no-discordant` - Only concordant pairs
- `-I 10 -X 700` - Insert size 10-700 bp
- For spike-in: add `--no-overlap --no-dovetail`
- Mapping quality filtering: controlled by `min_mapping_quality` in `config/config.yaml` (set to 0 to disable when converting SAM → BAM)

### SEACR Peak Calling

- **With control:** Uses IgG sample as background
- **Without control:** Top 0.01 (1%) threshold
- **Normalization:** `non` if spike-in normalized, `norm` otherwise

### Spike-in Calibration

```
Scaling factor S = C / (fragments mapped to E. coli)
Normalized coverage = primary_genome_coverage × S
```

Where C is typically 10,000.

## Downstream Analysis

### Differential Analysis (DESeq2)

Run the differential analysis script:

```bash
Rscript workflow/scripts/differential_analysis.R \
    results/peakCalling/SEACR \
    results/alignment/bam \
    "K27me3,K27me3,K4me3,K4me3" \
    "K27me3_rep1,K27me3_rep2,K4me3_rep1,K4me3_rep2" \
    results/differential/comparison
```

### Replicate Correlation

Automatically generated as `results/qc/replicate_correlation.pdf`. Re-run manually if needed:

```bash
Rscript workflow/scripts/replicate_correlation.R \
    "results/alignment/bed/K27me3_rep1_bowtie2.fragmentsCount.bin500.bed,..." \
    "K27me3_rep1,K27me3_rep2,..." \
    results/qc/replicate_correlation.pdf
```

### Peak QC and FRiP

Automatically generated under `results/qc/peak/`. Re-run manually if needed:

```bash
Rscript workflow/scripts/peak_qc.R \
    results/peakCalling/SEACR \
    results/alignment/bam \
    "K27me3_rep1,K27me3_rep2,K4me3_rep1,K4me3_rep2" \
    "K27me3,K27me3,K4me3,K4me3" \
    results/qc/peak/peak_qc
```

## Troubleshooting

### Per base sequence content fails FastQC

This is **normal** for CUT&Tag data due to Tn5 preference. The discordant sequence content at the beginning of reads does not affect alignment.

### High duplication rate in IgG samples

Expected behavior. IgG control samples have high duplication rates since reads derive from non-specific tagmentation.

### Low spike-in alignment rate

For abundant targets (e.g., H3K27me3), E. coli reads can be 0.01-10% of total. For IgG controls, this percentage is typically much higher.

## Citation

If you use this pipeline, please cite:

1. **CUT&Tag method:**
   > Kaya-Okur HS, et al. (2019). CUT&Tag for efficient epigenomic profiling of small samples and single cells. Nature Communications 10:1930.

2. **SEACR:**
   > Meers MP, et al. (2019). Peak calling by Sparse Enrichment Analysis for CUT&RUN chromatin profiling. Epigenetics & Chromatin 12:42.

3. **Original tutorial:**
   > Zheng Y, et al. (2020). CUT&Tag Data Processing and Analysis Tutorial. https://yezhengstat.github.io/CUTTag_tutorial/

## License

This pipeline implementation follows the original tutorial structure and is provided for educational and research purposes.

## Contact

For issues related to:
- **This Snakemake implementation:** Open an issue in this repository
- **Original tutorial:** yzheng23@fredhutch.org

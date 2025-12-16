#!/bin/bash
#' CUT&Tag Reference Genome Preparation Script
#' 
#' Downloads and prepares reference genomes and annotations for CUT&Tag analysis.
#' Based on: https://yezhengstat.github.io/CUTTag_tutorial/

set -e

REF_DIR="${1:-references}"
THREADS="${2:-8}"

mkdir -p "${REF_DIR}"
cd "${REF_DIR}"

echo "=== CUT&Tag Reference Preparation ==="
echo "Reference directory: ${REF_DIR}"
echo "Threads: ${THREADS}"
echo ""

# =============================================================================
# 1. Download hg38 reference genome
# =============================================================================
echo "=== Downloading hg38 reference genome ==="

if [ ! -f "hg38.fa" ]; then
    wget -O hg38.fa.gz \
        https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
    gunzip hg38.fa.gz
fi

# Build Bowtie2 index
if [ ! -f "bowtie2_index/hg38.1.bt2" ]; then
    echo "Building Bowtie2 index for hg38..."
    mkdir -p bowtie2_index
    bowtie2-build --threads ${THREADS} hg38.fa bowtie2_index/hg38
fi

# Generate chromosome sizes
if [ ! -f "hg38.chrom.sizes" ]; then
    echo "Generating chromosome sizes..."
    samtools faidx hg38.fa
    cut -f1,2 hg38.fa.fai > hg38.chrom.sizes
fi

# =============================================================================
# 2. Download E. coli spike-in reference
# =============================================================================
echo ""
echo "=== Downloading E. coli spike-in reference ==="

if [ ! -f "Ecoli_U00096.fa" ]; then
    # E. coli K-12 MG1655 (U00096.3)
    wget -O Ecoli_U00096.fa.gz \
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
    gunzip Ecoli_U00096.fa.gz
    mv GCF_000005845.2_ASM584v2_genomic.fna Ecoli_U00096.fa 2>/dev/null || true
fi

# Build Bowtie2 index for E. coli
if [ ! -f "bowtie2_index/Ecoli.1.bt2" ]; then
    echo "Building Bowtie2 index for E. coli..."
    bowtie2-build --threads ${THREADS} Ecoli_U00096.fa bowtie2_index/Ecoli
fi

# =============================================================================
# 3. Download gene annotation
# =============================================================================
echo ""
echo "=== Downloading gene annotations ==="

if [ ! -f "hg38_genes.bed" ]; then
    # Download RefSeq genes from UCSC
    wget -O hg38_refGene.txt.gz \
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz"
    gunzip hg38_refGene.txt.gz
    
    # Convert to BED format (chr, txStart, txEnd, name, score, strand)
    awk -v OFS="\t" '{print $3, $5, $6, $13, 0, $4}' hg38_refGene.txt | \
        sort -k1,1 -k2,2n | uniq > hg38_genes.bed
    
    rm hg38_refGene.txt
fi

# =============================================================================
# 4. Download SEACR
# =============================================================================
echo ""
echo "=== Downloading SEACR ==="

if [ ! -d "SEACR" ]; then
    git clone https://github.com/FredHutch/SEACR.git
fi

echo ""
echo "=== Reference Preparation Complete ==="
echo ""
echo "Reference files created:"
echo "  - hg38.fa (reference genome)"
echo "  - bowtie2_index/hg38 (Bowtie2 index)"
echo "  - hg38.chrom.sizes (chromosome sizes)"
echo "  - Ecoli_U00096.fa (E. coli spike-in genome)"
echo "  - bowtie2_index/Ecoli (E. coli Bowtie2 index)"
echo "  - hg38_genes.bed (gene annotations)"
echo "  - SEACR/ (SEACR peak caller)"
echo ""
echo "Update config.yaml with these paths:"
echo "  bowtie2_index: \"${REF_DIR}/bowtie2_index/hg38\""
echo "  chrom_sizes: \"${REF_DIR}/hg38.chrom.sizes\""
echo "  gene_bed: \"${REF_DIR}/hg38_genes.bed\""
echo "  spikein.bowtie2_index: \"${REF_DIR}/bowtie2_index/Ecoli\""
echo "  seacr_path: \"${REF_DIR}/SEACR/SEACR_1.3.sh\""

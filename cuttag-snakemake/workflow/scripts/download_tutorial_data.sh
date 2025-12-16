#!/bin/bash
#' CUT&Tag Tutorial Data Download Script
#' 
#' Downloads example data from GEO as described in the tutorial:
#' https://yezhengstat.github.io/CUTTag_tutorial/
#'
#' Data from Kaya-Okur et al. (2020)
#' GEO: GSE145187

set -e

# Project path
PROJ_PATH="${1:-.}"
DATA_DIR="${PROJ_PATH}/data"
FASTQ_DIR="${PROJ_PATH}/data/fastq"

mkdir -p "${DATA_DIR}"
mkdir -p "${FASTQ_DIR}"

echo "=== CUT&Tag Tutorial Data Download ==="
echo "Project path: ${PROJ_PATH}"
echo "Data will be downloaded to: ${DATA_DIR}"
echo ""

# H3K27me3 replicate 1 (SRX8754646)
echo "Downloading H3K27me3_rep1..."
mkdir -p "${DATA_DIR}/K27me3_rep1"
wget -O "${DATA_DIR}/K27me3_rep1/K27me3_rep1_R1_001.fastq.gz" \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR122/034/SRR12253034/SRR12253034_1.fastq.gz
wget -O "${DATA_DIR}/K27me3_rep1/K27me3_rep1_R2_001.fastq.gz" \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR122/034/SRR12253034/SRR12253034_2.fastq.gz

# H3K27me3 replicate 2 (SRX7713678)
echo "Downloading H3K27me3_rep2..."
mkdir -p "${DATA_DIR}/K27me3_rep2"
wget -O "${DATA_DIR}/K27me3_rep2/K27me3_rep2_R1_001.fastq.gz" \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/041/SRR11074241/SRR11074241_1.fastq.gz
wget -O "${DATA_DIR}/K27me3_rep2/K27me3_rep2_R2_001.fastq.gz" \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/041/SRR11074241/SRR11074241_2.fastq.gz

# H3K4me3 replicate 1 (SRX7713692)
echo "Downloading K4me3_rep1..."
mkdir -p "${DATA_DIR}/K4me3_rep1"
wget -O "${DATA_DIR}/K4me3_rep1/K4me3_rep1_R1_001.fastq.gz" \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/055/SRR11074255/SRR11074255_1.fastq.gz
wget -O "${DATA_DIR}/K4me3_rep1/K4me3_rep1_R2_001.fastq.gz" \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/055/SRR11074255/SRR11074255_2.fastq.gz

# H3K4me3 replicate 2 (SRX7713696)
echo "Downloading K4me3_rep2..."
mkdir -p "${DATA_DIR}/K4me3_rep2"
wget -O "${DATA_DIR}/K4me3_rep2/K4me3_rep2_R1_001.fastq.gz" \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/059/SRR11074259/SRR11074259_1.fastq.gz
wget -O "${DATA_DIR}/K4me3_rep2/K4me3_rep2_R2_001.fastq.gz" \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/059/SRR11074259/SRR11074259_2.fastq.gz

# IgG replicate 1 (SRX8468909)
echo "Downloading IgG_rep1..."
mkdir -p "${DATA_DIR}/IgG_rep1"
wget -O "${DATA_DIR}/IgG_rep1/IgG_rep1_R1_001.fastq.gz" \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/078/SRR11974578/SRR11974578_1.fastq.gz
wget -O "${DATA_DIR}/IgG_rep1/IgG_rep1_R2_001.fastq.gz" \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/078/SRR11974578/SRR11974578_2.fastq.gz

# IgG replicate 2 (SRX5545346)
echo "Downloading IgG_rep2..."
mkdir -p "${DATA_DIR}/IgG_rep2"
wget -O "${DATA_DIR}/IgG_rep2/IgG_rep2_R1_001.fastq.gz" \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR875/001/SRR8754611/SRR8754611_1.fastq.gz
wget -O "${DATA_DIR}/IgG_rep2/IgG_rep2_R2_001.fastq.gz" \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR875/001/SRR8754611/SRR8754611_2.fastq.gz
wget -O "${DATA_DIR}/IgG_rep2/IgG_rep2_R1_002.fastq.gz" \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR875/002/SRR8754612/SRR8754612_1.fastq.gz
wget -O "${DATA_DIR}/IgG_rep2/IgG_rep2_R2_002.fastq.gz" \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR875/002/SRR8754612/SRR8754612_2.fastq.gz

echo ""
echo "=== Merging technical replicates/lanes ==="

# Merge lanes for each sample
for sample in K27me3_rep1 K27me3_rep2 K4me3_rep1 K4me3_rep2 IgG_rep1 IgG_rep2; do
    echo "Merging ${sample}..."
    cat "${DATA_DIR}/${sample}"/*_R1_*.fastq.gz > "${FASTQ_DIR}/${sample}_R1.fastq.gz"
    cat "${DATA_DIR}/${sample}"/*_R2_*.fastq.gz > "${FASTQ_DIR}/${sample}_R2.fastq.gz"
done

echo ""
echo "=== Download Complete ==="
echo "FASTQ files are ready in: ${FASTQ_DIR}"
echo ""
echo "Sample files:"
ls -lh "${FASTQ_DIR}"

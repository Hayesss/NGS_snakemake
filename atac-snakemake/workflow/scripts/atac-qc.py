#!/usr/bin/env python
# coding: utf-8
"""
ATAC-seq QC Script
==================
This script generates quality control metrics and plots for ATAC-seq data.

Metrics computed:
  - Total raw reads
  - Overall alignment rate
  - Final mapped reads (after filtering)
  - Mitochondrial read percentage
  - Blacklist region overlap percentage
  - MAPQ filtered percentage
  - Duplicate percentage
  - TSS enrichment score and plot
  - Fragment size distribution plot

Usage:
    python atac-qc.py <qc_table> <frag_dist_fig> <tss_fig> <sample_name> \
                      <align_log> <mito_bam> <dedup_bam> <shift_bed> \
                      <blacklist> <picard_log> <tss_bed>

Reference:
    https://github.com/QuKunLab/ATAC-pipe
"""

import sys
import os
import re
import pysam
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from collections import Counter

def GetlineNum(inf):
    """Get line count for BAM or BED file."""
    if inf.endswith('.bam'):
        return int((os.popen(f"bedtools bamtobed -i {inf} | wc -l")).read().split()[0])
    elif inf.endswith('.bed'):
        return int((os.popen(f"wc -l {inf}")).read().split()[0])
    return 0

def GetlineNumOfAinB(inf1, inf2):
    """Get count of regions in inf1 that overlap with inf2."""
    return int(os.popen(f"bedtools intersect -a {inf1} -b {inf2} -u | wc -l").read().split()[0])

def GetlineNumOfMapping(inf):
    """Parse bowtie2 alignment log to get total reads and alignment rate."""
    total_reads_re = re.compile(r'(\d+) reads;')
    alignment_rate_re = re.compile(r'([0-9.]+)% overall')
    content = ''.join(open(inf, 'r').readlines())
    total_reads = int(total_reads_re.search(content).groups()[0])
    alignment_rate = float(alignment_rate_re.search(content).groups()[0])
    return total_reads, alignment_rate

def GetlineNumOfPicard(inf):
    """Parse Picard MarkDuplicates log to get duplicate count."""
    match = re.search(r'Marking (.*) records', open(inf, 'r').read())
    if match:
        return int(match.group(1))
    return 0

def FragmentDistPlot(inBAM, outFIG):
    """Generate fragment size distribution plot."""
    print(f"Generating fragment distribution plot...")
    sam = os.popen(f"samtools view -f 0x0002 {inBAM}")
    fragL = []
    for line in sam:
        items = line.rstrip("\n").split('\t')
        flag = int(items[1])
        flag_length = abs(int(items[8]))
        if (flag == 99 or flag == 163) and flag_length <= 600:
            fragL.append(flag_length)
    
    if not fragL:
        print("Warning: No valid fragments found for distribution plot")
        return
    
    count = Counter(fragL)
    plt.figure(figsize=(10.0, 4.0))
    count = dict(sorted(count.items(), key=lambda item: item[0]))
    plt.plot(np.array(list(count.keys())), 
             np.array(list(count.values())) / float(sum(count.values())), 'r')
    plt.xlabel('Fragment length (bp)')
    plt.ylabel('Fraction of reads')
    plt.title('Fragment Size Distribution')
    plt.axvline(x=120, color='gray', linestyle='--', alpha=0.5, label='NFR cutoff')
    plt.axvline(x=180, color='gray', linestyle='--', alpha=0.5)
    plt.axvline(x=247, color='gray', linestyle='--', alpha=0.5, label='Mono-nucleosome')
    plt.legend()
    plt.savefig(outFIG, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Fragment distribution saved: {outFIG}")

def TSSEnrichmentScore(inBAM, tssBED, outFIG):
    """Calculate TSS enrichment score and generate plot."""
    print(f"Calculating TSS enrichment score...")
    
    def asn_mat(val, mat, s_int, e_int, t, strand, weight):
        if float(val) > s_int and float(val) < e_int - 1 and t < 1000:
            if strand == '+':
                base = val - s_int
            else:
                base = e_int - val
            mat[int(t)][int(base)] += weight
        return mat
    
    def sub_Mat(inbed, inBAM):
        mat = np.zeros([1000, 4000])
        bedfile = np.loadtxt(inbed, 'str')
        bamfile = pysam.Samfile(inBAM, "rb")
        rows_bed = len(bedfile)
        
        for i in range(0, rows_bed):
            if len(bedfile[i]) < 6:
                continue
            chrom, start, end, name, score, strand = bedfile[i][0:6]
            tss = int(start) + (int(end) - int(start)) / 2
            tss_left = tss - 2000
            tss_right = tss + 2000
            
            try:
                for p2_rds in bamfile.fetch(str(chrom), max(0, int(tss_left - 2000)), int(tss_right + 2000)):
                    if p2_rds.mapq < 30:
                        continue
                    if p2_rds.is_reverse:
                        continue
                    else:
                        l_pos = p2_rds.pos + 4
                        ilen = abs(p2_rds.tlen) - 9
                        r_pos = l_pos + ilen
                    mat = asn_mat(l_pos, mat, tss_left, tss_right, ilen, strand, 1)
                    mat = asn_mat(r_pos, mat, tss_left, tss_right, ilen, strand, 1)
            except:
                continue
        return mat
    
    Mat = sub_Mat(tssBED, inBAM)
    mat0 = np.sum(Mat, 0)
    factor = np.mean(mat0[1:200])
    
    if factor == 0:
        print("Warning: Could not calculate TSS enrichment (no signal)")
        return 0
    
    plt.figure(figsize=(8.0, 5.0))
    plt.plot(range(-2000, 2000)[:-1], (mat0 / factor)[:-1], 'k.', alpha=0.3)
    smoothed = np.convolve(mat0, np.ones(20), 'same') / 20 / factor
    plt.plot(range(-2000, 2000)[:-1], smoothed[:-1], 'r', linewidth=2)
    plt.xlabel('Distance to TSS (bp)')
    plt.ylabel('Enrichment score')
    plt.title('TSS Enrichment')
    plt.axhline(y=1, color='gray', linestyle='--', alpha=0.5)
    plt.savefig(outFIG, dpi=150, bbox_inches='tight')
    plt.close()
    
    es = max(smoothed)
    print(f"TSS enrichment score: {es:.2f}")
    print(f"TSS enrichment plot saved: {outFIG}")
    return es

def main():
    if len(sys.argv) < 12:
        print("Usage: python atac-qc.py <qc_table> <frag_dist_fig> <tss_fig> <sample_name> "
              "<align_log> <mito_bam> <dedup_bam> <shift_bed> <blacklist> <picard_log> <tss_bed>")
        sys.exit(1)
    
    outqc, FragDistFig, tssEnrichFig, name, maplog, chrMbam, dechrMrmdupbam, dechrMrmdupbed, BL, picardlog, TSS = sys.argv[1:]
    
    print(f"Generating QC metrics for sample: {name}")
    
    # Create output directory if needed
    os.makedirs(os.path.dirname(outqc), exist_ok=True)
    
    # Generate plots
    FragmentDistPlot(dechrMrmdupbam, FragDistFig)
    TSSenrichScore = TSSEnrichmentScore(dechrMrmdupbam, TSS, tssEnrichFig)
    
    # Calculate QC metrics
    print("Calculating QC metrics...")
    
    totalReads, OverallAlignmentRate = GetlineNumOfMapping(maplog)
    print(f"  Total reads: {totalReads}")
    print(f"  Overall alignment rate: {OverallAlignmentRate}%")
    
    chrMCount = GetlineNum(chrMbam) / 2
    chrMPercent = chrMCount / float(totalReads) * 100
    print(f"  Mitochondrial reads: {chrMPercent:.2f}%")
    
    mappedReads = GetlineNum(dechrMrmdupbam) / 2
    mappedPercent = mappedReads / float(totalReads) * 100
    print(f"  Final mapped reads: {mappedReads:.0f} ({mappedPercent:.2f}%)")
    
    BLCount = GetlineNumOfAinB(dechrMrmdupbed, BL) / 2
    BLPercent = BLCount / float(totalReads) * 100
    print(f"  Blacklist overlap: {BLPercent:.2f}%")
    
    TSSCount = GetlineNumOfAinB(dechrMrmdupbed, TSS) / 2
    TSSPercent = TSSCount / float(totalReads) * 100
    print(f"  TSS overlap: {TSSPercent:.2f}%")
    
    DupCount = GetlineNumOfPicard(picardlog) / 2
    DupPercent = DupCount / float(totalReads) * 100
    print(f"  Duplicate rate: {DupPercent:.2f}%")
    
    qcPercent = OverallAlignmentRate - mappedPercent - chrMPercent - BLPercent - DupPercent
    
    # Write QC table
    qc = open(outqc, 'w')
    qc.write('\t'.join([
        "Sample", "TotalRawReads", "OverallAlignmentRate%", 
        "FinalMappedReads", "FinalMapped%", "chrM%", 
        "BlackListReads%", "MAPQFiltered%", "Duplicate%", 
        "TSS%", "TSSEnrichScore"
    ]) + '\n')
    qc.write(f"{name}\t")
    qc.write("%.0f\t%.2f\t%.0f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f" % (
        totalReads, OverallAlignmentRate, mappedReads, mappedPercent,
        chrMPercent, BLPercent, qcPercent, DupPercent, TSSPercent, TSSenrichScore
    ))
    qc.close()
    
    print(f"QC table saved: {outqc}")

if __name__ == "__main__":
    main()

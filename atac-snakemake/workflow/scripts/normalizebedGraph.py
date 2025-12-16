#!/usr/bin/env python
# coding: utf-8
"""
Normalize BedGraph Script for ATAC-seq Pipeline
================================================
This script normalizes sequencing depth to 10 million reads (RPM normalization).
This allows comparison of signal intensity across samples with different 
sequencing depths.

Usage:
    python normalizebedGraph.py <input.bedGraph> <output.bedGraph> [read_length]

Arguments:
    input.bedGraph:  Input bedGraph file
    output.bedGraph: Output normalized bedGraph file
    read_length:     Read length for estimating read count (default: 50)

Reference:
    https://github.com/QuKunLab/ATAC-pipe
"""

import sys
import os

def main():
    if len(sys.argv) < 3:
        print("Usage: python normalizebedGraph.py <input.bedGraph> <output.bedGraph> [read_length]")
        sys.exit(1)
    
    inBedGraph = sys.argv[1]
    outBedGraph = sys.argv[2]
    READLEN = 50 if len(sys.argv) <= 3 else int(sys.argv[3])
    FACTOR = 1e7  # Normalize to 10 million reads
    
    print(f"Input: {inBedGraph}")
    print(f"Output: {outBedGraph}")
    print(f"Read length: {READLEN}")
    
    # First pass: calculate total coverage
    total_coverage = 0
    with open(inBedGraph, 'r') as fin:
        for line in fin:
            items = line.rstrip('\n').split()
            if len(items) < 4:
                continue
            chrom, start, end, counts = items[:4]
            total_coverage += float(counts) * (int(end) - int(start))
    
    total_reads = total_coverage / READLEN
    print(f"Estimated total reads: {total_reads:.0f}")
    print(f"Normalization factor: {FACTOR / total_reads:.6f}")
    
    # Second pass: normalize and write output
    fo = open(outBedGraph, 'w')
    with open(inBedGraph, 'r') as fin:
        for line in fin:
            items = line.rstrip('\n').split()
            if len(items) < 4:
                continue
            chrom, start, end, counts = items[:4]
            value = round(float(counts) / total_reads * FACTOR, 4)
            fo.write('\t'.join([chrom, start, end, str(value)]) + '\n')
    fo.close()
    
    print(f"Normalization complete: {outBedGraph}")

if __name__ == "__main__":
    main()

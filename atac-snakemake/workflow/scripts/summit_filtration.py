#!/usr/bin/env python
# coding: utf-8
"""
Summit Filtration Script for ATAC-seq Pipeline
===============================================
This script filters summit positions based on quality-filtered peaks.
Only summits whose peak IDs are present in the filtered peak file are retained.
The output summit regions are extended by a specified bandwidth on each side.

Usage:
    python summit_filtration.py <summit_file> <peak_file> [bandwidth]

Arguments:
    summit_file: Path to the MACS2 summit file (*_summits.bed)
    peak_file:   Path to the filtered peak file (*.narrowPeak.Q0.05.rm_blacklist.bed)
    bandwidth:   Extend N bp on each side of summit (default: 250)

Output:
    Creates <summit_file>.filtered.bed with extended summit regions
"""

import sys
import os

def main():
    if len(sys.argv) < 3:
        print("Usage: python summit_filtration.py <summit_file> <peak_file> [bandwidth]")
        sys.exit(1)
    
    FILE_SUMMIT = sys.argv[1]
    FILE_PEAK = sys.argv[2]
    BANDWIDTH = int(sys.argv[3]) if len(sys.argv) > 3 else 250
    
    FILE_OUT = FILE_SUMMIT.replace(".bed", ".filtered.bed")
    
    # Put the filtered peak_ids into hash
    h0 = {}
    with open(FILE_PEAK) as fi:
        for line in fi:
            items = line.rstrip().split('\t')
            if len(items) >= 4:
                peak_id = items[3]
                h0[peak_id] = 1
    
    # Filter the summits and extend by bandwidth
    fo = open(FILE_OUT, 'w')
    with open(FILE_SUMMIT) as fi:
        for line in fi:
            items = line.rstrip().split('\t')
            if len(items) < 5:
                continue
            chrom = items[0]
            summit = int(items[1])
            peak_id = items[3]
            abs_log10_qval = items[4]
            start = max(0, summit - BANDWIDTH)
            end = summit + BANDWIDTH
            if peak_id in h0:
                results = '\t'.join([chrom, str(start), str(end), peak_id, abs_log10_qval])
                fo.write(results + '\n')
    fo.close()
    
    print(f"Filtered summits written to: {FILE_OUT}")
    print(f"Total filtered summits: {len(h0)}")

if __name__ == "__main__":
    main()


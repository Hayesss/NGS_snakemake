#!/usr/bin/env python
# coding: utf-8
"""
Summit Rename Script for ATAC-seq Pipeline
===========================================
This script renames merged summit regions with unique identifiers based on
their genomic coordinates (chr_start_end format).

Usage:
    python summit_rename.py <merged_summit_file>

Arguments:
    merged_summit_file: Path to the merged summit file from bedtools merge

Output:
    Creates two files:
    - <merged_summit_file>.rename.bed: BED file with renamed peak IDs
    - <merged_summit_file>.txt: Mapping of new IDs to original collapsed IDs
"""

import sys
import os

def main():
    if len(sys.argv) < 2:
        print("Usage: python summit_rename.py <merged_summit_file>")
        sys.exit(1)
    
    FILE_SUMMIT = sys.argv[1]
    FILE_OUT_BED = FILE_SUMMIT.replace('.bed', '.rename.bed')
    FILE_OUT_NAME = FILE_SUMMIT.replace('.bed', '.txt')
    
    fo1 = open(FILE_OUT_BED, 'w')
    fo2 = open(FILE_OUT_NAME, 'w')
    
    # Write header to name mapping file
    fo2.write("new_peak_id\toriginal_peak_ids\n")
    
    count = 0
    with open(FILE_SUMMIT) as fi:
        for line in fi:
            items = line.rstrip().split('\t')
            if len(items) < 4:
                continue
            chrom, start, end = items[0], items[1], items[2]
            original_ids = items[3]
            
            # Create new peak ID based on coordinates
            peak_id = '_'.join([chrom, start, end])
            
            # Write BED file with new peak ID
            new_line = '\t'.join([chrom, start, end, peak_id]) + '\n'
            fo1.write(new_line)
            
            # Write name mapping
            fo2.write(f'{peak_id}\t{original_ids}\n')
            count += 1
    
    fo1.close()
    fo2.close()
    
    print(f"Renamed summit file written to: {FILE_OUT_BED}")
    print(f"Name mapping file written to: {FILE_OUT_NAME}")
    print(f"Total merged summits: {count}")

if __name__ == "__main__":
    main()


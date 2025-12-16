#!/usr/bin/env python
# coding: utf-8
"""
Count Table Merge Script for ATAC-seq Pipeline
===============================================
This script merges individual sample count tables into a single matrix.

Usage:
    python count_table_merge.py <counts_directory>

Arguments:
    counts_directory: Directory containing *.readcount files

Output:
    Creates counts.tsv in the parent directory of counts_directory
    - Rows: Peak IDs
    - Columns: Sample names
    - Values: Read counts
"""

import sys
import os
import pandas as pd
from pathlib import Path

def main():
    if len(sys.argv) < 2:
        print("Usage: python count_table_merge.py <counts_directory>")
        sys.exit(1)
    
    DIR_COUNTS = sys.argv[1]
    FILE_OUT_COUNTS = os.path.join(DIR_COUNTS, "counts.tsv")
    
    # Find all readcount files
    files_counts = [
        os.path.join(DIR_COUNTS, fn) 
        for fn in os.listdir(DIR_COUNTS) 
        if fn.endswith('.readcount')
    ]
    
    if not files_counts:
        print(f"No .readcount files found in {DIR_COUNTS}")
        sys.exit(1)
    
    print(f"Found {len(files_counts)} count files to merge:")
    for f in files_counts:
        print(f"  - {os.path.basename(f)}")
    
    # Read and merge count tables
    list_df_counts = []
    for fn in files_counts:
        try:
            df_counts = pd.read_csv(fn, sep='\t', index_col=0)
            list_df_counts.append(df_counts)
            print(f"Loaded: {os.path.basename(fn)} ({len(df_counts)} peaks)")
        except Exception as e:
            print(f"Warning: Failed to read {fn}: {e}")
    
    if not list_df_counts:
        print("Error: No valid count files could be loaded")
        sys.exit(1)
    
    # Concatenate all count tables
    df_counts = pd.concat(list_df_counts, axis=1)
    
    # Fill missing values with 0 (peaks not found in some samples)
    df_counts = df_counts.fillna(0).astype(int)
    
    # Save merged count table
    df_counts.to_csv(FILE_OUT_COUNTS, sep='\t', index_label='Peak_ID')
    
    print(f"\nMerged count table written to: {FILE_OUT_COUNTS}")
    print(f"Total peaks: {len(df_counts)}")
    print(f"Total samples: {len(df_counts.columns)}")

if __name__ == "__main__":
    main()


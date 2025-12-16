#!/usr/bin/env python
# coding: utf-8
"""
Count Table Generation Script for ATAC-seq Pipeline
====================================================
This script counts reads overlapping each peak region.
Uses bedtools intersect to find overlaps and counts reads that overlap
by at least 25bp.

Usage:
    python countTable.py <shift.bed> <peaks.bed> <output_dir> <sample_name>

Arguments:
    shift.bed:    Input shift BED file (from Tn5 shift step)
    peaks.bed:    Peak file (narrowPeak or merged summits)
    output_dir:   Output directory
    sample_name:  Sample name for output file naming

Output:
    Creates <output_dir>/<sample_name>.readcount with columns:
    - PeakID: Peak identifier
    - <sample_name>: Read count for this sample

Reference:
    https://github.com/QuKunLab/ATAC-pipe
"""

import sys
import os

def main():
    if len(sys.argv) < 5:
        print("Usage: python countTable.py <shift.bed> <peaks.bed> <output_dir> <sample_name>")
        sys.exit(1)
    
    shiftBED, mergePEAK, outdir, sn = sys.argv[1:5]
    
    print(f"Shift BED: {shiftBED}")
    print(f"Peak file: {mergePEAK}")
    print(f"Sample: {sn}")
    
    # Create output directory if needed
    os.makedirs(outdir, exist_ok=True)
    
    # Run bedtools intersect
    # https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html
    tmp_file = os.path.join(outdir, f'{sn}.count')
    cmd = f"bedtools intersect -a {shiftBED} -b {mergePEAK} -wo > {tmp_file}"
    os.system(cmd)
    
    # Count reads per peak
    # A read is counted if it overlaps by >= 25bp
    peakcount = {}
    with open(tmp_file, 'r') as fi:
        for line in fi:
            items = line.rstrip('\n').split('\t')
            if len(items) < 10:
                continue
            # items[9] is the peak ID (4th column of peak file, 0-indexed as 9th in intersect output)
            # items[-1] is the overlap length
            peak_id = items[9]
            overlap = int(items[-1])
            # Count if overlap >= 25bp
            peakcount[peak_id] = peakcount.get(peak_id, 0) + (overlap >= 25)
    
    # Write output count table
    outPeakCount = os.path.join(outdir, f'{sn}.readcount')
    with open(outPeakCount, 'w') as fo:
        fo.write(f"PeakID\t{sn}\n")
        with open(mergePEAK) as fi:
            for line in fi:
                items = line.rstrip('\n').split('\t')
                if len(items) < 4:
                    continue
                peak_id = items[3]
                count = peakcount.get(peak_id, 0)
                fo.write(f"{peak_id}\t{count}\n")
    
    # Cleanup temp file
    os.system(f'rm -f {tmp_file}')
    
    print(f"Count table written to: {outPeakCount}")
    print(f"Total peaks: {sum(1 for _ in open(mergePEAK)) if os.path.exists(mergePEAK) else 0}")
    print(f"Peaks with reads: {len(peakcount)}")

if __name__ == "__main__":
    main()

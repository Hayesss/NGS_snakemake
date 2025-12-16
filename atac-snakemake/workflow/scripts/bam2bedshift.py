#!/usr/bin/env python
# coding: utf-8
"""
BAM to BED with Tn5 Shift Script for ATAC-seq Pipeline
=======================================================
This script performs +4/-5 bp shift to correct the Tn5-dimer binding sites.
See: https://www.biostars.org/p/428577/ for details about why this is needed.

The Tn5 transposase inserts adapters with a 9-bp stagger. To find the exact
location of Tn5 binding:
  - Forward strand reads: shift +4 bp
  - Reverse strand reads: shift -5 bp

Usage:
    python bam2bedshift.py <input.bam> <output.bed> [extend]

Arguments:
    input.bam:  Input BAM file (sorted, deduplicated)
    output.bed: Output BED file with shifted positions
    extend:     Extend N bp upstream and downstream relative to cut site
                (default: 25, creates 50bp windows)

Reference:
    https://github.com/QuKunLab/ATAC-pipe
"""

import sys
import os

def main():
    if len(sys.argv) < 3:
        print("Usage: python bam2bedshift.py <input.bam> <output.bed> [extend]")
        sys.exit(1)
    
    inBAM = sys.argv[1]
    outBED = sys.argv[2]
    EXTEND = 25 if len(sys.argv) <= 3 else int(sys.argv[3])
    
    print(f"Input BAM: {inBAM}")
    print(f"Output BED: {outBED}")
    print(f"Extend: {EXTEND} bp")
    
    bed = open(outBED, 'w')
    sam = os.popen(f"samtools view {inBAM}")  # no SAM header
    
    count = 0
    for line in sam:
        items = line.rstrip('\n').split('\t')
        # Reference:
        # 1) https://broadinstitute.github.io/picard/explain-flags.html
        # 2) https://ppotato.files.wordpress.com/2010/08/slide1.png?w=627
        # Sum of flags:
        #  99: read paired + proper pair + mate reverse strand + 1st in pair
        # 163: read paired + proper pair + mate reverse strand + 2nd in pair
        flags = int(items[1])
        chrom = items[2]
        if not chrom.startswith('chr'):
            continue
        
        # items[8]: insert fragment size (TLEN)
        insert_size = abs(int(items[8]))
        
        if flags == 99 or flags == 163:
            start = int(items[3])      # current reads (+)
            end = start + insert_size  # the paired reads (-)
            
            # Apply +4/-5 shift and extend
            start_extl = start + 4 - int(EXTEND)
            start_extr = start + 4 + int(EXTEND) if int(EXTEND) != 0 else start + 4 + 1
            end_extl = end - 5 - int(EXTEND)
            end_extr = end - 5 + int(EXTEND) if int(EXTEND) != 0 else end - 5 + 1
            
            # Ensure no negative coordinates
            if start_extl < 0:
                start_extl = 0
            if end_extl < 0:
                end_extl = 0
            
            # Write results
            read_name = items[0]
            if flags == 99:
                bed.write('\t'.join([chrom, str(start_extl), str(start_extr), 
                                    read_name, str(insert_size - 9), '+']) + '\n')
                bed.write('\t'.join([chrom, str(end_extl), str(end_extr), 
                                    read_name, str(insert_size - 9), '+']) + '\n')
            if flags == 163:
                bed.write('\t'.join([chrom, str(start_extl), str(start_extr), 
                                    read_name, str(insert_size - 9), '-']) + '\n')
                bed.write('\t'.join([chrom, str(end_extl), str(end_extr), 
                                    read_name, str(insert_size - 9), '-']) + '\n')
            count += 1
    
    bed.close()
    print(f"Processed {count} read pairs")

if __name__ == "__main__":
    main()

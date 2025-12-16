#!/usr/bin/env python
# coding: utf-8
"""
Generate Nucleosome Track Script for ATAC-seq Pipeline
=======================================================
This script separates ATAC-seq fragments into nucleosome-free and 
nucleosome-containing regions based on fragment length.

Fragment length interpretation:
  - < 120 bp: Nucleosome-free regions (NFR)
  - 180-247 bp: Mononucleosome
  - 315-473 bp: Dinucleosome

Usage:
    python generate_nucleosome_track.py <shift.bed> <nf.bed> <neu.bed> <chrom_sizes>

Arguments:
    shift.bed:   Input shift BED file
    nf.bed:      Output nucleosome-free BED file
    neu.bed:     Output nucleosome BED file (mono + di)
    chrom_sizes: Chromosome sizes file

Output:
    Creates BED files for NFR and nucleosome regions, plus BigWig tracks
"""

import sys
import os

def main():
    if len(sys.argv) < 5:
        print("Usage: python generate_nucleosome_track.py <shift.bed> <nf.bed> <neu.bed> <chrom_sizes>")
        sys.exit(1)
    
    inBED, nfBED, neuBED, genome_sizes = sys.argv[1:5]
    
    print(f"Input: {inBED}")
    print(f"NFR output: {nfBED}")
    print(f"Nucleosome output: {neuBED}")
    
    fo1 = open(nfBED, 'w')
    fo2 = open(neuBED, 'w')
    processed = {}
    
    nfr_count = 0
    mono_count = 0
    di_count = 0
    
    with open(inBED, 'r') as fi:
        for line in fi:
            items = line.rstrip().split()
            if len(items) < 6:
                continue
            chrom, start, end, reads_name, frag_len, strand = items[:6]
            start, end, frag_len = int(start), int(end), int(frag_len)
            
            # Skip if already processed (each fragment produces 2 lines)
            if reads_name in processed:
                continue
            else:
                processed[reads_name] = 1
            
            tn5bs_left = int((start + end) / 2)
            tn5bs_right = int(tn5bs_left + frag_len)
            
            # Nucleosome-free regions (< 120 bp)
            if frag_len < 120:
                tn5bs_center = int((tn5bs_left + tn5bs_right) / 2)
                tn5bs_left = max(0, tn5bs_center - 25)
                tn5bs_right = tn5bs_center + 25
                results = '\t'.join([chrom, str(tn5bs_left), str(tn5bs_right), 
                                    reads_name, '0', strand])
                fo1.write(results + '\n')
                nfr_count += 1
            
            # Mononucleosome (180-247 bp)
            if 180 < frag_len < 247:
                tn5bs_center = int((tn5bs_left + tn5bs_right) / 2)
                tn5bs_left = max(0, tn5bs_center - 25)
                tn5bs_right = tn5bs_center + 25
                results = '\t'.join([chrom, str(tn5bs_left), str(tn5bs_right), 
                                    reads_name, '0', strand])
                fo2.write(results + '\n')
                mono_count += 1
            
            # Dinucleosome (315-473 bp)
            if 315 < frag_len < 473:
                tn5bs_center1 = int((tn5bs_left + tn5bs_right) / 3)
                tn5bs_center2 = int((tn5bs_left + tn5bs_right) / 3 * 2)
                tn5bs_left1 = max(0, tn5bs_center1 - 25)
                tn5bs_right1 = tn5bs_center1 + 25
                tn5bs_left2 = max(0, tn5bs_center2 - 25)
                tn5bs_right2 = tn5bs_center2 + 25
                results = '\t'.join([chrom, str(tn5bs_left1), str(tn5bs_right1), 
                                    reads_name, '0', strand])
                fo2.write(results + '\n')
                results = '\t'.join([chrom, str(tn5bs_left2), str(tn5bs_right2), 
                                    reads_name, '0', strand])
                fo2.write(results + '\n')
                di_count += 1
    
    fo1.close()
    fo2.close()
    
    print(f"NFR fragments: {nfr_count}")
    print(f"Mononucleosome fragments: {mono_count}")
    print(f"Dinucleosome fragments: {di_count}")
    
    # Generate BigWig files
    nfBedGraph = nfBED.replace('.bed', '.bedGraph')
    neuBedGraph = neuBED.replace('.bed', '.bedGraph')
    
    # Sort and create bedGraph
    print("Generating bedGraph files...")
    os.system(f'sort -k1,1V -k2,2n {nfBED} -o {nfBED}.sorted')
    os.system(f'sort -k1,1V -k2,2n {neuBED} -o {neuBED}.sorted')
    os.system(f'genomeCoverageBed -bg -split -i {nfBED}.sorted -g {genome_sizes} > {nfBedGraph}')
    os.system(f'genomeCoverageBed -bg -split -i {neuBED}.sorted -g {genome_sizes} > {neuBedGraph}')
    
    # Convert bedGraph to BigWig
    nfBigWig = nfBED.replace('.bed', '.bw')
    neuBigWig = neuBED.replace('.bed', '.bw')
    print("Generating BigWig files...")
    os.system(f'sort -k1,1 -k2,2n {nfBedGraph} -o {nfBedGraph}')
    os.system(f'sort -k1,1 -k2,2n {neuBedGraph} -o {neuBedGraph}')
    os.system(f'bedGraphToBigWig {nfBedGraph} {genome_sizes} {nfBigWig}')
    os.system(f'bedGraphToBigWig {neuBedGraph} {genome_sizes} {neuBigWig}')
    
    # Cleanup
    os.system(f'rm -f {nfBED}.sorted {neuBED}.sorted {nfBedGraph} {neuBedGraph}')
    
    print(f"BigWig files created: {nfBigWig}, {neuBigWig}")

if __name__ == "__main__":
    main()

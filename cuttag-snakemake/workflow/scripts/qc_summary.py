#!/usr/bin/env python3
"""
CUT&Tag QC Summary Script
Generates alignment summary statistics from Bowtie2 and Picard outputs.

Based on: https://yezhengstat.github.io/CUTTag_tutorial/
"""

import os
import re
import pandas as pd


def parse_bowtie2_summary(filepath):
    """Parse Bowtie2 alignment summary file."""
    stats = {
        'sequencing_depth': 0,
        'mapped_fragments': 0,
        'alignment_rate': 0.0
    }
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        line = line.strip()
        # Parse total reads
        if 'reads; of these:' in line:
            stats['sequencing_depth'] = int(line.split()[0])
        # Parse concordantly mapped exactly 1 time
        elif 'aligned concordantly exactly 1 time' in line:
            stats['mapped_unique'] = int(line.split()[0])
        # Parse concordantly mapped >1 times
        elif 'aligned concordantly >1 times' in line:
            stats['mapped_multi'] = int(line.split()[0])
        # Parse overall alignment rate
        elif 'overall alignment rate' in line:
            rate_str = line.split('%')[0]
            stats['alignment_rate'] = float(rate_str)
    
    # Calculate total mapped fragments
    stats['mapped_fragments'] = stats.get('mapped_unique', 0) + stats.get('mapped_multi', 0)
    
    return stats


def parse_picard_metrics(filepath):
    """Parse Picard MarkDuplicates metrics file."""
    stats = {
        'read_pairs_examined': 0,
        'duplication_rate': 0.0,
        'estimated_library_size': 0
    }
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Find the header line and data line
    header_idx = None
    for i, line in enumerate(lines):
        if line.startswith('LIBRARY'):
            header_idx = i
            break
    
    if header_idx is not None and header_idx + 1 < len(lines):
        headers = lines[header_idx].strip().split('\t')
        values = lines[header_idx + 1].strip().split('\t')
        
        for h, v in zip(headers, values):
            if h == 'READ_PAIRS_EXAMINED':
                stats['read_pairs_examined'] = int(v) if v else 0
            elif h == 'PERCENT_DUPLICATION':
                stats['duplication_rate'] = float(v) * 100 if v else 0.0
            elif h == 'ESTIMATED_LIBRARY_SIZE':
                stats['estimated_library_size'] = int(v) if v else 0
    
    return stats


def parse_spikein_depth(filepath):
    """Parse spike-in sequencing depth file."""
    with open(filepath, 'r') as f:
        depth = int(f.read().strip())
    return depth


def main():
    # Get parameters from snakemake
    samples = snakemake.params.samples
    histones = snakemake.params.histones
    spikein_enabled = snakemake.params.spikein_enabled
    
    results = []
    
    for i, sample in enumerate(samples):
        # Parse alignment summary
        align_file = f"results/alignment/sam/bowtie2_summary/{sample}_bowtie2.txt"
        align_stats = parse_bowtie2_summary(align_file)
        
        # Parse duplicate metrics
        dup_file = f"results/alignment/removeDuplicate/picard_summary/{sample}_picard.dupMark.txt"
        dup_stats = parse_picard_metrics(dup_file)
        
        # Initialize record
        record = {
            'Sample': sample,
            'Histone': histones[i],
            'SequencingDepth': align_stats['sequencing_depth'],
            'MappedFragments_hg38': align_stats['mapped_fragments'],
            'AlignmentRate_hg38': f"{align_stats['alignment_rate']:.2f}%",
            'DuplicationRate': f"{dup_stats['duplication_rate']:.2f}%",
            'EstimatedLibrarySize': dup_stats['estimated_library_size'],
            'UniqueFragments': int(align_stats['mapped_fragments'] * (1 - dup_stats['duplication_rate'] / 100))
        }
        
        # Parse spike-in if enabled
        if spikein_enabled:
            spikein_file = f"results/alignment/sam/bowtie2_summary/{sample}_bowtie2_spikeIn.txt"
            spikein_depth_file = f"results/alignment/sam/bowtie2_summary/{sample}_bowtie2_spikeIn.seqDepth"
            
            if os.path.exists(spikein_file):
                spikein_stats = parse_bowtie2_summary(spikein_file)
                record['MappedFragments_spikeIn'] = spikein_stats['mapped_fragments']
                record['AlignmentRate_spikeIn'] = f"{spikein_stats['alignment_rate']:.2f}%"
            
            if os.path.exists(spikein_depth_file):
                spikein_depth = parse_spikein_depth(spikein_depth_file)
                record['SpikeInDepth'] = spikein_depth
                record['ScalingFactor'] = 10000 / spikein_depth if spikein_depth > 0 else 0
        
        results.append(record)
    
    # Create DataFrame and save
    df = pd.DataFrame(results)
    df.to_csv(snakemake.output.summary, sep='\t', index=False)
    
    print(f"QC summary saved to: {snakemake.output.summary}")


if __name__ == "__main__":
    main()

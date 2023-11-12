'''
Purpose: Filter variants if they don't have the same reads aligned to them on both haplotypes
Author: Mira Mastoras, mmastora@ucsc.edu
Usage: python3 check_haplotype_read_alignments.py
'''

import argparse
import pysam
import numpy as np

def main():
    parser = argparse.ArgumentParser(
        description='Filter variants if they don\'t have the same reads aligned to them on both haplotypes')
    parser.add_argument('--hap1Bam', type=str,
                        help='(Str) haplotype 1 bam file')
    parser.add_argument('--hap2Bam', type=str,
                        help='haplotype 2 bam file')
    parser.add_argument('--blocks', type=str,
                        help='tab separated matched bed file containing projectable variants and their projections on the other haplotype')

    # Fetch the arguments
    args = parser.parse_args()
    hap1BamPath = args.hap1Bam
    hap2BamPath = args.hap2Bam
    blocksPath = args.blocks

    # read in bams
    hap1Bam=pysam.AlignmentFile(hap1BamPath, "rb")
    hap2Bam = pysam.AlignmentFile(hap2BamPath, "rb")

    # read in variants
    with open(blocksPath, newline='') as file:
        result_list = list(csv.reader(file))
    print(result_list)
    result_2D = numpy.array(result_list)

    print(result_2D)
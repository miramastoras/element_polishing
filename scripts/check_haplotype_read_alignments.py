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
    parser.add_argument('--vcf', type=str,
                        help='Original vcf file')

    # Fetch the arguments
    args = parser.parse_args()
    hap1BamPath = args.hap1Bam
    hap2BamPath = args.hap2Bam
    blocksPath = args.blocks
    vcfPath = args.vcf

    # read in bams
    hap1Bam=pysam.AlignmentFile(hap1BamPath, "rb")
    hap2Bam = pysam.AlignmentFile(hap2BamPath, "rb")

    with open(blocksPath) as blocksFile:
        lines = blocksFile.readlines()
        lines = [line.rstrip().split('\t') for line in lines]

    for line in lines:
        print(line)
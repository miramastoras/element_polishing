'''
Purpose: Filter variants if they don't have the same reads aligned to them on both haplotypes
Author: Mira Mastoras, mmastora@ucsc.edu
Usage: python3 check_haplotype_read_alignments.py
'''

import argparse
import pysam
import numpy as np

def get_read_names(contig, start,end, alignmentFile):
    '''
    :param contig: contig name
    :param start: start coordinate
    :param end: end coordinate
    :param alignmentFile: pysam alignment file object
    :return: read names aligned to the specifed location
    '''

    readNames=[]
    for read in alignmentFile.fetch(contig=contig,start=start,end=end):
        readNames.append(read.query_name)
    return readNames

def filter_variants(projectableBam,projectionBam,blocksPath,outVcfPath):
    '''
    :param bam
    :param blocksPath:
    :return:
    '''

    # read blocks into 2d list
    with open(blocksPath) as blocksFile:
        lines = blocksFile.readlines()
        blocksList = [line.rstrip().split('\t') for line in lines]

    # for every variant & its projection to the other haplotype
    for line in blocksList:
        projectableContig = line[0]
        projectableStart = int(line[1])
        projectableEnd= int(line[2])

        projectionContig = line[12]
        projectionStart = int(line[13])
        projectionEnd = int(line[14])

        # get list of reads aligned to variant and its projected location on other haplotype
        projectableReadNames=get_read_names(projectableContig,projectableStart,projectableEnd,projectableBam)
        projectionReadNames=get_read_names(projectionContig,projectionStart,projectionEnd,projectionBam)

        projectableReadNames.sort()
        projectionReadNames.sort()

        if projectionReadNames==projectableReadNames:
            with open(outVcfPath, 'a') as outVcf:
                vcfEntry=[line[0]] + line[2:11]
                print(*vcfEntry,sep="\t",file = outVcf)

def main():
    parser = argparse.ArgumentParser(
        description='Filter variants if they don\'t have the same reads aligned to them on both haplotypes')
    parser.add_argument('--hap1Bam', type=str,
                        help='(Str) haplotype 1 bam file')
    parser.add_argument('--hap2Bam', type=str,
                        help='haplotype 2 bam file')
    parser.add_argument('--hap1Blocks', type=str,
                        help='paste -d "\\t" <outputProjectable> <outputProjection> > <hap1Blocks file>')
    parser.add_argument('--hap2Blocks', type=str,
                        help='paste -d "\\t" <outputProjectable> <outputProjection> > <hap2Blocks file>')
    parser.add_argument('--inVcf', type=str,
                        help='Original vcf file')
    parser.add_argument('--outVcf', type=str,
                        help='Filtered vcf file')

    # Fetch the arguments
    args = parser.parse_args()
    hap1BamPath = args.hap1Bam
    hap2BamPath = args.hap2Bam
    hap1BlocksPath = args.hap1Blocks
    hap2BlocksPath=args.hap2Blocks
    vcfPath = args.inVcf
    vcfOutPath = args.outVcf

    # print vcf header to output file
    outVcf = open(vcfOutPath, "a")
    with open(vcfPath) as inVcf:
        lines = inVcf.readlines()
        [outVcf.write(line) for line in lines if line.split("\t")[0][0]=="#" ]
    outVcf.close()

    # read in bams
    hap1Bam=pysam.AlignmentFile(hap1BamPath, "rb")
    hap2Bam = pysam.AlignmentFile(hap2BamPath, "rb")

    # filter hap1 variants
    filter_variants(hap1Bam,hap2Bam,hap1BlocksPath,vcfOutPath)
    # filter hap2 variants
    filter_variants(hap2Bam, hap1Bam, hap2BlocksPath, vcfOutPath)

if __name__ == '__main__':
    main()
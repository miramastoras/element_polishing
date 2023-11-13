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

        projectionContig = line[13]
        projectionStart = int(line[14])
        projectionEnd = int(line[15])

        # get list of reads aligned to variant and its projected location on other haplotype
        projectableReadNames=get_read_names(projectableContig,projectableStart,projectableEnd,projectableBam)
        projectionReadNames=get_read_names(projectionContig,projectionStart,projectionEnd,projectionBam)

        projectableReadNames.sort()
        projectionReadNames.sort()

        if projectionReadNames==projectableReadNames:
            with open(outVcfPath, 'a') as outVcf:
                vcfEntry=[line[0]] + [line[3]] + [line[4]] + [line[6]] + [line[7]] + [line[5]] + line[8:11] + [line[12]] + [line[11]]
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
    parser.add_argument('--hap1Vcf', type=str,
                        help='Original hap1 vcf file')
    parser.add_argument('--hap2Vcf', type=str,
                        help='Original hap1 vcf file')
    parser.add_argument('--outPrefix', type=str,
                        help='Filtered vcf file')

    # Fetch the arguments
    args = parser.parse_args()
    hap1BamPath = args.hap1Bam
    hap2BamPath = args.hap2Bam
    hap1BlocksPath = args.hap1Blocks
    hap2BlocksPath=args.hap2Blocks
    hap1VcfPath = args.hap1Vcf
    hap2VcfPath = args.hap2Vcf
    hap1OutVcfPath = args.outPrefix + ".hap1.vcf"
    hap2OutVcfPath = args.outPrefix + ".hap2.vcf"

    # print vcf header to output file
    outHap1Vcf = open(hap1OutVcfPath, "a")
    with open(hap1VcfPath) as inVcf:
        lines = inVcf.readlines()
        [outHap1Vcf.write(line) for line in lines if line.split("\t")[0][0]=="#" ]
    outHap1Vcf.close()

    outHap2Vcf = open(hap2OutVcfPath, "a")
    with open(hap2VcfPath) as inVcf:
        lines = inVcf.readlines()
        [outHap2Vcf.write(line) for line in lines if line.split("\t")[0][0] == "#"]
    outHap2Vcf.close()

    # read in bams
    hap1Bam=pysam.AlignmentFile(hap1BamPath, "rb")
    hap2Bam = pysam.AlignmentFile(hap2BamPath, "rb")

    # filter hap1 variants
    filter_variants(hap1Bam,hap2Bam,hap1BlocksPath,hap1OutVcfPath)
    # filter hap2 variants
    filter_variants(hap2Bam, hap1Bam, hap2BlocksPath, hap2OutVcfPath)

if __name__ == '__main__':
    main()
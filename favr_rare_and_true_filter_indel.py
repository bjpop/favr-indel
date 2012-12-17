#!/bin/env python

'''
Variant filter script.

Authors: Bernie Pope, Danny Park, Fabrice Odefrey, Tu Nguyen-Dumont.

Reads a list of variants from a TSV file for a sample and compares
them to the sequence reads in BAM files for other samples.

Revision history:

2 May 2011.    Initial incomplete version. Can read TSV file and BAMS and
               generate a table of which variants from the original sample
               are found in the other samples.

30 May 2011.   Filter out reads which are considered 'is_del'. These appear
               to come from deletions, and don't belong in the pileup for
               a given position.

31 May 2011.   Change the supported input format from a custom TSV file to
               a SIFT (csv) file.

25 Aug 2011.   Changed the binning rule to bin 35-length only reads, in
               addition to the previous binning rule.

29 Aug 2011.   Changed back to not bin 35-length only reads. This is now
               done in another tool.

19 Sep 2011.   Added classification for family samples.

20 Sep 2011.   Allowed the input to be tab separated, with coordinates comma separated.

14 Nov 2011.   Added classify arguments to the command line.

'''

import os
import pysam
import sys
import csv
import getopt
from favr_common_indel import (safeReadInt, getEvidence,
    makeSafeFilename, sortByCoord, parseVariantRowVCF, parseVariantRowAnnovar,
    initEvidence)
from favr_rare_and_true_classify_indel import classify
import argparse

parser = argparse.ArgumentParser(description='Filtering of variants based on their presence/absence and abundance in samples from comparators.')
parser.add_argument('bamFilenames',
                    nargs='+',
                    help='Comparator BAM files',
                    type=str)
parser.add_argument('--variants_annovar',
                    metavar='VARS',
                    type=str,
                    help='variant file in Annovar format')
parser.add_argument('--variants_vcf',
                    metavar='VARS',
                    type=str,
                    help='variant file in VCF format')
parser.add_argument('--bin',
                    metavar='BIN_FILE',
                    type=str,
                    required=True,
                    help='store binned variants in this file')
parser.add_argument('--keep',
                    metavar='KEEP_FILE',
                    type=str,
                    required=True,
                    help='store kept variants in this file')
parser.add_argument('--log',
                    metavar='LOG_FILE',
                    type=str,
                    required=True,
                    help='store log messages in this file')
parser.add_argument('--varLikePercent',
                    metavar='N',
                    type=int,
                    required=True,
                    help='percentage of reads for a sample same as the variant')
parser.add_argument('--samplesPercent',
                    metavar='N',
                    type=int,
                    required=True,
                    help='percent of total samples which pass the threshold')

def main():
    options = parser.parse_args()
    variants = []
    if options.variants_annovar and options.variants_vcf:
        exit('Please only specify variants in one format, either Annovar or VCF')
    elif options.variants_annovar:
        with open(options.variants_annovar, 'rU') as variantsFile:
            for row in csv.reader(variantsFile, delimiter=',', quotechar='"'):
                variant = parseVariantRowAnnovar(row)
                if variant:
                    variants.append(variant)
    elif options.variants_vcf:
        with open(options.variants_vcf, 'rU') as variantsFile:
            for row in csv.reader(variantsFile, delimiter='\t', quotechar='"'):
                variant = parseVariantRowVCF(row)
                if variant:
                    variants.append(variant)
    # compute the presence/absence of each variant in the bam files
    evidence = initEvidence(variants) 
    # update the evidence based on reads seen in the bam files
    getEvidence(evidence, variants, options.bamFilenames)
    # filter the variants
    filter(options, evidence)

def filter(options, evidence):
    '''Decide which variants to keep and which to bin.'''
    binFilename = options.bin
    keepFilename = options.keep
    logFilename = options.log
    with open (logFilename, 'w') as logFile:
        with open(binFilename, 'w') as binFile:
            with open(keepFilename, 'wb') as keepFile:
                csvWriter = csv.writer(keepFile, delimiter='\t', quotechar='|')
                # sort the variants by coordinate
                for key, info in sortByCoord(evidence):
                    chr, pos = key
                    pos += 1 # we use 0 based indexing internally and 1 based externally
                    classification = classify(options, info.counts)
                    # record the classification of this variant in the logfile
                    logFile.write("%s:%d: %s: %s\n" % (chr, pos, classification.action, classification.reason))
                    if classification.action == 'bin':
                        # bin the variant
                        binFile.write('%s:%d\n' % (chr, pos))
                        for readCount,depth in info.counts:
                            binFile.write('    <vars/coverage: %d/%d>\n' % (readCount,depth))
                    elif classification.action == 'keep':
                        # keep the variant
                        csvWriter.writerow(info.inputRow)

if __name__ == "__main__":
   main()
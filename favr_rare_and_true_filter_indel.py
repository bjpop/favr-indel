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
# XXX this windowvar is for both nearby clipping and nearby indels.
# should we have two options separately?
parser.add_argument('--windowvar',
                    metavar='N',
                    type=int,
                    default=0,
                    help='the window variance for finding clippings')
# XXX now allow only one feature at a time.
# should allow multiple simultaneous features.
parser.add_argument('--feature',
                    metavar='FEATURE',
                    type=str,
                    required=True,
                    help='the feature that would be used for filtering')
def main():
    options = parser.parse_args()
    variants = []
    if options.variants_annovar and options.variants_vcf:
        exit('Please only specify variants in one format, either Annovar or VCF')
    elif options.variants_annovar:
        variantsFilename = options.variants_annovar
        reader = lambda file: csv.reader(file, delimiter=',', quotechar='"')
        variantParser = parseVariantRowAnnovar 
    elif options.variants_vcf:
        variantsFilename = options.variants_vcf
        reader = lambda file: csv.reader(file, delimiter='\t', quotechar='"')
        variantParser = parseVariantRowVCF 
    else:
        exit('Please specify which variant format to use, either Annovar or VCF')
    with open(variantsFilename, 'rU') as variantsFile:
        for row in reader(variantsFile):
            variant = variantParser(row)
            if variant:
                variants.append(variant)
    # compute the presence/absence of each variant in the bam files
    evidence = initEvidence(variants) 
    # update the evidence based on reads seen in the bam files
    getEvidence(evidence, variants, options.bamFilenames, options.windowvar)
    # with respect to each type of variants (SVN, indel, indel with clipping)
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
                    classification = classify(options, info.features)
                    # record the classification of this variant in the logfile
                    logFile.write("%s:%d: %s: %s\n" % (chr, pos, classification.action, classification.reason))
                    if classification.action == 'bin':
                        # bin the variant
                        binFile.write('%s:%d\n' % (chr, pos))
                        # XXX should handle multiple features
                        # should comment out after finding better way
                        # to get counts for each feature
                        #for features in info.features:
                        #    binFile.write('    <vars/coverage: %d/%d>\n' % (readCount,depth))
                    elif classification.action == 'keep':
                        # keep the variant
                        csvWriter.writerow(info.inputRow)

if __name__ == "__main__":
   main()

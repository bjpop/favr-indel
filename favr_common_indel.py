'''
Common utilities for FAVR programs. Indel version.

Authors: Bernie Pope, Danny Park, Fabrice Odefrey, Tu Nguyen-Dumont.
'''

import os
import pysam
import sys

def safeReadInt(str):
    if str.isdigit():
        return int(str)
    else:
        raise Exception, 'not an integer: ' + str

# sort the variants by coordinate.
def sortByCoord(evidence):
    return sorted(evidence.items(), cmp=compareCoord)

# compare two chromosome coordinates for ordering.
def compareCoord(coord1, coord2):
    chr1,pos1 = coord1[0].split(':')
    chr2,pos2 = coord2[0].split(':')
    if chr1 == chr2:
        return cmp(int(pos1), int(pos2))
    else:
        return compareChrCode(chr1[3:], chr2[3:])

def compareChrCode(code1, code2):
    isNumerical1 = code1.isdigit()
    isNumerical2 = code2.isdigit()
    if isNumerical1:
        if isNumerical2:
            return cmp(int(code1), int(code2))
        else:
            return -1
    else:
        return cmp(code1, code2)

def getEvidence(variantList, bamFilenames):
    evidence, variants = initEvidence(variantList)
    # Iterate over sample BAM files.
    for bamFile in bamFilenames:
       with pysam.Samfile(bamFile, "rb") as bam:
           # Count how many samples have each particular variant.
           countVariants(evidence, variants, bam)
    return evidence, variants

class EvidenceInfo(object):
    def __init__(self, inputRow, counts):
        self.inputRow = inputRow
        self.counts = counts

def initEvidence(variantList):
    '''Parse the variant list, and initialised the evidence dictionary'''
    evidence = {}
    outputVariants = []
    for variant in variantList:
        info = parseVariantRow(variant)
        if info:
            evidence[info.id] = EvidenceInfo(inputRow = info.inputRow, counts = [])
            outputVariants.append(info)
    return evidence, outputVariants

def showEvidence(evidence):
    '''Print out the frequency counter for each variant.'''
    for chrPos,info in evidence.items():
        print("%s %s" % (info.inputRow,str(info.counts)))

def countVariants(evidence, variants, bam):
    '''For each variant in the list, check if it is evident in this particular sample BAM.'''
    for variant in variants:
        reads = lookupPileup(bam, variant.chromosome, variant.startPosition, variant.endPosition)
        sameAsVariant = 0
        coverage = len(reads)
        for read in reads:
            if variant.inAlignment(read):
                sameAsVariant += 1
        evidence[variant.id].counts.append((sameAsVariant,coverage))

class Variant(object):
    def __init__(self, id, chromosome, startPosition, endPosition, inputRow):
        self.id = id
        self.chromosome = chromosome
        self.startPosition = startPosition
        self.endPosition = endPosition
        self.inputRow = inputRow

class SNV(Variant):
    def __init__(self, id, chromosome, startPosition, inputRow, refBase, variantBase):
        super(SNV, self).__init__(id, chromosome, startPosition, startPosition, inputRow)
        self.refBase = refBase
        self.variantBase = variantBase

    def __str__(self):
        return "SNV({},{},{},{})".format(self.chromosome,
                self.startPosition, self.refBase, self.variantBase)

    # XXX untested
    # XXX fix bug with cigars past the position
    def inAlignment(self, alignment):
        offset = (self.startPosition - 1) - alignment.pos
        if offset < 0:
            return False
        cigar_coords = cigarToCoords(alignment)
        for cigar in cigar_coords:
            if cigar.code == 'I':
                offset += cigar.length
            elif cigar.code == 'D':
                offset -= cigar.length
        if offset > 0 and offset < alignment.qlen:
            return alignment[offset] == self.variantBase
        else:
            return False

class Deletion(Variant):
    def __init__(self, id, chromosome, startPosition, endPosition, inputRow):
        super(Deletion, self).__init__(id, chromosome, startPosition, endPosition, inputRow)

    def __str__(self):
        return "Deletion({},{},{})".format(self.chromosome,
                self.startPosition, self.endPosition)

    def inAlignment(self, alignment):
        cigar_coords = cigarToCoords(alignment)
        zeroStart = self.startPosition - 1
        zeroEnd = self.endPosition - 1
        for cigar in cigar_coords:
            if (cigar.code == 'D' and
                cigar.start == zeroStart and
                cigar.start + cigar.length - 1 == zeroEnd):
                return True
        return False

class Insertion(Variant):
    def __init__(self, id, chromosome, startPosition, inputRow, insertedBases):
        super(Insertion, self).__init__(id, chromosome, startPosition, startPosition, inputRow)
        self.insertedBases = insertedBases 

    def __str__(self):
        return "Insertion({},{},{},{})".format(self.chromosome,
                self.startPosition, self.endPosition, self.insertedBases)

    def inAlignment(self, alignment):
        cigar_coords = cigarToCoords(alignment)
        zeroStart = self.startPosition - 1
        zeroEnd = self.endPosition - 1
        for cigar in cigar_coords:
            # XXX should check that the inserted bases are the same
            if cigar.code == 'I' and cigar.start == zeroStart:
                return True
        return False

class CigarCoord(object):
    def __init__(self, code, start, length):
        self.code = code
        self.start = start # zero based
        self.length = length

def cigarToCoords(alignment):
    # pos is the reference position of the first aligned base in the read
    pos = alignment.pos # zero based
    cigar = alignment.cigar
    # the portion of the read which aligned to the reference, does not include
    # soft clipped bases on the flanks
    alignedBases = alignment.query
    if cigar == None:
        return [] 
    result = []
    for code,length in cigar:
        if code == 0:
            # Match
            result.append(CigarCoord('M', pos, length))
            pos += length
        elif code == 1:
            # Insert.
            # The read has bases which do not appear in the reference.
            # position of insert appears to be the position of the
            # base just prior to the insert.
            # pos should already point to the next position in
            # the reference so we don't increment it
            # We assume insertions cannot appear before the first matched base
            # (they would be in the soft-clipped region).
            # XXX should calculate the inserted bases and check them against
            # the variant
            previous_pos = pos - 1
            result.append(CigarCoord('I', previous_pos, length))
        elif code == 2:
            # Deletion. 
            # The read is missing bases which appear in the reference.
            #result.append(('D', pos, pos + length - 1))
            result.append(CigarCoord('D', pos, length))
            pos += length
        else:
            # 'N', 4 : 'S', 5 : 'H', 6 : 'P', 7 : '=', 8 : 'X' 
            # we ignore the other types of cigar codes for the moment
            # can can safely ignore soft and hard clipping because they
            # only appear at the ends of the read outside the mapped
            # section. We should probably handle N, which is a skipped
            # base. N is like a deletion. P is for padding which is
            # needed in multiple alignments but may not be relevent to us.
            pass
    return result


def parseVariantRow(row):
    '''Extract the interesting information about a variant from a SIFT TSV row.'''
    if len(row) >= 25 and row[0] != 'Func':
        chr, startPos, endPos, ref, var = row[21:26]
        if var == '-':
            # The variant is missing so it is a deletion.
            return Deletion(id = "%s:%s" % (chr,startPos), 
                            chromosome = chr,
                            startPosition = int(startPos),      
                            endPosition = int(endPos),
                            inputRow = row)
        elif ref == '-':
            # The reference is missing to it is an insertion.
            return Insertion(id = "%s:%s" % (chr,startPos), 
                             chromosome = chr,
                             startPosition = int(startPos),      
                             inputRow = row,
                             insertedBases = var)
        else:
            # Variant and reference are present, assume an SNV
            return SNV(id = "%s:%s" % (chr,startPos), 
                       chromosome = chr,
                       startPosition = int(startPos),
                       inputRow = row,
                       refBase = ref, 
                       variantBase = var)

    return None

# XXX fixme name
def lookupPileup(bam, chr, startCol, endCol):
    '''Retrieve the pileup for a particular chromosome:position coordinate.'''
    # assume argument is in 1-based numbering.
    # samtools gives back a range of pilups that cover the requested coordinates,
    # (no idea why) so we have to search for the particular one we want

    # XXX what's the best way to avoid duplicates?
    result = []
    read_ids = set()
    for read in bam.fetch(chr, startCol-1, endCol):
        if read.qname not in read_ids: 
            result.append(read)
            read_ids.add(read.qname)
    return result

def makeSafeFilename(name):
    if not os.path.exists(name):
        return name
    else:
        count = 2
        while (count < sys.maxint):
            newName = name + str(count)
            if not os.path.exists(newName):
                return newName
            count += 1
    # a safety condition in case we can't return a safe name
    return None

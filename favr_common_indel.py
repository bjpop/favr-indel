'''
Common utilities for FAVR programs. Indel version.

Authors: Bernie Pope, Danny Park, Fabrice Odefrey, Tu Nguyen-Dumont.

# We store all reference coordinates in zero-based indices, which
# is the way they are stored in BAM files, but not usually the way they
# are used in Variant Call formats, such as annovar and VCF, which are
# one-based.
# However, we always display coordinates in one-based format, to be
# consistent with external tools.
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
    #chr1,pos1 = coord1[0].split(':')
    #chr2,pos2 = coord2[0].split(':')
    chr1,pos1 = coord1[0]
    chr2,pos2 = coord2[0]
    if chr1 == chr2:
        #return cmp(int(pos1), int(pos2))
        return cmp(pos1, pos2)
    else:
        return compareChrCode(chr1, chr2)

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

def getEvidence(evidence, variants, bamFilenames, windowVar):
    #evidence, variants = initEvidence(variantList)
    # Iterate over sample BAM files.
    for bamFile in bamFilenames:
       with pysam.Samfile(bamFile, "rb") as bam:
           # Count how many samples have each particular variant.
           countVariants(evidence, variants, bam, windowVar)
    #return evidence 

# counts various kinds of signals we are looking for in the data
class Counter(object):
    def __init__(self):
        self.snv_exact = 0
        self.indel_exact = 0
        self.indel_with_clipping = 0

    def __str__(self):
        return "(snv_exact = {}, indel_exact = {}, indel_with_clipping = {})" \
               .format(self.snv_exact, self.indel_exact, self.indel_with_clipping)

class EvidenceInfo(object):
    def __init__(self, inputRow, counts):
        # original row from the input file
        self.inputRow = inputRow 
        # list of pairs: (reads same as variant, coverage)
        self.counts = counts
    def __str__(self):
        # XXX inputRow has too much information to just see the counts
        return "counts: {}, coverage:{} " \
                .format(' '.join([str(count) for count, coverage in self.counts]),
                        ' '.join([str(coverage) for count, coverage in self.counts]))

def initEvidence(variants):
    evidence = {}
    for variant in variants:
        evidence[variant.id] = EvidenceInfo(inputRow = variant.inputRow, counts = [])
    return evidence

def showEvidence(evidence):
    '''Print out the frequency counter for each variant.'''
    for chrPos,info in evidence.items():
        print("%s %s" % (info.inputRow, str(info.counts)))

def countVariants(evidence, variants, bam, windowVar):
    '''For each variant in the list, check if it is evident in this particular sample BAM.'''
    for variant in variants:
        # XXX need logger?
        print '='*30
        print 'variants start %d end %d (one-based)' % \
               (variant.startPosition+1, variant.endPosition+1)

        reads = lookupReads(bam, variant.chromosome, variant.startPosition, variant.endPosition)
        # sameAsVariant = 0
        counter = Counter() # initialise zero counter
        coverage = len(reads)
        for read in reads:
            # update the counter
            variant.inAlignment(read, counter)
        
        # count nearby clippings if variant is insertion or deletion
        if isinstance(variant, Insertion) or isinstance(variant, Deletion):
            # widen the window to get reads which
            # have clipping at variant position,
            # because pysam.fetch does not return the read that
            # has clipping on the query range.
            # look up reads that intersect window variance
            clipping_var_start = variant.startPosition - windowVar
            clipping_var_end = variant.endPosition + windowVar
            reads_clipping = lookupReads(bam, variant.chromosome,
                                         clipping_var_start, clipping_var_end)
            # XXX need to consider two different coverage,
            # one is for exact match, the other is for nearby clipping.
            coverage_clipping = len(reads_clipping)
            # while getting clipping for each reads,
            # increase counter by 1 if the clipping is on the variant position.
            for clipping in get_clipping_in_alignments(reads_clipping):
                if clipping.isNearby(variant):
                    counter.indel_with_clipping += 1
                    # XXX need logger?
                    print clipping, 'is on', variant

        evidence[variant.id].counts.append((counter, coverage))

class Clipping(object):
    def __init__(self, startPosition, endPosition):
        self.startPosition = startPosition
        self.endPosition = endPosition

    # Pretty print, showing start position in one-based form
    def __str__(self):
        return "Clipping({},{})".format(self.startPosition + 1,
                self.endPosition + 1)
        
    def isNearby(self, variant):
        '''Return True if this clipping intersects with or
           is within the variant.'''
        clipping_start = self.startPosition
        clipping_end = self.endPosition
        variant_start = variant.startPosition
        variant_end = variant.endPosition
        # intersect on the left
        if (clipping_start <= variant_start and
            clipping_end >= variant_start):
            return True
        # within variant
        elif (clipping_start > variant_start and
              clipping_end <= variant_end):
            return True
        # intersect on the right
        elif (clipping_start <= variant_end and
              clipping_end >= variant_end):
            return True
        return False

def get_clipping_in_alignments(alignments):
    for alignment in alignments:
        cigar_coords = cigarToCoords(alignment)
        for cigar in cigar_coords:
            if cigar.code == 'S' or cigar.code == 'H':
                yield Clipping(cigar.start, cigar.start + cigar.length - 1)

class Variant(object):
    def __init__(self, id, chromosome, startPosition, endPosition, inputRow):
        self.id = id # (chr,startpos)
        self.chromosome = chromosome
        self.startPosition = startPosition # zero based
        self.endPosition = endPosition     # zero based
        self.inputRow = inputRow

class SNV(Variant):
    def __init__(self, id, chromosome, startPosition, inputRow, refBase, variantBase):
        super(SNV, self).__init__(id, chromosome, startPosition, startPosition, inputRow)
        self.refBase = refBase
        self.variantBase = variantBase

    # Pretty print, showing start position in one-based form
    def __str__(self):
        return "SNV({},{},{},{})".format(self.chromosome,
                self.startPosition + 1, self.refBase, self.variantBase)

    # XXX untested
    def inAlignment(self, alignment, counter):
        target = self.startPosition - alignment.pos
        offset = 0
        for cigar in cigarToCoords(alignment):
            if target < 0:
                #return False
                # don't update the counter
                return
            # An SNV can only appear inside a Matched sequence
            elif cigar.code == 'M':
                if 0 <= target < cigar.length:
                    #return self.variantBase == alignment.query[offset + target]
                    if self.variantBase == alignment.query[offset + target]:
                        counter.snv_exact += 1
                        return
                else:
                    offset += cigar.length
                    target -= cigar.length
            elif cigar.code == 'I':
                offset += cigar.length
            elif cigar.code == 'D':
                target -= cigar.length
        #return False
        # don't update the counter
        return

class Deletion(Variant):
    def __init__(self, id, chromosome, startPosition, endPosition, inputRow):
        super(Deletion, self).__init__(id, chromosome, startPosition, endPosition, inputRow)

    # Pretty print, showing start and end position in one-based form
    def __str__(self):
        return "Deletion({},{},{})".format(self.chromosome,
                self.startPosition + 1, self.endPosition + 1)

    def inAlignment(self, alignment, counter):
        cigar_coords = cigarToCoords(alignment)
        for cigar in cigar_coords:
            if (cigar.code == 'D' and
                cigar.start == self.startPosition and
                cigar.start + cigar.length - 1 == self.endPosition):
                counter.indel_exact += 1
        # don't update the counter
        return

class Insertion(Variant):
    def __init__(self, id, chromosome, startPosition, inputRow, insertedBases):
        super(Insertion, self).__init__(id, chromosome, startPosition, startPosition, inputRow)
        self.insertedBases = insertedBases 

    # Pretty print, showing start and end position in one-based form
    def __str__(self):
        return "Insertion({},{},{},{})".format(self.chromosome,
                self.startPosition + 1, self.endPosition + 1, self.insertedBases)

    def inAlignment(self, alignment, counter):
        cigar_coords = cigarToCoords(alignment)
        for cigar in cigar_coords:
            # XXX should check that the inserted bases are the same
            if cigar.code == 'I' and cigar.start == self.startPosition:
                #return True
                counter.indel_exact += 1
        #return False
        # don't update the counter
        return

class CigarCoord(object):
    def __init__(self, code, start, length):
        # 'M' = match, 'I' = insert, 'D' = delete, 'S' = soft clipped
        # 'H' = hard clipped, 'P' = padded, 'N' = non-match, 'X' = ???
        self.code = code
        self.start = start # zero based
        self.length = length
    def __str__(self):
        return '(code %s start %d length %d)' % (self.code, self.start, self.length)

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
    for n, (code,length) in enumerate(cigar):
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
        elif code == 4:
            # Soft clipping
            # The sequence of read is not aligned from the first residue or
            # to the last residue.
            # If the clipped alignment is on the first residue,
            # the position should be alignment.pos - length,
            # since the alignment.pos starts from the first aligned base.
            # If the clipped alignmet is on the end residue,
            # the position should be pos.
            clipping_pos = 0
            if n == 0:  # the first residue
                clipping_pos = pos - length
            else:
                clipping_pos = pos
            result.append(CigarCoord('S', clipping_pos, length))
        elif code == 5:
            # Hard clipping
            # Similar to Soft clipping
            # The difference is that the hard clipped subsequence is not
            # present in the alignment record.
            # This means that there is no difference in handling
            # Hard Clipping for us.
            clipping_pos = 0
            if n == 0:  # the first residue
                clipping_pos = pos - length
            else:
                clipping_pos = pos
            result.append(CigarCoord('H', clipping_pos, length))
        else:
            # 3 : 'N', 6 : 'P', 7 : '=', 8 : 'X' 
            # we ignore the other types of cigar codes for the moment
            # We should probably handle N, which is a skipped
            # base. N is like a deletion. P is for padding which is
            # needed in multiple alignments but may not be relevent to us.
            pass
    return result

def parseVariantRowVCF(row):
    '''Extract the interesting information about a variant from VCF input file.'''
    if len(row) >= 5 and row[0].startswith('chr'):
        # starPos and endPos are input in one-based coordinates, but
        # we store them in zero-based coordinates
        chr, startPos, id, ref, var = row[:5]
        startPosZero = int(startPos) - 1
        #varId = "%s:%s" % (chr,str(startPosZero))
        varId = (chr, startPosZero)
        # XXX need to handle indels
        #if var == '-':
        #    return Deletion(id = varId, 
        #                    chromosome = chr,
        #                    startPosition = startPosZero,      
        #                    endPosition = endPosZero,
        #                    inputRow = row)
        #elif ref == '-':
        #    return Insertion(id = varId, 
        #                     chromosome = chr,
        #                     startPosition = startPosZero, 
        #                     inputRow = row,
        #                     insertedBases = var)
        #else:
            # Variant and reference are present, assume an SNV
        return SNV(id = varId, 
                   chromosome = chr,
                   startPosition = startPosZero,
                   inputRow = row,
                   refBase = ref, 
                   variantBase = var)
    return None

def parseVariantRowAnnovar(row):
    '''Extract the interesting information about a variant from annovar input file.'''
    if len(row) >= 25 and row[0] != 'Func':
        # starPos and endPos are input in one-based coordinates, but
        # we store them in zero-based coordinates
        chr, startPos, endPos, ref, var = row[21:26]
        startPosZero = int(startPos) - 1
        endPosZero = int(endPos) - 1
        #varId = "%s:%s" % (chr,str(startPosZero))
        varId = (chr, startPosZero)
        if var == '-':
            # The variant is missing; it is a deletion.
            return Deletion(id = varId, 
                            chromosome = chr,
                            startPosition = startPosZero,      
                            endPosition = endPosZero,
                            inputRow = row)
        elif ref == '-':
            # The reference is missing; it is an insertion.
            return Insertion(id = varId, 
                             chromosome = chr,
                             startPosition = startPosZero, 
                             inputRow = row,
                             insertedBases = var)
        else:
            # Variant and reference are present, assume an SNV
            return SNV(id = varId, 
                       chromosome = chr,
                       startPosition = startPosZero,
                       inputRow = row,
                       refBase = ref, 
                       variantBase = var)
    return None

def lookupReads(bam, chr, startCol, endCol):
    '''Retrieve the pileup for a particular chromosome:position coordinate.'''
    # arguments are in zero-based indices
    # samtools gives back a range of pilups that cover the requested coordinates,
    # (no idea why) so we have to search for the particular one we want

    # XXX what's the best way to avoid duplicates?
    result = []
    read_ids = set()
    for read in bam.fetch(chr, startCol, endCol+1):
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

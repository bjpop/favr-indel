'''
Classifier for the Variant filter script.

Authors: Bernie Pope, Danny Park, Fabrice Odefrey, Tu Nguyen-Dumont.
'''

# The classify function decides if a particular variant from the
# variant list should be kept or binned (discarded). It returns
# a classification of the variant which can be used to
# determine if it should be kept or discarded and a reason
# why.
#
# The decision to bin or keep is based on how many times
# we see the same variant in other samples.
#
# If you want to use a different binning criteria, then
# rewrite the classify function. The one provided below is just
# an example.

class Classify(object):
    def __init__(self, action, reason):
        self.action = action # 'bin' or 'keep'
        self.reason = reason # some text explaining why

# the argument variantInfo is a list of pairs:
#    (readCount, depth)
# one per sample, where
#    readCount is the number of reads in a sample which match the variant (at the same position)
#    depth is the number of reads in a sample covering the same position as the variant (aka coverage)

# XXX handle different kinds of counters
def classify(options, variantInfo):
    #varLikeThresh = options.varLikeThresh
    varLikePercent = options.varLikePercent
    samplesPercent = options.samplesPercent
    totalSamples = 0     # number of sample files
    binableSamples = 0   # number of samples which are considered binable

    # XXX should consider multiple features
    feature = options.feature

    # variantInfo is the list of Features
    for features in variantInfo:
        totalSamples += 1
        readCount, depth = getSingleCount(features, feature)
        if (depth > 0) and percentage(readCount, depth) >= varLikePercent:
            binableSamples += 1
        
    if totalSamples > 0:
        #if (binableSamples * 100 / totalSamples) >= samplesPercent:
        #if (float(binableSamples) / totalSamples * 100.0) >= samplesPercent:
        if percentage(binableSamples, totalSamples) >= samplesPercent:
            return Classify('bin', binThresholdMessage % (binableSamples, totalSamples, samplesPercent))
        else:
            return Classify('keep', keepMessage % (binableSamples, totalSamples, samplesPercent))
    else:
        return Classify('bin', binZeroSamplesMessage)

# assumes denominator is not zero
def percentage(numerator, denominator):
    return (float(numerator) / denominator) * 100.0

def getSingleCount(features, feature):
    # deal with single feature
    if feature == 'snv_exact':
        return (features.snv_exact.count, features.snv_exact.depth)
    if feature == 'indel_exact':
        return (features.indel_exact.count, features.indel_exact.depth)
    if feature == 'indel_similar':
        return (features.indel_similar.count, features.indel_similar.depth)
    if feature == 'indel_nearby':
        return (features.indel_nearby.count, features.indel_nearby.depth)
    if feature == 'indel_with_clipping':
        return (features.indel_with_clipping.count, features.indel_with_clipping.depth)


binThresholdMessage = '(binableSamples(=%d) * 100 / totalSamples(=%d)) >= samplesPercent(=%d)'
binZeroSamplesMessage = 'there were zero samples to compare with'
keepMessage = '(binableSamples(=%d) * 100 / totalSamples(=%d)) < samplesPercent(=%d)'

#!/usr/bin/env python


# CLASS annotation
#######################################################################################################################

# Getting CLASS annotation of a given variant
def getClassAnnotation(variant, transcript, protein, mutprotein, loc, ssrange):
    # Variants in UTR
    chkutr = checkUTR(transcript, variant)
    if chkutr == 'UTR5':
        if 'Ex' in loc and (variant.isDeletion() or variant.isComplex()): return 'IM'
        return '5PU'
    elif chkutr == 'UTR3':
        if 'Ex' in loc and (variant.isDeletion() or variant.isComplex()): return 'SL'
        return '3PU'

    # Variants crossing exon-intron boundaries
    if '-' in loc: return 'ESS'

    # Intronic, splice site and essential splice site variants
    if 'In' in loc:
        if transcript.isInEssentialSpliceSite(variant): return 'ESS'
        if transcript.intronLength(int(loc[loc.find('/') + 1:])) >= 9:
            if transcript.isIn_SS5_Site(variant): return 'SS5'
        if transcript.isInSplicingRegion(variant, ssrange): return 'SS'
        return 'INT'

    # Checking if variant affect the first or last three bases of an exon
    potSS = transcript.isInFirstOrLast3BaseOfExon(variant)

    protL = len(protein)
    mutprotL = len(mutprotein)

    # Frame-shift coding variants
    if not variant.isInFrame(): return 'FS'

    # Synonymous coding variants
    if protein == mutprotein:
        if potSS: return 'EE'
        return 'SY'

    # Variants affecting the initiation amino acid
    if protein[0] != mutprotein[0]: return 'IM'

    # Stop gain and stop lost variants
    while len(protein) > 0 and len(mutprotein) > 0:
        if protein[0] == mutprotein[0]:
            protein = protein[1:]
            mutprotein = mutprotein[1:]
        else:
            break

    if protein == '': return '3PU'

    if protein[0] == 'X' and len(mutprotein) == 0: return 'SL'
    if protein[0] == 'X' and mutprotein[0] != 'X': return 'SL'

    while len(protein) > 0 and len(mutprotein) > 0:
        if protein[-1] == mutprotein[-1]:
            protein = protein[:-1]
            mutprotein = mutprotein[:-1]
        else:
            break

    if 'X' in mutprotein: return 'SG'

    # Non-synonymous coding variants
    if protL == mutprotL:
        if potSS: return 'EE'
        return 'NSY'

    # In-frame variants
    if potSS: return 'EE'
    return 'IF'

#######################################################################################################################



# Sequence Ontology (SO)
#######################################################################################################################


def getSequenceOntologyAnnotation(variant, transcript, protein, mutprotein, loc):

    # Variants in UTR
    chkutr = checkUTR(transcript, variant)
    if chkutr == 'UTR5':
        if 'Ex' in loc and (variant.isDeletion() or variant.isComplex()): return 'initiator_codon_variant'
        return '5_prime_UTR_variant'
    elif chkutr == 'UTR3':
        if 'Ex' in loc and (variant.isDeletion() or variant.isComplex()): return 'stop_lost'
        return '3_prime_UTR_variant'

    where = transcript.whereIsThisVariant(variant)

    if '-' in where:
        first = where[:where.find('-')]
        if first.startswith('In'): return 'splice_acceptor_variant'
        if first.startswith('Ex'): return 'splice_donor_variant'
        if first.startswith('fsIn'): return 'splice_acceptor_variant'

    if 'In' in where:

        if isInSpliceDonor(transcript, variant): return 'splice_donor_variant'
        if isInSpliceAcceptor(transcript, variant): return 'splice_acceptor_variant'

        if transcript.intronLength(int(where[where.find('/') + 1:])) >= 9:
            if transcript.isIn_SS5_Site(variant): return 'splice_donor_5th_base_variant'

    out = []

    if isInSplicingRegion(transcript, variant):
        if where.startswith('In') or where.startswith('fsIn'):
            return 'intron_variant|splice_region_variant'
        else:
            out.append('splice_region_variant')

    if where.startswith('In') or where.startswith('fsIn'): return 'intron_variant'

    if variant.isInFrame():
        if variant.isDeletion(): out.append('inframe_deletion')
        if variant.isInsertion(): out.append('inframe_insertion')
        if variant.isComplex():    out.append('inframe_indel')
    else:
        out.append('frameshift_variant')
        return '|'.join(out)

    if protein == mutprotein:
        out.append('synonymous_variant')
        return '|'.join(out)

    if protein[0] != mutprotein[0]:
        out.append('initiator_codon_variant')
        return '|'.join(out)

    if (not protein == mutprotein) and len(protein) == len(mutprotein): out.append('missense_variant')

    while len(protein) > 0 and len(mutprotein) > 0:
        if protein[0] == mutprotein[0]:
            protein = protein[1:]
            mutprotein = mutprotein[1:]
        else:
            break

    if protein == '': return '3_prime_UTR_variant'

    if protein[0] == 'X' and len(mutprotein) == 0: out.append('stop_lost')
    else:
        if protein[0] == 'X' and mutprotein[0] != 'X': out.append('stop_lost')

    while len(protein) > 0 and len(mutprotein) > 0:
        if protein[-1] == mutprotein[-1]:
            protein = protein[:-1]
            mutprotein = mutprotein[:-1]
        else:
            break

    if 'X' in mutprotein: out.append('stop_gained')

    if ('stop_gained' in out) or ('stop_lost' in out):
        if 'missense_variant' in out: out.remove('missense_variant')

    return '|'.join(out)


def isInSpliceDonor(transcript, variant):
    if transcript.strand == 1:
        for exon in transcript.exons:
            isLastExon = (exon.index == len(transcript.exons))
            if not isLastExon and variant.overlap(exon.end + 1, exon.end + 2): return True
        return False
    else:
        for exon in transcript.exons:
            isLastExon = (exon.index == len(transcript.exons))
            if not isLastExon and variant.overlap(exon.start - 1, exon.start): return True
        return False


def isInSpliceAcceptor(transcript, variant):
    if transcript.strand == 1:
        for exon in transcript.exons:
            isFirstExon = (exon.index == 1)
            if not isFirstExon and variant.overlap(exon.start - 1, exon.start): return True
        return False
    else:
        for exon in transcript.exons:
            isFirstExon = (exon.index == 1)
            if not isFirstExon and variant.overlap(exon.end + 1, exon.end + 2): return True

        return False


def isInSplicingRegion(transcript, variant):
    if transcript.strand == 1:
        for exon in transcript.exons:
            isFirstExon = (exon.index == 1)
            isLastExon = (exon.index == len(transcript.exons))
            if not isLastExon and variant.overlap(exon.end - 2, exon.end + 8): return True
            if not isFirstExon and variant.overlap(exon.start - 7, exon.start + 3): return True
        return False
    else:
        for exon in transcript.exons:
            isFirstExon = (exon.index == 1)
            isLastExon = (exon.index == len(transcript.exons))
            if not isLastExon and variant.overlap(exon.start - 7, exon.start + 3): return True
            if not isFirstExon and variant.overlap(exon.end - 2, exon.end + 8): return True
        return False

#######################################################################################################################

def checkUTR(transcript, variant):
    if variant.isInsertion():
        x = variant.pos - 1
        y = variant.pos
    else:
        x = variant.pos
        y = variant.pos + len(variant.ref) - 1
    if transcript.isPositionOutsideCDS_5prime(x) or transcript.isPositionOutsideCDS_5prime(y): return 'UTR5'
    if transcript.isPositionOutsideCDS_3prime(x) or transcript.isPositionOutsideCDS_3prime(y): return 'UTR3'
    return ''
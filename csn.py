#!/usr/bin/env python


# CSN annotation
#######################################################################################################################

import core


# Class representing a CSN annotation
class CSNAnnot():
    # Constructor
    def __init__(self, coord1, intr1, coord2, intr2, dna, protein, coord1_ins, intr1_ins, coord2_ins, intr2_ins, dna_ins):
        self.coord1 = coord1
        self.intr1 = intr1
        self.coord2 = coord2
        self.intr2 = intr2
        self.dna = dna
        self.protein = protein
        self.coord1_ins = coord1_ins
        self.intr1_ins = intr1_ins
        self.coord2_ins = coord2_ins
        self.intr2_ins = intr2_ins
        self.dna_ins = dna_ins

    # Getting annotation as a single String
    def getAsString(self):
        # Adding the first part of the csn annotation (coordinates)
        ret = 'c.' + str(self.coord1)
        if self.intr1 != 0:
            if self.intr1 > 0:
                ret += '+' + str(self.intr1)
            else:
                ret += str(self.intr1)
        if self.coord2 is not None:
            ret += '_' + str(self.coord2)
            if self.intr2 != 0:
                if self.intr2 > 0:
                    ret += '+' + str(self.intr2)
                else:
                    ret += str(self.intr2)
        # Adding the second part of the csn annotation (DNA and protein level)
        ret += self.dna + self.protein
        return ret

    # Getting annotation as a list of separate values
    def getAsFields(self):
        if not self.coord1_ins == '':
            return self.coord1_ins, self.intr1_ins, self.coord2_ins, self.intr2_ins, self.dna_ins, self.protein
        else:
            return self.coord1, self.intr1, self.coord2, self.intr2, self.dna, self.protein


#######################################################################################################################

# Getting CSN annotation of a given variant
def getAnnotation(variant, transcript, reference, prot, mutprot):

    # Creating csn annotation coordinates
    coord1, intr1, coord2, intr2 = calculateCSNCoordinates(variant, transcript)

    # Creating DNA level annotation
    dna, dna_ins = makeDNAannotation(variant, transcript, reference)

    # Creating protein level annotation
    where = transcript.whereIsThisVariant(variant)
    if not '-' in where and where.startswith('Ex'):
        protein, protchange = makeProteinString(variant, transcript, reference, prot, mutprot, coord1)
    else:
        protein, protchange = '', ('.','.','.')

    # Transforming coordinates if the variant is a duplication
    if not dna_ins == '':
        coord1_ins, intr1_ins, coord2_ins, intr2_ins = coord1, intr1, coord2, intr2
        coord1, intr1, coord2, intr2 = duplicationCoordinates(variant, transcript)
    else:
        coord1_ins, intr1_ins, coord2_ins, intr2_ins = '', '', '', ''

    # Creating and returning csnAnnot object
    csn = CSNAnnot(coord1, intr1, coord2, intr2, dna, protein, coord1_ins, intr1_ins, coord2_ins, intr2_ins, dna_ins)
    return csn, protchange


# Calculating csn annotation coordinates
def calculateCSNCoordinates(variant, transcript):
    # Returning coordinates if variant is a base substitution
    if variant.isSubstitution():
        x, y = transformToCSNCoordinate(variant.pos, transcript)
        return x, y, None, None

    # Returning coordinates if variant is an insertion
    if variant.isInsertion():
        startx, starty = transformToCSNCoordinate(variant.pos - 1, transcript)
        endx, endy = transformToCSNCoordinate(variant.pos, transcript)
        if transcript.strand == 1:
            return startx, starty, endx, endy
        else:
            return endx, endy, startx, starty

    # Returning coordinates if variant is a deletion
    if variant.isDeletion():
        startx, starty = transformToCSNCoordinate(variant.pos, transcript)
        if len(variant.ref) == 1:
            return startx, starty, None, None
        endx, endy = transformToCSNCoordinate(variant.pos + len(variant.ref) - 1, transcript)
        if transcript.strand == 1:
            return startx, starty, endx, endy
        else:
            return endx, endy, startx, starty

    # Returning coordinates if variant is a complex indel
    if variant.isComplex():
        startx, starty = transformToCSNCoordinate(variant.pos, transcript)
        if len(variant.ref) == 1:
            return startx, starty, None, None
        endx, endy = transformToCSNCoordinate(variant.pos + len(variant.ref) - 1, transcript)
        if transcript.strand == 1:
            return startx, starty, endx, endy
        else:
            return endx, endy, startx, starty

    return None, None, None, None


# Calculating DNA level annotation of the variant
def makeDNAannotation(variant, transcript, reference):
    # Returning DNA level annotation if variant is a base substitution
    if variant.isSubstitution():
        if transcript.strand == 1:
            return variant.ref + '>' + variant.alt, ''
        else:
            return variant.ref.reverseComplement() + '>' + variant.alt.reverseComplement(), ''

    # Returning DNA level annotation if variant is an insertion
    if variant.isInsertion():
        insert = core.Sequence(variant.alt)
        if transcript.strand == 1:
            before = reference.getReference(variant.chrom, variant.pos - len(insert), variant.pos - 1)
            # Checking if variant is a duplication
            if insert == before and (variant.pos - len(insert) >= transcript.transcriptStart):
                if len(insert) > 4: return 'dup' + str(len(insert)), 'ins' + insert
                return 'dup' + insert, 'ins' + insert
            else:
                return 'ins' + insert, ''
        else:
            before = reference.getReference(variant.chrom, variant.pos, variant.pos + len(insert) - 1)
            # Checking if variant is a duplication
            if insert == before and (variant.pos + len(insert) - 1 <= transcript.transcriptEnd):
                if len(insert) > 4: return 'dup' + str(len(insert)), 'ins' + insert.reverseComplement()
                return 'dup' + insert.reverseComplement(), 'ins' + insert.reverseComplement()
            else:
                return 'ins' + insert.reverseComplement(), ''

    # Returning DNA level annotation if variant is a deletion
    if variant.isDeletion():
        if len(variant.ref) > 4:
            return 'del' + str(len(variant.ref)), ''
        if transcript.strand == 1:
            return 'del' + variant.ref, ''
        else:
            return 'del' + core.Sequence(variant.ref).reverseComplement(), ''

    # Returning DNA level annotation if variant is a complex indel
    if variant.isComplex():
        if transcript.strand == 1:
            return 'delins' + variant.alt, ''
        else:
            return 'delins' + variant.alt.reverseComplement(), ''


# Calculating protein level annotation of the variant
def makeProteinString(variant, transcript, reference, prot, mutprot, coord1):

    protcopy = prot[:]

    # Checking if there was no change in protein sequence
    if prot == mutprot:
        idx = coord1 / 3
        if coord1 % 3 > 0: idx += 1
        return '_p.=', (str(idx), prot[idx-1], prot[idx-1])

    # Checking if the variant affects the initiating amino acid
    if prot[0] != mutprot[0]: return '_p.' + changeTo3letters(prot[0]) + '1?', ('1', prot[0], mutprot[0])

    # Trimming common starting substring
    leftindex = 1
    rightindex = len(prot)
    while len(prot) > 0 and len(mutprot) > 0:
        if prot[0] == mutprot[0]:
            prot = prot[1:]
            mutprot = mutprot[1:]
            leftindex += 1
        else:
            break

    if prot == '': return '', ('.','.','.')

    # Frameshift mutations
    if (len(variant.alt) - len(variant.ref)) % 3 > 0:
        if len(prot) > 0: return '', (str(leftindex), prot[0], '.')
        else: return '', (str(leftindex), '.', '.')

    # Checking if variant results in a stop lost mutation
    if prot[0] == 'X' and len(mutprot) == 0:
        return '_p.X' + str(leftindex) + '?extX?', (str(leftindex), 'X', '?')

    if prot[0] == 'X' and mutprot[0] != 'X':
        nextstop = mutprot.find('X')
        if nextstop != -1:
            return '_p.X' + str(leftindex) + changeTo3letters(mutprot[0]) + 'extX' + str(nextstop), (str(leftindex), 'X', mutprot[0])
        else:
            return '_p.X' + str(leftindex) + changeTo3letters(mutprot[0]) + 'extX?', (str(leftindex), 'X', mutprot[0])

    # Trimming common ending substring
    while len(prot) > 0 and len(mutprot) > 0:
        if prot[-1] == mutprot[-1]:
            prot = prot[:-1]
            mutprot = mutprot[:-1]
            rightindex -= 1
        else:
            break

    # Checking if variant results in an amino acid change
    if len(prot) == 1 and len(mutprot) == 1: return '_p.' + changeTo3letters(prot) + str(leftindex) + changeTo3letters(
        mutprot), (str(leftindex), prot, mutprot)

    # Checking if variant results in a deletion
    if len(mutprot) == 0:
        if len(prot) == 1:
            return '_p.' + changeTo3letters(prot) + str(leftindex) + 'del', (str(leftindex), prot, '-')
        else:
            if leftindex == rightindex: protpos = str(leftindex)
            else: protpos = str(leftindex) + '-' + str(rightindex)
            return '_p.' + changeTo3letters(prot[0]) + str(leftindex) + '_' + changeTo3letters(prot[-1]) + str(
                rightindex) + 'del',  (protpos, prot, '-')

    # Checking if variant results in an insertion or duplication
    if len(prot) == 0:
        if protcopy[leftindex - len(mutprot) - 1:leftindex - 1] == mutprot:
            if len(mutprot) == 1:
                return '_p.' + changeTo3letters(protcopy[leftindex - len(mutprot) - 1]) + str(
                    leftindex - len(mutprot)) + 'dup',  (str(leftindex - 1)+'-'+str(rightindex + 1), '-', mutprot)
            else:
                return '_p.' + changeTo3letters(protcopy[leftindex - len(mutprot) - 1]) + str(
                    leftindex - len(mutprot)) + '_' + changeTo3letters(protcopy[leftindex - 2]) + str(
                    leftindex - 1) + 'dup',  (str(leftindex - 1)+'-'+str(rightindex + 1), '-', mutprot)
        if 'X' in mutprot: mutprot = mutprot[:mutprot.find('X') + 1]
        return '_p.' + changeTo3letters(protcopy[leftindex - 2]) + str(leftindex - 1) + '_' + changeTo3letters(
            protcopy[rightindex]) + str(rightindex + 1) + 'ins' + changeTo3letters(mutprot), (str(leftindex - 1)+'-'+str(rightindex + 1), '-', mutprot)

    # Checking if variant results in a complex change
    if len(prot) > 0 and len(mutprot) > 0:
        ret = '_p.'
        if len(prot) == 1:
            ret += changeTo3letters(prot) + str(leftindex)
        else:
            ret += changeTo3letters(prot[0]) + str(leftindex) + '_' + changeTo3letters(prot[-1]) + str(rightindex)
        ret += 'delins'
        if 'X' in mutprot: mutprot = mutprot[:mutprot.find('X') + 1]
        ret += changeTo3letters(mutprot)

        if leftindex == rightindex: protpos = str(leftindex)
        else: protpos = str(leftindex) + '-' + str(rightindex)
        return ret, (protpos, prot, mutprot)

    return "", ()


# Transforming a genomic position ot csn coordinate
def transformToCSNCoordinate(pos, transcript):

    prevExonEnd = 99999999

    # Checking if genomic position is not outside CDS
    if not transcript.isPositionOutsideCDS(pos):
        sumOfExonLengths = -transcript.codingStart + 1
        # Iterating through exons
        for i in range(len(transcript.exons)):
            exon = transcript.exons[i]
            if i > 0:
                if transcript.strand == 1:
                    # Checking if genomic position is within intron
                    if prevExonEnd < pos < exon.start + 1:
                        if pos <= (exon.start + 1 - prevExonEnd) / 2 + prevExonEnd:
                            x, y = transformToCSNCoordinate(prevExonEnd, transcript)
                            return x, pos - prevExonEnd
                        else:
                            x, y = transformToCSNCoordinate(exon.start + 1, transcript)
                            return x, pos - exon.start - 1
                else:
                    # Checking if genomic position is within intron
                    if exon.end < pos < prevExonEnd:
                        if pos >= (prevExonEnd - exon.end + 1) / 2 + exon.end:
                            x, y = transformToCSNCoordinate(prevExonEnd, transcript)
                            return x, prevExonEnd - pos
                        else:
                            x, y = transformToCSNCoordinate(exon.end, transcript)
                            return x, exon.end - pos
            # Checking if genomic position is within exon
            if exon.start + 1 <= pos <= exon.end:
                if transcript.strand == 1:
                    return sumOfExonLengths + pos - exon.start, 0
                else:
                    return sumOfExonLengths + exon.end - pos + 1, 0
            # Calculating sum of exon lengths up to this point
            sumOfExonLengths += exon.length
            if transcript.strand == 1:
                prevExonEnd = exon.end
            else:
                prevExonEnd = exon.start + 1

    # If genomic position is outside CDS
    else:

        sumpos = 0
        for i in range(len(transcript.exons)):
            exon = transcript.exons[i]

            if i > 0:
                if transcript.strand == 1:
                    # Checking if genomic position is within intron
                    if prevExonEnd < pos < exon.start + 1:
                        if pos <= (exon.start + 1 - prevExonEnd) / 2 + prevExonEnd:
                            x, y = transformToCSNCoordinate(prevExonEnd, transcript)
                            return x, pos - prevExonEnd
                        else:
                            x, y = transformToCSNCoordinate(exon.start + 1, transcript)
                            return x, pos - exon.start - 1
                else:
                    # Checking if genomic position is within intron
                    if exon.end < pos < prevExonEnd:
                        if pos >= (prevExonEnd - exon.end + 1) / 2 + exon.end:
                            x, y = transformToCSNCoordinate(prevExonEnd, transcript)
                            return x, prevExonEnd - pos
                        else:
                            x, y = transformToCSNCoordinate(exon.end, transcript)
                            return x, exon.end - pos

            if transcript.strand == 1:
                if pos > transcript.codingEndGenomic:
                    if transcript.codingEndGenomic < exon.start+1 and exon.end < pos: sumpos += exon.length
                    elif exon.contains(transcript.codingEndGenomic) and exon.end < pos: sumpos += exon.end - transcript.codingEndGenomic + 1
                    elif exon.contains(transcript.codingEndGenomic) and exon.contains(pos): sumpos += pos - transcript.codingEndGenomic
                    elif transcript.codingEndGenomic < exon.start+1 and exon.contains(pos): sumpos += pos - exon.start - 1
                if pos < transcript.codingStartGenomic:
                    if pos < exon.start+1 and exon.end < transcript.codingStartGenomic: sumpos += exon.length
                    elif pos < exon.start+1 and exon.contains(transcript.codingStartGenomic): sumpos += transcript.codingStartGenomic - exon.start
                    elif exon.contains(pos) and exon.contains(transcript.codingStartGenomic): sumpos += transcript.codingStartGenomic - pos
                    elif exon.contains(pos) and exon.end < transcript.codingStartGenomic: sumpos += exon.end - pos

            if transcript.strand == -1:
                if pos < transcript.codingEndGenomic:
                    if pos < exon.start+1 and exon.end < transcript.codingEndGenomic: sumpos += exon.length
                    elif pos < exon.start+1 and exon.contains(transcript.codingEndGenomic): sumpos += transcript.codingEndGenomic - exon.start
                    elif exon.contains(pos) and exon.contains(transcript.codingEndGenomic): sumpos += transcript.codingEndGenomic - pos
                    elif exon.contains(pos) and exon.end < transcript.codingEndGenomic: sumpos += exon.end - pos
                if transcript.codingStartGenomic < pos:
                    if transcript.codingStartGenomic < exon.start+1 and exon.end < pos: sumpos += exon.length
                    elif exon.contains(transcript.codingStartGenomic) and exon.end < pos: sumpos += exon.end - transcript.codingStartGenomic + 1
                    elif exon.contains(transcript.codingStartGenomic) and exon.contains(pos): sumpos += pos - transcript.codingStartGenomic
                    elif transcript.codingStartGenomic < exon.start+1 and exon.contains(pos): sumpos += pos - exon.start - 1

            if transcript.strand == 1: prevExonEnd = exon.end
            else: prevExonEnd = exon.start + 1


        if transcript.strand == 1:
            if pos > transcript.codingEndGenomic:  return '+'+str(sumpos), 0
            if pos < transcript.codingStartGenomic: return '-'+str(sumpos), 0
        else:
            if pos < transcript.codingEndGenomic: return '+'+str(sumpos), 0
            if pos > transcript.codingStartGenomic: return '-'+str(sumpos), 0

        return str(sumpos), 0




# Calculating csn coordinates for duplications
def duplicationCoordinates(variant, transcript):
    if transcript.strand == 1:
        coord1, intr1 = transformToCSNCoordinate(variant.pos - len(variant.alt), transcript)
        if len(variant.alt) == 1:
            coord2, intr2 = None, None
        else:
            coord2, intr2 = transformToCSNCoordinate(variant.pos - 1, transcript)
    else:
        coord1, intr1 = transformToCSNCoordinate(variant.pos + len(variant.alt) - 1, transcript)
        if len(variant.alt) == 1:
            coord2, intr2 = None, None
        else:
            coord2, intr2 = transformToCSNCoordinate(variant.pos, transcript)
    return coord1, intr1, coord2, intr2


# Changing protein sequence of 1-letter amino acid code to 3-letter code
def changeTo3letters(aas):
    ret = ''
    codes = {
        'I': 'Ile', 'M': 'Met', 'T': 'Thr', 'N': 'Asn',
        'K': 'Lys', 'S': 'Ser', 'R': 'Arg', 'L': 'Leu',
        'P': 'Pro', 'H': 'His', 'Q': 'Gln', 'V': 'Val',
        'A': 'Ala', 'D': 'Asp', 'E': 'Glu', 'G': 'Gly',
        'F': 'Phe', 'Y': 'Tyr', 'C': 'Cys', 'W': 'Trp',
        '*': 'X', 'X': 'X', 'x': 'X'}
    for aa in aas: ret += codes[aa]
    return ret

#######################################################################################################################

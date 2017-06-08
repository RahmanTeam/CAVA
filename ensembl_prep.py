#!/usr/bin/env python


# Transcript database preparation tool (ensembl_prep)
#######################################################################################################################

# Basic imports
from __future__ import division
import os
import sys
import gzip
import datetime
from optparse import OptionParser
from operator import itemgetter
import urllib

if sys.version_info[0] == 3:
    print '\nCAVA does not run on Python 3.\n'
    quit()

# Checking if installation is complete
if not os.path.isfile(os.path.dirname(os.path.realpath(__file__))+"/pysamdir/pysam/Pileup.py"):
    print '\nCAVA installation not complete.'
    print 'Please run install.sh in the cava directory.\n'
    quit()

# Pysam import
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)) + '/pysamdir')
import pysam


#######################################################################################################################

# Class representing a transcript
class Transcript(object):

    # Constructor
    def __init__(self):
        self.ENST = None
        self.GENE = None
        self.ENSG = None
        self.CHROM = None
        self.STRAND = None
        self.POS = None
        self.POSEND = None
        self.GENETYPE = None
        self.TRANSTYPE = None
        self.CODING_START = None
        self.CODING_END = None
        self.CODING_START_RELATIVE = None
        self.CCDS = None
        self.EXONS = []
        self.PROTL = None
        self.CDNAL = None
        self.isComplete = None


    # Get summary information about the transcript
    def getInfoString(self):
        if self.STRAND == '1': ret = '+/'
        else: ret = '-/'
        cdna = self.getcDNALength()
        return ret+str(round((self.POSEND-self.POS+1)/1000,1))+'kb/'+str(len(self.EXONS))+'/'+str(round(cdna/1000,1))+'kb/'+str(self.getProteinLength())


    # Get cDNA length of the transcript
    def getcDNALength(self):
        ret = 0
        for exon in self.EXONS: ret += exon.END - exon.START
        return ret


    # Get protein length of the transcript
    def getProteinLength(self):
        codingdna = 0
        if self.STRAND == '1':
            for exon in self.EXONS:
                if exon.END < self.CODING_START: continue
                if exon.START > self.CODING_END: continue
                if exon.START <= self.CODING_START <= exon.END: start = self.CODING_START
                else: start = exon.START + 1
                if exon.START <= self.CODING_END <= exon.END: end = self.CODING_END
                else: end = exon.END
                codingdna += end - start + 1
        else:
            for exon in self.EXONS:
                if exon.START > self.CODING_START: continue
                if exon.END < self.CODING_END: continue
                if exon.START <= self.CODING_START <= exon.END: end = self.CODING_START
                else: end = exon.END
                if exon.START <= self.CODING_END <= exon.END: start = self.CODING_END
                else: start = exon.START + 1
                codingdna += end - start + 1
        return int((codingdna - 3) / 3)


    # Check if it is a candidate transcript
    def isCandidate(self):
        if not (self.GENETYPE=='protein_coding' and self.TRANSTYPE=='protein_coding'): return False
        return (self.CODING_START is not None and self.CODING_END is not None) and self.isComplete


    # Output transcript
    def output(self, outfile, outfile_list):
        out = self.ENST + '\t' + self.GENE + '\t' + self.ENSG + '\t' + self.getInfoString() + '\t' + self.CHROM + '\t' + self.STRAND + '\t' + str(self.POS)
        out += '\t' + str(self.POSEND) + '\t' + str(self.CODING_START_RELATIVE) + '\t' + str(self.CODING_START)
        out += '\t' + str(self.CODING_END)
        for exondata in self.EXONS: out += '\t' + str(exondata.START) + '\t' + str(exondata.END)
        outfile.write(out + '\n')
        outfile_list.write(self.ENSG+'\t'+self.GENE+'\t'+self.ENST+'\n')


    # Finalize transcript
    def finalize(self):
        if self.STRAND == '1':
            self.POS = self.EXONS[0].START
            self.POSEND = self.EXONS[len(self.EXONS) - 1].END
            codingStartRelative = 0
            for exondata in self.EXONS:
                if exondata.START <= self.CODING_START <= exondata.END:
                    codingStartRelative += self.CODING_START - exondata.START
                    break
                else:
                    codingStartRelative += exondata.END - exondata.START
            self.CODING_START_RELATIVE = codingStartRelative
        else:
            self.POS = self.EXONS[len(self.EXONS) - 1].START
            self.POSEND = self.EXONS[0].END
            codingStartRelative = 0
            for exondata in self.EXONS:
                if exondata.START <= self.CODING_START <= exondata.END:
                    codingStartRelative += exondata.END - self.CODING_START + 1
                    break
                else:
                    codingStartRelative += exondata.END - exondata.START
            self.CODING_START_RELATIVE = codingStartRelative
        self.PROTL = self.getProteinLength()
        self.CDNAL = self.getcDNALength()


# Class representing an exon
class Exon(object):

    # Constructor
    def __init__(self, start, end):
        self.START = start
        self.END = end


# Class representing a gene
class Gene(object):

    # Constructor
    def __init__(self, symbol, ensg):
        self.SYMBOL = symbol
        self.ENSG = ensg
        self.TRANSCRIPTS = dict()


    # Select ICR transcript
    def selectTranscript(self):
        ccds_set = []
        nonccds_set = []
        for enst,transcript in self.TRANSCRIPTS.iteritems():
            if transcript.CCDS: ccds_set.append(transcript)
            else: nonccds_set.append(transcript)

        if len(ccds_set) > 0: candidates = ccds_set
        else: candidates = nonccds_set

        selected = Transcript()
        selected.PROTL = selected.CDNAL = -1
        for t in candidates:
            if t.PROTL > selected.PROTL: selected = t
            elif t.PROTL == selected.PROTL and t.CDNAL > selected.CDNAL: selected = t

        return selected


    # Output all or selected transcripts
    def output(self, outfile, outfile_list, select, mcg_transcripts):
        if select:
            if self.SYMBOL in mcg_transcripts.keys():
                ok = False
                for _,transcript in self.TRANSCRIPTS.iteritems():
                    if transcript.ENST in mcg_transcripts[self.SYMBOL]:
                        transcript.output(outfile,outfile_list)
                        ok = True
            if self.SYMBOL not in mcg_transcripts.keys() or not ok:
                transcript = self.selectTranscript()
                transcript.output(outfile,outfile_list)
        else:
            for _,transcript in self.TRANSCRIPTS.iteritems():
                transcript.output(outfile,outfile_list)


#######################################################################################################################

# Retrieve tag value
def getValue(tags, tag):
    ret=None
    for x in tags:
        x = x.strip()
        if x.startswith(tag):
            s = x[x.find('\"') + 1:]
            ret = s[:s.find('\"')]
            break
    return ret


# Retrieve boolean tag value
def getBooleanValue(tags, tag):
    for x in tags:
        x = x.strip()
        if x.startswith('tag'):
            s = x[x.find('\"') + 1:]
            value = s[:s.find('\"')]
            if value==tag: return True
    return False


# Read transcript IDs from file
def readTranscriptIDs(inputfn):
    ret = set()
    for line in open(inputfn): ret.add(line.strip())
    return ret


# Sort records in file
def sortRecords(records, idx1, idx2):
    ret = []
    chroms = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', 'MT', 'X', 'Y']
    for i in range(len(chroms)):
        chrom = chroms[i]
        if chrom in records.keys():
            records[chrom] = sorted(records[chrom], key=itemgetter(idx1,idx2))
    for i in range(len(chroms)):
        chrom = chroms[i]
        if chrom in records.keys():
            for record in records[chrom]: ret.append(record)
    return ret


# Write records to file
def writeToFile(sortedRecords, filename):
    outfile = open(filename, 'w')
    for record in sortedRecords:
        s = str(record[0]).rstrip()
        for i in range(1, len(record)): s += '\t' + str(record[i]).rstrip()
        outfile.write(s + '\n')
    outfile.close()


# Read records from file as a list
def readRecords(inputfn):
    ret = []
    for line in open(inputfn): ret.append(line.strip())
    return ret


# Process Ensembl data
def run(ver, options, genome_build):

    # Read manually selected MCG transcripts from file
    dir = os.path.dirname(os.path.realpath(sys.argv[0]))
    mcg_transcripts=dict()
    for line in open(dir+'/ensembl_prep/MCG_transcripts.txt'):
        line = line.strip()
        if line == '': continue
        cols = line.split('\t')
        if cols[0] not in mcg_transcripts.keys(): mcg_transcripts[cols[0]] = set()
        mcg_transcripts[cols[0]].add(cols[1])
    # Changing transcript for certain releases
    if int(options.ensembl)>=71: mcg_transcripts['BMPR1A'] = {'ENST00000372037'}
    if int(options.ensembl)>=69: mcg_transcripts['PRKAR1A'] = {'ENST00000392711'}
    if int(options.ensembl)>=76:
        mcg_transcripts['MEN1'] = {'ENST00000394374'}
        mcg_transcripts['RECQL4'] = {'ENST00000617875'}

    # Dictionary of Gene objects
    genesdata = dict()

    # Load custom transcript IDs
    transIDs = set()
    if options.input is not None:
        transIDs = readTranscriptIDs(options.input)
        print '\nOnly ' + str(len(transIDs)) + ' transcripts read from ' + options.input + ' are considered\n'
    else: print '\nAll transcripts from the Ensembl release are considered\n'

    # Print out info
    if options.select: print 'Transcript selection switched on\n'

    # Load candidate and CCDS data for Ensembl <75
    candidates = dict()
    if int(options.ensembl) < 75:
        for line in open(dir+'/ensembl_prep/info'+options.ensembl+'.txt'):
            line=line.strip()
            if line=='': continue
            cols = line.split('\t')
            if cols[0] not in candidates.keys(): candidates[cols[0]]=dict()
            candidates[cols[0]][cols[1]]=int(cols[2])

    # Download Ensembl data
    sys.stdout.write('Downloading Ensembl database... ')
    sys.stdout.flush()
    url = 'ftp://ftp.ensembl.org/pub/release-' + options.ensembl + '/gtf/homo_sapiens/Homo_sapiens.' + genome_build + '.' + options.ensembl + '.gtf.gz'
    try:
        urllib.urlretrieve(url, 'ensembl_data.gz')
    except:
        print '\n\nCannot connect to Ensembl FTP site. No internet connection?\n'
        quit()

    sys.stdout.write('OK\n')

    # Iterate through the lines in the ensembl data file
    sys.stdout.write('Extracting transcript data... ')
    sys.stdout.flush()
    first = True
    prevenst = ''
    transcript = None
    for line in gzip.open('ensembl_data.gz', 'r'):
        line = line.strip()
        if line.startswith('#'): continue
        cols = line.split('\t')

        # Only consider transcripts on the following chromosomes
        if cols[0] not in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', 'MT', 'X', 'Y']: continue

        # Consider only certain types of lines
        if cols[2] not in ['exon','transcript','start_codon','stop_codon']: continue

        # Annotation tags
        tags = cols[8].split(';')
        # Retrieve transcript ID
        enst = getValue(tags, 'transcript_id')

        # Do not consider transcript if it is not on the custom transcript list
        if options.input is not None and enst not in transIDs: continue

        # Finalize and output transcript object
        if not enst == prevenst:

            # Finalize transcript and add to Gene object if candidate
            if not first:
                transcript.finalize()
                if transcript.isCandidate():
                    if transcript.ENSG not in genesdata.keys(): genesdata[transcript.ENSG] = Gene(transcript.GENE, transcript.ENSG)
                    genesdata[transcript.ENSG].TRANSCRIPTS[transcript.ENST] = transcript

            # Initialize new Transcript object
            transcript = Transcript()
            transcript.ENST = enst
            transcript.GENE = getValue(tags, 'gene_name')
            transcript.ENSG = getValue(tags, 'gene_id')
            transcript.CHROM = cols[0]
            if cols[6] == '+': transcript.STRAND = '1'
            else: transcript.STRAND = '-1'

            # Retrieve gene biotype and transcript biotype
            transcript.GENETYPE = getValue(tags, 'gene_type')
            if transcript.GENETYPE is None: transcript.GENETYPE = getValue(tags, 'gene_biotype')
            transcript.TRANSTYPE = getValue(tags, 'transcript_type')
            if transcript.TRANSTYPE is None: transcript.TRANSTYPE = getValue(tags, 'transcript_biotype')
            if transcript.TRANSTYPE is None: transcript.TRANSTYPE = cols[1]

        # If line represents an exon
        if cols[2] == 'exon':
            idx = 0
            for x in tags:
                x = x.strip()
                if x.startswith('exon_number'):
                    s = x[x.find('\"') + 1:]
                    idx = int(s[:s.find('\"')]) - 1
                    break
            start = int(cols[3]) - 1
            end = int(cols[4])
            if idx >= len(transcript.EXONS):
                for _ in range(len(transcript.EXONS), idx + 1): transcript.EXONS.append(None)
            transcript.EXONS[idx] = Exon(start, end)

        if cols[2] == 'start_codon':
            if transcript.STRAND == '1':
                if transcript.CODING_START is None or int(cols[3]) < transcript.CODING_START: transcript.CODING_START = int(cols[3])
            else:
                if transcript.CODING_START is None or int(cols[4]) > transcript.CODING_START: transcript.CODING_START = int(cols[4])

        if cols[2] == 'stop_codon':
            if transcript.STRAND == '1':
                if transcript.CODING_END is None or int(cols[4]) > transcript.CODING_END: transcript.CODING_END = int(cols[4])
            else:
                if transcript.CODING_END is None or int(cols[3]) < transcript.CODING_END: transcript.CODING_END = int(cols[3])

        # Check if transcript is complete and is a CCDS transcript
        if transcript.isComplete is None:
            if int(options.ensembl) < 75:
                if transcript.ENST in candidates[transcript.CHROM].keys():
                    transcript.CCDS = (candidates[transcript.CHROM][transcript.ENST] == 1)
                    transcript.isComplete = True
                else:
                    transcript.isComplete = False
            else:
                transcript.isComplete = not (getBooleanValue(tags, 'cds_start_NF') or getBooleanValue(tags, 'cds_end_NF'))
                if getValue(tags, 'ccds_id') is not None: transcript.CCDS=True
                else: transcript.CCDS=False

        prevenst = enst
        if first: first = False

    # Finalize last transcript and add to Gene object if candidate
    if transcript is not None:
        transcript.finalize()
        if transcript.isCandidate():
            if transcript.ENSG not in genesdata.keys(): genesdata[transcript.ENSG] = Gene(transcript.GENE, transcript.ENSG)
            genesdata[transcript.ENSG].TRANSCRIPTS[transcript.ENST] = transcript

    # If no transcript ID from the input file was found in the Ensembl release
    if len(genesdata) == 0:
        print '\n\nNo transcripts from '+options.input+' found in Ensembl release.'
        print '\nNo transcript database created.'
        print "-----------------------------------------------------------------\n"
        os.remove('ensembl_data.gz')
        quit()

    # Initialize temporary output file
    outfile = open('temp.txt', 'w')

    # Initialize output list file if needed
    outfile_list = open(options.output+'.txt','w')
    outfile_list.write('# Created by ensembl_prep '+ver+' based on Ensembl release '+options.ensembl+' (genome build '+genome_build+')\n')
    outfile_list.write('ENSG\tGENE\tENST\n')

    # Output transcripts of each gene
    for ensg, gene in genesdata.iteritems(): gene.output(outfile,outfile_list,options.select,mcg_transcripts)

    # Close temporary output files
    outfile.close()
    outfile_list.close()

    # Sort temporary output file
    data = dict()
    counter = 0
    for line in open('temp.txt'):
        if not line.startswith('ENST'): continue
        counter += 1
        line.rstrip()
        record = line.split('\t')
        record[6] = int(record[6])
        if record[4] in data.keys():
            data[record[4]].append(record)
        else:
            data[record[4]] = []
            data[record[4]].append(record)

    sys.stdout.write('OK\n')
    sys.stdout.write('Sorting transcripts... ')
    sys.stdout.flush()
    sortedRecords = sortRecords(data, 6, 7)
    writeToFile(sortedRecords, options.output)

    # Remove temporary files
    sys.stdout.write('OK\n')
    sys.stdout.write('Removing temporary files... ')
    sys.stdout.flush()
    os.remove('temp.txt')
    os.remove('ensembl_data.gz')
    sys.stdout.write('OK\n')

    # Return sorted records
    return len(sortedRecords)


# Use Tabix to index output file     
def indexFile(options):
    sys.stdout.write('Compressing output file... ')
    sys.stdout.flush()
    pysam.tabix_compress(options.output, options.output + '.gz', force=True)
    sys.stdout.write('OK\n')
    sys.stdout.write('Indexing output file... ')
    sys.stdout.flush()
    pysam.tabix_index(options.output + '.gz', seq_col=4, start_col=6, end_col=7, meta_char='#', force=True)
    sys.stdout.write('OK\n')


# CHeck if string is a number (integer)
def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


#######################################################################################################################

if __name__ == '__main__':

    # Version number
    ver = 'v1.2.2'

    # Command line argument parsing
    descr = 'ensembl_prep '+ver+' is a simple tool for generating the local Ensembl transcript database file used by CAVA (via the @ensembl option flag).'
    epilog = '\nExample usage: ./ensembl_prep.py -i input.txt -e 70 -o out -s\n\n'
    OptionParser.format_epilog = lambda self, formatter: self.epilog
    parser = OptionParser(usage='python path/to/cava/ensembl_prep.py <options>', version=ver, description=descr,epilog=epilog)
    parser.add_option('-i', "--in", default=None, dest='input', action='store',help="Input filename (list of ENST IDs)")
    parser.add_option('-o', "--out", default=None, dest='output', action='store', help="Output filename prefix")
    parser.add_option('-e', "--ens", default=None, dest='ensembl', action='store', help="Ensembl release version")
    #parser.add_option('-g', "--genome", dest='genome', action='store', default='GRCh37',help="Human genome reference version (default: %default)")
    parser.add_option('-s', "--select", default=False, dest='select', action='store_true',help="Select transcript for each gene [default: %default]")

    (options, args) = parser.parse_args()

    # Checking if all required options specified
    if options.ensembl is None:
        print '\nError: no Ensembl release specified. Use option -h to get help!\n'
        quit()
    if not is_number(options.ensembl):
        print '\nError: Ensembl release specified is not an integer. Use option -h to get help!\n'
        quit()
    if options.output is None:
        print '\nError: no output file name specified. Use option -h to get help!\n'
        quit()

    # Must use Ensembl release >= 70
    if not (int(options.ensembl) >= 70 or int(options.ensembl) == 65) :
        print '\nError: This version works with Ensembl v65 or >= v70.\n'
        quit()

    # Genome build
    # genome_build = options.genome
    genome_build = 'GRCh37' if int(options.ensembl) <= 75 else 'GRCh38'

    # Printing out version information  
    print "\n---------------------------------------------------------------------------------------"
    print 'CAVA ' + ver + ' transcript database preparation tool (ensembl_prep) is now running'
    print 'Started: ', datetime.datetime.now(), '\n'

    # Print info
    print 'Ensembl version:  ' + options.ensembl
    print 'Reference genome: ' + genome_build

    # Creating compressed output file
    Nretrieved = run(ver, options, genome_build)
    print '\nA total of ' + str(Nretrieved) + ' transcripts have been retrieved\n'

    # Indexing output file with Tabix
    indexFile(options)

    # Removing uncompressed output file
    os.remove(options.output)

    # Printing out summary information
    print ''
    print '---------------------'
    print 'Output files created:'
    print '---------------------'
    print options.output + '.gz (transcript database)'
    print options.output + '.gz.tbi (index file)'
    print options.output + '.txt (list of transcripts)'

    print ''
    print 'CAVA ensembl_prep successfully finished: ', datetime.datetime.now()
    print "---------------------------------------------------------------------------------------\n"

#######################################################################################################################

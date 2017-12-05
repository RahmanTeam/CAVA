#!env/bin/python

from optparse import OptionParser
from ensembldb import main

# Version
version = '1.2.3'

# Command line argument parsing
descr = 'CAVA ensembl_db v' + version
epilog = '\nExample usage: CAVA-{}/ensembl_db -e 75 -o ensembl_db_75\n\n'.format(version)
OptionParser.format_epilog = lambda self, formatter: self.epilog
parser = OptionParser(usage='CAVA-{}/ensembl_db <options>'.format(version), version=version, description=descr, epilog=epilog)
parser.add_option('-i', "--input", default=None, dest='input', action='store',help="Input filename (list of ENST IDs)")
parser.add_option('-o', "--output", default=None, dest='output', action='store', help="Output filename prefix")
parser.add_option('-e', "--ensembl", default=None, dest='ensembl', action='store', help="Ensembl release version")
parser.add_option('-s', "--select", default=False, dest='select', action='store_true',help="Select transcript for each gene [default: %default]")
(options, args) = parser.parse_args()

main.run(options, version)


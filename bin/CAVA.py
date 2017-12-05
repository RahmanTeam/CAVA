#!env/bin/python

from optparse import OptionParser
from cava_ import main
from cava_ import helper

# Version
version = 'v1.2.3'

# Read default configuration file name from the default_config_path file
default_config_file = helper.defaultConfigPath()

descr = 'CAVA (Clinical Annotation of VAriants) is a lightweight, fast and flexible NGS variant annotation tool that provides consistent transcript-level annotation.'
epilog = '\nExample usage: CAVA-{}/cava -c config.txt -i input.vcf -o output\n\n'.format(version)
OptionParser.format_epilog = lambda self, formatter: self.epilog
parser = OptionParser(usage='CAVA-{}/cava <options>'.format(version), version=version, description=descr, epilog=epilog)
parser.add_option('-i', "--input", default='input.vcf', dest='input', action='store', help="Input file name [default value: %default]")
parser.add_option('-o', "--output", default='output', dest='output', action='store', help="Output file name prefix [default value: %default]")
parser.add_option('-c', "--config", default=default_config_file, dest='conf', action='store', help="Configuration file name [default value: %default]")
parser.add_option('-s', "--stdout", default=False, dest='stdout', action='store_true', help="Write output to standard output [default value: %default]")
parser.add_option('-t', "--threads", default=1, dest='threads', action='store', help="Number of threads [default value: %default]")
(copts, args) = parser.parse_args()

main.run(copts, version, default_config_file)


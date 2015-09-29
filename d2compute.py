#!/bin/python

import optparse

#Create a parser to read arguments
parser = optparse.OptionParser()

#Inputfilename
parser.add_option('-i', action='store', dest='inputfile')

#Outputprefix
parser.add_option('-o', action='store', dest='outputprefix', default='./')

options, reminder = parser.parse_args()

inputfile = options.inputfile
outputprefix = options.outputprefix

species = []
sequences = []

#Extract long sequence or NGS data from file
def analyzeSeq(file):
    f = open(file)
    lines = f.readlines()
    reads = []
    read = ''
    firstline = 1
    for line in lines:
       if line.startswith('>'):
           if not firstline:
               reads.append(read)
               read = ''
           else:
               firstline = 0
       else:
           read += line.strip()
    reads.append(read)
    f.close()
    return reads

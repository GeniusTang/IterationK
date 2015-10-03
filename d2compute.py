#!/bin/python

import errmsg
import optparse
import kcount

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

#Analyze species names and filenames from inputfile.
def analyzeInput(file):
    f = open(file)
    lines = f.readlines()
    f.close()
    species = []
    sequences = []
    for line in lines:
        species.append(line.strip().split()[1])
        sequences.append(analyzeSeq(line.strip().split()[0]))
    return species, sequences
  
#Create a class to save the species name and kmer
class Species:
    def __init__(self, name, readlist, k, inputtype):
        self.name = name
        self.kmer, self.kprob = kcount.create_readlist_kmer(readlist, k, inputtype)
        self.M = len(readlist)
        self.Beta = len(readlist[0])
        self.E = self.M * (self.Beta - k + 1) * self.kprob

'''
def d2compare(X, Y):
    D2 = sum(X.kmer * Y.kmer)
    #Read pairs
    M = X.M 
    #Read length
    Beta = X.Beta 
    #Expectaion of read appearance
    E = X.E 
    
    X2 = X.kmer - E
    Y2 = Y.kmer - E
    D2star = sum(X2 * Y2 / E)
    X
    XYsqrt = np.sqrt((X2**2+Y2**2))
    D2shepp = sum(X2 * Y2 / XYsqrt) 
    
    d2 = 
'''  
    
    


def main():
    #Create a parser to read arguments
    parser = optparse.OptionParser()
    
    #Inputfilename
    parser.add_option('-i', action='store', dest='inputfile')
    
    #Outputprefix
    parser.add_option('-o', action='store', dest='outputprefix', default='./')
    
    #Kvalue
    parser.add_option('-k', action='store', dest='kvalue')

    #Input type
    parser.add_option('-t', action='store', dest='inputtype',\
help='L for long sequences, A for fasta files')

    options, reminder = parser.parse_args()
    
    inputfile = options.inputfile
    outputprefix = options.outputprefix
    inputtype = options.inputtype

    kvalue = errmsg.check_k(options.kvalue)
    errmsg.check_inputfile(inputfile)
    errmsg.check_inputtype(inputtype)

    species, sequences = analyzeInput(inputfile)
    
    test = Species(species[0], sequences[0], kvalue, inputtype)
    print test.name
    print test.M
    print test.Beta
    print test.kprob
    print test.E
    for kmer in test.kmer:
        print kmer
    
     
if __name__ == '__main__':
    main()
     

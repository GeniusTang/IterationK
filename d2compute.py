#!/bin/python
import numpy as np
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
        self.E = sum(self.kmer) * self.kprob
        self.tilde = self.kmer - self.E

#Compare the cosinse similarity between x and y.
def cosine(x, y):
    return sum(x * y) / (np.sqrt(sum(x**2)) * np.sqrt(sum(y**2)))

#Compare d2 between X and Y.
def d2compare(X, Y):
    return (1-cosine(X.kmer, Y.kmer))
    
#Compare d2star between X and Y.
def d2starcompare(X, Y):
    return 0.5 * (1 - cosine(X.tilde/np.sqrt(X.E), Y.tilde/np.sqrt(Y.E)))

#Compare d2shepp between X and Y.
def d2sheppcompare(X, Y):
    denom = (X.tilde**2 + Y.tilde**2) ** 0.25
    return 0.5 * (1 - cosine(X.tilde/denom, Y.tilde/denom))

#Compare the pairwise distance by given method.
def pairwise(species, sequences, kvalue, inputtype, method):
    spelist = []
    spenumber = len(species)
    distmatrix = np.zeros([spenumber]*2)
    for i in range(spenumber):
        spelist.append(Species(species[i], sequences[i], kvalue, inputtype))
    for i in range(spenumber):
        for j in range(i+1, spenumber):
            distance = method(spelist[i], spelist[j])
            distmatrix[i][j] = distance
            distmatrix[j][i] = distance
    return distmatrix
    
#Print the distance matrix in given format.
def formatprint(species, distmatrix, outfile, method, kvalue):
    f = open(outfile+method+'_k'+str(kvalue), 'wt')
    spenumber = len(distmatrix)
    f.write(str(spenumber)+'\n')
    for i in range(spenumber):
        f.write('%-9s'%species[i])
        for j in range(spenumber):
            f.write(' '+'%.4f'%distmatrix[i][j])
        f.write('\n')
    f.close()
     


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

    #Method
    parser.add_option('-m', action='store', dest='d2method')

    options, reminder = parser.parse_args()
    
    inputfile = options.inputfile
    outputprefix = options.outputprefix
    inputtype = options.inputtype
    d2method = options.d2method

    kvalue = errmsg.check_k(options.kvalue)
    errmsg.check_inputfile(inputfile)
    errmsg.check_inputtype(inputtype)
    errmsg.check_outputdir(outputprefix)
    d2method = errmsg.check_d2method(d2method)
    
    species, sequences = analyzeInput(inputfile)
    for method in d2method:
        if method == 'd2':
            distmatrix = pairwise(species, sequences, kvalue, inputtype, d2compare)
        elif method == 'd2star':
            distmatrix = pairwise(species, sequences, kvalue, inputtype, d2starcompare)
        else:
            distmatrix = pairwise(species, sequences, kvalue, inputtype, d2sheppcompare)
        formatprint(species, distmatrix, outputprefix, method, kvalue)
    '''
    print distmatrix
    test1 = Species(species[0], sequences[0], kvalue, inputtype)
    test2 = Species(species[0], sequences[1], kvalue, inputtype)
    print test.name
    print 'prob: ', test.kprob
    print 'E: ', test.E
    print 'tilde: ', test.tilde
    print 'kmer: ', test.kmer 
    print d2compare(test1, test2)
    '''
     
if __name__ == '__main__':
    main()
     

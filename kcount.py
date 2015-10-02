#!/bin/python

import numpy as np


#Generate a quantary dictionary.
SEQDICT = dict(zip('ACGT', range(4)))

#Generate the kmer count array for one read
def create_read_kmer(read, k, inputtype):
    """
    If inputtype is A, count both forward and reverse strand,
    If inputtype is L, only count forward strand.
    """
    length = 4 ** k
    kcount = np.zeros(length)
    temp_f = 0
    for i in range(k):
        temp_f = (temp_f * 4) + SEQDICT[read[i]]
 
    kcount[temp_f] += 1

    for i in range(k, len(read)):
        temp_f = 4 * (temp_f - SEQDICT[read[i-k]] * (4**(k-1))) + SEQDICT[read[i]]
        kcount[temp_f] += 1

    #The hash value of forward kmer and reverse kmer add up to (4**k)
    if inputtype == 'A':
        return (kcount + kcount[::-1])
    else:
        return kcount
 
#Generate the kmer count array for all reads
def create_readlist_kmer(readlist, k, inputtype):
    length = 4 ** k
    kcount = np.zeros(length)
    for read in readlist:
        kcount += create_read_kmer(read, k, inputtype)
    return kcount

#!/bin/python

import numpy as np


#Generate a quantary dictionary.
SEQDICT = dict(zip('ATCG', range(4)))

#Generate the kmer count array for one read
def create_read_kmer(read, k):
    length = 4 ** k
    kcount = np.zeros(length)
    temp_f = 0
    temp_r = 0
    for i in range(k):
        temp_f = (temp_f * 4) + SEQDICT[read[i]]
        temp_r = (temp_r * 4) + SEQDICT[read[k-i-1]]
 
    kcount[temp_f] += 1
    kcount[temp_r] += 1

    for i in range(k, len(read)):
        temp_f = 4 * (temp_f - SEQDICT[read[i-k]] * (4**(k-1))) + SEQDICT[read[i]]
        temp_r = (temp_r - SEQDICT[read[i-k]]) / 4 + SEQDICT[read[i]] * (4**(k-1))
        kcount[temp_f] += 1
        kcount[temp_r] += 1
    return kcount
 
#Generate the kmer count array for all reads
def create_readlist_kmer(readlist, k):
    length = 4 ** k
    kcount = np.zeros(length)
    for read in readlist:
        kcount += create_read_kmer(read, k)
    return kcount

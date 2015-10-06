#!/bin/python

from d2compute import *
import numpy as np


#Get a new cluster given the old cluster and the clustering pair in this step.
def newcluster(pair, cluster):
    i, j = pair
    return (cluster[:i] + cluster[i+1:j] + cluster[j+1:] + [cluster[i]+cluster[j]])

#Calculate the average distance
def averagedist(distmatrix, cluster):
    L = len(cluster)
    newdist = np.zeros([L]*2)
    for i in range(L):
        for j in range(i+1, L):
            total = 0
            combination = [(spea, speb) for spea in cluster[i] for speb in cluster[j]] 
            for spea, speb in combination:
                total += distmatrix[spea][speb]
            average = total / len(combination)
            newdist[i][j] = average
            newdist[j][i] = average
    return newdist

#Choose the minimum distance pair to cluster
def iteration(newdist_smallk, newdist_bigk, cluster):
    i_smallk = np.where(newdist_smallk == newdist_smallk[newdist_smallk>0].min())[0][0]
    j_smallk = np.where(newdist_smallk == newdist_smallk[newdist_smallk>0].min())[1][0]
    i_bigk = np.where(newdist_bigk == newdist_bigk[newdist_bigk>0].min())[0][0]
    j_bigk = np.where(newdist_bigk == newdist_bigk[newdist_bigk>0].min())[1][0]
    if (i_smallk == i_bigk and j_smallk == j_bigk):
        return (False, newcluster((i_bigk, j_bigk), cluster))
    else:
        return (True, newcluster((i_smallk, j_smallk), cluster)) 

    
     

def main():
    #Create a parser to read arguments
    parser = optparse.OptionParser()

    #Iteration
    parser.add_option('--iter', action='store_true', dest='iteration', default=False)
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
    cluster = [[i] for i in range(len(species))]
    bigk = kvalue
    smallk = (bigk - 1) if options.iteration else bigk
    for method in d2method:
        if method == 'd2':
            compare = d2compare
        elif method == 'd2star':
            compare = d2starcompare
        else:
            compare = d2sheppcompare
        distmatrix_smallk = pairwise(species, sequences, smallk, inputtype, compare)
        if options.iteration: 
            distmatrix_bigk = pairwise(species, sequences, bigk, inputtype, compare) 
        else:
            distmatrix_bigk = distmatrix_smallk 
        while(len(cluster) > 1):
            print bigk, cluster
            newdist_smallk = averagedist(distmatrix_smallk, cluster)
            newdist_bigk = averagedist(distmatrix_bigk, cluster)
            decreasek, cluster = iteration(newdist_smallk, newdist_bigk, cluster) 
            if decreasek:
                bigk = smallk
                smallk = bigk - 1
                distmatrix_smallk = pairwise(species, sequences, smallk, inputtype, compare)
                distmatrix_bigk = pairwise(species, sequences, bigk, inputtype, compare)
        print bigk, cluster
             
if __name__ == '__main__':
    main()


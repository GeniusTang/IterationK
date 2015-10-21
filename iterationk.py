#!/bin/python

from d2compute import *
import numpy as np
import sys

#Print a list in Newick format
def listprint(cluster, seqdict):
   if isinstance(cluster, list):
       temp = '('
       for element in cluster:
           temp += listprint(element, seqdict)
       return temp+'),'
   else:
       return seqdict[cluster]+',' 
    
#Return a list containing all elements in this cluster.
def elements(cluster):
    if isinstance(cluster, list):
        temp = []
        for element in cluster:
            temp += elements(element)
        return temp
    else:
        return [cluster]
 
#Get a new cluster given the old cluster and the clustering pair in this step.
def newcluster(pair, cluster):
    i, j = pair
    return (cluster[:i] + cluster[i+1:j] + cluster[j+1:] + [[cluster[i], cluster[j]]])

#Calculate the average distance
def averagedist(distmatrix, cluster):
    L = len(cluster)
    newdist = np.zeros([L]*2)
    for i in range(L):
        for j in range(i+1, L):
            total = 0
            combination = [(spea, speb) for spea in elements(cluster[i]) for speb in elements(cluster[j])] 
            for spea, speb in combination:
                total += distmatrix[spea][speb]
            average = total / len(combination)
            newdist[i][j] = average
            newdist[j][i] = average
    return newdist

#Choose the minimum distance pair to cluster
def iteration(newdist_smallk, newdist_bigk, cluster, seqdict):
    i_smallk = np.where(newdist_smallk == newdist_smallk[newdist_smallk>0].min())[0][0]
    j_smallk = np.where(newdist_smallk == newdist_smallk[newdist_smallk>0].min())[1][0]
    i_bigk = np.where(newdist_bigk == newdist_bigk[newdist_bigk>0].min())[0][0]
    j_bigk = np.where(newdist_bigk == newdist_bigk[newdist_bigk>0].min())[1][0]
    print 'cluster', listprint(cluster[i_smallk], seqdict)[:-1].replace(',)', ')')\
,'and', listprint(cluster[j_smallk], seqdict)[:-1].replace(',)', ')')
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
 
    f = open(outputprefix + 'outfile', 'wt')
    sys.stdout = f
    species, sequences = analyzeInput(inputfile)
    seqdict = dict(zip(range(len(species)), species))

    cluster = [i for i in range(len(species))]
    bigk = kvalue
    smallk = max([(bigk - 1) if options.iteration else bigk, 3])
    for method in d2method:
        if method == 'd2':
            compare = d2compare
        elif method == 'd2star':
            compare = d2starcompare
        else:
            compare = d2sheppcompare
        distmatrix_smallk = pairwise(species, sequences, smallk, inputtype, compare)
        if (not options.iteration or smallk == bigk): 
            distmatrix_bigk = distmatrix_smallk 
        else:
            distmatrix_bigk = pairwise(species, sequences, bigk, inputtype, compare) 
        while(len(cluster) > 1):
            print "Cluster result:", listprint(cluster, seqdict).replace(',)', ')')[1:-2]+';'
            print '-'*30
            print 'Use k = '+str(bigk)+',',
            newdist_smallk = averagedist(distmatrix_smallk, cluster)
            newdist_bigk = averagedist(distmatrix_bigk, cluster)
            decreasek, cluster = iteration(newdist_smallk, newdist_bigk, cluster, seqdict) 
            if decreasek:
                #k is at least 3
                bigk = max([smallk, 3])
                smallk = max([bigk - 1, 3])
                distmatrix_smallk = pairwise(species, sequences, smallk, inputtype, compare)
                distmatrix_bigk = pairwise(species, sequences, bigk, inputtype, compare)
        print 'Final cluster:', listprint(cluster, seqdict).replace(',)', ')')[1:-2]+';'
        f.close()
        f = open(outputprefix + 'outtree', 'wt')
        sys.stdout = f
        print listprint(cluster, seqdict).replace(',)', ')')[1:-2]+';'
        f.close()
             
if __name__ == '__main__':
    main()


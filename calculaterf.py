#!/bin/python

RFPATH = "/home/rcf-40/kujin/panasas/IterationK/RF/"
RFVALUEPATH = "/home/rcf-40/kujin/panasas/IterationK/RFVALUE/"
METHOD = 'd2'


for i in range(1,7):
    rffile = open(RFVALUEPATH + 'T' + str(i) + '_' + METHOD, 'wt')
    rffile.write('k,iter,fixed')
    for j in range(3,12):
        f = open(RFPATH + 'T' + str(i) + '_iter_k' + str(j) + '_' + METHOD + '_rf') 
        lines = f.readlines()
        f.close()
        lines = map(lambda x:float(x.strip().split()[-1]), lines)
        iterrf = sum(lines) / len(lines) / 8

        f = open(RFPATH + 'T' + str(i) + '_k' + str(j) + '_' + METHOD + '_rf')
        lines = f.readlines()
        f.close()
        lines = map(lambda x:float(x.strip().split()[-1]), lines)
        rf = sum(lines) / len(lines) / 8 
        rffile.write('\n' + str(j) + ',' + str(iterrf) + ',' + str(rf))

rffile.close() 

        

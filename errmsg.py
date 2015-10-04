#!/bin/python
import sys
import os

INPUTTYPE = ['A', 'L']
METHOD = ['d2', 'd2shepp', 'd2star']

def check_k(k):
    try:
        k = int(k)
        if k<=0:
            raise
    except:
        sys.exit("k should be a positive integer!")
    else:
        return k
 
def check_inputfile(inputfile):
    if not os.path.exists(inputfile):
        sys.exit("Input file %s doesn't exist!"%inputfile)

def check_inputtype(inputtype):
    if not inputtype in INPUTTYPE:
        sys.exit("Please enter the correct input data type: only A or L.")

def check_d2method(d2method):
    method = map(lambda x:x.strip(), d2method.split(','))
    for i in method:
        if not i in METHOD:
            sys.exit("No method named %s!"%i)
    return method

def check_outputdir(outputprefix):
    dirname = os.path.dirname(outputprefix)
    if dirname and not os.path.isdir(dirname):
        sys.exit(dirname + " is not a directory!")

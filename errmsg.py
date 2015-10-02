#!/bin/python
import sys
import os

INPUTTYPE = ['A', 'L']

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


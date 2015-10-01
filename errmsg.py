#!/bin/python
import sys

INPUTTYPE = ['A', 'L']
def check_inputtype(inputtype):
    if not inputtype in INPUTTYPE:
        sys.exit("Please enter the correct input data type: only A or L.")


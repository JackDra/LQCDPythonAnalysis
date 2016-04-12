#!/usr/bin/env python

from Params import *
from SetLists import *

def InputGammaAndSet(inputparams):
    feedsetlist = DefSetList
    feedgammalist = ''
    for isys in inputparams:
        if isys[0] != '-':
            raise IOError("input arguments are specified with -")
        elif '-h' in isys:
            print 'commands are (with comma separated lists):'
            print '-g= specifies gamma matricies (see OppFuns.py)'
            print '-s= specifies set list to use (see SetLists.py)'
            exit()
        elif '-g' in isys:
            feedgammalist = isys.replace('-g=','').split(',')
        elif '-s' in isys:
            feedsetlist = isys.replace('-s=','').split(',')
        return feedgammalist,feedsetlist

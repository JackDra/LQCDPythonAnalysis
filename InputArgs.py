#!/usr/bin/env python

from Params import *
from SetLists import *

def InputGammaAndSet(inputparams):
    feedsetlist = DefSetList
    feedgammalist = CreateGammaList()
    for isys in inputparams:
        if isys[0] != '-':
            raise IOError("input arguments are specified with -")
        elif '-h' in isys:
            print 'commands are (with comma separated lists):'
            print '-gamma= or -g= specifies gamma matricies (see gammaopp file)'
            print '-set= or -s= specifies set list to use'
            exit()
        elif '-g' in isys:
            feedgammalist = isys.replace('-g=','').split(',')
        elif '-s' in isys:
            feedsetlist = isys.replace('-s=','').split(',')
        return feedgammalist,feedsetlist

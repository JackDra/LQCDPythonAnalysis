#!/usr/bin/env python

import xmltodict
from collections import OrderedDict


def RecFTDAS(dictin):
    if isinstance(dictin,dict) or isinstance(dictin,OrderedDict):
        for ikey,idata in dictin.iteritems():
            dictin[ikey] = RecFTDAS(idata)
        return dictin
    else:
        print dictin
        try:
            dictout = map(float,dictin)
        except:
            try:
                dictout = FormatToDictAvgStd(str(dictin))
            except:
                raise TypeError('final value in dictionary is not string')
        return dictout


def RecDictToKeys(dictin):
    if isinstance(dictin,dict) or isinstance(dictin,OrderedDict):
        RecKeys = RecDictToKeys(dictin[dictin.keys()[0]])
        return [dictin.keys()] + RecKeys
    else:
        return []


def RecDictToNp(dictin):
    if isinstance(dictin,dict) or isinstance(dictin,OrderedDict):
        dataout = []
        for idict in dictin.itervalues():
            dataout += [RecDictToNp(idict)]
        return dataout
    else:
        return dictin

    ## Keys going down the dictionary chane , total np array with correct dimension
    ## Caution, no checking to see if dimensions are same size.
def RecDictToNpWrap(dictin):
    return RecDictToKeys(dictin),RecDictToNp(dictin)

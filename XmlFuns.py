#!/usr/bin/env python

import xmltodict
from collections import OrderedDict
from XmlFormatting import *


def RecFTDAS(dictin):
    if isinstance(dictin,dict) or isinstance(dictin,OrderedDict):
        for ikey,idata in dictin.iteritems():
            dictin[ikey] = RecFTDAS(idata)
        return dictin
    else:
        try:
            dictout = map(float,dictin)
        except:
            try:
                if len(str(dictin).split()) == 3:
                    dictout = FormatToDictAvgStdChi(str(dictin))
                elif len(str(dictin).split()) == 2:
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

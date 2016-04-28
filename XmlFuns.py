#!/usr/bin/env python

import xmltodict
from collections import OrderedDict

def AvgStdToFormat(Avg,Std,frmtflag='f'):
    return ('{0:20.10'+frmtflag+'} {1:20.10'+frmtflag+'}').format(Avg,Std)

def DictAvgStdToFormat(Dict,frmtflag='f'):
    return AvgStdToFormat(Dict['Avg'],Dict['Std'],frmtflag=frmtflag)

def BootAvgStdToFormat(Dict,frmtflag='f'):
    return AvgStdToFormat(Dict.Avg,Dict.Std,frmtflag=frmtflag)

def FormatToAvgStd(String):
    return String.strip().split()

def FormatToDictAvgStd(String):
    return {zip(['Avg','Std'],FormatToAvgStd(String))}


def RecursiveFTDAS(dictin):
    if isinstance(dictin,dict) or isinstance(dictin,OrderedDict):
        for ikey,idata in dictin.iteritems():
            dictin[ikey] = RecursiveFTDAS(idata)
        return dictin
    else:
        if isinstance(dictin,str):
            return FormatToDictAvgStr(dictin)
        else:
            raise TypeError('final value in dictionary is not string')

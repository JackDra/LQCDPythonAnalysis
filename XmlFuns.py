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
    return map(float,String.strip().split())

def FormatToDictAvgStd(String):
    return dict(zip(['Avg','Std'],FormatToAvgStd(String)))


def RecursiveFTDAS(dictin):
    if isinstance(dictin,dict) or isinstance(dictin,OrderedDict):
        for ikey,idata in dictin.iteritems():
            dictin[ikey] = RecursiveFTDAS(idata)
        return dictin
    else:
        try:
            dictout = map(float,dictin)
        except:
            try:
                dictout = FormatToDictAvgStd(str(dictin))
            except:
                raise TypeError('final value in dictionary is not string')
        return dictout


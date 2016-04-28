#!/usr/bin/env python

import xmltodict
from collections import OrderedDict
from Params import *

def tstr(it):
    return 't'+str(it)

def untstr(it):
    return int(it.replace('t',''))

def xmlcut(icut):
    return 'cut'+str(icut)

def unxmlstr(icut):
    return int(icut.replace('cut',''))


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



def AvgStdChiToFormat(Avg,Std,Chi):
    return ('{0:20.10f} {1:20.10f} {2:20.10f}').format(Avg,Std,Chi)

def DictAvgStdChiToFormat(Dict,Chi):
    return AvgStdChiToFormat(Dict['Avg'],Dict['Std'],Chi)

def BootAvgStdChiToFormat(Dict,Chi):
    return AvgStdChiToFormat(Dict.Avg,Dict.Std,Chi)

def FormatToAvgStdChi(String):
    return map(float,String.strip().split())

def FormatToDictAvgStdChi(String):
    return dict(zip(['Avg','Std','Chi'],FormatToAvgStdChi(String)))


def LREVecToFormat(iLE,iRE,iEM,DoPoF):
    if DoPoF:
        SmPoFList = []
        for iPoF in ['PoF0','PoF1']:
            SmPoFList += [ism+'_'+iPoF for ism in DefSmList]
    else:
        SmPoFList = DefSmList
    dictout = OrderedDict()
    dictout['Emass'] = '{0:20.10f}'.format(iEMass)
    dictout['Left_Evec'] = OrderedDict()
    dictout['Right_Evec'] = OrderedDict()
    for eLE,thiskey in zip(iLE,SmPoFList):
        dictout['Left_Evec'][thiskey] = '{0:20.10f}'.format(eLE)
    for eRE,thiskey in zip(iRE,SmPoFList):
        dictout['Left_Evec'][thiskey] = '{0:20.10f}'.format(eRE)
    return dictout
        
    

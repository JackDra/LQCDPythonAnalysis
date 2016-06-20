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

def unxmlcut(icut):
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



def AvgStdChiToFormat(Avg,Std,Chi,frmtflag='f'):
    return ('{0:20.10'+frmtflag+'} {1:20.10'+frmtflag+'} {2:20.10f}').format(Avg,Std,Chi)

def DictAvgStdChiToFormat(Dict,Chi,frmtflag='f'):
    return AvgStdChiToFormat(Dict['Avg'],Dict['Std'],Chi,frmtflag=frmtflag)

def BootAvgStdChiToFormat(Dict,Chi,frmtflag='f'):
    return AvgStdChiToFormat(Dict.Avg,Dict.Std,Chi,frmtflag=frmtflag)

def FormatToAvgStdChi(String):
    return map(float,String.strip().split())

def FormatToDictAvgStdChi(String):
    return dict(zip(['Avg','Std','Chi'],FormatToAvgStdChi(String)))



def LREVecToFormat(iLE,iRE,iEM,DoPoF):
    if DoPoF:
        SmPoFList = []
        for iPoF in ['PoF'+str(ishift) for ishift in PoFShifts]:
            SmPoFList += [ism+'_'+iPoF for ism in DefSmList]
    else:
        SmPoFList = DefSmList
    dictout = OrderedDict()
    dictout['Emass'] = '{0:20.10f}'.format(iEM)
    dictout['Left_Evec'] = OrderedDict()
    dictout['Right_Evec'] = OrderedDict()
    for eLE,thiskey in zip(iLE,SmPoFList):
        dictout['Left_Evec'][thiskey] = '{0:20.10f}'.format(eLE)
    for eRE,thiskey in zip(iRE,SmPoFList):
        dictout['Right_Evec'][thiskey] = '{0:20.10f}'.format(eRE)
    return dictout
        
    
def xmlfitr(thefitr):
    return 'r'+'-'.join(map(str,thefitr))

def unxmlfitr(thefitr):
    return thefitr.replace('r','').split('-')
    
def xmlTSink(thestr):
    return 'tsink'+str(thestr)

def unxmlTSink(thestr):
    return int(thestr.replace('tsink',''))

def ParamsToFitFlag(theFitList):
    listout = []
    for icut,cutfitlist in enumerate(theFitList):    
        listout.append([])
        for ifit in cutfitlist:
            listout[icut].append(xmlfitr(ifit))
    return listout


def FitFlagXmlToOld(param,fitr):
    oldparam = param.replace('slope','sl').replace('constant','con')
    oldfitr = unxmlfitr(fitr)
    return 'fit '+oldparam+' '+oldfitr[0]+'-'+oldfitr[1]

def FitFlagXmlToOldSF(fitr):
    oldfitr = unxmlfitr(fitr)
    return oldfitr[0]+'-'+oldfitr[1]

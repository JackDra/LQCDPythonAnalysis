#!/usr/bin/env python

from Params import *
import xmltodict
from collections import OrderedDict
import numpy as np

def GetInfoFromFilelist(filedict):
    dictout = OrderedDict()
    dictout['nconfig'] = len(filedict.keys())
    lenths = []
    for ikey,cfglist in filedict.iteritems():
        lenths.append(len(cfglist))
    dictout['nsrcAvg'] = np.mean(lenths)
    dictout['nMeasure'] = np.sum(lenths)
    return dictout

def CreateDefaultInfodict():
    dictout = OrderedDict()
    dictout['nconfig'] = -1
    dictout['nsrcAvg'] = -1
    dictout['nMeasure'] = -1
    return dictout
    
    

def tflowstr(it):
    return 't_flow'+str(it)

def untflowstr(it):
    return float(it.replace('t_flow',''))


def tstr(it):
    return 't'+str(it)

def untstr(it):
    return int(it.replace('t',''))

def xmlcut(icut):
    if isinstance( icut, int ):
        return 'cut'+str(icut)
    else:
        return 'cut'+'-'.join(map(str,icut))
    
def unxmlcut(icut):
    if '-' in icut:
        return map(int,icut.replace('cut','').split('-'))
    else:
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



def MomToFormat(mom,frmtflag='f'):
    thisfmt = ' '.join(['{'+str(imom)+':20.10'+frmtflag+'} ' for imom in xrange(len(mom))])
    return (thisfmt).format(*mom)

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
        jSmPoFList = []
        iSmPoFList = []
        for iPoF in ['PoF'+str(ishift) for ishift in xrange(PoFShifts+1)]:
            iSmPoFList += [ism+'_'+iPoF for ism in DefiSmList]
            jSmPoFList += [jsm+'_'+iPoF for jsm in DefjSmList]
    else:
        iSmPoFList = DefiSmList
        jSmPoFList = DefjSmList
    dictout = OrderedDict()
    dictout['Emass'] = '{0:20.10f}'.format(iEM)
    dictout['Source_Evec'] = OrderedDict()
    dictout['Sink_Evec'] = OrderedDict()
    for eLE,thiskey in zip(iLE,jSmPoFList):
        dictout['Source_Evec'][thiskey] = '{0:20.10f}'.format(eLE)
    for eRE,thiskey in zip(iRE,iSmPoFList):
        dictout['Sink_Evec'][thiskey] = '{0:20.10f}'.format(eRE)
    return dictout
        
    
def xmlfitr(thefitr):
    return 'r'+'-'.join(map(str,thefitr))

def unxmlfitr(thefitr):
    if 'fitr' in thefitr:
        return thefitr.replace('fitr','').split('-')
    else:
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

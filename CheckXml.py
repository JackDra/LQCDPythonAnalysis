#!/usr/bin/env python

from Params import *
from Params import outputdir
from XmlFuns import *
from XmlFormatting import *
from FitParams import *
from OppFuns import *
from ReadXml import CheckMomFile,CheckNconfFile
from OutputXmlData import SetUpPDict
from SetLists import *
from MiscFuns import *
import time

def CheckNconf(thisGammaList,CheckSetList,thisMomList=RunMomList,CheckList=[''],cfuns=False):
    nconf = -2
    for CheckType in CheckList:
        if len(CheckType) > 0:
            CheckType += '/'
            if any([itype in CheckType for itype in ['SumMeth','TSF']]): CheckSetList = ReduceTsink(thisSetList)
            if cfuns: outputdir = outputdir + 'cfuns/'
        SFList = ['']
        if 'OSF' in CheckType:
            SFList = OneStateParList['C3']
        if 'TSF' in CheckType:
            SFList = TwoStateParList['C3']

        for iset in CheckSetList:
            for iSF in SFList:
                for igamma in thisGammaList:
                    gammadir = outputdir+CreateOppDir(igamma)+'/' + CheckType
                    for pstr in GetMomFromGamma(igamma,thisMomList=thisMomList):
                        ip = qstrTOqcond(pstr)
                        filename = iset+igamma+ iSF
                        dump,checkfile = SetUpPDict(ip,gammadir,filename)
                        nconf = min(nconf, CheckNconfFile(checkfile+'.xml'))
    return nconf
    

##Also works for cfuns##
##xmlinput = { Ratio_Factor , Boots/Values , thismomlist , tlist } 
##outputdict = { thismom , [tVals] / [Vals] / [Valserr] / [Boot] bs }
# def Check3ptFiles(thisGammaList,thisSetList,thisMomList,CheckType='',cfuns=False):
#     CheckSetList,thisdir = thisSetList,outputdir
#     outputbool = True
#     if len(CheckType) > 0:
#         CheckType += '/'
#         if any([itype in CheckType for itype in ['SumMeth','TSF']]): CheckSetList = ReduceTsink(thisSetList)
#         if cfuns: thisdir = outputdir + 'cfuns/'
#     SFList = ['']
#     if 'OSF' in CheckType:
#         SFList = OneStateParList['C3']
#     if 'TSF' in CheckType:
#         SFList = TwoStateParList['C3']

#     for iset in CheckSetList:
#         for iSF in SFList:
#             for igamma in thisGammaList:
#                 gammadir = thisdir+CreateOppDir(igamma)+'/' + CheckType
#                 for pstr in GetMomFromGamma(igamma,thisMomList=thisMomList):
#                     ip = qstrTOqcond(pstr)
#                     filename = iset+igamma+ iSF
#                     dump,checkfile = SetUpPDict(ip,gammadir,filename)
#                     outputbool = outputbool and CheckMomFile(checkfile+'.xml')
#     return outputbool
                

## list of booleans corresponding to what needs to be done relative to list thisMomList
def Check3ptArray(thisGammaList,thisSetList,thisMomList=RunMomList,CheckType='',cfuns=False,printout=True):
    CheckSetList,thisdir = thisSetList,outputdir
    totstart = time.time()
    outlist = {}
    if len(CheckType) > 0:
        CheckType += '/'
        if any([itype in CheckType for itype in ['SumMeth','TSF']]): CheckSetList = ReduceTsink(thisSetList)
        if cfuns: thisdir = outputdir + 'cfuns/'
    SFList = ['']
    if 'OSF' in CheckType:
        SFList = OneStateParList['C3']
    if 'TSF' in CheckType:
        SFList = TwoStateParList['C3']
    for icg,igamma in enumerate(thisGammaList):
        if printout: print 'Checking: ' , GetTimeForm(icg,len(thisGammaList),time.time()-totstart) , igamma ,  '          \r',
        outlist[igamma] = {}
        gammadir = thisdir+CreateOppDir(igamma)+'/' + CheckType
        for iset in CheckSetList:
            outlist[igamma][iset] = []
            for pstr in GetMomFromGamma(igamma,thisMomList=thisMomList):
                ip = qstrTOqcond(pstr)
                filename = iset+igamma
                if not all([CheckMomFile(SetUpPDict(ip,gammadir,filename+iSF)[1]+'.xml') for iSF in SFList]):
                    outlist[igamma][iset].append(pstr)
    if len(thisGammaList) < 5:
        if printout: print 'Checking complete, Total Time: ' , GetTimeStr(time.time()-totstart) +thisGammaList[0].replace('doub','').replace('sing','')+ ' '*40
    else:
        if printout: print 'Checking complete, Total Time: ' , GetTimeStr(time.time()-totstart) + ' '*40
    return outlist
                

## list of booleans corresponding to what needs to be done relative to list thisMomList
def Check3ptAllSets(thisGammaList,thisSetList,thisMomList=RunMomList,CheckType='',cfuns=False,printout=True):
    outlist = Check3ptArray(thisGammaList,thisSetList,thisMomList=thisMomList,CheckType=CheckType,cfuns=cfuns,printout = printout)
    CheckSetList = thisSetList
    outnoset = {}
    for igamma in thisGammaList:
        outnoset[igamma] = []
        for setlist in outlist[igamma].itervalues():
            for ip in setlist:
                outnoset[igamma].append(ip)
        outnoset[igamma] = OrderMomList(outnoset[igamma])
    return outnoset
            
        

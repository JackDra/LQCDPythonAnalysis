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
from collections import OrderedDict as OD
import time

def CheckNconf(thisGammaList,CheckSetList,thisMomList=RunMomList,CheckList=[''],cfuns=False,minmax='min'):
    nconf = 10e16
    thisdir = outputdir
    NconfDict = OD()
    NconfDict['Missing'] = []
    NconfDict['Dep'] = []
    for CheckType in CheckList:
        # print 'Checking' , CheckType
        thisSetList = CheckSetList
        if 'RF' == CheckType: CheckType = ''
        if len(CheckType) > 0:
            CheckType += '/'
            if any([itype in CheckType for itype in ['SumMeth','TSF']]): thisSetList = ReduceTsink(CheckSetList)
            if cfuns: thisdir = outputdir + 'cfuns/'
        SFList = ['']
        if 'OSF' in CheckType:
            SFList = OneStateParList['C3']
        if 'TSF' in CheckType:
            SFList = TwoStateParList['C3']

        for iset in thisSetList:
            # print '    Checking', iset
            for iSF in SFList:
                # if len(iSF) > 0: print '       Checking', iSF , ' '*50
                for igamma in thisGammaList:
                    if 'doub' not in igamma and 'sing' not in igamma and CheckType == '': continue
                    print 'Checking', iSF , igamma, ' '*50 ,' \r',
                    gammadir = thisdir+CreateOppDir(igamma)+'/' + CheckType
                    for pstr in GetMomFromGamma(igamma,thisMomList=thisMomList):
                        ip = qstrTOqcond(pstr)
                        filename = iset+igamma+ iSF
                        dump,checkfile = SetUpPDict(ip,gammadir,filename)
                        thisnconf = CheckNconfFile(checkfile+'.xml')
                        if 'File Missing' == thisnconf:
                            if not any(igamma in inconf for inconf in NconfDict['Missing']):
                                NconfDict['Missing'].append(igamma+' '+qstrTOqcond(pstr))
                            # return 'File Missing: ' + checkfile+'.xml' , NconfDict
                        elif 'depreciated' == thisnconf:
                            if not any(igamma in inconf for inconf in NconfDict['Dep']):
                                NconfDict['Dep'].append(igamma+' '+qstrTOqcond(pstr))
                        else:
                            thiskey = 'nconf'+str(thisnconf)
                            if thiskey not in NconfDict:
                                NconfDict[thiskey] = []
                            if not any(igamma in inconf for inconf in NconfDict[thiskey]):
                                NconfDict[thiskey].append(igamma+' '+qstrTOqcond(pstr))
                            # print ''
                            # print 'Changed nconfigs from ',nconf,' to ',thisnconf , ' in file:'
                            # print checkfile+'.xml'
                            # print ''
                            # elif thisnconf > nconf and thisnconf < 1e10:
                            #     print ''
                            #     print 'Larger nconfigs from ',nconf,' compared to ',thisnconf , ' in file:'
                            #     print checkfile+'.xml'
                            #     print ''
                            if nconf > 10e10 :
                                nconf = thisnconf
                            elif thisnconf > 0:
                                if minmax == 'min':
                                    nconf = min(nconf, thisnconf)
                                elif minmax == 'max':
                                    nconf = max(nconf, thisnconf)
    print ' '*50
    if nconf > 10e10: nconf = 'All files used depreciated code'
    return nconf,NconfDict
    

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
def Check3ptArray(thisGammaList,thisSetList,thisMomList=RunMomList,CheckType='',cfuns=False,printout=True,thisNconf = False):
    CheckSetList,thisdir = thisSetList,outputdir
    totstart = time.time()
    outlist = {}
    if len(CheckType) > 0:
        CheckType += '/'
        # if any([itype in CheckType for itype in ['SumMeth','TSF']]): CheckSetList = ReduceTsink(thisSetList)
        if cfuns: thisdir = outputdir + 'cfuns/'
    SFList = ['']
    if 'OSF' in CheckType:
        SFList = OneStateParList['C3']
    elif 'TSF' in CheckType:
        SFList = TwoStateParList['C3']
    for icg,igamma in enumerate(thisGammaList):
        if printout: print 'Checking: ' , GetTimeForm(icg,len(thisGammaList),time.time()-totstart) , igamma ,  '          \r',
        outlist[igamma] = {}
        for iset in CheckSetList:
            outlist[igamma][iset] = []
            for pstr in GetMomFromGamma(igamma,thisMomList=thisMomList):
                ip = qstrTOqcond(pstr)
                if ('doub' not in igamma) and ('sing' not in igamma):
                    CheckBool = True
                    for dsgamma in ['doub'+igamma,'sing'+igamma]:
                        filename = iset+dsgamma
                        gammadir = thisdir+CreateOppDir(dsgamma)+'/' + CheckType
                        CheckBool = CheckBool and any([CheckMomFile(SetUpPDict(ip,gammadir,filename+iSF)[1]+'.xml',nconftest=thisNconf) for iSF in SFList])
                else:
                    filename = iset+igamma
                    gammadir = thisdir+CreateOppDir(igamma)+'/' + CheckType
                    CheckBool = all([CheckMomFile(SetUpPDict(ip,gammadir,filename+iSF)[1]+'.xml',nconftest=thisNconf) for iSF in SFList])
                if not CheckBool:
                    outlist[igamma][iset].append(pstr)
    if len(thisGammaList) < 5:
        if printout: print 'Checking complete, Total Time: ' , GetTimeStr(time.time()-totstart) +thisGammaList[0].replace('doub','').replace('sing','')+ ' '*40
    else:
        if printout: print 'Checking complete, Total Time: ' , GetTimeStr(time.time()-totstart) + ' '*40
    return outlist
                

## list of booleans corresponding to what needs to be done relative to list thisMomList
def Check3ptAllSets(thisGammaList,thisSetList,thisMomList=RunMomList,CheckType='',cfuns=False,printout=True,thisNconf = False):
    outlist = Check3ptArray(thisGammaList,thisSetList,thisMomList=thisMomList,CheckType=CheckType,cfuns=cfuns,printout = printout,thisNconf=thisNconf)
    CheckSetList = thisSetList
    outnoset = {}
    for igamma in thisGammaList:
        outnoset[igamma] = []
        for setlist in outlist[igamma].itervalues():
            for ip in setlist:
                outnoset[igamma].append(ip)
        outnoset[igamma] = OrderMomList(outnoset[igamma])
    return outnoset
            
        

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

def CheckNconfMass(CheckSetList,thisMomList=RunMomList,CheckList=[''],cfuns=True,minmax='mim'):
    nconf = 10e16
    thisdir = outputdir
    NconfDict = OD()
    NconfDict['Missing'] = []
    NconfDict['Dep'] = []
    for CheckType in CheckList:
        # print 'Checking' , CheckType
        thisSetList = ReduceTsink(CheckSetList)
        if 'RF' == CheckType: CheckType = ''
        if len(CheckType) > 0: CheckType += '/'
        if cfuns: thisdir = outputdir + 'cfuns/'
        SFList = ['']
        if 'OSF' in CheckType:
            SFList = OneStateParList['C2']
        if 'TSF' in CheckType:
            SFList = TwoStateParList['C2']

        existsDep = False
        for iset in thisSetList:
            # print '    Checking', iset
            for iSF in SFList:
                # if len(iSF) > 0: print '       Checking', iSF , ' '*50
                print 'Checking', iSF , ' '*50 ,' \r',
                if cfuns:
                    twoptdir = thisdir+CreateOppDir('twopt')+'/' + CheckType
                else:
                    twoptdir = thisdir+CreateOppDir('Mass')+'/' + CheckType
                    
                for pstr in thisMomList:
                    ip = qstrTOqcond(pstr)
                    if cfuns:
                        if 'TSF' in CheckType and ('m0' in iSF or 'Dm' in iSF):
                            filename = 'twopt' + iSF
                        else:
                            filename = iset+ 'twopt' + iSF
                    else:
                        if 'TSF' in CheckType and ('m0' in iSF or 'Dm' in iSF):
                            filename = 'Mass' + iSF
                        else:
                            filename = iset+ 'Mass' + iSF
                    dump,checkfile = SetUpPDict(ip,twoptdir,filename)
                    thisnconf = CheckNconfFile(checkfile+'.xml')
                    if 'File Missing' == thisnconf:
                        NconfDict['Missing'].append('twopt '+ip)
                    elif 'Dep' == thisnconf:
                        existsDep = True
                        if not any(keygamma in inconf for inconf in NconfDict['Dep']):
                            NconfDict['Dep'].append('twopt '+ip)
                    else:
                        thiskey = 'nconf'+str(thisnconf)
                        if thiskey not in NconfDict: NconfDict[thiskey] = []
                        if not any('twopt' in inconf for inconf in NconfDict[thiskey]):
                            NconfDict[thiskey].append('twopt '+ip)
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
    if existsDep:
        nconf = 'Depreciated code results'
    elif nconf > 10e10:
        nconf = 'No Files Found'
    return nconf,NconfDict
    

def CheckNconf(inputGammaList,CheckSetList,thisMomList=RunMomList,CheckList=[''],cfuns=False,minmax='min'):
    nconf = 10e16
    thisdir = outputdir
    NconfDict = OD()
    NconfDict['Missing'] = []
    NconfDict['Dep'] = []

    if ('twopt' in inputGammaList or 'Mass' in inputGammaList) and not any(icheck in CheckList for icheck in ['RF','Fits','SumMeth']):
        massnconf,MassNconfDict = CheckNconfMass(CheckSetList,thisMomList=thisMomList,CheckList=CheckList,minmax=minmax)
    else:
        massnconf = False
    thisGammaList = list(inputGammaList)
    if 'twopt' in thisGammaList: thisGammaList.remove('twopt')
    if 'Mass' in thisGammaList: thisGammaList.remove('Mass')
        
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

        existsDep = False
        for iset in thisSetList:
            # print '    Checking', iset
            for iSF in SFList:
                # if len(iSF) > 0: print '       Checking', iSF , ' '*50
                for igamma in thisGammaList:
                    if ('doub' not in igamma) and ('sing' not in igamma) and len(CheckType) == 0: continue
                    keygamma = igamma
                    # keygamma = igamma.replace('doub','').replace('sing','')
                    print 'Checking', iSF , igamma, ' '*50 ,' \r',
                    gammadir = thisdir+CreateOppDir(igamma)+'/' + CheckType
                    for pstr in GetMomFromGamma(igamma,thisMomList=thisMomList):
                        ip = qstrTOqcond(pstr)
                        filename = iset+igamma+ iSF
                        dump,checkfile = SetUpPDict(ip,gammadir,filename)
                        thisnconf = CheckNconfFile(checkfile+'.xml')
                        if 'File Missing' == thisnconf:
                            if not any(keygamma in inconf for inconf in NconfDict['Missing']):
                                NconfDict['Missing'].append(keygamma+' '+ip)
                            # return 'File Missing: ' + checkfile+'.xml' , NconfDict
                        elif 'Dep' == thisnconf:
                            existsDep = True
                            if not any(keygamma in inconf for inconf in NconfDict['Dep']):
                                NconfDict['Dep'].append(keygamma+' '+ip)
                        else:
                            thiskey = 'nconf'+str(thisnconf)
                            if thiskey not in NconfDict:
                                NconfDict[thiskey] = []
                            if not any(keygamma in inconf for inconf in NconfDict[thiskey]):
                                NconfDict[thiskey].append(keygamma+' '+ip)
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
    if massnconf != False:
        if minmax == 'min':
            nconf = min(nconf, massnconf)
        elif minmax == 'max':
            if len(thisGammaList) == 0: nconf = -1
            nconf = max(nconf, massnconf)
        for imasskey in MassNconfDict.iterkeys():
            if imasskey in NconfDict.keys():
                NconfDict[imasskey] += MassNconfDict[imasskey]
            else:
                NconfDict[imasskey] = MassNconfDict[imasskey]
    if existsDep:
        nconf = 'Depreciated code results'
    elif nconf > 10e10:
        nconf = 'No Files Found'
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
def Check3ptArray(thisGammaList,thisSetList,thisMomList=RunMomList,CheckType='',cfuns=False,printout=True,thisNconf = RunNconfs):
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
    elif 'TSF' in CheckType:
        SFList = TwoStateParList['C3']
    for icg,igamma in enumerate(thisGammaList):
        keygamma = igamma.replace('doub','').replace('sing','')
        # keygamma = igamma
        if printout: print 'Checking: ' , GetTimeForm(icg,len(thisGammaList),time.time()-totstart) , igamma ,  '          \r',
        if keygamma not in outlist.keys(): outlist[keygamma] = {}
        for iset in CheckSetList:
            if iset not in outlist[keygamma].keys(): outlist[keygamma][iset] = []
            for pstr in GetMomFromGamma(igamma,thisMomList=thisMomList):
                ip = qstrTOqcond(pstr)
                filename = iset+igamma
                gammadir = thisdir+CreateOppDir(igamma)+'/' + CheckType
                CheckBool = all([CheckMomFile(SetUpPDict(ip,gammadir,filename+iSF)[1]+'.xml',nconftest=thisNconf) for iSF in SFList])
                if not CheckBool and pstr not in outlist[keygamma][iset]:
                    outlist[keygamma][iset].append(pstr)
    if len(thisGammaList) < 5:
        if printout: print 'Checking complete, Total Time: ' , GetTimeStr(time.time()-totstart) +thisGammaList[0].replace('doub','').replace('sing','')+ ' '*40
    else:
        if printout: print 'Checking complete, Total Time: ' , GetTimeStr(time.time()-totstart) + ' '*40
    return outlist
                

## list of booleans corresponding to what needs to be done relative to list thisMomList
def Check3ptAllSets(thisGammaList,thisSetList,thisMomList=RunMomList,CheckType='',cfuns=False,printout=True,thisNconf = RunNconfs):
    outlist = Check3ptArray(thisGammaList,thisSetList,thisMomList=thisMomList,CheckType=CheckType,cfuns=cfuns,printout = printout,thisNconf=thisNconf)
    CheckSetList = thisSetList
    outnoset = {}
    for igamma in outlist.keys():
        outnoset[igamma] = []
        for setlist in outlist[igamma].itervalues():
            for ip in setlist:
                outnoset[igamma].append(ip)
        outnoset[igamma] = OrderMomList(outnoset[igamma])

    return outnoset
            
        

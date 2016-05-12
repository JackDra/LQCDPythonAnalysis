#!/usr/bin/env python

import ReadXml
import xmltodict
from XmlFuns import *
from XmlFormatting import *
import os
from BootTest import BootStrap1
from Params import *
from FitParams import *
import cPickle as pickle
from OppFuns import *
from OutputXmlData import *
from SetLists import *


##Also works for cfuns##
##xmlinput = { Ratio_Factor , Boots/Values , thismomlist , tlist } 
##outputdict = { thismom , [tVals] / [Vals] / [Valserr] / [Boot] bs }
def Check3ptFiles(thisGammaList,thisSetList,thisMomList,CheckType='',cfuns=False):
    CheckSetList,thisdir = thisSetList,outputdir
    outputbool = True
    if len(CheckType) > 0:
        CheckType += '/'
        if any([itype in CheckType for itype in ['SumMeth','TSF']]): CheckSetList = ReduceTsink(thisSetList)
        if cfuns: thisdir = outputdir + 'cfuns/'
    SFList = ['']
    if 'OSF' in CheckType:
        SFList = OneStateParList['C3']
    if 'TSF' in CheckType:
        SFList = TwoStateParList['C3']

    for iset in CheckSetList:
        for iSF in SFList:
            for igamma in thisGammaList:
                gammadir = thisdir+CreateOppDir(igamma)+'/' + CheckType
                for pstr in GetMomFromGamma(igamma,thisMomList=thisMomList):
                    ip = qstrTOqcond(pstr)
                    filename = iset+igamma+ iSF
                    dump,checkfile = SetUpPDict(ip,gammadir,filename)
                    outputbool = outputbool and CheckMomFile(checkfile+'.xml')
    return outputbool
                

## list of booleans corresponding to what needs to be done relative to list thisMomList
def Check3ptArray(thisGammaList,thisSetList,thisMomList=RunMomList,CheckType='',cfuns=False):
    CheckSetList,thisdir = thisSetList,outputdir
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
    for igamma in thisGammaList:
        outlist[igamma] = {}
        gammadir = thisdir+CreateOppDir(igamma)+'/' + CheckType
        for iset in CheckSetList:
            outlist[igamma][iset] = []
            for pstr in GetMomFromGamma(igamma,thisMomList=thisMomList):
                ip = qstrTOqcond(pstr)
                filename = iset+igamma
                dump,checkfile = SetUpPDict(ip,gammadir,filename)
                print checkfile
                if not all([CheckMomFile(checkfile.replace(ip,iSF+ip)+'.xml') for iSF in SFList]):
                    outlist[igamma][iset].append(pstr)
    return outlist
                

## list of booleans corresponding to what needs to be done relative to list thisMomList
def Check3ptAllSets(thisGammaList,thisSetList,thisMomList=RunMomList,CheckType='',cfuns=False):
    outlist = Check3ptArray(thisGammaList,thisSetList,thisMomList=thisMomList,CheckType=CheckType,cfuns=cfuns)
    CheckSetList = thisSetList
    outnoset = {}
    for igamma in thisGammaList:
        outnoset[igamma] = []
        checkSet = outlist[igamma].keys()[0]
        for ip in outlist[igamma][checkSet]:
            if all([ip in outlist[igamma][iset] for iset in outlist[igamma].keys()]):
                outnoset[igamma].append(ip)
    return outnoset
            
        

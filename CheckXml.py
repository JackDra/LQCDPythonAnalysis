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
    outputbool = False
    if len(CheckType) > 0:
        CheckType += '/'
        if any([itype in CheckType for itype in ['SumMeth','TSF']]): CheckSetList = ReduceTsink(thisSetList)
        if cfuns: thisdir = outputdir + 'cfuns/'
    for igamma in thisGammaList:
        gammadir = thisdir+CreateOppDir(igamma)+'/' + CheckType
        for pstr in GetMomFromGamma(igamma):
            ip = qstrTOqcond(pstr)
            for iset in CheckSetList:
                SFList = ['']
                if 'OSF' in CheckType:
                    SFList = OneStateParList['C3']
                if 'TSF' in CheckType:
                    SFList = TwoStateParList['C3']
                for iSF in SFList:
                    filename = iset+igamma+ iSF
                    dump,checkfile = SetUpPDict(ip,gammadir,filename)
                    print checkfile
                    outputbool = outputbool or not CheckMomFile(checkfile+'.xml')
    return outputbool
                

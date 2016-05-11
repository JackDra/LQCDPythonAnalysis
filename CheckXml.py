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


##Also works for cfuns##
##xmlinput = { Ratio_Factor , Boots/Values , thismomlist , tlist } 
##outputdict = { thismom , [tVals] / [Vals] / [Valserr] / [Boot] bs }
def CheckFitFiles(thisGammaList,thisSetList,thisMomList):
    outputbool = False
    xmlMomList = map(qstrTOqcond,thisMomList)
    for igamma in thisGammaList:
        gammadir = outputdir+CreateOppDir(igamma)+'/Fits/'
        for ip in xmlMomList:
            for iset in thisSetList:
                filename = iset+igamma
                dump,checkfile = SetUpPDict(ip,gammadir,filename)
                print checkfile+'.xml'
                outputbool = outputbool or not CheckMomFile(checkfile+'.xml')
    return outputbool
                

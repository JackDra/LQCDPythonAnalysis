#!/usr/bin/env python

import numpy as np
from Params import *
from MiscFuns import *
from ReadTxt import ExtractValues
from SetLists import *
from FormFactors import CreateFF
from OutputData import PrintFFSet
from OppFuns import PrintOpps
from FFParams import *
import sys
from multiprocessing import Pool
from MultiWrap import *
import time
import datetime
# from guppy import hpy; h=hpy()
import resource



def PickMassSet(MassSet,theset):
    theMass = DefMassVal[DefMassVal.keys()[0]]
    thesetmass = 'Default'
    # if 'SF' in theset and RemoveTSink(theset) in MassSet.keys():
    #     thesetmass = RemoveTSink(theset)
    #     theMass = MassSet[thesetmass]['Boot'].Avg
    # elif 'SumMeth' in theset and 'sm32OSFTsink' in MassSet.keys():
    #     thesetmass = 'sm32OSFTsink'
    #     theMass = MassSet[thesetmass]['Boot'].Avg
    return theMass,thesetmass

def CreateFFWrap(thisMass,thesetmass,theset,setdict,thisCurr):
    # print 'Set:' + theset + ' MassSetPicked:'+thesetmass
## FF { { momsqrd } { Boot/Avg/Chi } }
    thisstart = time.time()
    FF = CreateFF(setdict,thisMass['Avg'],thisCurr)
    PrintFFSet(FF,theset,thisMass,thesetmass,thisCurr)
    # print 'Fit and Print for ' , theset , ' took: ',str(datetime.timedelta(seconds=time.time()-thisstart)) , ' h:m:s'    


#FitMasses later
def DoFF(thisMethodList,thisCurr,thisSetList,thisGammaList):
    # print 'Running : ' , thisCurr , ' ' , thisMethodList[0]
    # print 'SetList:\n' + '\n'.join(thisSetList)
    # print ''
    # PrintOpps(thisGammaList)
    # print 'Reading Methods:\n' + '\n'.join(thisMethodList)

    data,MassSet = ExtractValues(outputdir,thisGammaList,thisSetList,thisMethodList)
    if len(data.keys()) == 0:
        # print 'No Sets Found, returning'
        return

    # print 'All data Collected:'
    # for iCol in data.keys():
    #     print iCol
    # print ''
    # try:
    #     # print 'Mass data Collected:'
    #     for iCol in MassSet.keys():
    #         print iCol
    # except:
    #     print 'No Masses Found, using default'
    # data { { StateList } { gamma } { mom } { Fit(Boot/Avg/Std/Chi) } }
    # print ''
    # print 'Default Mass Set To: ' , DefMassVal[DefMassVal.keys()[0]]
    # print ''
    # print 'Creating Form Factors:' 
    inputparams = [PickMassSet(MassSet,theset)+(theset,setdict,thisCurr) for theset,setdict in data.iteritems()]
    start = time.time()
    for ipar in inputparams: CreateFFWrap(*ipar)
    print 'Fit and Print for ' , ' '.join(thisMethodList) , thisCurr , ' '.join(thisSetList) ,' in total took: ',str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s'    
    # print ''


# if kappa == 12090:
#     thisSmearList = DefSmearList
#     if OnlySelVar:
#         thisTvarList = [DefTvarPicked]
#     else:
#         thisTvarList = DefTvarList
#     thisTSinkList = AllTSinkList
#     thisStateList = ['1']
# elif kappa == 12104:
#     thisSmearList = []
#     if OnlySelVar:
#         thisTvarList = [DefTvarPicked]
#     else:
#         thisTvarList = DefTvarList
#     thisTSinkList = REvecTSinkList
#     thisStateList = ['1']
thisGammaList = []

# if len(sys.argv) == 1:
#     raise IOError('please select valid current type :\n' + '\n'.join(CurrTypes)+'\n or MethodType:\n'+'\n'.join(MethodList))

thisCurrList,thisMethodList = [],[]
for isys in sys.argv[1:]:
    if isys in MethodList:
        thisMethodList.append(isys)
    elif isys in CurrOpps.keys():
        thisCurrList.append(isys)
    else:
        print isys , ' Not valid current or method type, skipping'

if len(thisMethodList) == 0:
    print 'No Method Selected, doing all'
    thisMethodList = MethodList

if len(thisCurrList) == 0:
    print 'No Current Type Selected, doing all'
    thisCurrList = CurrTypes


    
print 'CurrList:\n' , '\n'.join(thisCurrList)
print ''
print 'MethodList:\n' , '\n'.join(thisMethodList)
print ''

inputparams = []
for thisCurr in thisCurrList:
    thisGammaList = CurrOpps[thisCurr] + [iC+'cmplx' for iC in CurrOpps[thisCurr]] + ['twopt']
    for imeth in thisMethodList:
        for iSet in DefSetList:
            print 'Adding to queue FF: ' , imeth , thisCurr , iSet
            inputparams.append(([imeth],thisCurr,[iSet],thisGammaList))

tottime = time.time()
if DoMulitcore:
    makeContextFunctions(DoFF)
    thisPool = Pool(min(AnaProc,len(inputparams)))
    thisPool.map(DoFF.mapper,inputparams)
    thisPool.close()
    thisPool.join()
else:
    for ip,iparam in enumerate(inputparams):
        print 'Total percent: ' ,GetPercent(ip,len(inputparams))
        print 'Time Left:' , GetTimeLeftStr(ip,len(inputparams),time.time()-tottime)
        DoFF(*iparam)
print 'Form Factor Creation Complete'

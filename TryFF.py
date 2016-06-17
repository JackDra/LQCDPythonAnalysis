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
from InputArgs import *



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
    # mprint( 'Set:' + theset + ' MassSetPicked:'+thesetmass)
## FF { { momsqrd } { Boot/Avg/Chi } }
    thisstart = time.time()
    FF,infodict = CreateFF(setdict,thisMass['Avg'],thisCurr)
    PrintFFSet(FF,theset,thisMass,thesetmass,thisCurr,infodict)
    if 'Vector' in thisCurr:
        NewFF = CombineVector(FF,thisMass)
        PrintFFSet(NewFF,theset,thisMass,thesetmass,'GeGm',infodict)
    mprint( 'Fit and Print for ' , theset , ' took: ',str(datetime.timedelta(seconds=time.time()-thisstart)) , ' h:m:s'    )


#FitMasses later
def DoFF(thisMethodList,thisCurr,thisSetList,thisGammaList):

    data,MassSet = ExtractValues(outputdir,thisGammaList,thisSetList,thisMethodList)
    if len(data.keys()) == 0:
        mprint( 'No Sets Found, returning')
        return

    mprint( 'All data Collected:')
    for iCol in data.keys():
        mprint( iCol)
    mprint( '')
    # data { { StateList } { gamma } { mom } { Fit(Boot/Avg/Std/Chi) } }
    mprint( '')
    mprint( 'Creating Form Factors:' )
    inputparams = [PickMassSet(MassSet,theset)+(theset,setdict,thisCurr) for theset,setdict in data.iteritems()]
    start = time.time()
    for ipar in inputparams: CreateFFWrap(*ipar)        
    print 'Fit and Print for ' , ' '.join(thisMethodList) , thisCurr , ' '.join(thisSetList) ,' in total took: ',str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s'
    mprint( '')



feedin = InputParams(sys.argv[1:])

print 'CurrList:\n' , '\n'.join(feedin['current'])
print ''
print 'MethodList:\n' , '\n'.join(feedin['method'])
print ''

inputparams = []
thisGammaList = []
if 'RF' in feedin['set']: del feedin['set']['RF']
for thisCurr in feedin['current']:
    thisGammaList = CurrOpps[thisCurr] + [iC+'cmplx' for iC in CurrOpps[thisCurr]] + ['twopt']
    for imeth in feedin['method']:
        if 'Fits' in imeth or 'OSF' in imeth:
            for iSet in feedin['set']:
                print 'Adding to queue FF: ' , imeth , thisCurr , iSet
                inputparams.append(([imeth],thisCurr,[iSet],thisGammaList))
        else:
            for iSet in GetTsinkSmLists(feedin['set'])[1]:
                print 'Adding to queue FF: ' , imeth , thisCurr , iSet
                inputparams.append(([imeth],thisCurr,[iSet],thisGammaList))
            
tottime = time.time()
if DoMulticore:
    makeContextFunctions(DoFF)
    thisPool = Pool(min(feedin['anaproc'],len(inputparams)))
    thisPool.map(DoFF.mapper,inputparams)
    thisPool.close()
    thisPool.join()
else:
    for ip,iparam in enumerate(inputparams):
        print 'Total percent: ' ,GetPercent(ip,len(inputparams))
        print 'Time Left:' , GetTimeLeftStr(ip,len(inputparams),time.time()-tottime)
        DoFF(*iparam)
print 'Form Factor Creation Complete'

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
    thisDS,baseCurr,dump = SplitDSCurr(thisCurr)
    if thisDS == '':
        thisCurrList = [thisCurr for ids in DefDSList]
        thisgflist = DefDSList
    else:
        thisCurrList = [thisCurr.replace(thisDS,'')]
        thisgflist = [thisDS]
    for iCurr,igf in zip(thisCurrList,thisgflist):
        combCurr = igf+iCurr
        FF,infodict = CreateFF(setdict,thisMass['Avg'],iCurr,gammaflag=igf)
        PrintFFSet(FF,theset,thisMass,thesetmass,combCurr,infodict)
        if 'Vector' in thisCurr and 'IsoVector' not in thisCurr and 'PsVector' not in thisCurr:
            NewFF = CombineVector(FF,thisMass)
            PrintFFSet(NewFF,theset,thisMass,thesetmass,combCurr.replace('Vector','GeGm'),infodict)

    mprint( 'Fit and Print for ' , theset , ' took: ',str(datetime.timedelta(seconds=time.time()-thisstart)) , ' h:m:s'    )


#FitMasses later
def DoFF(thisMethodList,thisCurr,thisSetList,thisGammaList,thisMomList):
    data,MassSet = ExtractValues(outputdir,thisGammaList,thisSetList,thisMethodList,thisMomList=thisMomList)
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

if 'RF' in feedin['method']: feedin['method'].remove('RF')
print 'MethodList:\n' , '\n'.join(feedin['method'])
print ''

inputparams = []
thisGammaList = []
if 'GeGm' in feedin['current']: feedin['current'].remove('GeGm')
for iDS in DefDSList:
    if iDS+'GeGm' in feedin['current']: feedin['current'].remove(iDS+'GeGm')
for thisCurr in feedin['current']:
    thisDS,baseCurr,dump = SplitDSCurr(thisCurr)
    if thisDS == '':
        thisGammaList = CurrOpps[baseCurr] + [iC+'cmplx' for iC in CurrOpps[baseCurr]]
        thisGammaList = ['doub'+ig for ig in thisGammaList] + ['sing'+ig for ig in thisGammaList]+ ['twopt']
    else:
        thisGammaList = CurrOpps[baseCurr] + [iC+'cmplx' for iC in CurrOpps[baseCurr]]
        thisGammaList = [thisDS+ig for ig in thisGammaList] + ['twopt']
        
    for imeth in feedin['method']:
        if 'Fits' in imeth or 'OSF' in imeth:
            for iSet in feedin['set']:
                print 'Adding to queue FF: ' , imeth , thisCurr , iSet
                inputparams.append(([imeth],thisCurr,[iSet],thisGammaList,feedin['mom']))
        else:
            for iSet in ReduceTsink(feedin['set']):
                print 'Adding to queue FF: ' , imeth , thisCurr , iSet
                inputparams.append(([imeth],thisCurr,[iSet],thisGammaList,feedin['mom']))
            
starttime = time.time()
feedin['anaproc'] = min(feedin['anaproc'],len(inputparams))
if DoMulticore and feedin['anaproc'] > 1:
    makeContextFunctions(DoFF)
    thisPool = Pool(feedin['anaproc'])
    thisPool.map(DoFF.mapper,inputparams)
    thisPool.close()
    thisPool.join()
else:
    for ip,iparam in enumerate(inputparams):
        print 'Total percent: ' ,GetPercent(ip,len(inputparams))
        print 'Time Left:' , GetTimeLeftStr(ip,len(inputparams),time.time()-tottime)
        DoFF(*iparam)
print 'Form Factor Creation Complete, time taken:', str(datetime.timedelta(seconds=time.time()-starttime)) , ' h:m:s '

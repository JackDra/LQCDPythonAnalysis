#!/usr/bin/env python

from Params import *
import numpy as np
from ReadTxt import *
from MiscFuns import *
from FitParams import *
from StringFix import *
from GraphDataNew import *
from SetLists import CutDupSet
from ReadDir import *
from FFParams import *
from FitParams import *
from GraphSummary import PlotFFSummary
import time,datetime
from InputArgs import *
from MultiWrap import *
from multiprocessing import Pool

feedin = InputParams(sys.argv[1:]+['-noprompt'])

def PlotMKFFWrap(a,b,c,d):
    PlotMKFFs(a,b,c,d,feedin['ForceTitle'])



def PickMKFFFewSets(currdata,thiscurr,thisSetList):
    PickedSetList = []
    # if kappa == 12090:
    if '12090' in currdata.keys():
        PickedSetList += ['k12090'+ifl for ifl in FlagList(thisSetList,'Fit','tsink29sm32',FitCutPicked['tsink29sm32'])]
        # PickedSetList += ['k12090'+ifl for ifl in FlagList(thisSetList,'Fit','tsink29sm64',FitCutPicked['tsink29sm64'])]
        # PickedSetList += ['k12090'+ifl for ifl in FlagList(thisSetList,'Fit','tsink29sm128',FitCutPicked['tsink29sm128'])]
        PickedSetList += ['k12090'+ifl for ifl in FlagList(thisSetList,'Fit','tsink29state1CMto20dt2',FitCutPicked['tsink29state1CM'])]
        PickedSetList += ['k12090'+ifl for ifl in FlagList(thisSetList,'TSFTsink','sm32',TSFCutPicked)]
        # elif kappa == 12104:
    if '12104' in currdata.keys():
        PickedSetList += ['k12104'+ifl for ifl in FlagList(thisSetList,'Fit','tsink29state1REvec',FitCutPicked['tsink29state1REvec'])]
        # PickedSetList += ['k12104'+ifl for ifl in FlagList(thisSetList,'OSFCM','tsink29state1REvec',OSFCutPicked)]
        
    # PickedSetList += FlagList(thisSetList,'OSFCM','tsink29sm32',OSFCutPicked)
    # PickedSetList += FlagList(thisSetList,'OSFCM','tsink29state1CMto20dt2',OSFCutPicked)
    # for Fitkey,FitCutVal in FitCutPicked.iteritems():
    # PickedSetList += FlagList(thisSetList,'Fit','tsink29'+PickedStateStr+'CM',FitCutPicked['tsink29'+PickedStateStr+'CM'])
    # PickedSetList += FlagList(thisSetList,'Fit','tsink26'+PickedStateStr+'PoF',FitCutPicked['tsink26'+PickedStateStr+'PoF'+str(PoFShifts)])
    # PickedSetList += FlagList(thisSetList,'Fit','tsink35sm32',FitCutPicked['tsink35sm32'])
    # PickedSetList += FlagList(thisSetList,'SumMeth',SumFitRPicked,SumCutPar)
    # PickedSetList += FlagList(thisSetList,'OSF',OSFCutPicked)
    # PickedSetList += FlagList(thisSetList,'TSFTsink',TSFCutPicked)
    return PickedSetList

def ReadAndPlotMKFF(thisCurrDict,DoList='All'):
    datadict = ReadMKFFDict(outputdir,thisCurrDict)
    start = time.time()
    currPSL = []
    for thiscurr,currdata in datadict.iteritems():
        if 'Few' in DoList or 'All' in DoList:
            print 'Plotting ' , thiscurr ,'6/6 Few            '
            PlotMKFFWrap(currdata,thiscurr,PickMKFFFewSets(currdata,thiscurr,thisCurrDict[thiscurr]),'')
        print 'Plotting ' , thiscurr ,'Complete, took: ', GetTimeStr(time.time()-start)
    



thisFFcomb = []
for icomb in feedin['comb']:
    if icomb in CombListForFFComb:
        thisFFcomb.append(icomb)
thisCurrcomb = []
for icomb in feedin['current']:
    if icomb in CurrListForFFComb:
        thisCurrcomb.append(icomb)

thisCurrDict = []
for iFFcomb in feedin['FFcomb']:
    if iFFcomb == '':
        for icurr in ElongateName(feedin['comb'],feedin['current']):
            print 'Looking in ', icurr
            thisCurrDict.append([GetCurrDict([icurr],klist=feedin['klist']),feedin['ffgraph']])
    else:
        for icurr in ElongateName(thisFFcomb,thisCurrcomb):
            print 'Looking in ', icurr+iFFcomb
            thisCurrDict.append([GetCurrDict([icurr+iFFcomb],klist=feedin['klist']),feedin['ffgraph']])

        
makeContextFunctions(ReadAndPlotMKFF)

starttime = time.time()
if len(thisCurrDict) > 0:
    if DoMulticore and feedin['anaproc'] > 1 and len(thisCurrDict) > 1:
        print 'Running Multicore'
        thisPool = Pool(min(len(thisCurrDict),feedin['anaproc']))
        thisPool.map(ReadAndPlotMKFF.mapper,thisCurrDict)
        thisPool.close()
        thisPool.join()
        print 'MKFF Plots Complete'
    else:
        print 'Running Single core'
        for icount,iin in enumerate(thisCurrDict):
            ReadAndPlotMKFF(*iin)
            print int((icount*100)/float(len(thisCurrDict))) , '% done' + ' '*50 + '\r',
        print 'MKFF Plots Complete'
    print 'Graphing Complete, time taken:', str(datetime.timedelta(seconds=time.time()-starttime)) , ' h:m:s '
                
else:
    print 'nothing to calculate'        


print 'All Plotting Complete'

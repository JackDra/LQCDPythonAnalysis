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

def PlotFFWrap(a,b,c,d):
    PlotFFs(a,b,c,d,feedin['ForceTitle'])


def FlagList(AllSetList,*flag):
    SLOut = []
    for iset in AllSetList:
        if all([str(iflag) in iset for iflag in flag]):
            SLOut.append(iset)
    return SLOut

def PlotTSFSets(currdata,thiscurr,thisSetList):
    for iTSF in TSFFileFlags:
        if iTSF == 'CM':
            for ism in DefSmList:
                PlotFFWrap(currdata,thiscurr,FlagList(thisSetList,'TSF',iTSF,ism),'TSF'+iTSF+ism+'CutComp')
        else:
            PlotFFWrap(currdata,thiscurr,FlagList(thisSetList,'TSF',iTSF),'TSF'+iTSF+'CutComp')
    for icut in TSFCutList:
        PlotFFWrap(currdata,thiscurr,FlagList(thisSetList,'TSF',icut),'TSFAll'+icut+'MethComp')


def PlotOSFSets(currdata,thiscurr,thisSetList):
    for iOSF in OSFFileFlags:
        for itsink in AllTSinkStrList:
            PlotFFWrap(currdata,thiscurr,FlagList(thisSetList,'OSF',iOSF,itsink),'OSF'+iOSF+itsink+'CutComp')
        for icut in OSFCutList:
            PlotFFWrap(currdata,thiscurr,FlagList(thisSetList,'OSF',iOSF,icut),'OSF'+iOSF+icut+'MethComp')
    for icut in OSFCutList:
        PlotFFWrap(currdata,thiscurr,FlagList(thisSetList,'OSF',icut),'OSFAll'+icut+'MethComp')
        
def PlotSumMethSets(currdata,thiscurr,thisSetList):
    for icut in SumCutList:
        PlotFFWrap(currdata,thiscurr,FlagList(thisSetList,'SumMeth',icut),'SumMeth'+icut+'FitRComp')
    for ifitr in SumFitRListFlags:
        PlotFFWrap(currdata,thiscurr,FlagList(thisSetList,'SumMeth',ifitr),'SumMeth'+ifitr+'CutComp')

def PlotFitMethSets(currdata,thiscurr,thisSetList):
    PlotFFWrap(currdata,thiscurr,FlagList(thisSetList,'Fit','tsink29','cut6')+FlagList(thisSetList,'Fit','PoF','cut6'),'FitMytsink29')
    PlotFFWrap(currdata,thiscurr,FlagList(thisSetList,'Fit','sm32','cut6'),'FitMysm32')
    for iset in DefSetList:
        PlotFFWrap(currdata,thiscurr,FlagList(thisSetList,'Fit',iset),'Fit'+iset+'CutComp')
    for icut in FitCutArgs:
        PlotFFWrap(currdata,thiscurr,FlagList(thisSetList,'Fit',icut),'Fit'+icut+'SetComp')
        

def PickFFAllSets(currdata,thiscurr,thisSetList):
    PickedSetList = []
    PickedSetList += FlagList(thisSetList,'Fit','tsink29sm32',FitCutPicked['tsink29sm32'])
    PickedSetList += FlagList(thisSetList,'Fit','tsink29state1CMto20dt2',FitCutPicked['tsink29state1CM'])
    PickedSetList += FlagList(thisSetList,'OSFCM','tsink29sm32',OSFCutPicked)
    PickedSetList += FlagList(thisSetList,'OSFCM','tsink29state1CMto20dt2',OSFCutPicked)
    PickedSetList += FlagList(thisSetList,'TSFTSink','sm32',TSFCutPicked)
    # for Fitkey,FitCutVal in FitCutPicked.iteritems():
    #     PickedSetList += FlagList(thisSetList,'Fit',Fitkey,FitCutVal)
    # for ifitr in SumFitRListFlags:
    #     for icut in SumCutList:
    #         PickedSetList += FlagList(thisSetList,'SumMeth',ifitr,icut)
    # for icut in OSFCutList:
    #     PickedSetList += FlagList(thisSetList,'OSF',icut)
    # for icut in TSFCutList:
    #     PickedSetList += FlagList(thisSetList,'TSFSmall',icut)
    # for icut in TSFCutList:
    #     PickedSetList += FlagList(thisSetList,'TSFTsink',icut)
    return PickedSetList

def PickFFFewSets(currdata,thiscurr,thisSetList):
    PickedSetList = []
    # for Fitkey,FitCutVal in FitCutPicked.iteritems():
    PickedSetList += FlagList(thisSetList,'Fit','tsink29'+PickedStateStr+'CM',FitCutPicked['tsink29'+PickedStateStr+'CM'])
    PickedSetList += FlagList(thisSetList,'Fit','tsink26'+PickedStateStr+'PoF',FitCutPicked['tsink26'+PickedStateStr+'PoF'+str(PoFShifts)])
    # PickedSetList += FlagList(thisSetList,'Fit','tsink35sm32',FitCutPicked['tsink35sm32'])
    # PickedSetList += FlagList(thisSetList,'SumMeth',SumFitRPicked,SumCutPar)
    # PickedSetList += FlagList(thisSetList,'OSF',OSFCutPicked)
    # PickedSetList += FlagList(thisSetList,'TSFTsink',TSFCutPicked)
    return PickedSetList

def ReadAndPlotFF(thisCurrDict,DoList='All'):
    datadict = ReadFFDict(outputdir,thisCurrDict)
    start = time.time()
    currPSL = []
    for thiscurr,currdata in datadict.iteritems():
        if 'TSF' in DoList or 'All' in DoList:
            print 'Plotting ' , thiscurr ,'1/6 TSF             '
            PlotTSFSets(currdata,thiscurr,thisCurrDict[thiscurr])
        if 'OSF' in DoList or 'All' in DoList:
            print 'Plotting ' , thiscurr ,'2/6 OSF             '
            PlotOSFSets(currdata,thiscurr,thisCurrDict[thiscurr])
        if 'Sum' in DoList or 'All' in DoList:
            print 'Plotting ' , thiscurr ,'3/6 Summation       '
            PlotSumMethSets(currdata,thiscurr,thisCurrDict[thiscurr])
        if 'Fits' in DoList or 'All' in DoList:
            print 'Plotting ' , thiscurr ,'4/6 Fits            '
            PlotFitMethSets(currdata,thiscurr,thisCurrDict[thiscurr])
        if 'Collect' in DoList or 'All' in DoList:
            print 'Collecting ' , thiscurr ,'5/6 for Summary            '
            currPSL.append(PickFFAllSets(currdata,thiscurr,thisCurrDict[thiscurr]))
        if 'Few' in DoList or 'All' in DoList:
            print 'Plotting ' , thiscurr ,'6/6 Few            '
            PlotFFWrap(currdata,thiscurr,PickFFFewSets(currdata,thiscurr,thisCurrDict[thiscurr]),'Summary')
        print 'Plotting ' , thiscurr ,'Complete, took: ', GetTimeStr(time.time()-start)
    return datadict,currPSL
    
def PlotFFqPick(datadict,thisPSL):
    start = time.time()
    for thisSL,(thiscurr,currdata) in zip(thisPSL,datadict.iteritems()):
        PlotFFSummary(thisSL,thiscurr,currdata,feedin['ForceTitle'])
        print 'Plotting Summary for qsqrd ' , thiscurr ,'Complete, took: ', GetTimeStr(time.time()-start)



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
            thisCurrDict.append([GetCurrDict([icurr]),feedin['ffgraph']])
    else:
        for icurr in ElongateName(thisFFcomb,thisCurrcomb):
            print 'Looking in ', icurr+iFFcomb
            thisCurrDict.append([GetCurrDict([icurr+iFFcomb]),feedin['ffgraph']])

        
makeContextFunctions(ReadAndPlotFF)
makeContextFunctions(PlotFFqPick)

starttime = time.time()
if len(thisCurrDict) > 0:
    if DoMulticore and feedin['anaproc'] > 1 and len(thisCurrDict) > 1:
        print 'Running Multicore'
        thisPool = Pool(min(len(thisCurrDict),feedin['anaproc']))
        output = thisPool.map(ReadAndPlotFF.mapper,thisCurrDict)
        thisPool.close()
        thisPool.join()
        print 'FF Plots Complete'
        if kappa == 12090:
            thisPool = Pool(min(len(output),feedin['anaproc']))
            output = thisPool.map(PlotFFqPick.mapper,output)
            thisPool.close()
            thisPool.join()
    else:
        print 'Running Single core'
        output = []
        for icount,iin in enumerate(thisCurrDict):
            output.append(ReadAndPlotFF(*iin))
            print int((icount*100)/float(len(thisCurrDict))) , '% done' + ' '*50 + '\r',
        print 'FF Plots Complete'
        if kappa == 12090:
            for iout in output:
                PlotFFqPick(*iout)
    print 'Graphing Complete, time taken:', str(datetime.timedelta(seconds=time.time()-starttime)) , ' h:m:s '
                
else:
    print 'nothing to calculate'        
    output = []

# datadict,currPSL = ReadAndPlotFF(thisCurrDict)
# # datadict,currPSL = ReadAndPlotFF(thisCurrDict)

print 'All Plotting Complete'

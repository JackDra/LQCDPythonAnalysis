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
                PlotFFs(currdata,thiscurr,FlagList(thisSetList,'TSF',iTSF,ism),'TSF'+iTSF+ism+'CutComp')
        else:
            PlotFFs(currdata,thiscurr,FlagList(thisSetList,'TSF',iTSF),'TSF'+iTSF+'CutComp')
    for icut in TSFCutList:
        PlotFFs(currdata,thiscurr,FlagList(thisSetList,'TSF',icut),'TSF'+icut+'MethComp')


def PlotOSFSets(currdata,thiscurr,thisSetList):
    for iOSF in OSFFileFlags:
        for itsink in AllTSinkStrList:
            PlotFFs(currdata,thiscurr,FlagList(thisSetList,'OSF',iOSF,itsink),'OSF'+iOSF+itsink+'CutComp')
        for icut in OSFCutList:
            PlotFFs(currdata,thiscurr,FlagList(thisSetList,'OSF',iOSF,icut),'OSF'+iOSF+icut+'MethComp')
    for icut in OSFCutList:
        PlotFFs(currdata,thiscurr,FlagList(thisSetList,'OSF',icut),'OSFAll'+icut+'MethComp')
        
def PlotSumMethSets(currdata,thiscurr,thisSetList):
    for icut in SumCutList:
        PlotFFs(currdata,thiscurr,FlagList(thisSetList,'SumMeth',icut),'SumMeth'+icut+'FitRComp')
    for ifitr in SumFitRListFlags:
        PlotFFs(currdata,thiscurr,FlagList(thisSetList,'SumMeth',ifitr),'SumMeth'+ifitr+'CutComp')

def PlotFitMethSets(currdata,thiscurr,thisSetList):
    PlotFFs(currdata,thiscurr,FlagList(thisSetList,'Fit','tsink29','cut6'),'FitMytsink29')
    PlotFFs(currdata,thiscurr,FlagList(thisSetList,'Fit','sm32','cut6'),'FitMysm32')
    for iset in DefSetList:
        PlotFFs(currdata,thiscurr,FlagList(thisSetList,'Fit',iset),'Fit'+iset+'CutComp')
    for icut in FitCutArgs:
        PlotFFs(currdata,thiscurr,FlagList(thisSetList,'Fit',icut),'Fit'+icut+'SetComp')
        

def PickFFAllSets(currdata,thiscurr,thisSetList):
    PickedSetList = []
    for Fitkey,FitCutVal in FitCutPicked.iteritems():
        PickedSetList += FlagList(thisSetList,'Fit',Fitkey,FitCutVal)
    for ifitr in SumFitRListFlags:
        for icut in SumCutList:
            PickedSetList += FlagList(thisSetList,'SumMeth',ifitr,icut)
    # for icut in OSFCutList:
    #     PickedSetList += FlagList(thisSetList,'OSF',icut)
    for icut in TSFCutList:
        PickedSetList += FlagList(thisSetList,'TSF',icut)
    return PickedSetList

def PickFFFewSets(currdata,thiscurr,thisSetList):
    PickedSetList = []
    # for Fitkey,FitCutVal in FitCutPicked.iteritems():
    PickedSetList += FlagList(thisSetList,'Fit','tsink29state1'+VarPref,FitCutPicked['tsink29state1'+VarPref])
    PickedSetList += FlagList(thisSetList,'Fit','tsink35sm32',FitCutPicked['tsink35sm32'])
    PickedSetList += FlagList(thisSetList,'SumMeth',SumFitRPicked,SumCutPar)
    # PickedSetList += FlagList(thisSetList,'OSF',OSFCutPicked)
    PickedSetList += FlagList(thisSetList,'TSFTsink',TSFCutPicked)
    return PickedSetList

def ReadAndPlotFFs(thisCurrDict):
    datadict = ReadFFDict(outputdir,thisCurrDict)
    start = time.time()
    currPSL = []
    for thiscurr,currdata in datadict.iteritems():
        print 'Plotting ' , thiscurr ,'1/5 TSF             \r',
        PlotTSFSets(currdata,thiscurr,thisCurrDict[thiscurr])
        print 'Plotting ' , thiscurr ,'2/5 OSF             \r',
        PlotOSFSets(currdata,thiscurr,thisCurrDict[thiscurr])
        print 'Plotting ' , thiscurr ,'3/5 Summation       \r',
        PlotSumMethSets(currdata,thiscurr,thisCurrDict[thiscurr])
        print 'Plotting ' , thiscurr ,'4/5 Fits            \r',
        PlotFitMethSets(currdata,thiscurr,thisCurrDict[thiscurr])
        print 'Plotting ' , thiscurr ,'5/5 Fits            \r',
        currPSL.append(PickFFAllSets(currdata,thiscurr,thisCurrDict[thiscurr]))
        PlotFFs(currdata,thiscurr,PickFFFewSets(currdata,thiscurr,thisCurrDict[thiscurr]),'Summary')
        print 'Plotting ' , thiscurr ,'Complete, took: ', GetTimeStr(time.time()-start)
    return datadict,currPSL
    
def PlotFFqPick(datadict,thisPSL):
    start = time.time()
    for thisSL,(thiscurr,currdata) in zip(currPSL,datadict.iteritems()):
        PlotFFSummary(thisSL,thiscurr,currdata)
        print 'Plotting Summary for qsqrd ' , thiscurr ,'Complete, took: ', GetTimeStr(time.time()-start)


feedin = InputParams(sys.argv[1:])

thisCurrDict = GetCurrDict(feedin['current'])

if DoMulticore:
    thisPool = Pool(min(len(thisCurrDict),AnaProc))
    makeContextFunctions(ReadAndPlotFFs)
    output = thisPool.map(ReadAndPlotFFs.mapper,thisCurrDict)
    thisPool.close()
    thisPool.join()
    # if kappa == 12090:
    #     thisPool = Pool(min(len(output),AnaProc))
    #     makeContextFunctions(PlotFFqPick)
    #     thisPool.map(PlotFFqPick.mapper,output)
    #     thisPool.close()
    #     thisPool.join()        
else:
    datadict,currPSL = ReadAndPlotFFs(thisCurrDict)
    if kappa == 12090: PlotFFqPick(datadict,currPSL)

print 'All Plotting Complete'

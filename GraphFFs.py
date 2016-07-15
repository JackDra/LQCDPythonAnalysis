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
        PlotFFs(currdata,thiscurr,FlagList(thisSetList,'TSF',icut),'TSFAll'+icut+'MethComp')


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
    PlotFFs(currdata,thiscurr,FlagList(thisSetList,'Fit','tsink29','cut6')+FlagList(thisSetList,'Fit','PoF','cut6'),'FitMytsink29')
    PlotFFs(currdata,thiscurr,FlagList(thisSetList,'Fit','sm32','cut6'),'FitMysm32')
    for iset in DefSetList:
        PlotFFs(currdata,thiscurr,FlagList(thisSetList,'Fit',iset),'Fit'+iset+'CutComp')
    for icut in FitCutArgs:
        PlotFFs(currdata,thiscurr,FlagList(thisSetList,'Fit',icut),'Fit'+icut+'SetComp')
        

def PickFFAllSets(currdata,thiscurr,thisSetList):
    PickedSetList = []
    for Fitkey,FitCutVal in FitCutPicked.iteritems():
        PickedSetList += FlagList(thisSetList,'Fit',Fitkey,FitCutVal)
    # for ifitr in SumFitRListFlags:
    #     for icut in SumCutList:
    #         PickedSetList += FlagList(thisSetList,'SumMeth',ifitr,icut)
    for icut in OSFCutList:
        PickedSetList += FlagList(thisSetList,'OSF',icut)
    for icut in TSFCutList:
        PickedSetList += FlagList(thisSetList,'TSFSmall',icut)
    for icut in TSFCutList:
        PickedSetList += FlagList(thisSetList,'TSFTsink',icut)
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
            print 'Plotting ' , thiscurr ,'1/5 TSF             '
            PlotTSFSets(currdata,thiscurr,thisCurrDict[thiscurr])
        if 'OSF' in DoList or 'All' in DoList:
            print 'Plotting ' , thiscurr ,'2/5 OSF             '
            PlotOSFSets(currdata,thiscurr,thisCurrDict[thiscurr])
        if 'Sum' in DoList or 'All' in DoList:
            print 'Plotting ' , thiscurr ,'3/5 Summation       '
            PlotSumMethSets(currdata,thiscurr,thisCurrDict[thiscurr])
        if 'Fits' in DoList or 'All' in DoList:
            print 'Plotting ' , thiscurr ,'4/5 Fits            '
            PlotFitMethSets(currdata,thiscurr,thisCurrDict[thiscurr])
        if 'Fits' in DoList or 'All' in DoList:
            print 'Plotting ' , thiscurr ,'5/5 Fits Summary            '
            currPSL.append(PickFFAllSets(currdata,thiscurr,thisCurrDict[thiscurr]))
        PlotFFs(currdata,thiscurr,PickFFFewSets(currdata,thiscurr,thisCurrDict[thiscurr]),'Summary')
        print 'Plotting ' , thiscurr ,'Complete, took: ', GetTimeStr(time.time()-start)
    return datadict,currPSL
    
def PlotFFqPick(datadict,thisPSL):
    start = time.time()
    for thisSL,(thiscurr,currdata) in zip(thisPSL,datadict.iteritems()):
        PlotFFSummary(thisSL,thiscurr,currdata)
        print 'Plotting Summary for qsqrd ' , thiscurr ,'Complete, took: ', GetTimeStr(time.time()-start)


feedin = InputParams(sys.argv[1:]+['-noprompt'])

DoList='Fits'
thisCurrDict = []
for icurr in ElongateName(ElongateName(feedin['comb'],feedin['current']),['/'+iCombFF for iCombFF in CombFFList]+['']):
    print 'Looking in ', icurr
    thisCurrDict.append([GetCurrDict([icurr]),DoList])
makeContextFunctions(ReadAndPlotFF)

if len(thisCurrDict) > 0:
    if DoMulticore and feedin['anaproc'] > 1 and len(thisCurrDict) > 1:
        print 'Running Multicore'
        thisPool = Pool(min(len(thisCurrDict),feedin['anaproc']))
        output = thisPool.map(ReadAndPlotFF.mapper,thisCurrDict)
        thisPool.close()
        thisPool.join()
    else:
        print 'Running Single core'
        output = []
        for icount,iin in enumerate(thisCurrDict):
            output.append(ReadAndPlotFF(*iin))
            print int((icount*100)/float(len(thisCurrDict))) , '% done' + ' '*50 + '\r',
else:
    print 'nothing to calculate'        
    output = []

# datadict,currPSL = ReadAndPlotFF(thisCurrDict)
# # datadict,currPSL = ReadAndPlotFF(thisCurrDict)
if kappa == 12090:
    for iout in output:
        PlotFFqPick(*iout)

print 'All Plotting Complete'

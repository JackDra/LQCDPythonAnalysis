#!/usr/bin/env python

from Params import *
import numpy as np
from ReadTxt import *
from MiscFuns import *
from FitParams import *
from StringFix import *
from GraphDataNew import *
from SetLists import CutDupSet
from OppFuns import CreateGammaList
import time,datetime
from MultiWrap import *
from multiprocessing import Pool
from InputArgs import InputGammaAndSet

##Hard Code lists:
thisMethodList = ['RF']
# thisMethodList = MethodList
# thisMomList = ['q = 0 0 0', 'q = -1 0 0', 'q = -1 -1 0','q = -1 -1 -1','q = -2 0 0' , 'q = -2 0 -1', 'q = -2 -1 -1','q = -2 -2 0','q = -2 -2 -1','q = -3 0 0'   ]
thisMomList = ['q = 0 0 0']
##


##datadict = { gamma } { mom } { method } { set }
##thisdatadict = { method } { set }
##thisGammaList format: D/S P4/3 Opp (e.g. doubP4g4)
##thisMomList in qstr format (e.g. 'q = 0 0 0')
##NB q = 0 0 0 MUST BE INCLUDED FOR TSF ETC
def progprint(numb,starttime,igamma):
    print 'Graphing Operator: ' + igamma , int(numb*100/float(9)) , '% time taken:' , str(datetime.timedelta(seconds=time.time()-starttime)) ,' h:m:s          \r',


def ReadAndPlotMass(thisMomList,thisSmearList,thisSetList,thisSetPoFLists,thisMethodList):
    thisAllSetList = thisSmearList+thisSetList
    for isetlist,dump in thisSetPoFLists[:len(thisSetPoFLists)/2]:
        thisAllSetList += isetlist
    datadict = ReadSetFitRFDict(outputdir,thisAllSetList,['twopt'],thisMethodList,thisMomList=thisMomList)
    thisMassdict = datadict['twopt']['q = 0 0 0']
    start = time.time()
    for imom in thisMomList[-1:]:
        momstart = time.time()
        # print 'Plotting ' , imom , 'Mass            \r',
        thistwoptdict = datadict['twopt'][imom]
        PlotMassData(thistwoptdict,thisSmearList+thisSetList,imom)
        # print 'Plotting ' , imom , 'State Fits      \r',
        PlotMassSFData(thistwoptdict,thisSmearList+thisSetList,imom)
        for ipof,(itpl,thistf) in enumerate(thisSetPoFLists):
            # print 'Plotting ' , imom , 'PoFSets ', str(int((ipof*100)/len(thisSetPoFLists))),'%        \r',

            PlotMassData(thistwoptdict,thisSmearList+itpl,imom,TitleFlag=thistf)
        print 'Plotting ' , imom, 'Took: ' , str(datetime.timedelta(seconds=(time.time()-momstart))) ,' h:m:s                      '
    
    
def ReadAndPlotDict(thisGammaList,thisMomList,thisSetList,thisMethodList):
    print thisSetList
    datadict = ReadSetFitRFDict(outputdir,thisSetList,thisGammaList,thisMethodList,thisMomList=thisMomList)
    thisGammaList = datadict.keys()
    thisMassdict = datadict['twopt']['q = 0 0 0']
    start = time.time()
    for imom in thisMomList:
        for icg,igamma in enumerate(thisGammaList):
            if any([idst in igamma for idst in ['doub','sing','twopt']]): continue
            gammastart = time.time()
            timeleft = GetTimeLeft(icg,len(thisGammaList),time.time()-start)
            if not CheckDict(datadict,igamma,imom): continue
            thisdatadict = datadict[igamma][imom]
            if kappa == 12090:
                progprint(0,time.time(),igamma)
                prevtime = time.time()
                PlotTSinkData(thisdatadict,thisSetList,igamma,imom)
                progprint(1,prevtime,igamma)
                prevtime = time.time()
                PlotTSinkData(thisdatadict,thisSetList,igamma,imom,thissm='state1'+REvecTvarList[0])
                progprint(2,prevtime,igamma)
                prevtime = time.time()
                # PlotTSinkSumData(thisdatadict,thisSetList,igamma,imom)
                # progprint(3,prevtime,igamma)
                # prevtime = time.time()
                # PlotTSinkSFData(thisdatadict,thisMassdict,thisSetList,igamma,imom,thisSF='TSFTsink')
                # progprint(4,prevtime,igamma)
                # prevtime = time.time()
                # PlotTSinkSFData(thisdatadict,thisMassdict,thisSetList,igamma,imom,thisSF='TSFtest32')
                # progprint(5,prevtime,igamma)
                # prevtime = time.time()
                # PlotTSinkSFData(thisdatadict,thisMassdict,thisSetList,igamma,imom,thisSF='TSFSmall')
                # progprint(6,prevtime,igamma)
                # prevtime = time.time()
                # PlotTSinkSFData(thisdatadict,thisMassdict,thisSetList,igamma,imom,thisSF='OSFTsink')
                # progprint(7,prevtime,igamma)
                prevtime = time.time()
                PlotCMData(thisdatadict,thisSetList,igamma,imom)
                progprint(8,prevtime,igamma)
                # PlotCMSFData(thisdatadict,thisMassdict,thisSetList,igamma,imom)
            elif kappa == 12104:
                PlotCMData(thisdatadict,thisSetList,igamma,imom)
                # PlotCMSFData(thisdatadict,thisMassdict,['tsink29state1to18dt2'],igamma,imom,thisSF='SFREvec')
            print 'Graphing Operator: ' + igamma + imom + ' took: ' , str(datetime.timedelta(seconds=(time.time()-gammastart))) ,' h:m:s                      '
                    


print sys.argv[1:]
feedgammalist,feedsetlist = InputGammaAndSet(sys.argv[1:])
            
thisGammaList = CreateGammaList(feedgammalist,twopt=True)

if thisGammaList == ['twopt']:
    # thisMethodList = ['RF','OSFCM','TSFCM']
    TvarPicked = ['tsink29'+str(istate) for istate in CreateMassSet([],['1'],[DefTvarPicked])]
    thisSmList = ['tsink29'+str(ism) for ism in CreateMassSet(DefSmearList,['1'],[])]
    # TvarLists = []
    TvarLists = [(['tsink29'+str(istate) for istate in CreateMassSet([],['1'],DefTvarDt1)],'PoFDt1'),
                 (['tsink29'+str(istate) for istate in CreateMassSet([],['1'],DefTvarDt2)],'PoFDt2'),
                 (['tsink29'+str(istate) for istate in CreateMassSet([],['1'],DefTvarDt3)],'PoFDt3'),
                 (['tsink29'+str(istate) for istate in CreateMassSet([],['1'],DefTvarDt4)],'PoFDt4'),
                 (['tsink29'+str(istate) for istate in CreateMassSet([],['1'],DefTvarto16)],'PoFto16'),
                 (['tsink29'+str(istate) for istate in CreateMassSet([],['1'],DefTvarto17)],'PoFto17'),
                 (['tsink29'+str(istate) for istate in CreateMassSet([],['1'],DefTvarto18)],'PoFto18'),
                 (['tsink29'+str(istate) for istate in CreateMassSet([],['1'],DefTvarto19)],'PoFto19'),
                 (['tsink29'+str(istate) for istate in CreateMassSet([],['1'],DefTvarto20)],'PoFto20')]
    thisAllSetList = thisSmList+TvarPicked
    print 'AllSetList:\n' + '\n'.join(thisAllSetList)
    print 'MethodList:\n' + '\n'.join(thisMethodList)
    if DoMulticore:
        inputparams = []
        for imom in thisMomList:
            if imom == 'q = 0 0 0':
                inputparams.append(([imom],thisSmList,TvarPicked,TvarLists,thisMethodList))
            else:
                inputparams.append((['q = 0 0 0',imom],thisSmList,TvarPicked,TvarLists,thisMethodList))                
        makeContextFunctions(ReadAndPlotMass)
        thisPool = Pool(min(len(inputparams),AnaProc))
        output = thisPool.map(ReadAndPlotMass.mapper,inputparams)
        thisPool.close()
        thisPool.join()
    else:
        ReadAndPlotMass(thisGammaList,thisMomList,feedsetlist,thisMethodList)
else:
    if thisMomList != ['q = 0 0 0']:
        thisMethodList = ['RF']
    print 'MethodList:\n' + '\n'.join(thisMethodList)
    print 'thisSetList:\n' + '\n'.join(feedsetlist)
    if DoMulticore:
        inputparams = []
        for igamma in thisGammaList:
            if any([idst in igamma for idst in ['doub','sing','twopt']]): continue
            inputparams.append((['doub'+igamma,'sing'+igamma,igamma,'twopt'],thisMomList,feedsetlist,thisMethodList))
        makeContextFunctions(ReadAndPlotDict)
        thisPool = Pool(min(len(inputparams),AnaProc))
        output = thisPool.map(ReadAndPlotDict.mapper,inputparams)
        thisPool.close()
        thisPool.join()
    else:
        ReadAndPlotDict(thisGammaList,thisMomList,feedsetlist,thisMethodList)
        
print 'Graphing all complete'
    

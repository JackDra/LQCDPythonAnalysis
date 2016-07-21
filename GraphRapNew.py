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
from InputArgs import *
from CreateCombs import CreateDictOldCombs
from CombParams import *

DoDS = True


##datadict = { gamma } { mom } { method } { set }
##thisdatadict = { method } { set }
##thisGammaList format: D/S P4/3 Opp (e.g. doubP4g4)
##thisMomList in qstr format (e.g. 'q = 0 0 0')
##NB q = 0 0 0 MUST BE INCLUDED FOR TSF ETC
def progprint(numb,starttime,igamma):
    print 'Graphing Operator: ' + igamma , int(numb*100/float(9)) , '% time taken:' , str(datetime.timedelta(seconds=time.time()-starttime)) ,' h:m:s          '


def ReadAndPlotMass(thisMomList,thisSmearList,thisSetList,thisSetPoFLists,thisMethodList):
    thisAllSetList = thisSmearList+thisSetList
    for isetlist,dump in thisSetPoFLists[:len(thisSetPoFLists)/2]:
        thisAllSetList += isetlist
    datadict = ReadSetFitRFDict(outputdir,thisAllSetList,['twopt'],thisMethodList,thisMomList=thisMomList)
    thisMassdict = datadict['twopt']['q = 0 0 0']
    start = time.time()
    for imom in thisMomList[-1:]:
        momstart = time.time()
        mprint( 'Plotting ' , imom , 'Mass            \r',)
        thistwoptdict = datadict['twopt'][imom]
        PlotMassData(thistwoptdict,thisSmearList+thisSetList,imom)
        mprint( 'Plotting ' , imom , 'State Fits      \r',)
        PlotMassSFData(thistwoptdict,thisSmearList+thisSetList,imom)
        for ipof,(itpl,thistf) in enumerate(thisSetPoFLists):
            mprint( 'Plotting ' , imom , 'PoFSets ', str(int((ipof*100)/len(thisSetPoFLists))),'%        \r',)
            PlotMassData(thistwoptdict,thisSmearList+itpl,imom,TitleFlag=thistf)
        print 'Plotting ' , imom, 'Took: ' , str(datetime.timedelta(seconds=(time.time()-momstart))) ,' h:m:s                      '
    
    
def ReadAndPlotDict(thisGammaList,thisMomList,thisSetList,thisMethodList,thisCombList):
    datadict = ReadSetFitRFDict(outputdir,thisSetList,thisGammaList,thisMethodList,thisMomList=thisMomList)
    combdatadict = CreateDictOldCombs(datadict,thisCombList)
    thisGammaList = datadict.keys()
    if not CheckDict(datadict,'twopt','q = 0 0 0'):
        # raise IOError('Mass data dict not found')
        # print 'WARNING: mass dictionary not found, is needed for OSF and TSF'
    else:
        thisMassdict = datadict['twopt']['q = 0 0 0']
    start = time.time()
    for imom in thisMomList:
        # if imom == 'q = 0 0 0' and len(thisMomList) > 1 and DoMulticore: continue
        for icomb in combdatadict.keys():
            for icg,igamma in enumerate(combdatadict[icomb].keys()):
                combgamma = icomb+igamma
                if any([idst in igamma for idst in ['twopt']]): continue
                if any([idst in igamma for idst in ['doub','sing']]) and DoDS==False: continue
                gammastart = time.time()
                timeleft = GetTimeLeft(icg,len(thisGammaList),time.time()-start)
                # if not CheckDict(datadict,igamma,imom): continue
                # thisdatadict = datadict[igamma][imom]
                if not CheckDict(combdatadict,icomb,igamma,imom): continue
                thisdatadict = combdatadict[icomb][igamma][imom]
                if kappa == 12090:
                    progprint(0,time.time(),combgamma)
                    prevtime = time.time()
                    PlotTSinkData(thisdatadict,thisSetList,combgamma,imom)
                    progprint(1,prevtime,combgamma)
                    prevtime = time.time()
                    PlotTSinkData(thisdatadict,thisSetList,combgamma,imom,thissm='state1'+REvecTvarList[0])
                    progprint(2,prevtime,combgamma)
                    prevtime = time.time()
                    PlotTSinkData(thisdatadict,thisSetList,combgamma,imom,thissm='state1'+PoFTvarList[0])
                    progprint(2,prevtime,combgamma)
                    if 'SumMeth' in thisMethodList:
                        prevtime = time.time()
                        PlotTSinkSumData(thisdatadict,thisSetList,combgamma,imom)
                        progprint(3,prevtime,combgamma)
                    if 'TSFTsink' in thisMethodList:
                        prevtime = time.time()
                        PlotTSinkSFData(thisdatadict,thisMassdict,thisSetList,combgamma,imom,thisSF='TSFTsink')
                        progprint(4,prevtime,combgamma)
                    if 'TSFtest32' in thisMethodList:
                        prevtime = time.time()
                        PlotTSinkSFData(thisdatadict,thisMassdict,thisSetList,combgamma,imom,thisSF='TSFtest32')
                        progprint(5,prevtime,combgamma)
                    if 'TSFSmall' in thisMethodList:
                        prevtime = time.time()
                        PlotTSinkSFData(thisdatadict,thisMassdict,thisSetList,combgamma,imom,thisSF='TSFSmall')
                        progprint(6,prevtime,combgamma)
                    if 'OSFTsink' in thisMethodList:
                        prevtime = time.time()
                        PlotTSinkSFData(thisdatadict,thisMassdict,thisSetList,combgamma,imom,thisSF='OSFTsink')
                        progprint(7,prevtime,combgamma)

                    prevtime = time.time()
                    PlotCMData(thisdatadict,thisSetList,combgamma,imom)
                    progprint(8,prevtime,combgamma)
                    if 'OSFCM' in thisMethodList:
                        PlotCMOSFData(thisdatadict,thisMassdict,thisSetList,combgamma,imom)
                    if 'TSFCM' in thisMethodList:
                        PlotCMTSFData(thisdatadict,thisMassdict,thisSetList,combgamma,imom)
                elif kappa == 12104:
                    PlotCMData(thisdatadict,thisSetList,combgamma,imom)
                    # PlotCMSFData(thisdatadict,thisMassdict,['tsink29state1to18dt2'],combgamma,imom,thisSF='SFREvec')
                print 'Graphing Operator: ' +icomb + igamma + imom + ' took: ' , str(datetime.timedelta(seconds=(time.time()-gammastart))) ,' h:m:s                      '
                    


# if all('-m' not in iin for iin in sys.argv[1:]):
#     feedin = InputParams(sys.argv[1:] + ['-noprompt'] + ['-m=Fits,TSFTsink,TSFtest32,OSF'])
# else:
feedin = InputParams(sys.argv[1:] + ['-noprompt'])
    
thisGammaList = CreateGammaList(feedin['gamma'],twopt=True)

if thisGammaList == ['twopt']:
    # feedin['method'] = ['RF','OSFCM','TSFCM']
    TvarPicked = ['tsink29'+str(istate) for istate in CreateMassSet([],['1'],DefTvarPicked)]
    # TvarPicked = ['tsink29'+str(istate) for istate in CreateMassSet([],['1'],DefTvarList)]
    # TvarPicked = ['tsink29'+str(istate) for istate in CreateMassSet([],['1'],DefTvarDt1)]
    # TvarPicked = ['tsink29'+str(istate) for istate in CreateMassSet([],['1'],DefTvarDt2)]
    # TvarPicked = ['tsink29'+str(istate) for istate in CreateMassSet([],['1'],DefTvarDt3)]
    # TvarPicked = ['tsink29'+str(istate) for istate in CreateMassSet([],['1'],DefTvarDt4)]
    # TvarPicked = ['tsink29'+str(istate) for istate in CreateMassSet([],['1'],DefTvarto17)]
    # TvarPicked = ['tsink29'+str(istate) for istate in CreateMassSet([],['1'],DefTvarto18)]
    # TvarPicked = ['tsink29'+str(istate) for istate in CreateMassSet([],['1'],DefTvarto19)]
    # TvarPicked = ['tsink29'+str(istate) for istate in CreateMassSet([],['1'],DefTvarto20)]
    TvarPicked += ['tsink26'+str(istate) for istate in CreateMassSet([],['1'],PoFTvarList)]
    thisSmList = ['tsink29'+str(ism) for ism in CreateMassSet(DefSmearList,['1'],[])]
    TvarLists = []
    # TvarLists = [(['tsink29'+str(istate) for istate in CreateMassSet([],['1'],DefTvarDt1)],'PoFDt1'),
    #              (['tsink29'+str(istate) for istate in CreateMassSet([],['1'],DefTvarDt2)],'PoFDt2'),
    #              (['tsink29'+str(istate) for istate in CreateMassSet([],['1'],DefTvarDt3)],'PoFDt3'),
    #              (['tsink29'+str(istate) for istate in CreateMassSet([],['1'],DefTvarDt4)],'PoFDt4'),
    #              (['tsink29'+str(istate) for istate in CreateMassSet([],['1'],DefTvarto16)],'PoFto16'),
    #              (['tsink29'+str(istate) for istate in CreateMassSet([],['1'],DefTvarto17)],'PoFto17'),
    #              (['tsink29'+str(istate) for istate in CreateMassSet([],['1'],DefTvarto18)],'PoFto18'),
    #              (['tsink29'+str(istate) for istate in CreateMassSet([],['1'],DefTvarto19)],'PoFto19'),
    #              (['tsink29'+str(istate) for istate in CreateMassSet([],['1'],DefTvarto20)],'PoFto20')]
    thisAllSetList = thisSmList+TvarPicked
    print 'AllSetList:\n' + '\n'.join(thisAllSetList)
    print ''
    print 'MethodList:\n' + '\n'.join(feedin['method'])
    if DoMulticore and len(feedin['mom']) > 1:
        inputparams = []
        for imom in feedin['mom']:
            if imom == 'q = 0 0 0':
                inputparams.append(([imom],thisSmList,TvarPicked,TvarLists,feedin['method']))
            else:
                inputparams.append((['q = 0 0 0',imom],thisSmList,TvarPicked,TvarLists,feedin['method']))                
        makeContextFunctions(ReadAndPlotMass)
        thisPool = Pool(min(len(inputparams),feedin['anaproc']))
        output = thisPool.map(ReadAndPlotMass.mapper,inputparams)
        thisPool.close()
        thisPool.join()
    else:
        ReadAndPlotMass(feedin['mom'],thisSmList,TvarPicked,TvarLists,feedin['method'])
else:
    # if any([imom != 'q = 0 0 0' for imom in feedin['mom']]):
    #     feedin['method'] = ['RF']
    print 'MethodList:\n' + '\n'.join(feedin['method'])
    print ''
    print 'thisSetList:\n' + '\n'.join(feedin['set'])
    if DoMulticore and len(thisGammaList) > 1 and feedin['anaproc'] > 1:
        inputparams = []
        for igamma in thisGammaList:
            if any([idst in igamma for idst in ['doub','sing','twopt']]): continue
            for imom in feedin['mom']:
                if imom == 'q = 0 0 0' :
                    inputparams.append((['doub'+igamma,'sing'+igamma,'twopt'],[imom],feedin['set'],feedin['method'],feedin['comb'])) 
                else:
                    inputparams.append((['doub'+igamma,'sing'+igamma,'twopt'],['q = 0 0 0',imom],feedin['set'],feedin['method'],feedin['comb']))
        makeContextFunctions(ReadAndPlotDict)
        thisPool = Pool(min(len(inputparams),feedin['anaproc']))
        output = thisPool.map(ReadAndPlotDict.mapper,inputparams)
        thisPool.close()
        thisPool.join()
    else:
        passgammalist = ['twopt']
        for igamma in thisGammaList:
            if any([idst in igamma for idst in ['doub','sing','twopt']]): continue
            passgammalist += ['doub'+igamma,'sing'+igamma]
        ReadAndPlotDict(passgammalist,feedin['mom'],feedin['set'],feedin['method'],feedin['comb'])
        
print 'Graphing all complete'
    

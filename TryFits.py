#!/usr/bin/env python


from OppFuns import CreateGammaList,WipeSet
from Params import *
import numpy as np
from BootTest import BootStrap1
from ReadTxt import *
from SetLists import *
from Fitting import FitRFSet
from OutputData import PrintFitSetToFile
from FitParams import *
from multiprocessing import Pool
from MultiWrap import *
from FFParams import *
import itertools as it
import sys
import time
import datetime
from InputArgs import *
from CheckXml import *



def TryFitsFun(thisGammaList,thisSetList,thisReadMomList,thisTSinkList,thischunk):
    # for ig,gammadata in enumerate(dataRF):
    #     for im,momdata in enumerate(gammadata):
    #         for iset,setdata in enumerate(momdata):
    #             print thisGammaList[ig] , ' ' , thisGammaMomList[thisGammaList[ig]] , ' ' , SetList[iset]
    #             print Pullflag(setdata,'Avg')
    #             print ''
    # dataRF = [ gamma , mom , set , iflow , it ] bs IF reading flow results
    # dataRF = [ gamma , mom , set , it ] bs
    [dataRF,data2pt,thisGammaMomList,BorA,infolistRF,infolist2pt,flowlist] = ReadRFnp(thisGammaList,thisSetList,thisMomList=thisReadMomList)
    start = time.time()
    FitDataBoot,FitDataAvg,FitDataChi = [],[],[]
    for igamma,(thisgamma,thismomlist) in enumerate(thisGammaMomList.iteritems()):
        if Debug: print 'thisgamma', thisgamma
        FitDataBoot.append([])
        FitDataAvg.append([])
        FitDataChi.append([])
        for imom,thismom in enumerate(thismomlist):
            # print 'Fitting ' , thisgamma , thismom , '       \r',
            if Debug: print '     thismom', thismom
            start = time.time()
            FitDataBoot[igamma].append([])
            FitDataAvg[igamma].append([])
            FitDataChi[igamma].append([])
            if len(flowlist) > 0:
                for icf,iflow in enumerate(flowlist):
                    FitDataBoot[igamma][imom].append([])
                    FitDataAvg[igamma][imom].append([])
                    FitDataChi[igamma][imom].append([])
                    for icut in FitCutList:
                        momdata = np.array(dataRF[igamma][imom])[:,icf,:]
                        dataoutBoot,dataoutAvg,dataoutChi = FitRFSet(momdata,thisTSinkList,icut)
                        FitDataBoot[igamma][imom][icf].append(dataoutBoot)
                        FitDataAvg[igamma][imom][icf].append(dataoutAvg)
                        FitDataChi[igamma][imom][icf].append(dataoutChi)
            else:
                for icut in FitCutList:
                    momdata = np.array(dataRF[igamma][imom])
                    dataoutBoot,dataoutAvg,dataoutChi = FitRFSet(momdata,thisTSinkList,icut)
                    FitDataBoot[igamma][imom].append(dataoutBoot)
                    FitDataAvg[igamma][imom].append(dataoutAvg)
                    FitDataChi[igamma][imom].append(dataoutChi)
            print thisgamma , thismom , ' Took: ',GetTimeStr(time.time()-start)
    #FitData = [ igamma , ip , icut , iset ]
    #FitData = [ igamma , ip , iflow , icut , iset ]
    return FitDataBoot,FitDataChi,thisGammaMomList,thisSetList,FitCutList,infolistRF,flowlist



feedin = InputParams(sys.argv[1:])
thisGammaList = CreateGammaList(feedin['gamma'])

ShowSetLists(feedin['set'])

totstart = time.time()
inputparams = []
# RunGammaList = []
for igamma in thisGammaList:
    print 'adding to que: ' , igamma
    if 'doub' not in igamma and 'sing' not in igamma and 'twopt' not in igamma:    
        parsegammalist = ['doub'+igamma,'sing'+igamma]
        for iChunk,(iSet,iTS) in enumerate(zip(feedin['set'],DefTSinkSetList)):
            if DefWipe:
                thisMomList = GetMomFromGamma(igamma,thisMomList=feedin['mom'])
            else:
                thisMomList = Check3ptArray(parsegammalist,[iSet],thisMomList=feedin['mom'],CheckType='Fits',printout=False)
                thisMomList = thisMomList[igamma][iSet]
            for imom in thisMomList:
                # RunGammaList.append(iDS+igamma)
                inputparams.append((parsegammalist,[iSet],[imom],[iTS],(iChunk*100)/float(len(feedin['set']))))
    elif igamma.replace('sing','').replace('doub','') not in thisGammaList and 'twopt' not in igamma:
        parsegammalist = [igamma]
        for iChunk,(iSet,iTS) in enumerate(zip(feedin['set'],DefTSinkSetList)):
            if DefWipe:
                thisMomList = GetMomFromGamma(igamma,thisMomList=feedin['mom'])
            else:
                thisMomList = Check3ptArray(parsegammalist,[iSet],thisMomList=feedin['mom'],CheckType='Fits',printout=False)
                thisMomList = thisMomList[igamma.replace('sing','').replace('doub','')][iSet]
            for imom in thisMomList:
                # RunGammaList.append(iDS+igamma)
                inputparams.append((parsegammalist,[iSet],[imom],[iTS],(iChunk*100)/float(len(feedin['set']))))

if len(inputparams) > 0:
    if DoMulticore and feedin['anaproc'] > 1 and len(inputparams) > 1:
        print 'Running Multicore'
        makeContextFunctions(TryFitsFun)
        thisPool = Pool(min(len(inputparams),feedin['anaproc']))
        output = thisPool.map(TryFitsFun.mapper,inputparams)
        thisPool.close()
        thisPool.join()
    else:
        print 'Running Single core'
        output = []
        for icount,iin in enumerate(inputparams):
            output.append(TryFitsFun(*iin))
            print int((icount*100)/float(len(inputparams))) , '% done' + ' '*50 + '\r',
else:
    print 'nothing to calculate'        
    output = []
    
# WipeSet(outputdir[0],RunGammaList,feedin['set'],filepref='Fits/')
print 'Done Fits ' + ' '*50
for iout in output:
    PrintFitSetToFile(*iout)

print 'Total fit time took: ' , str(datetime.timedelta(seconds=time.time()-totstart)) , ' h:m:s '

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
    # dataRF = [ gamma , mom , set , it ] bs
    # for ig,gammadata in enumerate(dataRF):
    #     for im,momdata in enumerate(gammadata):
    #         for iset,setdata in enumerate(momdata):
    #             print thisGammaList[ig] , ' ' , thisGammaMomList[thisGammaList[ig]] , ' ' , SetList[iset]
    #             print Pullflag(setdata,'Avg')
    #             print ''
    [dataRF,data2pt,thisGammaMomList,BorA,infolistRF,infolist2pt] = ReadRFnp(thisGammaList,thisSetList,thisMomList=thisReadMomList)
    start = time.time()
    FitDataBoot,FitDataAvg,FitDataChi = [],[],[]
    for igamma,(thisgamma,thismomlist) in enumerate(thisGammaMomList.iteritems()):
        FitDataBoot.append([])
        FitDataAvg.append([])
        FitDataChi.append([])
        for imom,thismom in enumerate(thismomlist):
            # print 'Fitting ' , thisgamma , thismom , '       \r',
            FitDataBoot[igamma].append([])
            FitDataAvg[igamma].append([])
            FitDataChi[igamma].append([])
            for icut in FitCutList:
                momdata = dataRF[igamma][imom]
                dataoutBoot,dataoutAvg,dataoutChi = FitRFSet(momdata,thisTSinkList,icut)
                FitDataBoot[igamma][imom].append(dataoutBoot)
                FitDataAvg[igamma][imom].append(dataoutAvg)
                FitDataChi[igamma][imom].append(dataoutChi)
    #FitData = [ igamma , ip , icut , iset ]
    # print ' '.join(thisGammaMomList.keys()) , ' at ' , thischunk,'% took: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s '
    return FitDataBoot,FitDataChi,thisGammaMomList,thisSetList,FitCutList,infolistRF



feedin = InputParams(sys.argv[1:])
thisGammaList = CreateGammaList(feedin['gamma'])

ShowSetLists(feedin['set'])

totstart = time.time()
inputparams = []
# RunGammaList = []
for igamma in thisGammaList:
    if 'doub' not in igamma and 'sing' not in igamma:
        print 'adding to que: ' , igamma
        for iChunk,(iSet,iTS) in enumerate(zip(feedin['set'],DefTSinkSetList)):
            if DefWipe:
                thisMomList = feedin['mom']
            else:
                thisMomList = Check3ptArray(['doub'+igamma,'sing'+igamma,igamma],[iSet],thisMomList=feedin['mom'],CheckType='Fits',printout=False)
                thisMomList = thisMomList[igamma][iSet]
            for imom in thisMomList:
                # RunGammaList.append(igamma)
                inputparams.append((['doub'+igamma,'sing'+igamma,igamma],[iSet],[imom],[iTS],(iChunk*100)/float(len(feedin['set']))))


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
    
# WipeSet(outputdir,RunGammaList,feedin['set'],filepref='Fits/')
print 'Done Fits ' + ' '*50
for iout in output:
    PrintFitSetToFile(*iout)

print 'Total fit time took: ' , str(datetime.timedelta(seconds=time.time()-totstart)) , ' h:m:s '

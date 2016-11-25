#!/usr/bin/env python
import numpy as np
import sys
from Params import *
from FitParams import *
from ReadTxt import ReadRFnp
from SumMeth import *
from SetLists import *
from OppFuns import CreateGammaList,WipeSet
from OutputData import PrintSumSetToFile
from FFParams import *
import time,datetime
import copy
from InputArgs import *
from XmlFormatting import *
from CheckXml import *


def FitSumWrap(thisGammaList,thisReadSetList,thisTSinkList,this2ptSetList,thisReadMomList):
    thisstart = time.time()
    [dataRF,data2pt,thisGammaMomList,BorA,infoRF,info2pt] = ReadRFnp(thisGammaList,thisReadSetList,thisMomList=thisReadMomList)
    # dataRF = [ gamma , mom , set , it ] bs
    if 'twopt' in thisGammaMomList.keys(): del thisGammaMomList['twopt']    
    SummedRF,SumFitBoot,SumFitAvg,SumFitChi,SumFitList = [],[],[],[],[]
    for (thisgamma,thismomlist),gammadata in zip(thisGammaMomList.iteritems(),dataRF):
        thisout = FitSum(SumMeth3ptCuts,thisTSinkList,thisgamma,thismomlist,gammadata)
        SummedRF.append(thisout[0])
        SumFitBoot.append(thisout[1])
        SumFitAvg.append(thisout[2])
        SumFitChi.append(thisout[3])
        SumFitList.append(thisout[4])
    WipeSet(outputdir[0],thisGammaMomList.keys(),this2ptSetList,filepref='SumMeth/')
    PrintSumSetToFile(SummedRF,SumFitBoot,SumFitChi,SumFitList,str(this2ptSetList[0]),thisGammaMomList,thisTSinkList,SumMeth3ptCuts,infoRF)
    print 'Fitting ' , ' '.join(thisGammaList) , 'complete, took: ' + GetTimeStr(time.time()-thisstart)



# if len(sys.argv) < 2: raise IOError("Enter gamma matrix stuff (See OppFuns.py : CreateGammaList)")
feedin = InputParams(sys.argv[1:])



thisGammaList = CreateGammaList(feedin['gamma'])

ParsTSinkList = []
for ict,its in enumerate(AllTSinkList):
    if any(['tsink'+str(its) in rts for rts in feedin['set']]):
        ParsTSinkList.append(its)

# thisReadSetList = feedin['set']
# thisTSinkList,thisReadSet2pt = GetTsinkSmLists(feedin['set'])
# thisTSinkList = map(unxmlTSink,thisTSinkList)

[thisReadSetList,thisReadSet2pt,thisTSinkList] = CreateSet(thisSingSmearL=SingSmearList,thisTSinkL=ParsTSinkList,
                                                           thisSmearL=[],thisREvecTvarL=[],thisPoFTvarL=[],thisTvarL=[])

ShowSetLists(thisReadSetList)

start = time.time()
inputparams = []
for igamma in thisGammaList:
    if 'doub' not in igamma and 'sing' not in igamma and 'twopt' not in igamma:
        parsegammalist = ['doub'+igamma,'sing'+igamma]
        if DefWipe:
            QueMomList = GetMomFromGamma(igamma,thisMomList=feedin['mom'])
        else:
            QueMomList = Check3ptAllSets(parsegammalist,thisReadSetList,thisMomList=feedin['mom'],CheckType='SumMeth')
            QueMomList = QueMomList[igamma]
        for imom in QueMomList:
            print 'adding to que: ' , igamma , imom 
            inputparams.append((parsegammalist,thisReadSetList,thisTSinkList,thisReadSet2pt,[imom]))
    elif igamma.replace('sing','').replace('doub','') not in thisGammaList and 'twopt' not in igamma:
        parsegammalist = [igamma]
        if DefWipe:
            QueMomList = GetMomFromGamma(igamma,thisMomList=feedin['mom'])
        else:
            QueMomList = Check3ptAllSets([igamma],thisReadSetList,thisMomList=feedin['mom'],CheckType='SumMeth')
            QueMomList = QueMomList[igamma.replace('sing','').replace('doub','')]
        for imom in QueMomList:
            print 'adding to que: ' , igamma , imom
            inputparams.append((parsegammalist,thisReadSetList,thisTSinkList,thisReadSet2pt,[imom]))
makeContextFunctions(FitSumWrap)

if len(inputparams) > 0:
    if DoMulticore and feedin['anaproc'] > 1 and len(inputparams) > 1:
        print 'Running multiprocessor fits'
        thisPool = Pool(processes=min(len(inputparams),feedin['anaproc']))
        thisPool.map(FitSumWrap.mapper,inputparams)
        thisPool.close()
        thisPool.join()
    else:
        print 'Running single processor fits'
        for ipar in inputparams : FitSumWrap(*ipar)
else:
    print 'nothing to calculate'
print 'All Done, time taken:' , GetTimeStr(time.time()-start)

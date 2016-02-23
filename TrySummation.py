#!/usr/bin/env python
import numpy as np
import sys
from Params import *
from FitParams import *
from ReadTxt import ReadRFnp
from SumMeth import *
from SetLists import *
from OutputData import PrintTSFMassToFile,PrintTSFToFile,PrintOSFMassToFile,PrintOSFToFile
from CreateCombs import MakeUmD
from OppFuns import CreateGammaList,WipeSet
from OutputData import PrintSumSetToFile
from FFParams import *
import time,datetime
import copy




def FitSumWrap(thisGammaList,thisReadSetList,thisTSinkList,this2ptSetList,thisReadMomList):
    thisstart = time.time()
    [dataRF,data2pt,thisGammaMomList,BorA] = ReadRFnp(thisGammaList,thisReadSetList,thisMomList=thisReadMomList)
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
    WipeSet(outputdir,thisGammaMomList.keys(),this2ptSetList,filepref='SumMeth/')
    PrintSumSetToFile(SummedRF,SumFitBoot,SumFitChi,SumFitList,'sm'+str(this2ptSetList[0]),thisGammaMomList,thisTSinkList,SumMeth3ptCuts)
    print 'Fitting ' , ' '.join(thisGammaList) , 'complete, took: ' + GetTimeStr(time.time()-thisstart)



# if len(sys.argv) < 2: raise IOError("Enter gamma matrix stuff (See OppFuns.py : CreateGammaList)")
ReadSmearList = ['32']
ReadTSinkList = [26,29,32,35,38]
thisGammaList = CreateGammaList(sys.argv[1:])
# thisReadMomList = ['q = 0 0 0']
thisReadMomList = qvecSet
[thisReadSetList,thisReadSet2pt,thisTSinkList] = CreateSet(thisSmearL=ReadSmearList,thisTvarL=[],thisTSinkL=ReadTSinkList)

start = time.time()
print 'Running multiprocessor fits'
inputparams = []
for igamma in thisGammaList:
    if 'doub' not in igamma and 'sing' not in igamma and 'twopt' not in igamma:
        inputparams.append((['doub'+igamma,'sing'+igamma,igamma],thisReadSetList,thisTSinkList,thisReadSet2pt,thisReadMomList))
    elif igamma.replace('sing','').replace('doub','') not in thisGammaList and 'twopt' not in igamma:
        inputparams.append(([igamma],thisReadSetList,thisTSinkList,thisReadSet2pt,thisReadMomList))        
makeContextFunctions(FitSumWrap)

if DoMulticore:
    thisPool = Pool(processes=min(len(inputparams),AnaProc))
    thisPool.map(FitSumWrap.mapper,inputparams)
    thisPool.close()
    thisPool.join()
else:
    for ipar in inputparams : FitSumWrap(*ipar)
print 'All Done, time taken:' , GetTimeStr(time.time()-start)

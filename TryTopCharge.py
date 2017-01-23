#!/usr/bin/env python

from array import array
import os
import numpy as np
import sys
from Params import *
from ReadCMCfuns import ReadSetAlpha,ReadListAlpha
from CMSTech import CreateCM2ptCfuns,CreatePoF2ptCfuns,PreptwoptCorr
from OutputData import PrintSetToFile,PrintCfunToFile,PrintAlphaSetToFile
from OutputXmlData import PrintLREvecMassToFile
from CreateCombs import CreategiDi
from SetLists import CreateMassSet
from MiscFuns import *
import time,datetime
from MultiWrap import *
from multiprocessing import Pool

def CreateTwoPtTop(thisMomList,thisSmearList,feedin= {'anaproc':AnaProc}):
    logfile = logdir+'LogTopCharge.log'
    errfile = logdir+'LogTopCharge.log'
    touch(logfile)
    touch(errfile)
    sys.stdout = open(logfile,'a',0)
    sys.stderr = open(errfile,'a',0)
    # print 'Running ' + ipTOqstr(thisMomList[0]) + ' ' +  str(int((thisMomList[0]*100)/float(len(qvecSet))))+'%' 

    if 'ReadList' in ListOrSet:
        [data2pt,dataTop,thisTopList,filelist] = ReadListAlpha(thisSmearList,thisMomList,conflist,Interps=DefInterpList,thistsourceList=PoFtsourceList)
    elif 'ReadSet' in ListOrSet:
        [data2pt,dataTop,thisTopList,filelist] = ReadSetAlpha(thisSmearList,thisMomList,dirread,Interps=DefInterpList,thistsourceList=PoFtsourceList)

    thisMomList = GetAvgMomListip(thisMomList)
    data2pt = np.array(PreptwoptCorr(np.array(data2pt)))
    ncon = np.size(filelist)
    InfoDict = {'nconfig':ncon}
    print 'ncon = ' + str(ncon)
    # print 'nboot = ' + str(nboot)
    ## data2pt = [ t_src, ism , jsm , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
    ## dataTop = [ t_flow, t_src, ism , jsm , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
    C2out = DiagSmearWithTsrc(data2pt).tolist()
    C2outTop = []
    dataTopHold = []
    for topdata in dataTop:
        dataTopHold.append(np.array(PreptwoptCorr(np.array(topdata))))
    dataTop = np.array(dataTopHold)
    for topdata in dataTop:
        C2outTop.append(DiagSmearWithTsrc(topdata).tolist())
    ## C2out = [ t_src*ism*ism , ip , it ] 
    ## C2outTop = [ tflow , t_src*ism*ism , ip , it ] 
    
    start = time.time()
    CMinputparams,PoFinputparams = [],[]
    makeContextFunctions(CreatePoF2ptCfuns)
    makeContextFunctions(CreateCM2ptCfuns)
    for icount,(itodt,itodtPoF) in enumerate(zip(DeftodtPicked,DefPoFVarList)):
        CMinputparams.append((data2pt[0],itodt,thisMomList))
        PoFinputparams.append((data2pt,itodtPoF,thisMomList))


    if DoMulticore and feedin['anaproc'] > 1 and len(PoFinputparams) > 1 and len(CMinputparams) > 1:
        thisPool = Pool(min(len(CMinputparams),feedin['anaproc']))
        outputPoF = thisPool.map(CreatePoF2ptCfuns.mapper,PoFinputparams)
        if len(thisSmearList) > 1 and DoCM:
            outputCM = thisPool.map(CreateCM2ptCfuns.mapper,CMinputparams)
        thisPool.close()
        thisPool.join()
    else:
        outputPoF,outputCM = [],[]
        for iin in PoFinputparams: outputPoF.append(CreatePoF2ptCfuns.mapper(iin))
        if len(thisSmearList) > 1 and DoCM:
            for iin in CMinputparams: outputCM.append(CreateCM2ptCfuns.mapper(iin))

    thisPoFTvarList = PoFTvarPicked
    if len(thisSmearList) > 1 and DoCM:
        thisCMTvarList = DefTvarPicked
        for iout,iTvar in zip(outputCM,thisCMTvarList):
            [CMdata2pt,LEvec,REvec,Emass] = iout
            ## CMdata2pt [ istate , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
            C2out += CMdata2pt.tolist()
            PrintLREvecMassToFile(LEvec,REvec,Emass,thisMomList,iTvar,AddDict=InfoDict,DoPoF=False)
    else:
        thisCMTvarList = []
        
    for iout,iTvar in zip(outputPoF,thisPoFTvarList):
        [CMdata2pt,LEvec,REvec,Emass] = iout
        ## CMdata2pt [ istate , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
        C2out += CMdata2pt.tolist()
        PrintLREvecMassToFile(LEvec,REvec,Emass,thisMomList,iTvar,AddDict=InfoDict,DoPoF=True)

    SetList = []
    SetList += CreateMassSet(thisSmearList,StateSet,[],tsrclist=PoFtsourceList,flipord=True)
    if len(thisSmearList) > 1 and DoCM: SetList += CreateMassSet([],CMStateSet,thisCMTvarList,flipord=True)
    SetList += CreateMassSet([],StateSet,thisPoFTvarList,flipord=True)
    PrintCfunToFile([C2out],SetList,thisMomList,['twopt'],AddDict=InfoDict,Top=True)


    CMinputparams,PoFinputparams = [],[]    
    for itop,topdata in enumerate(dataTop):
        # print 'CM Technique for Tflow=',itop,' \r'
        # for icount,itodt in enumerate(DeftodtList):
        for icount,(itodt,itodtPoF) in enumerate(zip(DeftodtPicked,DefPoFVarList)):
            CMinputparams.append((topdata[0],itodt,thisMomList))
            PoFinputparams.append((topdata,itodtPoF,thisMomList))

            
    if DoMulticore and feedin['anaproc'] > 1:
        thisPool = Pool(min(len(CMinputparams),feedin['anaproc']))
        outputPoFTop = thisPool.map(CreatePoF2ptCfuns.mapper,PoFinputparams)
        if len(thisSmearList) > 1 and DoCM:
            outputCMTop = thisPool.map(CreateCM2ptCfuns.mapper,CMinputparams)
        thisPool.close()
        thisPool.join()
    else:
        outputPoFTop,outputCMTop = [],[]
        for iin in PoFinputparams: outputPoFTop.append(CreatePoF2ptCfuns.mapper(iin))
        if len(thisSmearList) > 1 and DoCM:
            for iin in CMinputparams: outputCMTop.append(CreateCM2ptCfuns.mapper(iin))
            
    
    if len(thisSmearList) > 1 and DoCM:
        thisCMTvarList = DefTvarPicked
        for itop,CMTop in enumerate(outputCMTop):
            [CMdata2pt,LEvec,REvec,Emass] = CMTop
            ## CMdata2pt [ istate , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
            C2outTop[itop] += CMdata2pt.tolist()
    else:
        thisCMTvarList = []

        
    thisPoFTvarList = PoFTvarPicked
    for itop,PoFTop in enumerate(outputPoFTop):
        [PoFdata2pt,LEvec,REvec,Emass] = PoFTop
        ## CMdata2pt [ istate , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
        C2outTop[itop] += PoFdata2pt.tolist()

    
    print 'CMTech Total Time Taken: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s                                         '
    start = time.time()
    print 'Printing to file \r',

    SetList = []
    SetList += CreateMassSet(thisSmearList,StateSet,[],tsrclist=PoFtsourceList,flipord=True)
    if len(thisSmearList) > 1 and DoCM: SetList += CreateMassSet([],CMStateSet,thisCMTvarList,flipord=True)
    SetList += CreateMassSet([],StateSet,thisPoFTvarList,flipord=True)
    PrintAlphaSetToFile(np.swapaxes(np.array(C2outTop),0,1),C2out,SetList,thisMomList,thisTopList,AddDict=InfoDict)

    print 'Printing took ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s  '
    # print 'Completed ' + ipTOqstr(thisMomList[0])

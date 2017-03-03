#!/usr/bin/env python

from array import array
import os
import numpy as np
import sys
from Params import *
from ReadCMCfuns import ReadSet,ReadList
from CMSTech import CreateCM2ptCfuns,CreatePoF2ptCfuns,PreptwoptCorr
from OverCMTech import *
from OutputData import PrintSetToFile,PrintCfunToFile
from OutputXmlData import PrintLREvecMassToFile
from CreateCombs import CreategiDi
from SetLists import CreateMassSet,CreateOverDetTvarSet
from MiscFuns import *
import time,datetime
from MultiWrap import *
from multiprocessing import Pool

def CreateTwoPt(thisMomList,thisiSmearList,thisjSmearList,feedin= {'anaproc':AnaProc}):
    logfile = logdir+'LogTwoPt.log'
    errfile = logdir+'LogTwoPt.log'
    touch(logfile)
    touch(errfile)
    sys.stdout = open(logfile,'a',0)
    sys.stderr = open(errfile,'a',0)
    # print 'Running ' + ipTOqstr(thisMomList[0]) + ' ' +  str(int((thisMomList[0]*100)/float(len(qvecSet))))+'%' 

    if 'ReadList' in ListOrSet:
        [data2pt,data3pt,filelist] = ReadList(thisiSmearList,thisjSmearList,thisMomList,{},{},[],[],
                                              conflist,[],Interps=DefInterpList,thistsourceList=PoFtsourceList)
    elif 'ReadSet' in ListOrSet:
        [data2pt,data3pt,filelist] = ReadSet(thisiSmearList,thisjSmearList,thisMomList,{},{},[],[],
                                             dirread,[],Interps=DefInterpList,thistsourceList=PoFtsourceList)

    thisMomList = GetAvgMomListip(thisMomList)
    data2pt = np.array(PreptwoptCorr(np.array(data2pt)))
    ncon = np.size(filelist)
    InfoDict = {'nconfig':ncon}
    print 'ncon = ' + str(ncon)
    # print 'nboot = ' + str(nboot)
    ## data2pt = [ t_src, ism , jsm , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
    C2out = FlattenSmearWithTsrc(data2pt).tolist()
    ## C2out = [ t_src*ism*ism , ip , it ] 
    
    if len(thisiSmearList) != len(thisjSmearList): print 'Warning: source and sink smearing lists are different sizes, CM analysis skipped'
    start = time.time()
    CMinputparams,PoFinputparams = [],[]
    makeContextFunctions(CreatePoF2ptCfuns)
    makeContextFunctions(CreateCM2ptCfuns)
    for icount,itodt in enumerate(DeftodtList):
        CMinputparams.append((data2pt[0],itodt,thisMomList))
        PoFinputparams.append((data2pt,itodt,thisMomList,DefPoFTvarRef))


    if DoMulticore and feedin['anaproc'] > 1:
        thisPool = Pool(min(len(CMinputparams),feedin['anaproc']))
        outputPoF = thisPool.map(CreatePoF2ptCfuns.mapper,PoFinputparams)
        if len(thisiSmearList) == len(thisjSmearList) and DoCM and len(thisiSmearList) > 1:
            outputCM = thisPool.map(CreateCM2ptCfuns.mapper,CMinputparams)
        thisPool.close()
        thisPool.join()
    else:
        outputPoF,outputCM = [],[]
        for iin in PoFinputparams: outputPoF.append(CreatePoF2ptCfuns.mapper(iin))
        if len(thisiSmearList) == len(thisjSmearList) and DoCM and len(thisiSmearList) > 1:
            for iin in CMinputparams: outputCM.append(CreateCM2ptCfuns.mapper(iin))
    
    
    thisPoFTvarList = ['PoF'+str(PoFShifts)+iTvar for iTvar in TwoPtDefTvarList]
    if len(thisiSmearList) == len(thisjSmearList) and DoCM and len(thisiSmearList) > 1:
        thisCMTvarList = ['CM'+iTvar for iTvar in TwoPtDefTvarList]
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


    print 'CMTech Total Time Taken: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s  '
    start = time.time()
    print 'Printing to file \r',

    SetList = []
    SetList += CreateMassSet(thisiSmearList,thisjSmearList,StateSet,[],tsrclist=PoFtsourceList,flipord=True)
    if len(thisiSmearList) == len(thisjSmearList) and DoCM and len(thisiSmearList) > 1: SetList += CreateMassSet([],[],CMStateSet,thisCMTvarList,flipord=True)
    SetList += CreateMassSet([],[],StateSet,thisPoFTvarList,flipord=True)
    PrintCfunToFile([C2out],SetList,thisMomList,['twopt'],AddDict=InfoDict)
    PrintSetToFile([C2out],SetList,thisMomList,['Mass'],0,AddDict=InfoDict)
    print 'Printing took ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s  '
    # print 'Completed ' + ipTOqstr(thisMomList[0])


def CreateTwoPtOverDet(thisMomList,thisiSmearList,thisjSmearList,feedin= {'anaproc':AnaProc}):
    logfile = logdir+'LogTwoPt.log'
    errfile = logdir+'LogTwoPt.log'
    touch(logfile)
    touch(errfile)
    sys.stdout = open(logfile,'a',0)
    sys.stderr = open(errfile,'a',0)
    # print 'Running ' + ipTOqstr(thisMomList[0]) + ' ' +  str(int((thisMomList[0]*100)/float(len(qvecSet))))+'%' 

    if 'ReadList' in ListOrSet:
        [data2pt,data3pt,filelist] = ReadList(thisiSmearList,thisjSmearList,thisMomList,{},{},[],[],
                                              conflist,[],Interps=DefInterpList,thistsourceList=PoFtsourceList)
    elif 'ReadSet' in ListOrSet:
        [data2pt,data3pt,filelist] = ReadSet(thisiSmearList,thisjSmearList,thisMomList,{},{},[],[],
                                             dirread,[],Interps=DefInterpList,thistsourceList=PoFtsourceList)

    thisMomList = GetAvgMomListip(thisMomList)
    data2pt = np.array(PreptwoptCorr(np.array(data2pt)))
    ncon = np.size(filelist)
    InfoDict = {'nconfig':ncon}
    print 'ncon = ' + str(ncon)
    # print 'nboot = ' + str(nboot)
    ## data2pt = [ t_src, ism , jsm , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
    C2out = FlattenSmearWithTsrc(data2pt).tolist()
    ## C2out = [ t_src*ism*jsm , ip , it ] 

    makeContextFunctions(CreatePoF2ptCfuns)
    for thistvar,(itomin,itomax,idt) in zip(DefTvarPicked,OverDettodtlist):
        PoFinputparams.append(data2pt,(itomin,itomax),idt,thisMomList)

    ##TO DO WORKING FROM HERE DEBUG BLAHHHH##
    if DoMulticore and feedin['anaproc'] > 1:
        thisPool = Pool(min(len(PoFinputparams),feedin['anaproc']))
        outputPoF = thisPool.map(CreatePoF2ptCfuns.mapper,PoFinputparams)
        thisPool.close()
        thisPool.join()
    else:
        outputPoF,outputCM = [],[]
        for iin in PoFinputparams: outputPoF.append(CreatePoF2ptCfuns.mapper(iin))
    
    
    thisPoFTvarList = ['PoF'+str(PoFShifts)+iTvar for iTvar in TwoPtDefTvarList]
        
    for iout,iTvar in zip(outputPoF,thisPoFTvarList):
        [CMdata2pt,LEvec,REvec,Emass] = iout
        ## CMdata2pt [ istate , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
        C2out += CMdata2pt.tolist()
        PrintLREvecMassToFile(LEvec,REvec,Emass,thisMomList,iTvar,AddDict=InfoDict,DoPoF=True)


    
    if len(thisiSmearList) != len(thisjSmearList): print 'Warning: source and sink smearing lists are different sizes, CM analysis skipped'
    start = time.time()
    for thistvar,(itomin,itomax,idt) in zip(DefTvarPicked,OverDettodtlist):
        print 'Starting Overdetumined Eigenvalue problem for '+ thistvar
        if Debug: print 'Check actual todt params: to',itomin,'-',itomax,'dt',idt 
        [CMdata2pt,LEvec,REvec,Emass] = CreateOvdet2ptCfun(data2pt,(itomin,itomax),idt,thisMomList)
        ## CMdata2pt [ istate , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
        C2out += CMdata2pt.tolist()
        PrintLREvecMassToFile(LEvec,REvec,Emass,thisMomList,thistvar,AddDict=InfoDict,DoPoF=False)


    print 'CMTech Total Time Taken: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s  '
    start = time.time()
    print 'Printing to file \r',

    SetList = []
    SetList += CreateMassSet(thisiSmearList,thisjSmearList,StateSet,[],tsrclist=PoFtsourceList,flipord=True)
    SetList += CreateOverDetTvarSet(StateSet,DefTvarPicked)
    PrintCfunToFile([C2out],SetList,thisMomList,['twopt'],AddDict=InfoDict)
    PrintSetToFile([C2out],SetList,thisMomList,['Mass'],0,AddDict=InfoDict)
    print 'Printing took ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s  '
    # print 'Completed ' + ipTOqstr(thisMomList[0])

    

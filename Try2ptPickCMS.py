#!/usr/bin/env python

from array import array
from os import walk
import os
import numpy as np
import sys
from Params import *
from ReadCMCfuns import ReadSet,ReadList
from CMSTech import CreateCM2ptCfuns,CreatePoF2ptCfuns,PreptwoptCorr
from OutputData import PrintSetToFile,PrintLREvecMassToFile,PrintCfunToFile
from CreateCombs import CreategiDi
from SetLists import CreateMassSet
from MiscFuns import *
import time,datetime
from MultiWrap import *
from multiprocessing import Pool

def CreateTwoPt(thisMomList,thisSmearList):
    logfile = logdir+'LogTwoPt.log'
    errfile = logdir+'LogTwoPt.log'
    touch(logfile)
    touch(errfile)
    sys.stdout = open(logfile,'a',0)
    sys.stderr = open(errfile,'a',0)
    # print 'Running ' + ipTOqstr(thisMomList[0]) + ' ' +  str(int((thisMomList[0]*100)/float(len(qvecSet))))+'%' 

    if ListOrSet == 'ReadList':
        [data2pt,data3pt,filelist] = ReadList(thisSmearList,thisMomList,{},{},[],[],
                                              conflist,[],Interps=DefInterpList)
    elif ListOrSet == 'ReadSet':
        [data2pt,data3pt,filelist] = ReadSet(thisSmearList,thisMomList,{},{},[],[],
                                             dirread,[],Interps=DefInterpList)

    data2pt = np.array(PreptwoptCorr(np.array(data2pt)))
    ncon = np.size(filelist)
    print 'ncon = ' + str(ncon)
    # print 'nboot = ' + str(nboot)
## data2pt = [ ism , jsm , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
    C2out = DiagSmear(data2pt).tolist()

    start = time.time()
    inputparams = []
    makeContextFunctions(CreatePoF2ptCfuns)
    makeContextFunctions(CreateCM2ptCfuns)
    for icount,itodt in enumerate(DeftodtList):
        inputparams.append((data2pt,itodt,thisMomList))


    if DoMulticore:
        thisPool = Pool(min(len(inputparams),AnaProc))
        output = thisPool.map(CreatePoF2ptCfuns.mapper,inputparams)
        output += thisPool.map(CreateCM2ptCfuns.mapper,inputparams)
        thisPool.close()
        thisPool.join()
    else:
        output = []
        for iin in inputparams: output.append(CreatePoF2ptCfuns.mapper(iin))
        for iin in inputparams: output.append(CreateCM2ptCfuns.mapper(iin))
    
    
    thisTvarList = ['PoF'+str(PoFShifts)+iTvar for iTvar in TwoPtDefTvarList]
    thisTvarList += ['CM'+iTvar for iTvar in TwoPtDefTvarList]
    for iout,iTvar in zip(output,thisTvarList):
        [CMdata2pt,LEvec,REvec,Emass] = iout
## CMdata2pt [ istate , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
        C2out += CMdata2pt.tolist()
        PrintLREvecMassToFile(LEvec,REvec,Emass,thisMomList,iTvar,DoPoF=DoPoF)
    print 'CMTech Total Time Taken: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s  '
    start = time.time()
    print 'Printing to file \r',
    SetList = CreateMassSet(thisSmearList,StateSet,thisTvarList,flipord=True)
    PrintCfunToFile([C2out],SetList,thisMomList,['twopt'])
    PrintSetToFile([C2out],SetList,thisMomList,['Mass'],0)
    print 'Printing took ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s  '
    # print 'Completed ' + ipTOqstr(thisMomList[0])

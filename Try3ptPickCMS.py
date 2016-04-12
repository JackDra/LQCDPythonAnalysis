#!/usr/bin/env python
from array import array
from os import walk
import os
import numpy as np
import sys
from Params import *
from RFCalc import CalcRatioFactor
from ReadCMCfuns import ReadSet,ReadList
from CMSTech import CreateCMCfuns, CreateREvecCfuns,CreateREPoFCfuns,PreptwoptCorr
from OutputData import PrintSetToFile,PrintCfunToFile
from CreateCombs import CreategiDi, CreateDS
from Fitting import FitRFSet
from SetLists import CreateDataTsinkSet,CreateDataSet,CreateREvecSet,CreateMassSet
from MiscFuns import touch
import copy


def CreateRF(RunType,thisTSinkList,thisSmearList,thisPrefList,thisMomList,thisPGList={},thisPDList={},thisDSList=DefDSList,giDi=False,DontWriteZero=False):
    thisGammaList = []
    for iDS in thisDSList:
        for Proj,GL in thisPGList.iteritems():
            thisGammaList += [iDS+'P'+Proj[3]+iGL for iGL in GL]
            if not giDi: iglog = GL[0]
        for Proj,GL in thisPDList.iteritems():
            thisGammaList += [iDS+'P'+Proj[3]+iGL for iGL in GL]
    if giDi: iglog = 'giDi'
    logfilestart = logdir+'LogAll.log.start'
    logfilemid = logdir+'Log'+RunType+str(thisTSinkList[0])+iglog+'.log'
    logfileend = logdir+'LogAll.log.end'
    errfile = logdir+'Log'+RunType+str(thisTSinkList[0])+iglog+'.log'
    touch(logfilestart)
    touch(logfilemid)
    touch(logfileend)
    touch(errfile)
    sys.stdout = open(logfilestart,'a',0)
    print 'Running  '+ RunType + ' tsink='+str(thisTSinkList[0])+iglog
    sys.stdout = open(logfilemid,'a',0)
    sys.stderr = open(errfile,'a',0)

    # print thisPGList
    # print thisPDList
    # print thisDSList
    # print thisSmearList
    # print thisPrefList
    # print thisMomList
    # print thisTSinkList

    if len(thisTSinkList) > 1 and RunType != 'PoF':
        raise IOError('Not supporting multiple tsinks, handle externally')
    
    if ListOrSet == 'ReadList':
        [data2pt,data3pt,filelist] = ReadList(thisSmearList,thisMomList,thisPGList,thisPDList,
                                              thisDSList,thisTSinkList,conflist,thisPrefList)
    elif ListOrSet == 'ReadSet':
        [data2pt,data3pt,filelist] = ReadSet(thisSmearList,thisMomList,thisPGList,thisPDList,
                                             thisDSList,thisTSinkList,dirread,thisPrefList)
    print 'Read Complete'
    print 'ncon=',len(filelist)
    ## data2pt = [ ism , jsm , ip ,  it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
    ##data3pt = [ ism , jsm , igamma , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)


    if 'CM' == RunType:
## CMdata2pt [ istate , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
## CMdata3pt  [ istate , igamma , ip , it] = bootstrap1 class (.Avg, .Std, .values, .nboot)
        if giDi: data3pt,thisGammaList = CreategiDi(data3pt,thisGammaList,thisDSList)
        data2pt = np.array(PreptwoptCorr(data2pt))
        data3pt = np.array(data3pt)[0,:,:,:,:,:]
        data2ptset = DiagSmear(data2pt).tolist()
        data3ptset = DiagSmear(data3pt).tolist()
        start = time.time()
        for icount,(itodt,iTvar) in enumerate(zip(DeftodtList,DefTvarList)):
            thisstart = time.time()
            timeleft = GetTimeLeft(icount,len(DeftodtList),(time.time()-start))
            print 'CMTech ' , iTvar , ' Time Left: ' , str(datetime.timedelta(seconds=timeleft)) , ' h:m:s   \r',
            [CMdata2pt,CMdata3pt] = CreateCMCfuns(data3pt,data2pt,itodt,thisMomList)
            data2ptset += CMdata2pt.tolist()
            data3ptset += CMdata3pt.tolist()
            print 'CMTech ' , iTvar , ' Time Taken: ' , str(datetime.timedelta(seconds=time.time()-thisstart)) , ' h:m:s  '
        data3ptset = np.array(data3ptset)
        data2ptset = np.array(data2ptset)
        SetList = CreateMassSet(thisSmearList,StateSet,DefTvarList,flipord=True)
        SetList = ['tsink'+str(thisTSinkList[0])+iS for iS in SetList]
        print 'CMTech Total Time Taken: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s  '
    elif 'REvec' == RunType:
## CMdata2pt [ istate , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
## CMdata3pt  [ istate , igamma , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
        if giDi: data3pt,thisGammaList = CreategiDi(data3pt,thisGammaList,thisDSList)
        data2pt = np.array(PreptwoptCorr(data2pt))
        data3pt = np.array(data3pt)[0,:,:,:,:,:]
        print 'Creating REvec CM Tech ' , REvecTvarList[0]
        [data2ptset,data3ptset] = CreateREvecCfuns(data3pt,data2pt,DefREvecVarList,thisMomList)
        SetList,dump = CreateREvecSet(thisTSinkList,StateSet,REvecTvarList)
    elif 'PoF' == RunType:
## CMdata2pt [ istate , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
## CMdata3pt  [ istate , igamma , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
        if giDi:
            for idata,idata3pt in enumerate(data3pt):
                dumpdata3pt,dumpGammaList = CreategiDi(idata3pt,thisGammaList,thisDSList)
                data3pt[idata] = dumpdata3pt
            thisGammaList = dumpGammaList
        if len(data3pt) < 2: raise IOError("PoF needs atleast two tsinks")
        data2pt = np.array(PreptwoptCorr(data2pt))
        print 'Creating PoF CM Tech ' , PoFTvarList[0]
        [data2ptset,data3ptset] = CreateREPoFCfuns(np.array(data3pt[0]),np.array(data3pt[1]),data2pt,DefPoFVarList,thisMomList)
        SetList,dump = CreateREvecSet(thisTSinkList,StateSet,PoFTvarList)
        # for it in range(15,30):
        #     print it, data3ptset[0][0][0][it].Avg, data2ptset[0][0][it].Avg
            
    elif 'TSink' == RunType:
        if giDi: data3pt,thisGammaList = CreategiDi(data3pt,thisGammaList,thisDSList)
        data2pt = np.array(PreptwoptCorr(data2pt))
        data3pt = np.array(data3pt)[0,:,:,:,:,:]
        data2ptset,data3ptset = DiagSmear(data2pt),DiagSmear(data3pt)
        SetList = ['tsink'+str(thisTSinkList[0])+'sm'+ism for ism in thisSmearList]
    print 'Analysis Complete'

    [RFr,SqrtFac] = CalcRatioFactor(data2ptset,data3ptset,str(thisTSinkList[0]))
    print 'RF Construction Complete'
    
    if DontWriteZero:
        thisMomList = thisMomList[1:]
        data3ptset = data3ptset[:,:,1:,:]
        RFr = RFr[:,:,1:,:]
    ## data2ptset [ iset , ip , it ]
    ## data3ptset [ iset , igamma , ip , it ] bs1
    PrintCfunToFile(np.rollaxis(data3ptset,1),SetList,thisMomList,thisGammaList)

    ## RFr = [  iset , igamma , ip , it ] bs1
    PrintSetToFile(np.rollaxis(RFr,1),SetList,thisMomList,thisGammaList,thisTSinkList[0])

    print 'Finished '+ RunType + ' tsink='+str(thisTSinkList[0]) + ' '+iglog
    sys.stdout = open(logfileend,'a',0)
    print 'Finished '+ RunType + ' tsink='+str(thisTSinkList[0]) + ' '+iglog

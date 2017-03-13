#!/usr/bin/env python
from array import array
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
from SetLists import CreateREvecSet,CreateMassSet
from MiscFuns import touch
import copy


def CreateRF(RunType,thisTSinkList,thisiSmearList,thisjSmearList,thisPrefList,thisMomList,thisPGList={},thisPDList={},thisDSList=DefDSList,giDi=False,DontWriteZero=False):
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
    
    if 'ReadList' in ListOrSet:
        if 'PoF' in RunType:
            [data2pt,data3pt,filelist] = ReadList(thisiSmearList,thisjSmearList,thisMomList,thisPGList,thisPDList,
                                                  thisDSList,thisTSinkList,conflist,thisPrefList,thistsourceList=PoFtsourceList)
        else:
            [data2pt,data3pt,filelist] = ReadList(thisiSmearList,thisjSmearList,thisMomList,thisPGList,thisPDList,
                                                  thisDSList,thisTSinkList,conflist,thisPrefList)
    elif 'ReadSet' in ListOrSet:
        if 'PoF' in RunType:
            [data2pt,data3pt,filelist] = ReadSet(thisiSmearList,thisjSmearList,thisMomList,thisPGList,thisPDList,
                                                 thisDSList,thisTSinkList,dirread,thisPrefList,thistsourceList=PoFtsourceList)
        else:
            [data2pt,data3pt,filelist] = ReadSet(thisiSmearList,thisjSmearList,thisMomList,thisPGList,thisPDList,
                                                 thisDSList,thisTSinkList,dirread,thisPrefList)
    print 'Read Complete'
    InfoDict = GetInfoFromFilelist(filelist)
    ## data2pt = [ tsource, ism , jsm , ip ,  it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
    ##data3pt = [ tsink, tsource, ism , jsm , igamma , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
    # if Debug:
    # for it in xrange(0,15):
    #     # RF = data3pt[0][0][0][0][0][0][it]/data2pt[0][0][-1][0][15]
    #     # RF.Stats()
    #     print
    #     # for iBoot,(b3,b2,br) in enumerate(zip(data3pt[0][0][0][0][0][0][it].values,data2pt[0][0][0][0][15].values,RF.values)):
    #     #     print it, iBoot,b3,b2,br 
    #     for iBoot,b3 in enumerate(data3pt[0][0][0][0][0][0][it].values):
    #         print it, iBoot,b3 
        # print it, data3pt[0][0][0][0][0][0][it].Avg, data2pt[0][0][0][0][it].Avg, data2pt[0][0][0][0][it].Std, RF.Avg, RF.Std
        
    if 'CM' == RunType:
## CMdata2pt [ istate , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
## CMdata3pt  [ istate , igamma , ip , it] = bootstrap1 class (.Avg, .Std, .values, .nboot)
        data2pt = np.array(PreptwoptCorr(data2pt[0]))
        data3pt = np.array(data3pt)[0,0,:,:,:,:,:]
        if giDi: data3pt,thisGammaList = CreategiDi(data3pt,thisGammaList,thisDSList)
        data2ptset = FlattenSmear(data2pt).tolist()
        data3ptset = FlattenSmear(data3pt).tolist()
        start = time.time()
        for icount,(itodt,iTvar) in enumerate(zip(AnatodtList,AnaTvarList)):
            thisstart = time.time()
            timeleft = GetTimeLeft(icount,len(AnatodtList),(time.time()-start))
            print 'CMTech ' , iTvar , ' Time Left: ' , str(datetime.timedelta(seconds=timeleft)) , ' h:m:s   \r',
            [CMdata2pt,CMdata3pt] = CreateCMCfuns(data3pt,data2pt,itodt,thisMomList)
            data2ptset += CMdata2pt.tolist()
            data3ptset += CMdata3pt.tolist()
            print 'CMTech ' , iTvar , ' Time Taken: ' , str(datetime.timedelta(seconds=time.time()-thisstart)) , ' h:m:s  '
        data3ptset = np.array(data3ptset)
        data2ptset = np.array(data2ptset)
        SetList = CreateMassSet(thisiSmearList,thisjSmearList,CMStateSet,AnaTvarList,flipord=True)
        SetList = ['tsink'+str(thisTSinkList[0])+iS for iS in SetList]
        print 'CMTech Total Time Taken: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s  '
    elif 'REvec' == RunType:
## CMdata2pt [ istate , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
## CMdata3pt  [ istate , ip , igamma , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
        data2pt = np.array(PreptwoptCorr(data2pt[0]))
        data3pt = np.array(data3pt)[0,0,:,:,:,:,:]
        if giDi: data3pt,thisGammaList = CreategiDi(data3pt,thisGammaList,thisDSList)
        print 'Creating REvec CM Tech ' , REvecTvarList[0]
        [data2ptset,data3ptset] = CreateREvecCfuns(data3pt,data2pt,DefREvecVarList,thisMomList)
        SetList,dump = CreateREvecSet(thisTSinkList,[PickedState],REvecTvarList)
        MassSetList = CreateMassSet([],[],[PickedState],REvecTvarList,flipord=True)
        PrintCfunToFile([data2ptset],MassSetList,thisMomList,['twopt'],AddDict=InfoDict)
        PrintSetToFile([data2ptset],MassSetList,thisMomList,['Mass'],0,AddDict=InfoDict)
    elif 'PoF' == RunType:
## CMdata2pt [ istate , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
## CMdata3pt  [ istate , ip , igamma , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
        # if giDi:
        #     for idata,idata3pt in enumerate(data3pt):
        #         dumpdata3pt,dumpGammaList = CreategiDi(idata3pt,thisGammaList,thisDSList)
        #         data3pt[idata] = dumpdata3pt
        #     thisGammaList = dumpGammaList
        if len(data3pt) < 2 and TimeInv: raise IOError("PoF needs atleast two tsinks with time invariance")
        data2pt = np.array(PreptwoptCorr(np.array(data2pt)))
        data3ptset,data2ptset = [],[]
        start = time.time()
        for itodt,iTvar in zip(DefPoFVarList,PoFTvarList):
            print 'Creating PoF CM Tech ' , iTvar
            thisstart = time.time()
            [CMdata2pt,CMdata3pt] = CreateREPoFCfuns(np.array(data3pt),data2pt,itodt,thisMomList,todtvalsLeft = DefPoFTvarRef)
            # if Debug:
            # for it in xrange(0,15):
            #     print
            #     print 'it',it, 'c3ptAvg',data3ptset[0][0][0][it].Avg, 'c2ptAvg',data2ptset[0][0][it].Avg
            #     for iboot in xrange(0,10):
            #         print 'it',it, 'c3pt',data3ptset[0][0][0][it].values[iboot], 'c2pt',data2ptset[0][0][it].values[iboot]
            data2ptset += CMdata2pt.tolist()
            data3ptset += CMdata3pt.tolist()
            print 'CMTech ' , iTvar , ' Time Taken: ' , str(datetime.timedelta(seconds=time.time()-thisstart)) , ' h:m:s  '
        data3ptset = np.array(data3ptset)
        data2ptset = np.array(data2ptset)
        SetList,dump = CreateREvecSet(thisTSinkList,StateSet,PoFTvarList,fliptodt=True)
        print 'CMTech Total Time Taken: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s  '
    elif 'TSink' == RunType:
        data2pt = np.array(PreptwoptCorr(data2pt[0]))
        data3pt = np.array(data3pt)[0,0,:,:,:,:,:]
        if giDi: data3pt,thisGammaList = CreategiDi(data3pt,thisGammaList,thisDSList)
        data2ptset,data3ptset = FlattenSmear(data2pt),FlattenSmear(data3pt)
        SetList = []
        for ism in thisiSmearList:
            SetList += ['tsink'+str(thisTSinkList[0])+'ism'+ism+'jsm'+jsm for jsm in thisjSmearList]
    print 'Analysis Complete'

    [RFr,SqrtFac] = CalcRatioFactor(data2ptset,data3ptset,str(thisTSinkList[0]),thisMomList)
    print 'RF Construction Complete'
    
    
    if DontWriteZero:
        thisMomList = thisMomList[1:]
        data3ptset = data3ptset[:,:,1:,:]
        RFr = RFr[:,:,1:,:]
    ## data2ptset [ iset , ip , it ]
    ## data3ptset [ iset , igamma , ip , it ] bs1
    PrintCfunToFile(np.rollaxis(data3ptset,1),SetList,thisMomList,thisGammaList,AddDict=InfoDict)

    ## RFr = [  iset , igamma , ip , it ] bs1
    PrintSetToFile(np.rollaxis(RFr,1),SetList,thisMomList,thisGammaList,thisTSinkList[0],AddDict=InfoDict)

    print 'Finished '+ RunType + ' tsink='+str(thisTSinkList[0]) + ' '+iglog + '                       '
    sys.stdout = open(logfileend,'a',0)
    print 'Finished '+ RunType + ' tsink='+str(thisTSinkList[0]) + ' '+iglog

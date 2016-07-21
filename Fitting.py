#!/usr/bin/env python

import numpy as np
from MultiWrap import *
from FitFunctions import *
from LLSBoot import * 
from FitParams import *
from BootTest import BootStrap1
# from numpy.polynomial.polynomial import polyfit
from SetLists import GetTsinkSmLists
from Params import *
from collections import OrderedDict
import re
import time
import datetime

def FitSMList(data,fitr,nsm):
    tdata = np.arange(fitr[0]-tsource,fitr[1]-tsource+1)
    tsmdata = [[],[]]
    dataout = []
    for it in tdata:
        for ism in range(nsm):
            tsmdata[0].append(it)
            tsmdata[1].append(ism)
            dataout.append(data[ism][it+tsource-1])
    return dataout,tsmdata

def PickBoot2pt(data,thisnsm,ism):
    return [data[ism],data[thisnsm+ism],data[-2],data[-1]]

def Trange(ydata,tdata):
    return np.array(ydata)[tdata]


# def TMassrange(ydata,tdata):
#     return np.abs(np.log(ydata[np.array(tdata)]/ydata[np.array(tdata)+int(massdt)]))/massdt

def FitRFWrap(ydata,tdata):
    boot,Avg,Err = FitBoots(Trange(ydata,tdata),tdata,ConstantFitFun)
    return [boot[0],Avg[0],Err[0]]

# RFin = [ ism/istate , it ] bs1
#fitAvgOut = [ ism/istate ]
#fitBoot = [ ism/istate ] bs1

def FitRFSet(RFin,thisTSinkList,icut):
    fitBoot,fitAvg,fitChi = [],[],[]
    # print icut , len(RFin) , len(thisTSinkList)
    for its,itsink in enumerate(thisTSinkList):
        iRF = RFin[its]
        tdata = np.arange(icut, int(itsink)-tsource-icut)
        if tdata.size < 1: 
            fitBoot.append(iRF[(int(itsink)-tsource)/2])
            fitAvg.append(iRF[(int(itsink)-tsource)/2].Avg)
            fitChi.append(float('NaN'))
    # if tdata.size < 1: raise IndexError("Icut reduced a set to 0 time slices (do separatly please)")
        else:
            [fitBoothold,fitAvghold,fitChihold] = FitRFWrap(iRF,tdata)
            fitBoot.append(fitBoothold)
            fitAvg.append(fitAvghold)
            fitChi.append(fitChihold)
    return [fitBoot,fitAvg,fitChi]

def FitMassSet(Massin,tmin,tmax):
    Massin = np.array(Massin)
    tdata = np.arange(tsource-1+tmin,tsource-1+tmax)
    MassRange = NDimOpp(Massin,1,Trange,tdata)
    [fitBoot,fitAvg,fitChi] = NDimOpp(np.rollaxis(MassRange,0,len(MassRange.shape)),1,FitRFWrap,tdata)
    return [fitBoot,fitAvg,fitChi]

## C2pt = [ ip , istate/ism  , it ] bs1
## C3pt = [ igamma , ip , istate/ism , it ]
#___2pt = [ ip , params ]
#___3pt = [ igamma , ip , icut , params ]


def MomTSSetFit(TSF2ptarray,C3pt,this3ptCutList,thisSetList,thisGammaMomList,this2ptFitRvec):
    thisDoMulticore = False
    def smfitwrap(thisBoot2ptmom,thisBoot2ptZ,C3mom,this3ptCutList,thisTSinkList,thisSmList):
        def GetTsinkInSm(C3,funsm,funSetList):
            C3out = []
            for iS,iSet in enumerate(funSetList):
                if funsm in iSet:
                    C3out.append(C3[iS])
            return C3out
        Boot3pt,Avg3pt,Chi3pt = [],[],[]
        for ism,thissm in enumerate(thisSmList):
            Params2pt,Params2ptZero = PickBoot2pt(thisBoot2ptmom,thisnsm,ism),PickBoot2pt(thisBoot2ptZ,thisnsm,ism)
            isC3 = GetTsinkInSm(C3mom,thissm,thisSetList)
            thisod = TwoStateFitMom3pt(Params2ptZero+Params2pt,isC3,
                                       this3ptCutList,thisTSinkList)
            Boot3pt.append(thisod[0])
            Avg3pt.append(thisod[1])
            Chi3pt.append(thisod[2])
        return Boot3pt,Avg3pt,Chi3pt

    start = time.time()
    thisTSinkList,thisSmList = GetTsinkSmLists(thisSetList)
    thisTSinkList = [int(its.replace('tsink','')) for its in thisTSinkList]
    thisnsm = len(thisSmList)
    this2ptFitR,perdone = this2ptFitRvec
    Boot3pt,Avg3pt,Chi3pt = [],[],[]
    Boot2pt,Avg2pt,Chi2pt = TSF2ptarray
    
    inputparams = []
    Boot2ptZ = Boot2pt[0]
    thisigamma = -1
    for igamma,(thisgamma,thismomlist) in enumerate(thisGammaMomList.iteritems()):
        if thisgamma == 'twopt': continue
        thisigamma += 1
        # print thisgamma , ' thismomlistlen=',len(thismomlist) , ' C3momlen=' , len(C3pt[thisigamma])
        for imom,thismom in enumerate(thismomlist):
            imom2pt = thisGammaMomList['twopt'].index(thismom)
            C3mom = C3pt[thisigamma][imom]
            Boot2ptMom = Boot2pt[imom2pt]
            inputparams.append((Boot2ptMom,Boot2ptZ,C3mom,this3ptCutList,thisTSinkList,thisSmList))
    if thisDoMulticore:
        output3pt = thisPool.map(smfitwrap.mapper,inputparams)
    else:
        output3pt = []
        for icn,iin in enumerate(inputparams): output3pt.append(smfitwrap(*iin))
    totcount = 0
    thisigamma = -1
    for igamma,(thisgamma,thismomlist) in enumerate(thisGammaMomList.iteritems()):
        if thisgamma == 'twopt': continue
        thisigamma += 1
        Boot3pt.append([])
        Avg3pt.append([])
        Chi3pt.append([])
        # print thisgamma , ' thismomlistlen=',len(thismomlist) , ' C3momlen=' , len(C3pt[thisigamma])
        for imom,thismom in enumerate(thismomlist):
            Boot3pt[thisigamma].append(output3pt[totcount][0])
            Avg3pt[thisigamma].append(output3pt[totcount][1])
            Chi3pt[thisigamma].append(output3pt[totcount][2])
            totcount += 1
        # print 'fit range ' , this2ptFitR , ' ThreePt Fit At: ' ,int((igamma*100) / float(len(thisGammaMomList.keys()))), '%                             \r',
    if thisDoMulticore:
        thisPool.close()
        thisPool.join()
    mprint( 'fit range ' , this2ptFitR , ' ' , perdone, '% took:  ',str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s                    ')
    return [Boot3pt,Avg3pt,Chi3pt]


def MomTSSetFit2pt(C2pt,thisSetList,thisGammaMomList,this2ptFitRvec):
    this2ptFitR,perdone = this2ptFitRvec
    Boot2pt,Avg2pt,Chi2pt = [],[],[]
    start = time.time()
    thisTSinkList,thisSmList = GetTsinkSmLists(thisSetList)
    thisnsm = len(thisSmList)
    thisTSinkList = [int(its.replace('tsink','')) for its in thisTSinkList]
    inputparams = [FitSMList(C2pt[imom],this2ptFitR,thisnsm) + (thisnsm,) for imom in range(len(thisGammaMomList['twopt']))]
    if DoMulticore:
        makeContextFunctions(TwoStateFit2pt)
        thisPool = Pool(min(len(thisGammaMomList['twopt']),AnaProc))
        output = thisPool.map(TwoStateFit2pt.mapper,inputparams)
    else:
        output = []
        for iin in inputparams: output.append(TwoStateFit2pt(*iin))
    for imom,thismom in enumerate(thisGammaMomList['twopt']):
        Boot2pt.append(output[imom][0])
        Avg2pt.append(output[imom][1])
        Chi2pt.append(output[imom][2])
    if DoMulticore:
        thisPool.close()
        thisPool.join()
    print 'fit range ' , this2ptFitR , ' ' , perdone, '% took:  ',str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s                    '
    return [Boot2pt,Avg2pt,Chi2pt]


# C2pt = [ it ] bs1
# C3pt = [ itsink , it ] bs1
#fitBoot2pt = [ params ]bs1
#fitAvg2pt =  [ params ]
#fitChi2pt =  #
#fitBoot3pt = [ icut , params ] bs1
#fitAvg3pt =  [ icut , params ]
#fitChi3pt =  [ icut ]

def TwoStateFit2pt(fitdata2pt,tdata2pt,parl):
    fitBoot2pt,fitAvg2pt,fitAvg2ptChi = FitBoots(fitdata2pt,tdata2pt,C2TwoStateFitFunCM,parlen=parl)
    return [fitBoot2pt,fitAvg2pt,fitAvg2ptChi[0]]

def TwoStateFitMom3pt(fitBoot2pt,C3pt,this3ptCutList,thisTSinkList):
    tsvar = True
    fitBoot3pt = []
    fitAvg3pt = []
    fitAvg3ptChi = []
    BootPars = np.rollaxis(Pullflag(fitBoot2pt,'values'),1)
    AvgPars = Pullflag(fitBoot2pt,'Avg')
    for icut,cutval in enumerate(this3ptCutList):
        tdata3pt = []
        fitdata3pt = []
        for its,thistsink in enumerate(thisTSinkList):
            tsinkdata = C3pt[its]
            tsinkhalf = (thistsink-tsource)/2
            if cutval > tsinkhalf:
                tdata3pt = [[tsinkhalf, thistsink-tsource]]
                fitdata3pt = [C3pt[tsinkhalf-1]]
            else:
                for it in np.arange(cutval,thistsink-tsource-cutval+1):
                    tdata3pt.append([it,thistsink-tsource])
                    # fitdata3pt.append(tsinkdata[it+tsource-1]/C2pt[thistsink-1])
                    fitdata3pt.append(tsinkdata[it+tsource-1])
        tdata3pt = np.rollaxis(np.array(tdata3pt),1)
        tdata3ptout = []
        for iboot,ipar in enumerate([AvgPars]+BootPars.tolist()):            
            tdata3ptout.append([])
            tdata3ptout[iboot].append(tdata3pt[0])
            tdata3ptout[iboot].append(tdata3pt[1])
            for ipp in ipar:
                tdata3ptout[iboot].append([ipp]*len(tdata3pt[0]))
        # if igamma == 0 and icut == 0:
        #     print tdata3pt
        #     print Pullflag(fitdata3pt,'Avg')
        #     print Pullflag(fitBoot2pt,'Avg')
        #     print C2pt[thistsink-1].Avg
        boothold,Avghold,Chihold = FitBoots(fitdata3pt,tdata3ptout,C3MomTwoStateFitFun,tBooted=True)
        fitAvg3pt.append(Avghold)
        fitAvg3ptChi.append(Chihold[0])
        fitBoot3pt.append(boothold)
    return [fitBoot3pt,fitAvg3pt,fitAvg3ptChi]


##ONE STATE FIT##


## C2pt = [ ip , istate/ism , it ] bs1
## C3pt = [ igamma , ip , iset , it ]
#___2pt = [ ip , istate/ism  , params ]
#___3pt = [ igamma , ip , iset , i3cut , params ]

def OneStateSet2pt(C2pt,thisSetList,thisGammaMomList,this2ptFitRvec):
    def sm2ptwrap(C2ptmom,thisSmList,this2ptFitR):
        Bootthis2pt,Avgthis2pt,Chithis2pt = [],[],[]
        for ism,thissm in enumerate(thisSmList):
            isC2 = C2ptmom[ism]
            thisod2 = OneStateFit2pt(isC2,this2ptFitR)
            Bootthis2pt.append(thisod2[0])
            Avgthis2pt.append(thisod2[1])
            Chithis2pt.append(thisod2[2])
        return Bootthis2pt,Avgthis2pt,Chithis2pt

    this2ptFitR,perdone = this2ptFitRvec
    thisTSinkList,thisSmList = GetTsinkSmLists(thisSetList,NoREvec=True)
    Boot2pt,Avg2pt,Chi2pt = [],[],[]
    start = time.time()
    inputparams = [(C2pt[imom],thisSmList,this2ptFitR) for imom in range(len(thisGammaMomList['twopt']))]
    makeContextFunctions(sm2ptwrap)
    if DoMulticore:
        thisPool = Pool(min(len(thisGammaMomList['twopt']),AnaProc))
        output = thisPool.map(sm2ptwrap.mapper,inputparams)
    else:
        output = []
        for iin in inputparams: output.append(sm2ptwrap(*iin))
    for imom,thismom in enumerate(thisGammaMomList['twopt']):
        Boot2pt.append(output[imom][0])
        Avg2pt.append(output[imom][1])
        Chi2pt.append(output[imom][2])
    if DoMulticore:
        thisPool.close()
        thisPool.join()
    print 'fit range ' , this2ptFitR , ' ' , perdone, '% took:  ',str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s                    '
    return [Boot2pt,Avg2pt,Chi2pt]


def OneStateSetFit(OSF2ptarray,C3pt,this3ptCutList,thisSetList,thisGammaMomList,this2ptFitRvec):
    thisDoMulticore = False
    def sm3ptwrap(thisBoot2ptmom,thisBoot2ptZ,C3mom,this3ptCutList,thisSetList,thisSML):
        def SplitIset(thisiset,thisSML):
            # thisiset = thisiset.replace('REvec','CM')
            sm = re.sub('tsink..','',thisiset)
            itsink = thisiset.replace(sm,'').replace('tsink','')
            print thisSML
            ism = thisSML.index(sm)
            return ism,int(itsink)

        Boot3pt,Avg3pt,Chi3pt = [],[],[]
        # print len(C3mom) , thisSetList
        for iset,thisset in enumerate(thisSetList):
            isC3 = C3mom[iset]
            thisism,thistsink = SplitIset(thisset,thisSML)
            thisod3 = OneStateFit3pt(thisBoot2ptZ[thisism]+thisBoot2ptmom[thisism],
                                     isC3,this3ptCutList,thistsink)
            Boot3pt.append(thisod3[0])
            Avg3pt.append(thisod3[1])
            Chi3pt.append(thisod3[2])
        return Boot3pt,Avg3pt,Chi3pt

    start = time.time()
    this2ptFitR,perdone = this2ptFitRvec
    thisTSinkList,thisSmList = GetTsinkSmLists(thisSetList,NoREvec=True)
    Boot2pt,Avg2pt,Chi2pt = OSF2ptarray
    Boot3pt,Avg3pt,Chi3pt = [],[],[]
    Boot2ptZ = Boot2pt[0]
    inputparams = []
    thisigamma = -1
    for igamma,(thisgamma,thismomlist) in enumerate(thisGammaMomList.iteritems()):
        if thisgamma == 'twopt': continue
        thisigamma += 1
        for imom,thismom in enumerate(thismomlist):
            imom2pt = thisGammaMomList['twopt'].index(thismom)
            C3mom =C3pt[thisigamma][imom]
            Boot2ptMom = Boot2pt[imom2pt]
            inputparams.append((Boot2ptMom,Boot2ptZ,C3mom,this3ptCutList,thisSetList,thisSmList))
    if thisDoMulticore:
        output3pt = thisPool.map(sm3ptwrap.mapper,inputparams)
    else:
        output3pt = []
        for icin,iin in enumerate(inputparams): output3pt.append(sm3ptwrap(*iin))
    totcount = 0
    thisigamma = -1
    for igamma,(thisgamma,thismomlist) in enumerate(thisGammaMomList.iteritems()):
        if thisgamma == 'twopt': continue
        thisigamma += 1
        Boot3pt.append([])
        Avg3pt.append([])
        Chi3pt.append([])
        for imom,thismom in enumerate(thismomlist):
            Boot3pt[thisigamma].append(output3pt[totcount][0])
            Avg3pt[thisigamma].append(output3pt[totcount][1])
            Chi3pt[thisigamma].append(output3pt[totcount][2])
            totcount += 1
        # print 'fit range ' , this2ptFitR , ' ThreePt Fit At: ' ,int((igamma*100) / float(len(thisGammaMomList.keys()))), '%                             \r',
    if thisDoMulticore:
        thisPool.close()
        thisPool.join()
    mprint( 'fit range ' , this2ptFitR , ' ' , perdone, '% took:  ',str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s                    ')
    return [Boot3pt,Avg3pt,Chi3pt]


# C2pt = [ it ] bs1
# C3pt = [ it ] bs1
#fitBoot2pt = [ params ]bs1
#fitAvg2pt =  [ params ]
#fitChi2pt =   #
#fitBoot3pt = [ icut  , params ] bs1
#fitAvg3pt =  [ icut  , params ]
#fitChi3pt =  [ icut  ]

## NB, this3ptCutList is how much to fit over from the centre out

def OneStateFit2pt(data2pt,fitr):
    tdata2pt = np.arange(fitr[0]-tsource,fitr[1]-tsource+1)
    fitdata2pt = np.array(data2pt)[tdata2pt+tsource-1]
    fitBoot2pt,fitAvg2pt,fitAvg2ptChi = FitBoots(fitdata2pt,tdata2pt,C2OneStateFitFun)
    return [fitBoot2pt,fitAvg2pt,fitAvg2ptChi[0]]


def OneStateFit3pt(fitBoot2pt,C3pt,this3ptCutList,thistsink):
    fitBoot3pt = []
    fitAvg3pt = []
    fitAvg3ptChi = []
    AvgPars  = Pullflag(fitBoot2pt,'Avg') 
    BootPars = np.rollaxis(Pullflag(fitBoot2pt,'values'),1)
    for icut,cutval in enumerate(this3ptCutList):
        tdata3pt = []
        fitdata3pt = []
        tsinkhalf = (thistsink-tsource)/2
        if cutval > tsinkhalf:
            tdata3pt = [[tsinkhalf, thistsink-tsource]]
            fitdata3pt = [C3pt[tsinkhalf-1]]
        else:
            for it in np.arange(cutval,thistsink-tsource-cutval+1):
                tdata3pt.append([it,thistsink-tsource])
                fitdata3pt.append(C3pt[it+tsource-1])
                # tdata3pt.append(it)
                # fitdata3pt.append(tsinkdata[it+tsource-1]/C2pt[thistsink-1])
        tdata3pt = np.rollaxis(np.array(tdata3pt),1)
        tdata3ptout = []
        for iboot,ipar in enumerate([AvgPars]+BootPars.tolist()):            
            tdata3ptout.append([])
            tdata3ptout[iboot].append(tdata3pt[0])
            tdata3ptout[iboot].append(tdata3pt[1])
            for ipp in ipar:
                tdata3ptout[iboot].append([ipp]*len(tdata3pt[0]))

        # if igamma == 0 and icut == 0:
        #     print tdata3pt
        #     print Pullflag(fitdata3pt,'Avg')
        #     print Pullflag(fitBoot2pt,'Avg')
        #     print C2pt[thistsink-1].Avg
        # boothold,Avghold,Chihold = FitBoots(fitdata3pt,tdata3pt,ConstantFitFun)
        boothold,Avghold,Chihold = FitBoots(fitdata3pt,tdata3ptout,C3OneStateFitFun,tBooted=True)
        fitAvg3pt.append(Avghold)
        fitAvg3ptChi.append(Chihold[0])
        fitBoot3pt.append(boothold)
    return [fitBoot3pt,fitAvg3pt,fitAvg3ptChi]



















# ## C2pt = [ istate/ism , ip , it ] bs1
# ## C3pt = [ itsink , ism/istate , igamma , ip , it ]
# #fit___ = [ istate/ism , "below" ]

# def TwoStateSetFit(C2pt,C3pt,this2ptFitR,this3ptCutList,thisTSinkList):
#     C2ptin = np.rollaxis(np.array(C2pt)[:,ZeroMomIndex,:],1)
#     C3ptin = np.rollaxis(np.array(C3pt)[:,:,:,ZeroMomIndex,:],1)
#     Boot3pt = []
#     Avg3pt = []
#     Chi3pt = []
#     thisnsm = C2ptin.shape[1]
#     C2ptfit,thistsmlist = FitSMList(C2ptin,this2ptFitR,thisnsm)
#     [Boot2pt,Avg2pt,Chi2pt] = TwoStateFit2pt(C2ptfit,thistsmlist,thisnsm)
#     for ism,(isC2pt,isC3pt) in enumerate(zip(np.swapaxes(C2ptin,0,1),C3ptin)):
#         Params2pt = PickBoot2pt(Boot2pt,thisnsm,ism)
#         thisod = TwoStateFit3pt(Params2pt,isC2pt,isC3pt,this3ptCutList,thisTSinkList)
#         Boot3pt.append(thisod[0])
#         Avg3pt.append(thisod[1])
#         Chi3pt.append(thisod[2])

#     return [np.array(Boot2pt),np.swapaxes(np.array(Boot3pt),0,1),np.array(Avg2pt),np.swapaxes(np.array(Avg3pt),0,1),np.array(Chi2pt),np.swapaxes(np.array(Chi3pt),0,1)]

# # C2pt = [ it ] bs1
# # C3pt = [ itsink , igamma , it ] bs1
# #fitAvg2pt = [ params ]
# #fitChi2pt = #
# #fitBoot2pt = [ params ]bs1
# #fitAvg3pt = [ igamma , params ]
# #fitChi3pt = [ igamma ]
# #fitBoot3pt = [ igamma , params ] bs1

# ## NB, this3ptCutList is how much to fit over from the centre out



# def TwoStateFit3pt(fitBoot2pt,C2pt,C3pt,this3ptCutList,thisTSinkList):
#     tsvar = True
#     if len(set(thisTSinkList)) < 2 : tsvar = False
#     fitBoot3pt = []
#     fitAvg3pt = []
#     fitAvg3ptChi = []
#     for icut,cutval in enumerate(this3ptCutList):
#         fitBoot3pt.append([])
#         fitAvg3pt.append([])
#         fitAvg3ptChi.append([])
#         for igamma,gammadata in enumerate(np.rollaxis(C3pt,1)):
#             tdata3pt = []
#             fitdata3pt = []
#             for thistsink,tsinkdata in zip(thisTSinkList,gammadata):
#                 for it in np.arange(cutval,thistsink-tsource-cutval+1):
#                     tdata3pt.append([it,thistsink-tsource])
#                     # fitdata3pt.append(tsinkdata[it+tsource-1]/C2pt[thistsink-1])
#                     fitdata3pt.append(tsinkdata[it+tsource-1])
#             tdata3pt = np.rollaxis(np.array(tdata3pt),1)
#             boothold,Avghold,Chihold = FitVarFunBoots(fitdata3pt,tdata3pt,CreateC3TSFitFun,fitBoot2pt,tsvar)
#             fitAvg3pt[icut].append(Avghold)
#             fitAvg3ptChi[icut].append(Chihold[0])
#             fitBoot3pt[icut].append(boothold)

#     return [fitBoot3pt,fitAvg3pt,fitAvg3ptChi]



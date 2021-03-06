#!/usr/bin/env python

import numpy as np
from Params import *
from MiscFuns import *
from FitFunctions import *
from LLSBoot import *
import time,datetime

def SumRF(thisRF,thisCut):
    SummedData = BootStrap1(nboot,0)
    for dataRF in thisRF[thisCut:len(thisRF)-thisCut]:
        SummedData = SummedData + dataRF
    SummedData.Stats()
    return SummedData


#Rin = [ itsink , it ] bs1
#SumOut = [ icut , itsink ] bs1

def CreateSumMeth(Rin,Cuts):
    SumOut = []
    for ic,icut in enumerate(Cuts):
        SumOut.append([])
        for tsinkR in Rin:
            SumOut[ic].append(SumRF(tsinkR,icut))
    return SumOut
# def CreateSumMeth(Rin,Cuts):
#     SumOut = []
#     for ic,icut in enumerate(Cuts):
#         SumOut.append([])
#         for tsinkR in Rin:
#             SumOut[ic].append(NDimOpp(tsinkR,1,SumRF,icut))
#     return SumOut

#Rin = [ imom , itsink , it ] 
#SumRin = [ imom , icut , itsink ] bs1
#data = [ imom , icut , fitr , par ] bs1

def FitSum(Cuts,thisTSinkList,thisGamma,thisMomList,Rin):
    bootdata,Avgdata,Chidata,FitList = [],[],[],[]
    TSL = np.array(thisTSinkList)
    SumRin = []
    start = time.time()
    for imom,thismom in enumerate(thisMomList):
        SumRin.append(CreateSumMeth(Rin[imom],Cuts))
        Chidata.append([])
        Avgdata.append([])
        bootdata.append([])
        FitList.append([])
        for icut,(thiscut,cRin) in enumerate(zip(Cuts,SumRin[-1])):
            Chidata[imom].append([])
            Avgdata[imom].append([])
            bootdata[imom].append([])
            FitList[imom].append([])
            for tmin,tmax in itertools.product(range(len(cRin)),range(len(cRin))):
                if tmin+1 < tmax:
                    datafit = np.array(cRin)[range(tmin,tmax+1)]
                    bootdatahold,Avgdatahold,Chidatahold = FitBoots(datafit,np.array(map(int,TSL[range(tmin,tmax+1)])),LinearFitFun)
                    bootdata[imom][icut].append(bootdatahold)
                    Avgdata[imom][icut].append(Avgdatahold)
                    Chidata[imom][icut].append(Chidatahold[0])
                    FitList[imom][icut].append((tmin,tmax))
                # print thisGamma , ' at ' , int((imom*100)/float(len(thisMomList))) , '% '
    # print thisGamma , ' complete, took:' , str(datetime.timedelta(seconds=(time.time()-start))) ,' h:m:s '
    return SumRin,bootdata,Avgdata,Chidata,FitList

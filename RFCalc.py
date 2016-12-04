#!/usr/bin/env python
from array import array
import os.path
import numpy as np
from Params import *
from MiscFuns import *
import operator as op
from copy import deepcopy

# def CalcRatioFactor(data2pt,data3pt,ZeroMomIndex):
#     RFOut = []
#     SqrtMomFac = []
#     for ip,data2ptp in enumerate(data2pt):
#         SqrtMomFac.append(CalcSqrtFac(data2pt[ZeroMomIndex],data2pt[ip],tsink))
#     for igamma,datagamma in enumerate(data3pt):
#         RFOut.append([])
#         for ip,datap in enumerate(datagamma):
#             RFOut[igamma].append([])
#             for it in range(tsource-1,tsink):
#                 RFOut[igamma][ip].append(datap[it]*SqrtMomFac[ip][it-tsource+1])
#                 RFOut[igamma][ip][it-tsource+1].Stats()
#     return [RFOut,SqrtMomFac]

def CalcSqrtFac(data2ptz,data2ptp,thistsink):
    SqrtFac = []
    if CHROMA:
        inttsink = int(thistsink)+1
        thistsrc = 0
    else:
        inttsink = int(thistsink)
        thistsrc = tsource-1
    two = data2ptz[inttsink-1]*data2ptp[inttsink-1]
    NegSQRTRatFac = False
    for ict in range(len(two.values)):
        if two.values[ict] < 0.0:
            NegSQRTRatFac = True
            two.values[ict] = abs(two.values[ict])
    two = two**(0.5)
    for it in range(thistsrc,inttsink): 
        tflip = inttsink-it+tsource-2
        one = data2ptz[it]/data2ptp[it]
        three = data2ptp[tflip]/data2ptz[tflip]
        inter = one*three
        for ict in range(len(inter.values)):
            if inter.values[ict] < 0.0:
                NegSQRTRatFac = True
                inter.values[ict] = abs(inter.values[ict])
        ott = ((inter)**(0.5))/two
        # ott = ott**(0.5)
        # ott.Stats()
        SqrtFac.append(ott)
    return [SqrtFac,NegSQRTRatFac]

    
# def DivSqrtFac(data2pt,thistsink):
#     SqrtFac = []
#     inttsink = int(thistsink)
#     for it in range(tsource-1,inttsink): 
#         SqrtFac.append(data2pt[inttsink-1])
#     return [SqrtFac,False]


#data3pt [ igamma , iset , it  ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
#data2pt [ iset , it ]  = bootstrap1 class (.Avg, .Std, .values, .nboot)
#RFout [ igamma , iset , it ] bs1

# def CalcRatioFactor(data2pt,data3pt,itsink,thisMomList,thisSetList):
#     RFOut = []
#     SqrtMomFac = []
#     for icset,(theset,setdata,setZPdata) in enumerate(zip(thisSetList,datap,data2pt[ZeroMomIndex])):
#         [CSFhold,Error] = CalcSqrtFac(setZPdata,setdata,itsink)
#         SqrtMomFac[ip].append(CSFhold)
#     # if Error == True: print 'Neg SqrtFactor for ' +qvecSet[pindex]+ ' ' + theset
#     for igamma,datagamma in enumerate(data3pt):
#         RFOut.append(np.rollaxis(NDimOpp2(datagamma,SqrtMomFac,1,RFMultFun,itsink),0,3))
#     return [np.array(RFOut),np.array(SqrtMomFac)]
                        
#data3pt [ iset , igamma , ip , it  ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
#data2pt [ iset , ip , it ]  = bootstrap1 class (.Avg, .Std, .values, .nboot)
#RFout [ iset , igamma , ip , it ] bs1

# def CalcRatioFactor(data2pt,data3pt,itsink):
#     RFOut = []
#     SqrtMomFac = []
#     for icset,setdata2pt in enumerate(data2pt):
#         SqrtMomFac.append([])
#         for ip,pdata2pt in enumerate(setdata2pt):
#             [CSFhold,Error] = CalcSqrtFac(deepcopy(setdata2pt[0]),deepcopy(pdata2pt),itsink)
#             SqrtMomFac[icset].append(CSFhold)
#             # if Error == True: print 'Neg SqrtFactor for ip: ' ,ip, ' iset: ' , icset
#     for icset,(setdata3pt,setsqrt) in enumerate(zip(data3pt,SqrtMomFac)):
#         RFOut.append([])
#         for igamma,gammadata3pt in enumerate(setdata3pt):
#             RFOut[icset].append([])
#             for ip,(pdata3pt,psqrt) in enumerate(zip(gammadata3pt,setsqrt)):
#                 RFOut[icset][igamma].append(RFMultFun(pdata3pt,psqrt,itsink))
#     return [np.array(RFOut),np.array(SqrtMomFac)]

def CalcRatioFactor(data2pt,data3pt,itsink,thisMomList):
    RFOut = []
    SqrtMomFac = []
    twoptMomList= GetAvgMomListip(thisMomList)
    for icset,setdata2pt in enumerate(data2pt):
        SqrtMomFac.append([])
        for ip,pdata2pt in enumerate(setdata2pt):
            passitsink = itsink
            # passitsink = len(pdata2pt)-1
            [CSFhold,Error] = CalcSqrtFac(deepcopy(setdata2pt[0]),deepcopy(pdata2pt),passitsink)
            SqrtMomFac[icset].append(CSFhold)
            # if Error == True: print 'Neg SqrtFactor for ip: ' ,ip, ' iset: ' , icset
    for icset,(setdata3pt,setsqrt) in enumerate(zip(data3pt,SqrtMomFac)):
        RFOut.append([])
        for igamma,gammadata3pt in enumerate(setdata3pt):
            RFOut[icset].append([])
            for ipc,(pdata3pt,ip) in enumerate(zip(gammadata3pt,thisMomList)):
                psqrt = setsqrt[GetAvgMomip(ip)]
                # if Debug:
                # print
                # # for it,(i3pt,i2pt) in enumerate(zip(pdata3pt,psqrt)):
                # for it2pt,d2pt in enumerate(data2pt[icset][GetAvgMomip(ip)]):
                #     print
                #     for it,i3pt in enumerate(pdata3pt):
                #     # print 'iGamma',igamma,'iset',icset,'pstr',ipTOqstr(ip),'pstr2pt',ipTOqstr(GetAvgMomip(ip),Avg=True),'it',it,'3pt',i3pt.Avg,'2pt',1/data2pt[icset][GetAvgMomip(ip)][-2].Avg,'Rfac',(i3pt/data2pt[icset][GetAvgMomip(ip)][-2]).Avg
                #         rfac = i3pt/d2pt
                #         rfac.Stats()
                #         print 'it2pt',it2pt,'it',it,'Rfac',rfac.Avg,rfac.Std

                RFOut[icset][igamma].append(RFMultFun(pdata3pt,psqrt,itsink))
                
    return [np.array(RFOut),np.array(SqrtMomFac)]



def RFMultFun(data,SQRTdata,itsink):
    dataout = []
    if CHROMA:
        datar = data[0:int(itsink)]
    else:
        datar = data[tsource-1:int(itsink)]
    for itz,(itdata,itSQRTdata) in enumerate(zip(datar,SQRTdata)):
        dataout.append(itdata*itSQRTdata)
        dataout[-1].Stats()
        # print 
        # print 'Average'
        # print 'c3pt',itdata.Avg,'c2pt',itSQRTdata.Avg,'Rfac',dataout[-1].Avg
        # for iboot,(i2pt,i3pt,iRfac) in enumerate(zip(itSQRTdata.values,itdata.values,dataout[-1].values)):
        #     print itz, iboot, 'c3pt',i3pt,'c2pt',i2pt,'Rfac',iRfac
    return np.array(dataout)


class SQRTRatFacError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

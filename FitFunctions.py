#!/usr/bin/env python

import numpy as np
from MiscFuns import VecDelta

#def CalcChiSqrdPDF(funct,params,datax,datayAvg,datayStd):
#     if len(datayAvg)-len(params) == 0:
#         return 0
#     else:
#         summedChi = 0        
#         for ix,iy,iyerr in zip(np.rollaxis(np.array(datax),1),datayAvg,datayStd):
#             summedChi = summedChi + ((iy-funct([ix],params))/iyerr)**2
#         return summedChi/(len(datayAvg)-len(params))
def makexunit(xin):
    try:
        return np.array(len(xin)*[1.0])
    except:
        return 1.0
                

def FormFactorO1(x,p):return p[0]*x[0]
def FormFactorO1Der(x,p):return [x[0]]

def FormFactorO2(x,p):return p[0]*x[0]+p[1]*x[1]
def FormFactorO2Der(x,p): return [x[0],x[1]]

def FormFactorO3(x,p):return p[0]*x[0]+p[1]*x[1]+p[2]*x[2]
def FormFactorO3Der(x,p):return [x[0],x[1],x[2]]


def ConstantFitFun(x,p):
    return p[0]*makexunit(x)

def ConstFFDer(x,p):
    return makexunit(x)

def LinearFitFun(x,p):
    return p[0]*np.array(x)+p[1]

def LinFFDer(x,p):
    return [x,makexunit(x)]

def TestTwoVarFitFun(x,p):
    return p[0]*x[0]+ p[1]*x[1]**2

def TestTwoVarFFDer(x,p):
    return [x[0],x[1]**2]

##ONE STATE FIT FUNCTIONS###

def C2OneStateFitFun(t,p):    
    A0,Ep = p[0],np.exp(p[1])
    return A0*np.exp(-Ep*t)

def C2OneStateFitFunNoExp(t,p):    
    A0,Ep = p[0],p[1]
    return A0*np.exp(Ep*(-t) )



def C2OSFFDer(t,p):
    A0,Ep = p[0],np.exp(p[1])
    return [np.exp(-Ep*t),-t*Ep*A0*np.exp(-Ep*t)]

def C3OneStateFitFun(tf,p):
    FitA0 = tf[2]
    FitAp0 = tf[4]
    Fitm0 = np.exp(tf[3])
    FitEp0 = np.exp(tf[5])
    return np.sqrt(np.abs(FitA0*FitAp0))*np.exp(-Fitm0*tf[1])*np.exp(-(FitEp0-Fitm0)*tf[0])*p[0]

def C3OneStateFitFunNoExp(tf,p):
    FitA0 = tf[2]
    FitAp0 = tf[4]
    Fitm0 = tf[3]
    FitEp0 = tf[5]
    return np.sqrt(np.abs(FitA0*FitAp0))*np.exp(-Fitm0*tf[1])*np.exp(-(FitEp0-Fitm0)*tf[0])*p[0]

def C3OSFFDer(tf,p):
    FitA0 = tf[2]
    FitAp0 = tf[4]
    Fitm0 = np.exp(tf[3])
    FitEp0 = np.exp(tf[5])
    return [np.sqrt(np.abs(FitA0*FitAp0))*np.exp(-Fitm0*tf[1])*np.exp(-(FitEp0-Fitm0)*tf[0])]


##Momenta TwoStateFit##

def C2TwoStateFitFunCM(t,p):
    Alen = (len(p)-2)/2
    Am,Amp,m,dm = np.array(p[:Alen]),np.array(p[Alen:-2]),np.exp(p[-2]),np.exp(p[-1])
    smpick = map(int,t[1])
    return Am[smpick]*(np.exp(-m*t[0]) + Amp[smpick]*np.exp(-(m+dm)*t[0]))

def C2TwoStateFitFunCMNoExp(t,p):
    Alen = (len(p)-2)/2
    Am,Amp,m,dm = np.array(p[:Alen]),np.array(p[Alen:-2]),p[-2],p[-1]
    smpick = map(int,t[1])
    return Am[smpick]*(np.exp(-m*t[0]) + Amp[smpick]*np.exp(-(m+dm)*t[0]))


def C2TSFLineFun(t,p):
    Am,Amp,m,dm = p[0],p[1],p[2],p[3] 
    try:
        return Am*(np.exp(-m*t) + Amp*np.exp(-(m+dm)*t))
    except:
        return Am*((m*(-t)).exp(1) + Amp*((m+dm)*(-t)).exp(1))
        

def C2TSFFCMDer(t,p):
    Alen = (len(p)-2)/2
    Am,Amp,m,dm = np.array(p[:Alen]),np.array(p[Alen:-2]),np.exp(p[-2]),np.exp(p[-1])
    smpick = map(int,t[1])
    pder = []
    for ism,Ami in enumerate(Am): # Am derivatives
        delta = VecDelta(smpick,ism)
        pder.append(delta*(np.exp(-m*t[0]) + Amp[smpick]*np.exp(-(m+dm)*t[0])))
    for ism,Ampi in enumerate(Amp): # Amp derivatives
        delta = VecDelta(smpick,ism)
        pder.append(delta*Am[smpick]*np.exp(-(m+dm)*t[0]))
    pder.append(-t[0]*m*Am[smpick]*(np.exp(-m*t[0]) + Amp[smpick]*np.exp(-(m+dm)*t[0]))) #mass
    pder.append(-t[0]*dm*Am[smpick]*Amp[smpick]*np.exp(-(m+dm)*t[0])) #Dm
    return pder

def C3TSFLineFun(Vals,thistcurr,thistsink):
    FitA0 = Vals[0]
    FitA1 = Vals[1]
    Fitm0 = Vals[2]
    FitDm = Vals[3]
    FitB00 = Vals[4]
    FitB10 = Vals[5]
    FitB01 = Vals[6]
    FitB11 = Vals[7]
    output =  (FitA0*np.exp(-Fitm0*thistsink) *
               (FitB00 + FitB10*(np.exp(-FitDm*thistcurr)+ np.exp(-FitDm*(thistsink-thistcurr))) +
                FitB11 *np.exp(-FitDm*thistsink)))
    return output



def C3MomTwoStateFitFun(tf,p):
    FitA0 = tf[2]
    FitA1 = tf[3]
    FitAp0 = tf[6]
    FitAp1 = tf[7]
    Fitm0 = np.exp(tf[4])
    FitDm = np.exp(tf[5])
    FitEp = np.exp(tf[8])
    FitDEp = np.exp(tf[9])
    if all(FitEp == Fitm0) and all(FitA0 == FitAp0):
        if all([tf[1][0]==itf for itf in tf[1]]):
            
            output = (FitA0*np.exp(-Fitm0*tf[1])*
                      (p[0] + p[1]*(np.exp(-FitDm*tf[0])+np.exp(-FitDm*(tf[1]-tf[0])))))
        else:
            output =  (FitA0*np.exp(-Fitm0*tf[1]) *
                       (p[0] + p[1]*(np.exp(-FitDm*tf[0])+ np.exp(-FitDm*(tf[1]-tf[0]))) +
                        p[3] *np.exp(-FitDm*tf[1])))
    else:
        if all([tf[1][0]==itf for itf in tf[1]]):
            output =  (np.sqrt(np.abs(FitA0*FitAp0))*np.exp(-Fitm0*tf[1]) * np.exp(-(FitEp-Fitm0)*tf[0])*
                       (p[0] + p[1]*np.exp(-FitDEp*tf[0])+ p[2]*np.exp(-FitDm*(tf[1]-tf[0]))))
        else:
            output =  (np.sqrt(np.abs(FitA0*FitAp0))*np.exp(-Fitm0*tf[1]) * np.exp(-(FitEp-Fitm0)*tf[0])*
                       (p[0] + p[1]*np.exp(-FitDEp*tf[0])+ p[2]*np.exp(-FitDm*(tf[1]-tf[0])) +
                        p[3] *np.exp(-FitDEp*tf[0])*np.exp(-FitDm*(tf[1]-tf[0]))))
    return output


def C3MomTwoStateFitFunNoExp(tf,p):
    FitA0 = tf[2]
    FitA1 = tf[3]
    FitAp0 = tf[6]
    FitAp1 = tf[7]
    Fitm0 = tf[4]
    FitDm = tf[5]
    FitEp = tf[8]
    FitDEp = tf[9]
    if all(FitEp == Fitm0) and all(FitA0 == FitAp0):
        if all([tf[1][0]==itf for itf in tf[1]]):
            
            output = (FitA0*np.exp(-Fitm0*tf[1])*
                      (p[0] + p[1]*(np.exp(-FitDm*tf[0])+np.exp(-FitDm*(tf[1]-tf[0])))))
        else:
            output =  (FitA0*np.exp(-Fitm0*tf[1]) *
                       (p[0] + p[1]*(np.exp(-FitDm*tf[0])+ np.exp(-FitDm*(tf[1]-tf[0]))) +
                        p[3] *np.exp(-FitDm*tf[1])))
    else:
        if all([tf[1][0]==itf for itf in tf[1]]):
            output =  (np.sqrt(np.abs(FitA0*FitAp0))*np.exp(-Fitm0*tf[1]) * np.exp(-(FitEp-Fitm0)*tf[0])*
                       (p[0] + p[1]*np.exp(-FitDEp*tf[0])+ p[2]*np.exp(-FitDm*(tf[1]-tf[0]))))
        else:
            output =  (np.sqrt(np.abs(FitA0*FitAp0))*np.exp(-Fitm0*tf[1]) * np.exp(-(FitEp-Fitm0)*tf[0])*
                       (p[0] + p[1]*np.exp(-FitDEp*tf[0])+ p[2]*np.exp(-FitDm*(tf[1]-tf[0])) +
                        p[3] *np.exp(-FitDEp*tf[0])*np.exp(-FitDm*(tf[1]-tf[0]))))
    return output


def C3MomTSFFDer(tf,p):
    FitA0 = tf[2]
    FitA1 = tf[3]
    FitAp0 = tf[6]
    FitAp1 = tf[7]
    Fitm0 = np.exp(tf[4])
    FitDm = np.exp(tf[5])
    FitEp = np.exp(tf[8])
    FitDEp = np.exp(tf[9])
    unitt = [0]*len(tf[1])
    if all(FitEp == Fitm0) and all(FitA0 == FitAp0):
        if all([tf[1][0]==itf for itf in tf[1]]):
            pder = [FitA0*np.exp(-Fitm0*tf[1])]
            pder.append(pder[0]*(np.exp(-FitDm*tf[0])+np.exp(-FitDm*(tf[1]-tf[0]))))
            pder.append(unitt)
            pder.append(unitt)
        else:
            pder = [FitA0*np.exp(-Fitm0*tf[1])]
            pder.append(pder[0]*(np.exp(-FitDm*tf[0])+np.exp(-FitDm*(tf[1]-tf[0]))))
            pder.append(unitt)
            pder.append(pder[0]*np.exp(-FitDm*tf[1]))
    else:
        if all([tf[1][0]==itf for itf in tf[1]]):
            pder = [np.sqrt(np.abs(FitA0*FitAp0))*np.exp(-Fitm0*tf[1]) * np.exp(-(FitEp-Fitm0)*tf[0])]
            pder.append(pder[0]*np.exp(-FitDEp*tf[0]))
            pder.append(pder[0]*np.exp(-FitDm*(tf[1]-tf[0])))
            pder.append(unitt)
        else:
            pder = [np.sqrt(np.abs(FitA0*FitAp0))*np.exp(-Fitm0*tf[1]) * np.exp(-(FitEp-Fitm0)*tf[0])]
            pder.append(pder[0]*np.exp(-FitDEp*tf[0]))
            pder.append(pder[0]*np.exp(-FitDm*(tf[1]-tf[0])))
            pder.append(pder[0]*np.exp(-FitDEp*tf[0])*np.exp(-FitDm*(tf[1]-tf[0])))    
    return pder



#aprox B2 <<<< B0




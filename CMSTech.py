#!/usr/bin/env python

import numpy as np
from scipy.linalg import eigvals,inv,sqrtm
from scipy.linalg import eigh,eig
# from numpy.linalg import eigh,eig
from Params import *
from copy import deepcopy
import operator as opp
from ReadTxt import ReadLREM
import time

# ## Cfun[t]
# def ProjectCorrJustPoF2pt(LEvec,Cfun,REvec,Booted=False,thisPoFShifts=1):
#     CMCfun = []
#     CfunShift = np.roll(Cfun,-1)
#     CfunShift2 = np.roll(CfunShift,-1)
#     CfunShift3 = np.roll(CfunShift2,-1)
#     CfunShift4 = np.roll(CfunShift3,-1)
#     if PoFShifts==1:
#         CfunExt = np.array([[Cfun,CfunShift],[CfunShift,CfunShift2]])
#     elif PoFShifts==2:
#         CfunExt = np.array([[Cfun,CfunShift,CfunShift2],[CfunShift,CfunShift2,CfunShift3],[CfunShift2,CfunShift3,CfunShift4]])
#     for istate,(stateRE,stateLE) in enumerate(zip(REvec,LEvec)):
#         CMCfun.append(np.dot(stateRE,np.dot(stateLE,CfunExt)))
#         if Booted:
#             for it,itCM in enumerate(CMCfun[istate]):
#                 CMCfun[istate][it].Stats()
#     return CMCfun


# #Cfun [ t ]
# def GetJustPoF(Cfunin,thistodtvals,Booted=False,thisPoFShifts=1):
#     thisto,thisdt = thistodtvals
#     Cfun = np.array(Cfunin)
#     if Booted:
#         Cfunto = Pullflag(Cfun[thisto],'Avg')
#         Cfuntop1 = Pullflag(Cfun[thisto+1],'Avg')
#         Cfuntop2 = Pullflag(Cfun[thisto+2],'Avg')
#         Cfuntop3 = Pullflag(Cfun[thisto+3],'Avg')
#         Cfuntop4 = Pullflag(Cfun[thisto+4],'Avg')
#         Cfuntodt = Pullflag(Cfun[thisto+thisdt],'Avg')
#         Cfuntodtp1 = Pullflag(Cfun[thisto+1+thisdt],'Avg')
#         Cfuntodtp2 = Pullflag(Cfun[thisto+2+thisdt],'Avg')
#         Cfuntodtp3 = Pullflag(Cfun[thisto+3+thisdt],'Avg')
#         Cfuntodtp4 = Pullflag(Cfun[thisto+4+thisdt],'Avg')
#     else:
#         Cfunto = Cfun[thisto]
#         Cfuntop1 = Cfun[thisto+1]
#         Cfuntop2 = Cfun[thisto+2]
#         Cfuntop3 = Cfun[thisto+3]
#         Cfuntop4 = Cfun[thisto+4]
#         Cfuntodt = Cfun[thisto+thisdt]
#         Cfuntodtp1 = Cfun[thisto+1+thisdt]
#         Cfuntodtp2 = Cfun[thisto+2+thisdt]
#         Cfuntodtp3 = Cfun[thisto+3+thisdt]
#         Cfuntodtp4 = Cfun[thisto+4+thisdt]
#     if thisPoFShifts==1:
#         Cfuntoout = np.array([[Cfunto,Cfuntop1],[Cfuntop1,Cfuntop2]])
#         Cfuntodtout = np.array([[Cfuntodt,Cfuntodtp1],[Cfuntodtp1,Cfuntodtp2]])
#     elif thisPoFShifts==2:
#         Cfuntoout = np.array([[Cfunto,Cfuntop1,Cfuntop2],[Cfuntop1,Cfuntop2,Cfuntop3],[Cfuntop2,Cfuntop3,Cfuntop4]])
#         Cfuntodtout = np.array([[Cfuntodt,Cfuntodtp1,Cfuntodtp2],[Cfuntodtp1,Cfuntodtp2,Cfuntodtp3],[Cfuntodtp2,Cfuntodtp3,Cfuntodtp4]])
#     [Emass,LEvec,REvec] = CreateLREves(Cfuntoout,Cfuntodtout,thisdt)
#     return Emass,LEvec,REvec


def CreatePoFMatrixtodt(Cfunto,Cfuntodt,thisPoFShifts=PoFShifts):
    if thisPoFShifts==0:
        Cfuntoout = Cfunto[0][0]
        Cfuntodtout = Cfuntodt[0][0]
    else:
        if thisPoFShifts==1:
            Cfuntoout = np.concatenate((np.concatenate((Cfunto[0][0],Cfunto[1][0]),1),
                                        np.concatenate((Cfunto[0][1],Cfunto[1][1]),1)))
            Cfuntodtout = np.concatenate((np.concatenate((Cfuntodt[0][0],Cfuntodt[1][0]),1),
                                          np.concatenate((Cfuntodt[0][1],Cfuntodt[1][1]),1)))
        elif thisPoFShifts==2:
            Cfuntoout = np.concatenate((np.concatenate((Cfunto[0][0],Cfunto[1][0],Cfunto[2][0]),1),
                                        np.concatenate((Cfunto[0][1],Cfunto[1][1],Cfunto[2][1]),1),
                                        np.concatenate((Cfunto[0][2],Cfunto[1][2],Cfunto[2][2]),1)))
            Cfuntodtout = np.concatenate((np.concatenate((Cfuntodt[0][0],Cfuntodt[1][0],Cfuntodt[2][0]),1),
                                          np.concatenate((Cfuntodt[0][1],Cfuntodt[1][1],Cfuntodt[2][1]),1),
                                          np.concatenate((Cfuntodt[0][2],Cfuntodt[1][2],Cfuntodt[2][2]),1)))
    return Cfuntoout,Cfuntodtout


## Cfun[tsource,ism,jsm,tcurr]

def CreatePoFMatrix(thisCfun,thisPoFShifts=PoFShifts):
    if thisPoFShifts==0:
        thisCfunExt = np.array(thisCfun[0])
    else:
        thisCfunShift = np.roll(thisCfun,-PoFDelta,axis=3)
        # if DEBUGPoF: thisCfunShift = np.roll(thisCfun,-PoFDelta/2,axis=3)
        if thisPoFShifts==1:
            if TimeInv:
                thisCfunShift2 = np.roll(thisCfunShift,-PoFDelta,axis=3)
                thisCfunExt = np.concatenate((np.concatenate((thisCfun[0],thisCfunShift[0]),1),
                                              np.concatenate((thisCfunShift[0],thisCfunShift2[0]),1)))
            else:
                thisCfunExt = np.concatenate((np.concatenate((thisCfun[0],thisCfunShift[0]),1),
                                              np.concatenate((thisCfun[1],thisCfunShift[1]),1)))
        elif thisPoFShifts==2:
            thisCfunShift2 = np.roll(thisCfunShift,-PoFDelta,axis=3)
            thisCfunExt = np.concatenate((np.concatenate((thisCfun[0],thisCfunShift[0],thisCfunShift2[0]),1),
                                          np.concatenate((thisCfun[1],thisCfunShift[1],thisCfunShift2[1]),1),
                                          np.concatenate((thisCfun[2],thisCfunShift[2],thisCfunShift2[2]),1)))
    return thisCfunExt

# ##DEBUGGING POF DECORRELATION
# def CreatePoFMatrixDECORRELATION(thisCfun,thisPoFShifts=PoFShifts):
#     if thisPoFShifts==0:
#         thisCfunExt = thisCfun
#     else:
#         if thisPoFShifts==1:
#             thisCfunShift = np.roll(thisCfun,-1,axis=2)
#             bottomCfunShift = DeCorrBoot(thisCfunShift)
#             bottomCfunShift2 = DeCorrBoot(np.roll(thisCfunShift,-1,axis=2))
#             thisCfunExt = np.concatenate((np.concatenate((thisCfun,thisCfunShift),1),
#                                           np.concatenate((bottomCfunShift,bottomCfunShift2),1)))
#         elif thisPoFShifts==2:
#             thisCfunShift3 = np.roll(thisCfunShift2,-1,axis=2)
#             thisCfunShift4 = np.roll(thisCfunShift3,-1,axis=2)
#             thisCfunExt = np.concatenate((np.concatenate((thisCfun,thisCfunShift,thisCfunShift2),1),
#                                           np.concatenate((thisCfunShift,thisCfunShift2,thisCfunShift3),1),
#                                           np.concatenate((thisCfunShift2,thisCfunShift3,thisCfunShift4),1)))
#     return thisCfunExt


## Cfun[tsink,tsource,ism,1,tcurr]
def CreateREvecProjPoFMatrix(thisCfun):
    thisCfunjsm = np.array(thisCfun)[0,0,:,0,:]
    if PoFShifts==0:
        thisCfunExt = thisCfunjsm
    else:
        if PoFShifts==1:
            if TimeInv:
                thisCfunjsmShift = np.array(thisCfun)[1,0,:,0,:]
                thisCfunjsmShift = np.roll(thisCfunjsmShift,-1,axis=2)
            else:
                thisCfunjsmShift = np.array(thisCfun)[0,1,:,0,:]
                # thisCfunjsmShift = np.roll(thisCfunjsmShift,-1,axis=2)
            thisCfunExt = np.concatenate((thisCfunjsm,thisCfunjsmShift))
        else:
            raise EnvironmentError('thisPoFShifts > 1 is not support yet for 3pt correlators, combine correlators here if you want to implement')
    return thisCfunExt

## Emass = 101 is invalid eigenvalue out
## Emass = 202 is b matrix not positive definite, removing higher smearings.

## Evec = [weights][state]
## SEmass = [state]

def AddNullState(Em,Le,Re,buffindex,thisdt=1):
    for ibuff in buffindex:
    # for ibuff in reversed(buffindex):
        Em = np.insert(Em,ibuff,np.exp(-202.0*thisdt))
        Le = [np.insert(iLe,ibuff,0.0) for iLe in Le]
        Re = [np.insert(iRe,ibuff,0.0) for iRe in Re]
        Le = np.insert(Le,ibuff,[0.0]*len(Em),axis=0)
        Re = np.insert(Re,ibuff,[0.0]*len(Em),axis=0)
    return Em,Le,Re
    

#sourts Evec and Emass with respect to Emass lowest to highest
def sortEvec(Evals,LEvec,REvec,thisdt):
    # for ice,iev in enumerate(Evals):
    #     if abs(iev.imag) > myeps:
    #         Evals[ice] = np.exp(101)+0j
    #     elif iev.real < 0.0:
    #         Evals[ice] = -iev
    Emass = -np.log(abs(Evals))/float(thisdt)
    # for ice,iem in enumerate(Emass):
    #     if iem < VarMassCutoff:
    #         Emass[ice] = 101
    
    sortindex = Emass.argsort()
    Emass = Emass[sortindex]
    LEvec = np.swapaxes(LEvec,0,1)[sortindex,:]
    REvec = np.swapaxes(REvec,0,1)[sortindex,:]
    return Emass,LEvec,REvec



# Mat*Evec = Evals * Evec
# Mat = C(to)^-1 * C(to+dt)
def CreateLREves(Cfunto,Cfuntodt,thisdt,masscutoff):
    # Ctoinv = inv(Cfunto)
    # Mat = np.dot(Ctoinv,Cfuntodt)
    # [Evals,REvec] = eig(Mat,right=True,left=False)
    # LEvec = REvec
    # [Evals,LEvec,REvec] = eig(Cfuntodt,b=Cfunto,left=True,right=True)
    Simto = np.array(Cfunto)
    Simtodt = np.array(Cfuntodt)
    buffindex = []
    ci = np.array([0])
    Ctolen = len(Cfunto)
    for cutindex in xrange(1,Ctolen):        
        ci = np.append(ci,cutindex)
        if 'Symmetric' in VarMethodMethod:
            ShalfInv = inv(sqrtm(Simto[ci[:,None],ci]))
            ThisMat = ShalfInv.dot(Simtodt[ci[:,None],ci].dot(ShalfInv))
            # for iindex in xrange(len(ThisMat)):
            #     for jindex in xrange(len(ThisMat)):
            #         print iindex, jindex, ThisMat[iindex][jindex], ThisMat[jindex][iindex]
            if 'eigh' in VarMethodMethod:
                thiseig,thisevec = eigh(ThisMat)
            else:
                thiseig,thisevec = eig(ThisMat)
        elif VarMethodMethod == 'Regular':
            ThisMat = Simtodt[ci[:,None],ci].dot(inv(Simto[ci[:,None],ci]))
            thiseig,thisevec = eig(ThisMat)
        elif VarMethodMethod == 'AxBxlSolve':
            thiseig,thisevec = eig(Simtodt[ci[:,None],ci],b=Simto[ci[:,None],ci])
        else:
            raise LookupError('VarMethodMethod not recognised (Params.py) : ' + VarMethodMethod)
        evecreal,evecimag = SplitCmplxReal(thisevec.flatten())
        eigreal,eigimag = SplitCmplxReal(thiseig)
        if any(-np.log(np.abs(eigreal))/float(thisdt) < VarMassCutoff) or any(np.array(eigreal) < 0) or any(np.abs(eigimag) > 0)  or any(np.abs(evecimag) > 0):
            # ibad = [ie < 0 for ie in thiseig].index(True)
            ci = np.delete(ci,ci.tolist().index(cutindex))
            buffindex.append(cutindex)
    # Simto = Cfunto
    # Simtodt = Cfuntodt
    # buffindex = []
    # for cutindex in xrange(len(Cfunto)-1):
    #     thiseig = eigvals(Simtodt,b=Simto)
    #     posdef = eigvals(Simto)
    #     if all(-np.log(abs(thiseig))/float(thisdt) > masscutoff) and all(posdef > 0):
    #         break
    #     else:
    #         # ibad = [ie < 0 for ie in thiseig].index(True)
    #         ibad = cutcmList[cutindex]
    #         buffindex.append(ibad)
    #         Simto = np.delete(np.delete(Simto,ibad,0),ibad,1)
    #         Simtodt = np.delete(np.delete(Simtodt,ibad,0),ibad,1)

    if len(ci) < 2:
        Evals,LEvec,REvec = AddNullState(np.ones(1),np.array([[1]]),np.array([[1]]),buffindex,thisdt=thisdt)
    else:
        if 'Symmetric' in VarMethodMethod:
            Simto = np.array(Cfunto)[ci[:,None],ci]
            Simtodt = np.array(Cfuntodt)[ci[:,None],ci]
            ShalfInv = inv(sqrtm(Simto))
            ThisMat = ShalfInv.dot(Simtodt.dot(ShalfInv))
            if 'eigh' in VarMethodMethod:
                [Evals,REvec] = eigh(ThisMat)
                REvec = ShalfInv.dot(REvec)
                LEvec = REvec
            else:
                [Evals,LEvec,REvec] = eig(ThisMat,left=True)
                REvec = ShalfInv.dot(REvec)
                LEvec = np.swapaxes(np.swapaxes(LEvec,0,1).dot(ShalfInv),0,1)
        elif VarMethodMethod == 'Regular':
            Simto = np.array(Cfunto)[ci[:,None],ci]
            Simtodt = np.array(Cfuntodt)[ci[:,None],ci]
            ThisMat = Simtodt.dot(inv(Simto))
            [Evals,LEvec,dump] = eig(ThisMat,left=True)
            ThisMat = inv(Simto).dot(Simtodt)
            [Evals,REvec] = eig(ThisMat)
        elif VarMethodMethod == 'AxBxlSolve':
            Simto = np.array(Cfunto)[ci[:,None],ci]
            Simtodt = np.array(Cfuntodt)[ci[:,None],ci]
            [Evals,LEvec,REvec] = eig(Simtodt,b=Simto,left=True)
        else:
            raise LookupError('VarMethodMethod not recognised (Params.py) : ' + VarMethodMethod)
        ## w = G^-1/2 u
        Evals,LEvec,REvec = AddNullState(Evals,LEvec,REvec,buffindex,thisdt=thisdt)
    return sortEvec(Evals,LEvec,REvec,thisdt)


# Mat*Evec = Evals * Evec
# Mat = C(to)^-1* C(to+dt)
# def CreateLREvesSum(Mat):
#     [Evals,LEvec] = eig(Mat)
#     REvec = LEvec
#     return sortEvec(Evals,LEvec,REvec,thisdt)

# ## Cfun[ism,jsm,tcurr]
# def ProjectCorr(LEvec,Cfun,REvec):
#     CMCfun = []
#     for istate,(stateRE,stateLE) in enumerate(zip(REvec,LEvec)):
#         CMCfun.append(np.dot(stateRE,np.dot(stateLE,Cfun)))
#         for it,itCM in enumerate(CMCfun[istate]):
#             CMCfun[istate][it].Stats()
#     return CMCfun

## Cfun[tsource,ism,jsm,tcurr]
def ProjectCorrPoF2pt(LEvec,Cfun,REvec,thisPoFShifts=PoFShifts):
    if DeCorrPoF:
        CfunExt = CreatePoFMatrixDECORRELATION(Cfun,thisPoFShifts=thisPoFShifts)
    else:
        CfunExt = CreatePoFMatrix(Cfun,thisPoFShifts=thisPoFShifts)
    # CfunShift = np.roll(Cfun,-1,axis=2)
    # CfunShift2 = np.roll(CfunShift,-1,axis=2)
    # if thisPoFShifts==0:
    #     CfunExt = Cfun
    # elif thisPoFShifts==1:
    #     CfunExt = np.concatenate((np.concatenate((Cfun,CfunShift),1),np.concatenate((CfunShift,CfunShift2),1)))
    # elif thisPoFShifts==2:
    #     CfunShift3 = np.roll(CfunShift2,-1,axis=2)
    #     CfunShift4 = np.roll(CfunShift3,-1,axis=2)
    #     CfunExt = np.concatenate((np.concatenate((Cfun,CfunShift,CfunShift2),1),
    #                               np.concatenate((CfunShift,CfunShift2,CfunShift3),1),
    #                               np.concatenate((CfunShift2,CfunShift3,CfunShift4),1)))

    ##DEBUG##
    CMCfun = []
    for istate,(stateRE,stateLE) in enumerate(zip(REvec,LEvec)):
        CMCfun.append(np.dot(stateRE,np.dot(stateLE,CfunExt)))
        for it,itCM in enumerate(CMCfun[istate]):
            CMCfun[istate][it].Stats()
            # print 'istate',istate, 'it',it, 'values',CMCfun[istate][it].Avg, CMCfun[istate][it].Std
    if Debug:
        print 'TwoPoint Run:'
        for ic,(iRE,iCfun) in enumerate(zip(REvec[0],np.dot(LEvec[0],CfunExt))):
            iCfun[5].Stats()
            print '5',ic,iRE,iCfun[5].Avg, CMCfun[ic][5].Avg
        print ''
    return CMCfun


# ## Cfun[ism,1,tcurr]
# def ProjectREvecCorr(Cfun,REvec):
#     CMCfun = []
#     Cfunjsm = np.array(Cfun)[:,0,:]
#     for istate,stateRE in enumerate(REvec):
#         CMCfun.append(np.dot(stateRE,Cfunjsm))
#         for it,itCM in enumerate(CMCfun[istate]):
#             CMCfun[istate][it].Stats()
#     return CMCfun


## Cfun[tsink,tsource,ism,1,tcurr]
def ProjectREvecCorrPoF(Cfun,REvec):
    CMCfun = []
    CfunExt = CreateREvecProjPoFMatrix(Cfun)
    for istate,stateRE in enumerate(REvec):
        CMCfun.append(np.dot(stateRE,CfunExt))
        for it,itCM in enumerate(CMCfun[istate]):
            CMCfun[istate][it].Stats()
            # print istate, it, CMCfun[istate][it].Avg, CMCfun[istate][it].Std
    ##DEBUG##
    # if Debug:
    #     print 'ThreePoint Run:'
    #     for ic,(iRE,iCfun) in enumerate(zip(REvec[0],CfunExt)):
    #         print ic,iRE,iCfun[5].Avg, CMCfun[ic][5].Avg
    #     print ''
    return CMCfun

def ProjectREvecCorr(Cfun,REvec):
    CMCfun = []
    CfunExt = np.array(Cfun)[:,0,:]
    ##DEBUG##
    if Debug:
        print 'ThreePoint Run:'
        # for ic,(iRE,iCfun) in enumerate(zip(REvec[0],CfunExt)):
        #     print ic,iRE,iCfun[24].Avg
        print ''
    for istate,stateRE in enumerate(REvec):
        CMCfun.append(np.dot(stateRE,CfunExt))
        for it,itCM in enumerate(CMCfun[istate]):
            CMCfun[istate][it].Stats()
    return CMCfun

def GetTvarREves(Cfunin,thistodtvals,masscut):
    thisto,thisdt = thistodtvals
    Cfun = np.array(Cfunin)
    Cfunto = Pullflag(Cfun[:,:,thisto-1],'Avg')
    Cfuntodt = Pullflag(Cfun[:,:,thisto-1+thisdt],'Avg')
    [Emass,LEvec,REvec] = CreateLREves(Cfunto,Cfuntodt,thisdt,masscut)
    return Emass,LEvec,REvec



#Cfun [ t_src, t , ism , jsm ]
def GetTvarREvesPoF(Cfunin,thistodtvals,masscut,thisPoFShifts=PoFShifts):
    thisto,thisdt = thistodtvals
    thisto = thisto-1 #python arrays start from 0, not 1
    Cfun = np.array(Cfunin)
    Cfuntomat,Cfuntodtmat = [[]],[[]]
    
    if CHROMA:
        Cfuntomat[0].append(Pullflag(Cfun[-1,:,:,thisto-PoFDelta],'Avg'))
        Cfuntodtmat[0].append(Pullflag(Cfun[-1,:,:,thisto+thisdt-PoFDelta],'Avg'))
    else:
        Cfuntomat[0].append(Pullflag(Cfun[-1,:,:,thisto],'Avg'))
        Cfuntodtmat[0].append(Pullflag(Cfun[-1,:,:,thisto+thisdt],'Avg'))
    if thisPoFShifts > 0:
        if TimeInv:
            Cfuntomat[0].append(Pullflag(Cfun[-1,:,:,thisto+PoFDelta],'Avg'))
            Cfuntomat.append([])
            Cfuntomat[1].append(Pullflag(Cfun[-1,:,:,thisto+PoFDelta],'Avg'))
            Cfuntomat[1].append(Pullflag(Cfun[-1,:,:,thisto+(2*PoFDelta)],'Avg'))
            Cfuntodtmat[0].append(Pullflag(Cfun[-1,:,:,thisto+thisdt+PoFDelta],'Avg'))
            Cfuntodtmat.append([])
            Cfuntodtmat[1].append(Pullflag(Cfun[-1,:,:,thisto+thisdt+PoFDelta],'Avg'))
            Cfuntodtmat[1].append(Pullflag(Cfun[-1,:,:,thisto+thisdt+(2*PoFDelta)],'Avg'))
        else:
            # if DEBUGPoF:
            #     Cfuntomat[0].append(Pullflag(Cfun[-1,:,:,thisto+PoFDelta],'Avg'))
            #     Cfuntomat.append([])
            #     Cfuntomat[1].append(Pullflag(Cfun[-2,:,:,thisto+PoFDelta/2],'Avg'))
            #     Cfuntomat[1].append(Pullflag(Cfun[-2,:,:,thisto+(3*PoFDelta)/2],'Avg'))
            #     Cfuntodtmat[0].append(Pullflag(Cfun[-1,:,:,thisto+thisdt+PoFDelta],'Avg'))
            #     Cfuntodtmat.append([])
            #     Cfuntodtmat[1].append(Pullflag(Cfun[-2,:,:,thisto+thisdt+PoFDelta/2],'Avg'))
            #     Cfuntodtmat[1].append(Pullflag(Cfun[-2,:,:,thisto+thisdt+(3*PoFDelta)/2],'Avg'))
            # else:
            if CHROMA:
                Cfuntomat[0].append(Pullflag(Cfun[-1,:,:,thisto],'Avg'))
                Cfuntomat.append([])
                Cfuntomat[1].append(Pullflag(Cfun[-2,:,:,thisto],'Avg'))
                Cfuntomat[1].append(Pullflag(Cfun[-2,:,:,thisto+PoFDelta],'Avg'))
                Cfuntodtmat[0].append(Pullflag(Cfun[-1,:,:,thisto+thisdt],'Avg'))
                Cfuntodtmat.append([])
                Cfuntodtmat[1].append(Pullflag(Cfun[-2,:,:,thisto+thisdt],'Avg'))
                Cfuntodtmat[1].append(Pullflag(Cfun[-2,:,:,thisto+thisdt+PoFDelta],'Avg'))
            else:
                Cfuntomat[0].append(Pullflag(Cfun[-1,:,:,thisto+PoFDelta],'Avg'))
                Cfuntomat.append([])
                Cfuntomat[1].append(Pullflag(Cfun[-2,:,:,thisto],'Avg'))
                Cfuntomat[1].append(Pullflag(Cfun[-2,:,:,thisto+PoFDelta],'Avg'))
                Cfuntodtmat[0].append(Pullflag(Cfun[-1,:,:,thisto+thisdt+PoFDelta],'Avg'))
                Cfuntodtmat.append([])
                Cfuntodtmat[1].append(Pullflag(Cfun[-2,:,:,thisto+thisdt],'Avg'))
                Cfuntodtmat[1].append(Pullflag(Cfun[-2,:,:,thisto+thisdt+PoFDelta],'Avg'))
        
        if thisPoFShifts > 1:
            Cfuntomat[0].append(Pullflag(Cfun[-1,:,:,thisto+(2*PoFDelta)],'Avg'))
            Cfuntomat[1].append(Pullflag(Cfun[-2,:,:,thisto+(2*PoFDelta)],'Avg'))
            Cfuntomat.append([])
            Cfuntomat[2].append(Pullflag(Cfun[-3,:,:,thisto],'Avg'))
            Cfuntomat[2].append(Pullflag(Cfun[-3,:,:,thisto+PoFDelta],'Avg'))
            Cfuntomat[2].append(Pullflag(Cfun[-3,:,:,thisto+(2*PoFDelta)],'Avg'))
            Cfuntodtmat[0].append(Pullflag(Cfun[-1,:,:,thisto+thisdt+(2*PoFDelta)],'Avg'))
            Cfuntodtmat[1].append(Pullflag(Cfun[-2,:,:,thisto+thisdt+(2*PoFDelta)],'Avg'))
            Cfuntodtmat.append([])
            Cfuntodtmat[2].append(Pullflag(Cfun[-3,:,:,thisto+thisdt],'Avg'))
            Cfuntodtmat[2].append(Pullflag(Cfun[-3,:,:,thisto+thisdt+PoFDelta],'Avg'))
            Cfuntodtmat[2].append(Pullflag(Cfun[-3,:,:,thisto+thisdt+(2*PoFDelta)],'Avg'))
    Cfuntoout,Cfuntodtout = CreatePoFMatrixtodt(Cfuntomat,Cfuntodtmat,thisPoFShifts=PoFShifts)
    [Emass,LEvec,REvec] = CreateLREves(Cfuntoout,Cfuntodtout,thisdt,masscut)
    return Emass,LEvec,REvec

# CM
# #Cfuns2pt [ ism , jsm , ip , it ]  = bootstrap1 class (.Avg, .Std, .values, .nboot)
# PoF
# #Cfuns2pt [ tsource, ism , jsm , ip , it ]  = bootstrap1 class (.Avg, .Std, .values, .nboot)

# #CMCfun2pt [ istate , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
# R/L Evecs [ ip , istate , ival ]
# Emass [ ip , istate ]

def CreatePoF2ptCfuns(Cfuns2pt,todtvals,thisMomList,printout=True,DoPoF=True,todtvalsLeft=False):
    if todtvalsLeft == False: todtvalsLeft = todtvals
    start = time.time()
    Cfuns2pt = np.array(Cfuns2pt)
    CMCfun2pt = []
    ##if ip**2 > pcutoff, make cfuns exclude sm128
    # LEvec,REvec,Emass = ReadLREM(todtvals,thisMomList,'PoF'+str(PoFShifts),NoWar=True)
    Dump,REvec,Emass = ReadLREM(todtvals,thisMomList,'CM',NoWar=True)
    LEvec,Dump,EmassLeft = ReadLREM(todtvalsLeft,thisMomList,'CM',NoWar=True)
    if (not DoPoF) or (not ReadPoF2pt) or (LEvec == None):
        Emass,LEvec,REvec = [],[],[]
        for ip,thisp in enumerate(thisMomList):
            if DoPoF:
                # Emasshold,LEvechold,REvechold = GetTvarREvesPoF(Cfuns2pt[:,:,ip],todtvals,ipTOE(thisp,VarMassCutoff))
                Emasshold,LEvechold,REvechold = GetTvarREvesPoF(Cfuns2pt[:,:,:,ip],todtvals,ipTOE(thisp,VarMassCutoff))
                # for istate in xrange(len(Emasshold)):
                #     print ' '.join(map(str,LEvechold[istate])) , Emasshold[istate]
            else:
                Emasshold,LEvechold,REvechold = GetTvarREves(Cfuns2pt[:,:,ip],todtvals,ipTOE(thisp,VarMassCutoff))
            Emass.append(Emasshold)
            LEvec.append(LEvechold)
            REvec.append(REvechold)
        LEvec,REvec = SignEvec(np.array(LEvec),np.array(REvec))
    for ip,thisp in enumerate(thisMomList):
        if DoPoF:
            CMCfun2pt.append(ProjectCorrPoF2pt(LEvec[ip],Cfuns2pt[:,:,:,ip],REvec[ip]))
            # CMCfun2pt.append(ProjectCorrPoF2pt(np.array([[1,0,0],[0,1,0],[0,0,1]]),Cfuns2pt[:,:,:,ip],np.array([[1,0,0],[0,1,0],[0,0,1]])))
        else:
            CMCfun2pt.append(ProjectCorrPoF2pt(LEvec[ip],np.array([Cfuns2pt[:,:,ip]]),REvec[ip],thisPoFShifts=0))
            # CMCfun2pt.append(ProjectCorrPoF2pt(np.array([[1,0,0],[0,1,0],[0,0,1]]),[Cfuns2pt[:,:,ip]],np.array([[1,0,0],[0,1,0],[0,0,1]]),thisPoFShifts=0))
    if printout:
        if DoPoF:
            print 'CM PoF Creation shift'+str(PoFShifts)+' to'+str(todtvals[0])+' dt'+str(todtvals[1])+ ' took: ' , GetTimeStr(time.time()-start)
        else:
            print 'CM Creation to'+str(todtvals[0])+' dt'+str(todtvals[1])+ ' took: ' , GetTimeStr(time.time()-start)
    return [np.rollaxis(np.array(CMCfun2pt),1),LEvec,REvec,Emass]


def CreateCM2ptCfuns(Cfuns2pt,todtvals,thisMomList,printout=True):
    return CreatePoF2ptCfuns(Cfuns2pt,todtvals,thisMomList,printout=printout,DoPoF=False)

# #Cfuns3pt [ tsink , ism , 1 , igamma , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
# #CMCfun3pt  [ istate , igamma , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
def CreateREvecCfuns(Cfuns3pt,Cfuns2pt,todtvals,thisMomList):
    LEvec,REvec,Emass = ReadLREM(todtvals,thisMomList,'CM')
    if REvec == None:
        CMCfun2pt,LEvec,REvec,Emass = CreatePoF2ptCfuns(Cfuns2pt,todtvals,thisMomList,DoPoF=False,printout=False)
        LEvec,REvec = SignEvec(np.array(LEvec),np.array(REvec))
    else:
        CMCfun2pt = []
        for ip,thisp in enumerate(thisMomList):
            CMCfun2pt.append(ProjectCorrPoF2pt(LEvec[ip],np.array([Cfuns2pt[:,:,ip]]),REvec[ip],thisPoFShifts=0))
    CMCfun3pt = []
    for igamma,gammaCfun in enumerate(np.rollaxis(np.array(Cfuns3pt),2)):
        CMCfun3pt.append([])
        for ip,pCfun in enumerate(np.rollaxis(gammaCfun,2)): 
            CMCfun3pt[igamma].append(ProjectREvecCorr(pCfun,None,REvec[ip],thisPoFShifts=0))
    return [np.rollaxis(np.array(CMCfun2pt),1),np.rollaxis(np.array(CMCfun3pt),2)]


# #Cfuns3pt [ tsink, tsource , ism , 1 , igamma , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
# #CMCfun3pt  [ istate , igamma , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
def CreateREPoFCfuns(Cfuns3pt,Cfuns2pt,todtvals,thisMomList,todtvalsLeft=False):
    if todtvalsLeft == False: todtvalsLeft = todtvals
    twoptMomList = GetAvgMomListip(thisMomList)
    if PoFShifts == 0:
        Dump,REvec,Emass = ReadLREM(todtvals,twoptMomList,'CM')
        LEvec,Dump,EmassLeft = ReadLREM(todtvalsLeft,twoptMomList,'CM')
    else:
        Dump,REvec,Emass = ReadLREM(todtvals,twoptMomList,'PoF'+str(PoFShifts))
        LEvec,Dump,EmassLeft = ReadLREM(todtvalsLeft,twoptMomList,'PoF'+str(PoFShifts))
        
    if REvec == None:
        CMCfun2pt,LEvec,REvec,Emass = CreatePoF2ptCfuns(Cfuns2pt,todtvals,twoptMomList,DoPoF=True,printout=False)
        LEvec,REvec = SignEvec(np.array(LEvec),np.array(REvec))
    else:
        CMCfun2pt = []
        for ip,thisp in enumerate(twoptMomList):
            CMCfun2pt.append(ProjectCorrPoF2pt(LEvec[ip],Cfuns2pt[:,:,:,ip],REvec[ip]))
        CMCfun2pt = np.rollaxis(np.array(CMCfun2pt),1)
    CMCfun3pt = []
    for igamma,gammaCfun in enumerate(np.rollaxis(np.array(Cfuns3pt),4)):
        CMCfun3pt.append([])
        for ipc,(pCfun,ip) in enumerate(zip(np.rollaxis(gammaCfun,4),thisMomList)): 
            # if Debug: print 'Three Point :(igamma,ip) ',igamma , ip, 'with twopt mom' , GetAvgMomip(ip)
            CMCfun3pt[igamma].append(ProjectREvecCorrPoF(pCfun,REvec[twoptMomList.index(GetAvgMomip(ip))]))
    return [np.array(CMCfun2pt),np.rollaxis(np.array(CMCfun3pt),2)]


# #Cfuns3pt [ ism , jsm , igamma , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
# #CMCfun3pt  [ istate , igamma , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)

def CreateCMCfuns(Cfuns3pt,Cfuns2pt,todtvals,thisMomList):
    CMCfun2pt,LEvec,REvec,Emass = CreatePoF2ptCfuns(Cfuns2pt,todtvals,thisMomList,DoPoF=False,printout=False)
    CMCfun3pt = []
    LEvec,REvec = SignEvec(LEvec,REvec)
    for igamma,gammaCfun in enumerate(np.rollaxis(np.array(Cfuns3pt),2)):
        CMCfun3pt.append([])
        for ip,pCfun in enumerate(np.rollaxis(gammaCfun,2)):
            CMCfun3pt[igamma].append(ProjectCorrPoF2pt(LEvec[0],np.array([pCfun]),REvec[ip],thisPoFShifts=0))
            
    return [np.array(CMCfun2pt),np.rollaxis(np.array(CMCfun3pt),2)]


def SignEvec(LEvec,REvec):
    LEvecOut = LEvec
    REvecOut = REvec
    for ip,(pLE,pRE) in enumerate(zip(LEvec,REvec)):
        for istate,(stateLE,stateRE) in enumerate(zip(pLE,pRE)):
            stateLEZ,stateREZ = LEvec[0][istate].tolist(),REvec[0][istate].tolist()
            LMaxI = stateLE.tolist().index(max(stateLE,key=abs))
            RMaxI = stateRE.tolist().index(max(stateRE,key=abs))
            LMaxIZ = stateLEZ.index(max(np.array(stateLEZ),key=abs))
            RMaxIZ = stateREZ.index(max(np.array(stateREZ),key=abs))
            ZMLsign = np.sign(stateREZ[LMaxIZ]*stateLE[LMaxI])
            ZMRsign = np.sign(stateREZ[RMaxIZ]*stateRE[RMaxI])
            # normL,normR = np.sum(stateLE),np.sum(stateRE)
            for ival,(LEvecVal,REvecVal) in enumerate(zip(stateLE,stateRE)):
                LEvecOut[ip][istate][ival] = ZMLsign*LEvecVal
                REvecOut[ip][istate][ival] = ZMRsign*REvecVal
    return LEvecOut,REvecOut

def Normalise(data2pt,tval):
    dataout = deepcopy(data2pt)
    for itsrc,datatsrc in enumerate(data2pt):
        for ism,dataism in enumerate(datatsrc):
            for jsm,datajsm in enumerate(dataism):
                for ip,datap in enumerate(datajsm):
                    norm = np.sqrt(data2pt[itsrc][ism][ism][ip][tval].Avg*data2pt[itsrc][jsm][jsm][ip][tval].Avg)
                    for it,datat in enumerate(datap):
                        dataout[itsrc][ism][jsm][ip][it] = datat/norm
                        dataout[itsrc][ism][jsm][ip][it].Stats()
    return dataout

def NormaliseNoP(data2pt,tval):
    dataout = deepcopy(data2pt)
    for itsrc,datatsrc in enumerate(data2pt):
        for ism,dataism in enumerate(datatsrc):
            for jsm,datajsm in enumerate(dataism):
                norm = np.sqrt(data2pt[itsrc][ism][ism][tval].Avg*data2pt[itsrc][jsm][jsm][tval].Avg)
                for it,datat in enumerate(datajsm):
                    dataout[itsrc][ism][jsm][it] = datat/norm
                    dataout[itsrc][ism][jsm][it].Stats()
    return dataout

def Symmetrize(data2pt):
    dataout = deepcopy(data2pt)
    for itsrc,datatsrc in enumerate(data2pt):
        for ism,dataism in enumerate(datatsrc):
            for jsm,datajsm in enumerate(dataism):
                for ip,datap in enumerate(datajsm):
                    dataout[itsrc][ism][jsm][ip] = (datap+datatsrc[jsm][ism][ip])/2.0
    GetBootStats(dataout)
    return dataout

def SymmetrizeNoPAvg(data2pt):
    dataout = deepcopy(data2pt)
    for itsrc,datatsrc in enumerate(data2pt):
        for ism,dataism in enumerate(datatsrc):
            for jsm,datajsm in enumerate(dataism):
                dataout[itsrc][ism][jsm] = (datajsm+datatsrc[jsm][ism])/2.0
    return dataout

def PreptwoptCorr(data2pt):
    if DoNorm and DoSym: return Symmetrize(Normalise(np.array(data2pt),tsource))
    elif DoNorm: return  Normalise(np.array(data2pt),tsource)
    elif DoSym: return Symmetrize(np.array(data2pt))
    # else: return data2pt
    else: return np.array(data2pt)
    # else: return np.array(data2pt)




#  LocalWords:  thisMomList

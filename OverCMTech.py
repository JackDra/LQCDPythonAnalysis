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
from CMSTech import SignEvec



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
def CreateLREves(Cfunto,Cfuntodt,thisdt):
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
    for cutindex in range(1,Ctolen):
        
        ci = np.append(ci,cutindex)
        if 'Symmetric' in VarMethodMethod:
            ShalfInv = inv(sqrtm(Simto[ci[:,None],ci]))
            ThisMat = ShalfInv.dot(Simtodt[ci[:,None],ci].dot(ShalfInv))
            # for iindex in range(len(ThisMat)):
            #     for jindex in range(len(ThisMat)):
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



## Cfun[tsource,ism,jsm,tcurr]
def ProjectCorrPoF2pt(LEvec,Cfun,REvec,thisPoFShifts=PoFShifts):
    if DeCorrPoF:
        CfunExt = CreatePoFMatrixDECORRELATION(Cfun,thisPoFShifts=thisPoFShifts)
    else:
        CfunExt = CreatePoFMatrix(Cfun,thisPoFShifts=thisPoFShifts)

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






#Cfun [ t_src , ism , jsm , t ]
def GetTvarOvdet(Cfunin,thistodtvals,thisPoFShifts=PoFShifts):
    (thistomin,thisdtmin),(thistomax,thisdtmax) = thistodtvals
    thistomin = thistomin-1 #python arrays start from 0, not 1
    thistomax = thistomax-1 #python arrays start from 0, not 1
    Cfun = np.array(Cfunin)
    Cfuntomat,Cfuntodtmat = False,False
    
    for thisto in range(thistomin,thistomax+1):
        for thisdt in range(thisdtmin,thisdtmax+1):
        if thisPoFShifts == 0:
            if Cfuntomat == False:
                Cfuntomat = np.matrix(Pullflag(Cfun[-1,:,:,thisto-PoFDelta],'Avg'))
                Cfuntodtmat = np.matrix(Pullflag(Cfun[-1,:,:,thisto+thisdt-PoFDelta],'Avg'))
            else:
                Cfuntomat = np.bmat([Cfuntomat,Pullflag(Cfun[-1,:,:,thisto-PoFDelta],'Avg')])
                Cfuntodtmat = np.bmat([Cfuntomat,Pullflag(Cfun[-1,:,:,thisto+thisdt-PoFDelta],'Avg')])
                # Cfuntomat = np.bmat([[Cfuntomat],Pullflag(Cfun[-1,:,:,thisto-PoFDelta],'Avg')])
                # Cfuntodtmat = np.bmat([[Cfuntomat],Pullflag(Cfun[-1,:,:,thisto+thisdt-PoFDelta],'Avg')])
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
    [Emass,LEvec,REvec] = CreateLREves(Cfuntoout,Cfuntodtout,thisdt)
    return Emass,LEvec,REvec

# CM
# #Cfuns2pt [ ism , jsm , ip , it ]  = bootstrap1 class (.Avg, .Std, .values, .nboot)
# PoF
# #Cfuns2pt [ tsource, ism , jsm , ip , it ]  = bootstrap1 class (.Avg, .Std, .values, .nboot)

# #CMCfun2pt [ istate , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
# R/L Evecs [ ip , istate , ival ]
# Emass [ ip , istate ]

def CreateOvdet2ptCfun(Cfuns2pt,todtrange,thisMomList):
    start = time.time()
    Cfuns2pt = np.array(Cfuns2pt)
    CMCfun2pt = []
    Emass,LEvec,REvec = [],[],[]
    for ip,thisp in enumerate(thisMomList):
        Emasshold,LEvechold,REvechold = GetTvarOvdet(Cfuns2pt[:,:,:,ip],todtvals)
        Emass.append(Emasshold)
        LEvec.append(LEvechold)
        REvec.append(REvechold)
    LEvec,REvec = SignEvec(np.array(LEvec),np.array(REvec))
    return [np.rollaxis(np.array(CMCfun2pt),1),LEvec,REvec,Emass]




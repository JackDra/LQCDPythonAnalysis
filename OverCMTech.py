#!/usr/bin/env python

import numpy as np
from scipy.linalg import eigvals,inv,sqrtm
from scipy.linalg import eigh,eig
# from numpy.linalg import eigh,eig
from Params import *
from copy import deepcopy
import operator as opp
import time
from CMSTech import *
from OvDetEigen import OverdetEigen

## Some main function call are taken from CMSTech

## Emass = 101 is invalid eigenvalue out
## Emass = 202 is b matrix not positive definite, removing higher smearings.



# Mat*Evec = Mat2 * Evals * Evec
def CreateLREvesOvDet(Cfunto,Cfuntodt,thisdt):    
    if VarMethodMethod == 'AxBxlSolve' or VarMethodMethod == 'Regular':
        Evals,REvec,Dump = OverdetEigen(Cfuntodt,Cfunto,OverDetIter)
        # Evals2,REvec2 = eig(Cfuntodt[:Cfuntodt.shape[1],:],b=Cfunto[:Cfunto.shape[1],:])
        # print 'Overdet'
        # print Cfunto
        # print Cfuntodt
        # print -np.log(Evals)/thisdt
        # print REvec
        # print 'reg'
        # print Cfunto[:Cfunto.shape[1],:]
        # print Cfuntodt[:Cfuntodt.shape[1],:]
        # print -np.log(Evals2)/thisdt
        # print REvec2
        LEvec = REvec ## Assume both left and right eigenvalues are equal (maybe think about turning Symmetrise on?)
    else:
        raise LookupError('VarMethodMethod Must be AxBxlSolve or Regular for Over Detumined eigenvalue problem: ' + VarMethodMethod)
    # print 'Evals:'
    # print Evals
    return sortEvec(Evals,LEvec,REvec,thisdt)



#Cfun [ t_src , ism , jsm , t ]
def GetTvarOvdet(Cfunin,thistovals,thisdt,thisPoFShifts=PoFShifts):
    (thistomin,thistomax) = thistovals
    thistomin = thistomin-1 #python arrays start from 0, not 1
    thistomax = thistomax-1 #python arrays start from 0, not 1
    Cfun = np.array(Cfunin)
    Cfuntomat,Cfuntodtmat = [],[]
    
    if thisPoFShifts > 0:  raise IOError('Code not implemented for thisPoFShift > 0 yet')
    for thisto in xrange(thistomin,thistomax+1):
        thistomat = np.matrix(Pullflag(Cfun[0,:,:,thisto-PoFDelta],'Avg'))
        thistodtmat = np.matrix(Pullflag(Cfun[0,:,:,thisto+thisdt-PoFDelta],'Avg'))
        if len(Cfuntomat) == 0:
            if VarMethodMethod == 'AxBxlSolve':
                Cfuntomat = thistomat
                Cfuntodtmat = thistodtmat
            elif VarMethodMethod == 'Regular':
                Cfuntomat = thistomat* inv(thistodtmat)
                Cfuntodtmat = np.eye(*thistomat.shape)
        else:
            if VarMethodMethod == 'AxBxlSolve':
                Cfuntomat = np.bmat([[Cfuntomat],[thistomat]])
                Cfuntodtmat = np.bmat([[Cfuntodtmat],[thistodtmat]])
            elif VarMethodMethod == 'Regular':
                Cfuntomat = np.bmat([[Cfuntomat],[thistomat* inv(thistodtmat)]])
                Cfuntodtmat = np.bmat([[Cfuntodtmat],[np.eye(*thistomat.shape)]])
    # Cfuntoout,Cfuntodtout = CreatePoFMatrixtodt(Cfuntomat,Cfuntodtmat,thisPoFShifts=PoFShifts)
    [Emass,LEvec,REvec] = CreateLREvesOvDet(Cfuntomat,Cfuntodtmat,thisdt)
    return Emass,LEvec,REvec

# CM
# #Cfuns2pt [ ism , jsm , ip , it ]  = bootstrap1 class (.Avg, .Std, .values, .nboot)
# PoF
# #Cfuns2pt [ tsource, ism , jsm , ip , it ]  = bootstrap1 class (.Avg, .Std, .values, .nboot)

# #CMCfun2pt [ istate , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
# R/L Evecs [ ip , istate , ival ]
# Emass [ ip , istate ]

def CreateOvdet2ptCfun(Cfuns2pt,torange,thisdt,thisMomList,DoPoF=True):
    Cfuns2pt = np.array(Cfuns2pt)
    CMCfun2pt = []
    Emass,LEvec,REvec = [],[],[]
    for ip,thisp in enumerate(thisMomList):
        Emasshold,LEvechold,REvechold = GetTvarOvdet(Cfuns2pt[:,:,:,ip],torange,thisdt)
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
    return [np.rollaxis(np.array(CMCfun2pt),1),LEvec,REvec,Emass]




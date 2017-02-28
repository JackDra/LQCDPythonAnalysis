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
from CMSTech import *
from OvDetEigen import OverdetEigen

## Some main function call are taken from CMSTech

## Emass = 101 is invalid eigenvalue out
## Emass = 202 is b matrix not positive definite, removing higher smearings.



# Mat*Evec = Mat2 * Evals * Evec
def CreateLREvesOvDet(Cfunto,Cfuntodt,thisdt):    
    if VarMethodMethod == 'AxBxlSolve':
        ## got dimensions wrong, take transpose
        # Evals,REvec,Dump = OverdetEigen(Cfunto,Cfuntodt,OverDetIter)
        print np.matrix(Cfunto.T)
        print np.matrix(Cfuntodt.T)
        Evals,REvec,Dump = OverdetEigen(Cfunto.T,Cfuntodt.T,OverDetIter)
        LEvec = REvec ## Assume both left and right eigenvalues are equal (maybe think about turning Symmetrise on?)
    else:
        raise LookupError('VarMethodMethod Must be AxBxlSolve for Over Detumined eigenvalue problem: ' + VarMethodMethod)
    print 'Evals:'
    print Evals
    return sortEvec(Evals,LEvec,REvec,thisdt)



#Cfun [ t_src , ism , jsm , t ]
def GetTvarOvdet(Cfunin,thistovals,thisdt,thisPoFShifts=PoFShifts):
    (thistomin,thistomax) = thistovals
    thistomin = thistomin-1 #python arrays start from 0, not 1
    thistomax = thistomax-1 #python arrays start from 0, not 1
    Cfun = np.array(Cfunin)
    Cfuntomat,Cfuntodtmat = [],[]
    
    if thisPoFShifts > 0:  raise IOError('Code not implemented for thisPoFShift > 0 yet')
    for thisto in range(thistomin,thistomax+1):
        if len(Cfuntomat) == 0:
            Cfuntomat = np.matrix(Pullflag(Cfun[0,:,:,thisto-PoFDelta],'Avg'))
            Cfuntodtmat = np.matrix(Pullflag(Cfun[0,:,:,thisto+thisdt-PoFDelta],'Avg'))
        else:
            Cfuntomat = np.bmat([Cfuntomat,Pullflag(Cfun[0,:,:,thisto-PoFDelta],'Avg')])
            Cfuntodtmat = np.bmat([Cfuntodtmat,Pullflag(Cfun[0,:,:,thisto+thisdt-PoFDelta],'Avg')])
            # Cfuntomat = np.bmat([[Cfuntomat],Pullflag(Cfun[-1,:,:,thisto-PoFDelta],'Avg')])
            # Cfuntodtmat = np.bmat([[Cfuntomat],Pullflag(Cfun[-1,:,:,thisto+thisdt-PoFDelta],'Avg')])
                
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




#!/usr/bin/env python

import numpy as np
from scipy.linalg import eigvals,inv,sqrtm
from scipy.linalg import eigh,eig
from scipy.linalg import qr,qz
from scipy.linalg import solve

# from numpy.linalg import eigh,eig
from Params import *
from copy import deepcopy
import operator as opp
from ReadTxt import ReadLREM
import time

## qr factorisation of (B,A)
def QRAppend(A,B):
    return qr(np.bmat([B,A]),mode='r')[0]


## R = [ R11 R12 ]
##     [  0  R22 ]
def GetRs(A,B):
    Rtild = QRAppend(A,B)
    print 'Rtild'
    print Rtild
    Rtop,Rbot = np.split(Rtild,2)
    print 
    print 'Rtop'
    print Rtop
    print 
    print 'Rbot'
    print Rbot
    R11,R12 = np.split(Rtop,2,1)
    R21,R22 = np.split(Rbot,2,1)
    print 'R11'
    print R11
    print 
    print 'R12'
    print R12
    print 
    print 'R21'
    print R21
    print 
    print 'R22'
    print R22
    print 

    if np.any(R21!=0.):
        print 'warning, did not QR factorise properly'
    return np.matrix(R11),np.matrix(R12),np.matrix(R22)
    

def GetQZ(R11,R12):
    R,R0,Q,Z = qz(R12,R11)
    print 'R'
    print R
    print 
    print 'R0'
    print R0
    print 
    print 'Q'
    print Q
    print 
    print 'Z'
    print Z
    print 
    return np.matrix(R0),np.matrix(R),np.matrix(Z)


## Ax = Bx l
def OverdetEigen(A,B,Niter):
    R11,R12,R22 = GetRs(np.matrix(A),np.matrix(B))
    R,R0,Z = GetQZ(R11,R12)
    E = R22*Z
    C = E.getH()*E
    Eval,Evec = eig(R0,b=R)
    EvalOut,EvecOut = Eval,Evec
    print 'Initial eval/evec'
    iterlist = []
    for istate,(sEval,sEvec) in enumerate(zip(Eval,Evec)):        
        smEvec = np.matrix(sEvec)
        print 'eval: ', sEval
        print 'evec: ', sEvec
        iiter = 0
        for iiter in range(Niter):
            Rmin1 = R0-(sEval*R)
            ydata = np.bmat([[-C*smEvec.T],[smEvec*smEvec.H]])
            xdata = np.bmat([[R.H*R,smEvec.T],[smEvec.H.T,np.matrix([[0]])]])
            try:
                vtildk = solve(xdata,ydata)[:,0]
            except:
                iterlist.append(iiter)
                continue
            vtild,kappa = np.array(vtildk[:-1]),vtildk[-1]
            EvecOut[istate] = vtild/np.sqrt(vtild.dot(vtild))
            vi = np.matrix(EvecOut[istate])
            nexteval = EvalOut[istate] + (((R*vi.T).H *(Rmin1*vi.T))/((R*vi.T).H *(R*vi.T)))
            EvalOut[istate] = nexteval[0,0]
            # EvalOut[istate] = EvalOut[istate] + (((R*vi.T).H *(Rmin1*vi.T))/complex(((R*vi.T).H *(R*vi.T))))
        iterlist.append(iiter)
    return EvalOut,EvecOut,iterlist

            
            

            

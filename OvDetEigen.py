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


rollindex=1

## qr factorisation of (B,A)
def QRAppend(A,B):
    return qr(np.bmat([B,A]),mode='r')[0]


## R = [ R11 R12 ]
##     [  0  R22 ]
def GetRs(A,B):
    mdim,ndim = A.shape
    Rtild = QRAppend(A,B)
    # print 'Rtild'
    # print Rtild
    Rtop,Rbot = Rtild[:ndim],Rtild[ndim:]
    # print 
    # print 'Rtop'
    # print Rtop
    # print 
    # print 'Rbot'
    # print Rbot
    R11,R12 = Rtop[:,:ndim],Rtop[:,ndim:]
    R21,R22 = Rbot[:,:ndim],Rbot[:,ndim:]
    # print
    # print 'R11'
    # print R11
    # print 
    # print 'R12'
    # print R12
    # print 
    # print 'R21'
    # print R21
    # print 
    # print 'R22'
    # print R22
    # print 

    if np.any(R21!=0.):
        print 'warning, did not QR factorise properly'
    return np.matrix(R11),np.matrix(R12),np.matrix(R22)
    

def GetQZ(R11,R12):
    R0,R,Q,Z = qz(R12,R11)
    # print 'R'
    # print R
    # print 
    # print 'R0'
    # print R0
    # print 
    # print 'Q'
    # print Q
    # print 
    # print 'Z'
    # print Z
    # print 
    return np.matrix(R0),np.matrix(R),np.matrix(Z)


## Ax = Bx l
def OverdetEigen(A,B,Niter,Tol=OverdetTol):
    R11,R12,R22 = GetRs(np.matrix(A),np.matrix(B))
    R0,R,Z = GetQZ(R11,R12)
    E = R22*Z
    C = E.getH()*E
    # Eval,Evec = eig(R12,b=R11)
    Eval,Evec = eig(R0,b=R)
    Evec = np.rollaxis(Evec,rollindex)
    # print
    # print 'Z*Evec'
    # print Z*np.matrix(Evec)[:,0]
    # EvalOut,EvecOut = Eval,Evec
    # print
    # print 'Initial eval/evec'
    # print Eval
    # print Evec
    iterlist = []
    # for istate,(sEval,sEvec) in enumerate(zip(Eval,Evec)):        
    EvalOut = []
    EvecOut = []
    for istate,(sEval,sEvec) in enumerate(zip(Eval,Evec)): 
        smEvec = np.matrix(sEvec)
        iiter = 0
        # print 'eval #'+str(istate+1)+': ', sEval
        # print 'evec #'+str(istate+1)+': ', sEvec
        # print 
        vi = deepcopy(smEvec)
        nexteval = sEval
        EvalOut.append(sEval)
        for iiter in xrange(Niter):
            Rmin1 = R0-(sEval*R)
            # minloc = np.argmin(Rmin1)
            # minloc = minloc/(Rmin1.shape[0]-1),minloc%(Rmin1.shape[1]-1)
            # Rmin1[minloc] = Rmin1[minloc] + np.max(Rmin1) ##create S matrix to remove singularities
            # Rmin1 = R12-(sEval*R11)
            # print
            # print 'C'
            # print C
            ydata = np.bmat([[-C*smEvec.T],[smEvec*smEvec.H]])
            xdata = np.bmat([[Rmin1.H*Rmin1,smEvec.T],[smEvec.H.T,np.matrix([[0]])]])
            try:
                # print 
                # print 'solve data'
                # print 'x'
                # print xdata[0]
                # print xdata[1]
                # print xdata[2]
                # print xdata[3]
                # print xdata[4]
                # print 'y'
                # print ydata
                vtildk = solve(xdata,ydata)[:,0]
                # print
                # print 'eh?'
                # print smEvec
                # print vtildk
                # print
            except:
                print
                print 'istate',istate,' iters =',iiter , ' Failed inversion'
                iterlist.append(iiter)
                break
            vtild,kappa = np.array(vtildk[:-1]),vtildk[-1]
            vi = np.matrix(vtild/np.linalg.norm(vtild))
            # nexteval = EvalOut[istate] + (((R*vi.T).H *(Rmin1*vi.T))/((R*vi.T).H *(R*vi.T)))
            preveval = deepcopy(nexteval)
            nexteval = nexteval + (((R*vi.T).H *(Rmin1*vi.T))/((R*vi.T).H *(R*vi.T)))
            nexteval = nexteval[0,0]
            # print 'Nexteval',nexteval
            EvalOut[istate] = nexteval
            # print 'Evals'
            # print preveval
            # print nexteval
            # print np.abs(preveval-nexteval)
            # print Tol
            # print 
            if np.abs(preveval-nexteval) < Tol:
                # print
                # print 'istate',istate,' iters =',iiter 
                break
        # EvalOut[istate] = EvalOut[istate] + (((R*vi.T).H *(Rmin1*vi.T))/complex(((R*vi.T).H *(R*vi.T))))
        # print 
        # print 'BLAH'
        # print Z
        # print vi.T
        # print np.array(Z*vi.T)[:,0]
        EvecOut.append(np.array(Z*vi.T)[:,0])
        iterlist.append(iiter)
    # for ival,ivec in zip(EvalOut,EvecOut):
        
    return np.array(EvalOut),np.rollaxis(np.array(EvecOut),1),iterlist

            
            

            

#!/usr/bin/env python
from array import array
import numpy as np
from Params import *
from MiscFuns import *
from collections import OrderedDict
from copy import deepcopy
import operator

giDiVecSet = ['P4g1D1','P4g2D2','P4g3D3']
##Proton: doublet is up quark, singlet is down quark
##Neutron: doublet is down quark, singlet is up quark
upCharge = 2.0/3.0
downCharge = -1.0/3.0
DSCombs = [('P4I',upCharge,downCharge),('P4g4',upCharge,downCharge),('P4g3g5',-1,1)]
ops = { "+": operator.add, "-": operator.sub } 

## data = [ itsink , ism/istate , igamma , ip , it ]
def MakeUmD(data,gammalist):
    dataout = deepcopy(data)
    newgammalist,pairindex = UDIndex(gammalist)
    for tsink,tsinkdata in enumerate(data):
        for state,statedata in enumerate(tsinkdata):
            for idu,index in enumerate(pairindex):
                dataout[tsink][state].append([])
                dgamma,ugamma = statedata[index[0]],statedata[index[1]]
                for ip,(pd,pu) in enumerate(zip(dgamma,ugamma)):
                    dataout[tsink][state][len(data[0][0])+idu].append([])
                    for it,(td,tu) in enumerate(zip(pd,pu)):
                        dataout[tsink][state][len(data[0][0])+idu][ip].append(td-tu)
                        dataout[tsink][state][len(data[0][0])+idu][ip][it].Stats()
    return dataout,newgammalist

def CreateDS(datadoub,datasing,thisGammaList):
    dataout = Swap3ptSS(datadoub)
    datahold = Swap3ptSS(datasing)
    for igamma, Lcoef,Rcoef in DSCombs:
        if igamma in thisGammaList:
            igloc = thisGammaList.index(igamma)
            dataout[igloc] = Lcoef*dataout[igloc]+Rcoef*datahold[igloc]
    return SwapBack3ptSS(dataout)

def CreategiDi(data3pt,thisGammaList,thisDSList):
    data3ptout = Swap3ptSS(data3pt)
    if len(data3ptout) != 4:
        ##DEBUG##
        print 'dimensions of correlator (ism, jsm, igamma)'
        for ic,idata in enumerate(data3pt):
            for jc,jdata in enumerate(idata):
                print ic , jc , len(jdata)
        raise IOError("Need giDi with i =1,2,3,4 for giDi combination")
    thisGammaListOut = np.array(thisGammaList)
    for iDS in thisDSList:
        g4D4i = thisGammaListOut.tolist().index(iDS+'P4g4D4')
        thisGammaListOut = np.where(thisGammaListOut==iDS+'P4g4D4',iDS+'P4giDi',thisGammaListOut)
        giDii = []
        for ivec,igamma in enumerate(giDiVecSet):
            giDii.append(thisGammaListOut.tolist().index(iDS+igamma))
            data3ptout[g4D4i] -= data3ptout[giDii[ivec]]/3.0
            
            data3ptout = np.delete(data3ptout,giDii,axis=0)
            thisGammaListOut = np.delete(thisGammaListOut,giDii)
    return [SwapBack3ptSS(data3ptout),thisGammaListOut.tolist()]
    

#cfunin [  ism , jsm , igamma , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
#cfunout [ igamma  , ism , jsm , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)

def Swap3ptSS(cfunin):
    return np.rollaxis(np.array(cfunin),2)


#cfunin [ igamma  , ism , jsm , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
#cfunout [  ism , jsm , igamma , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)

def SwapBack3ptSS(cfunin):
    return np.rollaxis(np.array(cfunin),0,3)

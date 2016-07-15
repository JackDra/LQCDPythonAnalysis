#!/usr/bin/env python

import operator

giDiVecSet = ['P4g1D1','P4g2D2','P4g3D3']
##Proton: doublet is up quark, singlet is down quark
##Neutron: doublet is down quark, singlet is up quark
upCharge = 2.0/3.0
downCharge = -1.0/3.0
DSCombs = [('P4I',upCharge,downCharge),('P4g4',upCharge,downCharge),('P4g3g5',-1,1)]
ops = { "+": operator.add, "-": operator.sub } 


def IsoVector(val1,val2):
    return val1 - val2

def Vector(val1,val2):
    return val1 + val2

def FFProton(val1,val2):
    return upCharge*val1 + downCharge*val2

def FFNeutron(val1,val2):
    return downCharge*val1 + upCharge*val2

def F1overF2(F1,F2):
    print F1
    print F2
    if F1 != False and F2 != False:
        return F1/F2
    else:
        return -1337.0

CombFFFunsDict = {'F1divF2':F1overF2}

CombFunsDict = {'IsoVector':IsoVector,
                'Vector':Vector,
                'Proton':FFProton,
                'Neutron':FFNeutron}

CombList = CombFunsDict.keys()
CombFFList = CombFFFunsDict.keys()

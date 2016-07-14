#!/usr/bin/env python


def IsoVector(val1,val2):
    return val1 - val2

def Vector(val1,val2):
    return val1 + val2

def FFProton(val1,val2):
    return upCharge*val1 + downCharge*val2

def FFNeutron(val1,val2):
    return downCharge*val1 + upCharge*val2

CombFunsDict = {'IsoVector':IsoVector,
                'Vector':Vector,
                'Proton':FFProton,
                'Neutron':FFNeutron}
CombList = CombFunsDict.keys()

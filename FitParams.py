#!/usr/bin/env python

import numpy as np
import os, errno
from Params import *
import itertools

def CheckCut(thisset,thisdict):
    for ikey,icut in thisdict.iteritems():
        if ikey in thisset and icut in thisset:
            return True
    return False

def GetCut(thisset,thisdict):
    for ikey,icut in thisdict.iteritems():
        if ikey in thisset:
            return icut
    return False


MassDtList = range(1,7)
# MassDtList = [2]
FitCutMin,FitCutMax = 3,9
FitCutList = range(FitCutMin,FitCutMax+1)
FitCutArgs = ['cut'+str(ic) for ic in FitCutList]
##TryFits.py
FitCutPicked =  {'tsink26sm32':'cut6',
                 'tsink26'+PickedStateStr+'REvec':'cut4',
                 'tsink29sm32':'cut6',
                 'tsink29sm64':'cut5',
                 'tsink29sm128':'cut5',
                 'tsink29'+PickedStateStr+'CM':'cut5',
                 'tsink32sm32':'cut6',
                 'tsink32'+PickedStateStr+'REvec':'cut6',
                 'tsink35sm32':'cut6',
                 'tsink38sm32':'cut6',
                 'tsink26'+PickedStateStr+'PoF'+str(PoFShifts):'cut4',
                 'tsink27'+PickedStateStr+'PoF'+str(PoFShifts):'cut4'}
##



SumMeth3ptCuts = range(4)
##picked for plotting
SumCutList = ['cut'+str(isum) for isum in SumMeth3ptCuts]
SumCutPar = 'cut3'
SumFitRList = ['fit sl 0-4','fit sl 1-4','fit sl 2-4']
SumFitRPicked = 'fitr1-4'

SumFitRListFlags = [ifitr.replace(' sl ','r') for ifitr in SumFitRList]
SumFitRSrip = [ifitr.replace('fit sl ','') for ifitr in SumFitRList]

def CreateFitList(thisTwoMin,thisTwoMinMax,thisTwoMax,thisTwoMaxMax,this3ptCutMin,this3ptCutMax):
    this2ptList = []
    for fmin in range(thisTwoMin,thisTwoMinMax+1):
        for fmax in range(thisTwoMax,thisTwoMaxMax+1):
            if fmin+4 < fmax:
                this2ptList.append((fmin,fmax))
    return this2ptList

OSF3ptCutMin,OSF3ptCutMax = 3,7
OSF2ptMinStart,OSF2ptMinEnd,OSF2ptMaxStart,OSF2ptMaxEnd = 19,29,32,35
OSF3ptCutList = range(OSF3ptCutMin,OSF3ptCutMax+1)
OneStateParList = {'C2':['Am','m0'] , 'C3':['B00']}
#picked for plotting
OSFCutList = ['cut4','cut5','cut6']
OSFCutPicked = 'cut5'
OSFfitvals = {'sm32': [28,35], 'sm64': [28,35] , 'sm128':[28,35],
              PickedStateStr+'REvec':[25,35],PickedStateStr+'PoF'+str(PoFShifts):[19,32]}
    
OSFfitr = {'sm32':str(OSFfitvals['sm32'][0])+'-'+str(OSFfitvals['sm32'][1]),
           'sm64':str(OSFfitvals['sm64'][0])+'-'+str(OSFfitvals['sm64'][1]),
           'sm128':str(OSFfitvals['sm128'][0])+'-'+str(OSFfitvals['sm128'][1]),
           PickedStateStr+'REvec':str(OSFfitvals[PickedStateStr+'REvec'][0])+'-'+str(OSFfitvals[PickedStateStr+'REvec'][1]),
           PickedStateStr+'PoF'+str(PoFShifts):str(OSFfitvals[PickedStateStr+'PoF'+str(PoFShifts)][0])+'-'+str(OSFfitvals[PickedStateStr+'PoF'+str(PoFShifts)][1])}
for itvar in DefTvarPicked:
    OSFfitvals[PickedStateStr+itvar] = [24,35]
    OSFfitr[PickedStateStr+itvar] = str(OSFfitvals[PickedStateStr+itvar][0])+'-'+str(OSFfitvals[PickedStateStr+itvar][1])

TSF3ptCutMin,TSF3ptCutMax = 2,5
TSF2ptMinStart,TSF2ptMinEnd,TSF2ptMaxStart,TSF2ptMaxEnd = 18,25,29,35
TSF3ptCutList = range(TSF3ptCutMin,TSF3ptCutMax+1)
TwoStateParList = {'C2':['Am','Amp','m0','Dm'] , 'C3':['B00','B10','B01','B11']}
## picked for plotting
TSFCutList = ['cut3','cut4','cut5']
TSFCutPicked = 'cut4'
##maybe make set depenant like above
TSFfitvals = [20,35]
TSFfitr = str(TSFfitvals[0])+'-'+str(TSFfitvals[1])


StateParList = {'Two':TwoStateParList,'One':OneStateParList}

MaxIters = 10000
def FitDefGuess(Fun,Len=1):
    if Fun.__name__ == 'FormFactorO1':
        return [1]
    if Fun.__name__ == 'FormFactorO2':
        return [1,1]
    if Fun.__name__ == 'FormFactorO3':
        return [1,1,1]
    if Fun.__name__ == 'FormFactorO':
        return [1]*Len
    if Fun.__name__ == 'ConstantFitFun':
        return [1]
    elif Fun.__name__ == 'LinearFitFun':
        return [1,1]
    elif Fun.__name__ == 'C2TwoStateFitFun':
        return [3.1861104305e-06,np.log(.45),4,np.log(0.3385347963)]
    elif Fun.__name__ == 'C2OneStateFitFun':
        return [1,np.log(.45)]
        # return [2.7174048910e-07,np.log(.45)]
    elif Fun.__name__ == 'C3OneStateFitFun':
        return [0]
    elif Fun.__name__ == 'C2TwoStateFitFunCM':
        output = []
        for i in range(Len):
            output.append(3.1861104305e-06)
            output.append(4)
        return output+[np.log(.45),np.log(0.3385347963)]
    elif Fun.__name__ == 'C3MomTwoStateFitFun':
        return [0,0,0,0]
    elif Fun.__name__ == 'C3TwoStateFitFun2par':
        return [0,0]
    elif Fun.__name__ == 'C3TwoStateFitFun3par':
        return [0,0,0]
    elif Fun.__name__ == 'TestTwoVarFitFun':
        return [1,1]
    else:
        return [1]*100

    




##GraphData.py
FitMassPicked = {'sm32':(15,18),'sm64':(14,18),'sm128':(14,18)} ## GraphData.py
for itvar in DefTvarPicked:
    FitMassPicked[PickedStateStr+itvar] = (10,19)
PickedTwoMax = 29

def WritePars(val,fileflag):
    fout = open('./'+fileflag+'.par','w')
    for ival in val:
        fout.write(str(ival)+'\n')
    fout.close()

def ReadPars(fileflag):
    fin = open('./'+fileflag+'.par','r')
    valout = []
    for line in fin:
        valout.append(float(line.replace('\n','')))
    return valout
##
##TryMassFits.py
FitMassMin = range(10,16)
FitMassMax = range(15,21)
FitMassList = []
for fmax,fmin in itertools.product(FitMassMax,FitMassMin):
    if fmin+2 < fmax:
        FitMassList.append((fmin,fmax))
##


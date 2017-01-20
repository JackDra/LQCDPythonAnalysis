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


MassDtList = range(1,5)
DoFList = range(3,17)
ChiThreshold = 0.001
# MassDtList = [2]
FitCutMin,FitCutMax = 1,9
FitCutList,FitCutArgs = [],[]
for icut in range(FitCutMin,FitCutMax+1):
    for jcut in range(FitCutMin,FitCutMax+1):
        if icut + jcut >= 10: continue
        FitCutList.append([icut,jcut])
        FitCutArgs.append('cut'+str(icut)+'-'+str(jcut))
##TryFits.py
FitCutPicked =  {'tsink26sm32':'cut6-6',
                 'tsink26'+PickedStateStr+'REvec':'cut4-4',
                 'tsink29sm32':'cut5-5',
                 'tsink29sm64':'cut5-5',
                 'tsink29sm128':'cut5-5',
                 'tsink29'+PickedStateStr+'CM':'cut5-5',
                 'tsink29'+PickedStateStr+'REvec':'cut5-5',
                 'tsink32sm32':'cut6-6',
                 'tsink32'+PickedStateStr+'REvec':'cut6-6',
                 'tsink35sm32':'cut6-6',
                 'tsink38sm32':'cut6-6',
                 'tsink13'+PickedStateStr+'PoF'+str(PoFShifts):'cut1-5',
                 'tsink26'+PickedStateStr+'PoF'+str(PoFShifts):'cut3-3',
                 'tsink27'+PickedStateStr+'PoF'+str(PoFShifts):'cut3-3'}
##



SumMeth3ptCuts = range(4)
##picked for plotting
SumCutList = ['cut'+str(isum) for isum in SumMeth3ptCuts]
SumCutPar = 'cut3'
# SumFitRList = ['fit sl 0-4','fit sl 1-4','fit sl 2-4']
SumFitRList = ['fit sl 0-4','fit sl 2-4']
SumFitRPicked = 'fitr0-4'

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
# OSF2ptMinStart,OSF2ptMinEnd,OSF2ptMaxStart,OSF2ptMaxEnd = 1,11,8,20
OSF2ptMinStart,OSF2ptMinEnd,OSF2ptMaxStart,OSF2ptMaxEnd = 2,10,11,21
FitMaxList = range(11,21)
OSF3ptCutList = range(OSF3ptCutMin,OSF3ptCutMax+1)
OneStateParList = {'C2':['Am','m0'] , 'C3':['B00']}
#picked for plotting
OSFCutList = ['cut4','cut5','cut6']
OSFCutPicked = 'cut5'
OSFfitvals = {'sm16': [9,21],'sm32': [8,19], 'sm64': [7,21] , 'sm128':[28,35],
              PickedStateStr+'REvec':[25,35],PickedStateStr+'PoF'+str(PoFShifts):[4,9],PickedStateStr+'CM':[5,21]}

OSFfitvals100 = {'sm16': [9,21],'sm32': [8,19], 'sm64': [7,21] , 'sm128':[28,35],
                 PickedStateStr+'REvec':[25,35],PickedStateStr+'PoF'+str(PoFShifts):[4,9],PickedStateStr+'CM':[5,21]}

OSFfitvals110 = {'sm16': [9,21],'sm32': [8,19], 'sm64': [7,21] , 'sm128':[28,35],
                 PickedStateStr+'REvec':[25,35],PickedStateStr+'PoF'+str(PoFShifts):[4,9],PickedStateStr+'CM':[5,21]}

OSFfitvals111 = {'sm16': [9,21],'sm32': [8,19], 'sm64': [7,21] , 'sm128':[28,35],
                 PickedStateStr+'REvec':[25,35],PickedStateStr+'PoF'+str(PoFShifts):[4,9],PickedStateStr+'CM':[5,15]}

OSFfitvals200 = {'sm16': [9,21],'sm32': [8,19], 'sm64': [7,21] , 'sm128':[28,35],
                 PickedStateStr+'REvec':[25,35],PickedStateStr+'PoF'+str(PoFShifts):[4,9],PickedStateStr+'CM':[5,15]}


OSFfitr = {'sm32':str(OSFfitvals['sm32'][0])+'-'+str(OSFfitvals['sm32'][1]),
           'sm64':str(OSFfitvals['sm64'][0])+'-'+str(OSFfitvals['sm64'][1]),
           'sm16':str(OSFfitvals['sm16'][0])+'-'+str(OSFfitvals['sm16'][1]),
           PickedStateStr+'REvec':str(OSFfitvals[PickedStateStr+'REvec'][0])+'-'+str(OSFfitvals[PickedStateStr+'REvec'][1]),
           PickedStateStr+'CM':str(OSFfitvals[PickedStateStr+'CM'][0])+'-'+str(OSFfitvals[PickedStateStr+'CM'][1]),
           PickedStateStr+'PoF'+str(PoFShifts):str(OSFfitvals[PickedStateStr+'PoF'+str(PoFShifts)][0])+'-'+str(OSFfitvals[PickedStateStr+'PoF'+str(PoFShifts)][1])}

OSFfitr100 = {'sm32':str(OSFfitvals100['sm32'][0])+'-'+str(OSFfitvals100['sm32'][1]),
              'sm64':str(OSFfitvals100['sm64'][0])+'-'+str(OSFfitvals100['sm64'][1]),
              'sm16':str(OSFfitvals100['sm16'][0])+'-'+str(OSFfitvals100['sm16'][1]),
              PickedStateStr+'REvec':str(OSFfitvals100[PickedStateStr+'REvec'][0])+'-'+str(OSFfitvals100[PickedStateStr+'REvec'][1]),
              PickedStateStr+'CM':str(OSFfitvals100[PickedStateStr+'CM'][0])+'-'+str(OSFfitvals100[PickedStateStr+'CM'][1]),
              PickedStateStr+'PoF'+str(PoFShifts):str(OSFfitvals100[PickedStateStr+'PoF'+str(PoFShifts)][0])+'-'+str(OSFfitvals100[PickedStateStr+'PoF'+str(PoFShifts)][1])}

OSFfitr110 = {'sm32':str(OSFfitvals110['sm32'][0])+'-'+str(OSFfitvals110['sm32'][1]),
              'sm64':str(OSFfitvals110['sm64'][0])+'-'+str(OSFfitvals110['sm64'][1]),
              'sm16':str(OSFfitvals110['sm16'][0])+'-'+str(OSFfitvals110['sm16'][1]),
              PickedStateStr+'REvec':str(OSFfitvals110[PickedStateStr+'REvec'][0])+'-'+str(OSFfitvals110[PickedStateStr+'REvec'][1]),
              PickedStateStr+'CM':str(OSFfitvals110[PickedStateStr+'CM'][0])+'-'+str(OSFfitvals110[PickedStateStr+'CM'][1]),
              PickedStateStr+'PoF'+str(PoFShifts):str(OSFfitvals110[PickedStateStr+'PoF'+str(PoFShifts)][0])+'-'+str(OSFfitvals110[PickedStateStr+'PoF'+str(PoFShifts)][1])}

OSFfitr111 = {'sm32':str(OSFfitvals111['sm32'][0])+'-'+str(OSFfitvals111['sm32'][1]),
              'sm64':str(OSFfitvals111['sm64'][0])+'-'+str(OSFfitvals111['sm64'][1]),
              'sm16':str(OSFfitvals111['sm16'][0])+'-'+str(OSFfitvals111['sm16'][1]),
              PickedStateStr+'REvec':str(OSFfitvals111[PickedStateStr+'REvec'][0])+'-'+str(OSFfitvals111[PickedStateStr+'REvec'][1]),
              PickedStateStr+'CM':str(OSFfitvals111[PickedStateStr+'CM'][0])+'-'+str(OSFfitvals111[PickedStateStr+'CM'][1]),
              PickedStateStr+'PoF'+str(PoFShifts):str(OSFfitvals111[PickedStateStr+'PoF'+str(PoFShifts)][0])+'-'+str(OSFfitvals111[PickedStateStr+'PoF'+str(PoFShifts)][1])}

OSFfitr200 = {'sm32':str(OSFfitvals200['sm32'][0])+'-'+str(OSFfitvals200['sm32'][1]),
              'sm64':str(OSFfitvals200['sm64'][0])+'-'+str(OSFfitvals200['sm64'][1]),
              'sm16':str(OSFfitvals200['sm16'][0])+'-'+str(OSFfitvals200['sm16'][1]),
              PickedStateStr+'REvec':str(OSFfitvals200[PickedStateStr+'REvec'][0])+'-'+str(OSFfitvals200[PickedStateStr+'REvec'][1]),
              PickedStateStr+'CM':str(OSFfitvals200[PickedStateStr+'CM'][0])+'-'+str(OSFfitvals200[PickedStateStr+'CM'][1]),
              PickedStateStr+'PoF'+str(PoFShifts):str(OSFfitvals200[PickedStateStr+'PoF'+str(PoFShifts)][0])+'-'+str(OSFfitvals200[PickedStateStr+'PoF'+str(PoFShifts)][1])}

OSFfitrMom = {'q = 0 0 0':OSFfitr,
              'q = 1 0 0':OSFfitr100,
              'q = 1 1 0':OSFfitr110,
              'q = 1 1 1':OSFfitr111,
              'q = 2 0 0':OSFfitr200}
    

# for itvar in DefTvarPicked[0]:
#     OSFfitvals[PickedStateStr+itvar] = [25,35]
#     OSFfitr[PickedStateStr+itvar] = str(OSFfitvals[PickedStateStr+itvar][0])+'-'+str(OSFfitvals[PickedStateStr+itvar][1])

TSF3ptCutMin,TSF3ptCutMax = 2,5
TSF2ptMinStart,TSF2ptMinEnd,TSF2ptMaxStart,TSF2ptMaxEnd = 1,7,7,20
# TSF2ptMinStart,TSF2ptMinEnd,TSF2ptMaxStart,TSF2ptMaxEnd = 4,5,10,11
TSF3ptCutList = range(TSF3ptCutMin,TSF3ptCutMax+1)
TwoStateParList = {'C2':['Am','Amp','m0','Dm'] , 'C3':['B00','B10','B01','B11']}
## picked for plotting
TSFCutList = ['cut3','cut4','cut5']
TSFCutPicked = 'cut4'
##maybe make set depenant like above
TSFfitvals = [3,20]
TSFfitr = str(TSFfitvals[0])+'-'+str(TSFfitvals[1])


StateParList = {'Two':TwoStateParList,'One':OneStateParList}

MaxIters = 10
def FitDefGuess(Fun,Len=1):
    if Fun.__name__ == 'DPfitfun':
        return [2.7,1]
        # return [-1,-1.6]
    if Fun.__name__ == 'DPfitfunOnePar':
        return [1.6]
        # return [-1,-1.6]
    if Fun.__name__ == 'DPfitfun2':
        return [2.7,1]
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
        return [3.1861104305e-06,np.log(.6),4,np.log(0.3385347963)]
    elif Fun.__name__ == 'C2OneStateFitFun':
        return [1.3422466805e-10,np.log(.6)]
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
    elif Fun.__name__ == 'OneOnRootNFitFun':
        return [1]
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


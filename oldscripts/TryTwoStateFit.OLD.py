#!/usr/bin/env python
from array import array
from os import walk
import os
import numpy as np
import sys
from Params import *
from RFCalc import CalcRatioFactorSS
from ReadCMCfuns import ReadSet,ReadList
from CMTech import CreateCMCfuns
from OutputData import PrintC2ToFile, PrintC3ToFile
from CreateCombs import CreategiDi, CreateDS
from Fitting import FitRFSet,TwoStateSetFit
dir = "/raid/jdragos/scratch/cfun/2ndk12090/"
logfile = '/home/accounts/jdragos/scripts/PythonAnalysis/TryTwoStateFit.log'
sys.stdout = open(logfile,'a',0)
sys.stderr = sys.stdout

print '----------------------------------------------------------------------------------'







ReadTSinkList = [29]
ReadSmearList = DefSmearList
FlagList = ['cm']
ReadTSinkList = [29,32,35,38]
# ReadSmearList = ['32']
# FlagList = ['cm','tsink','tsink','tsink']
ReadMomList = [iqTOip(0)]



if ListOrSet == 'ReadList':
    [data2pt,data3pt,filelist] = ReadList(ReadSmearList,ReadMomList,DefProjGammaList,ReadProjDerList,
                                          DefDSList,ReadTSinkList,conflist,FlagList)
elif ListOrSet == 'ReadSet':
    [data2pt,data3pt,filelist] = ReadSet(ReadSmearList,ReadMomList,DefProjGammaList,ReadProjDerList,
                                         DefDSList,ReadTSinkList,dir,FlagList)


## data2pt = [ ism , jsm , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
## data3pt = [ itsink , ism , jsm , igamma , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)

print 'Creating giDi'
data3ptgiDi = []
for tsinkdata in data3pt:
    [data3pthold,CombGammaList] = CreategiDi(tsinkdata,ReadGammaList,DefDSList)
    data3ptgiDi.append(data3pthold)

ncon = np.size(filelist)
print ''
print 'ncon = ' + str(ncon)
print 'nboot = ' + str(nboot)
print 'All Opperators: \n'+'\n'.join(CombGammaList)
print 'All Momenta: '+', '.join(qvecSet[ReadMomList].tolist())
print 'All T Sinks: '+ ', '.join(map(str,ReadTSinkList))
print 'All Smearings: ' + ', '.join(ReadSmearList)







print 'Doing Two State Fits'
data2ptDiag = DiagSmear(data2pt)
data3ptDiag = Diag3ptSmear(data3pt)


print 'Writing Cfuns to file'
PrintC2ToFile(data2ptDiag,ReadSmearList,ReadMomList)
for itsink,tsinkdata in zip(ReadTSinkList,data3ptDiag):
    PrintC3ToFile(tsinkdata,ReadSmearList,ReadMomList,CombGammaList,itsink)

[TSF2pt,TSF3pt,TSF2ptAvg,TSF3ptAvg] = TwoStateSetFit(data2ptDiag,data3ptDiag,[FitMin2pt,FitMax2pt],Fit3ptCut,ReadTSinkList)

for sm,smTF,smTFAvg in zip(ReadSmearList,TSF2pt,TSF2ptAvg):
    print ''
    for thispar,parTF,parTFAvg in zip(TwoStateParList['C2'],smTF,smTFAvg):
        print 'sm{0:21}{1:2}{2:20.10f}{3:20.10f}{4:20.10f}'.format(sm,thispar,parTF.Avg , parTF.Std,parTFAvg[1] )
        

for sm,smTF,smTFAvg in zip(ReadSmearList,TSF3pt,TSF3ptAvg):
    for thistsink,tsinkTF,tsinkTFAvg in zip(ReadTSinkList,smTF,smTFAvg):
        print ''
        for thispar,parTF,parTFAvg in zip(TwoStateParList['C3'],
                                          np.swapaxes(np.array(tsinkTF),0,1),
                                          np.swapaxes(np.array(tsinkTFAvg),0,1)):
            print ''
            for thisgamma, gammaTF,gammaTFAvg in zip(CombGammaList,parTF,parTFAvg):
                print ('sm{0:3}ts{1:2}{2:14}{3:2}{4:20.10f}{5:20.10f}{6:20.10f}'
                       .format(sm,thistsink,thisgamma,thispar, gammaTF.Avg , gammaTF.Std , gammaTFAvg[1]))








# print ''
# print 'Calculating Ratios of Diagonally Smeared Correlators'

# RFrff = []
# for itsink,tsinkdata3pt in zip(ReadTSinkList,data3ptgiDi):
#     [RFrffhold,SqrtFac] = CalcRatioFactorSS(data2pt,tsinkdata3pt,str(itsink),ReadMomList,ReadSmearList)
#     RFrff.append(RFrffhold)







# ## RFrff = [ itsink , igamma , ip , ism , it ] bs1


# print ''
# for itsink , tsinkRF in zip(ReadTSinkList,RFrff):
#     print 'Writing RF data tsink' + str(itsink)
#     PrintSSSetToFile(tsinkRF,ReadSmearList,ReadMomList,CombGammaList,itsink)


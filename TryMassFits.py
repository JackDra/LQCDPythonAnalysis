#!/usr/bin/env python

from Params import *
import numpy as np
from BootTest import BootStrap1
from ReadTxt import ReadSetnp
from SetLists import CreateMassSet,SplitSet,PickTSinkSet
from Fitting import FitMassSet
from OutputData import PrintFitMassSetToFile
from FitParams import *
import itertools as it
import sys

##thisReadList must have form tsink*sm* or tsink*state*to*dt*

thisStateList = ['1']
thisTvarList = ['to18dt2']
thisGammaList = ['twopt']
thisSmearList = DefInterpSmearList
# tolist = ['17','18','19']
# dtlist = ['1','2','3']
# thisTvarlist = []
# for ito in tolist:
#     for idt in dtlist:
#         thisTvarList.append('to'+ito+'dt'+idt)


# thisFitMinList = [19,20,21]
# thisFitMaxList = [23,24,25]

# thisFitList = list(it.product(thisFitMinList,thisFitMaxList))
# for ifit,thisfit in enumerate(thisFitList):
#     if thisfit[0] >= thisfit[1]: del thisFitList[ifit]

SetList = CreateMassSet(thisSmearList,thisStateList,thisTvarList)
print 'SetList:\n' + '\n'.join(SetList)+'\n'


[SetData,thisMomList] = ReadSetnp(thisGammaList,SetList)

print 'ReadMomList:\n' + '\n'.join(thisMomList)+'\n'
SetData = SetData[0]
## SetData = [ ip , istate , it ]



print 'Data Read in from textfiles, Finding Fits.'
FitData,FitDataAvg,FitDataChi =  [],[],[]
sys.stdout.write('fitrange: ')
sys.stdout.flush()
for ifit in FitMassList:
    sys.stdout.write(str(ifit)+' ')
    sys.stdout.flush()
    [FDhold,FDAhold,FDChold] = FitMassSet(SetData,ifit[0],ifit[1])
    FitData.append(FDhold)
    FitDataAvg.append(FDAhold)
    FitDataChi.append(FDChold)
#FitData = [ ifit , ip , istate ] bs1
#FitDataAve = [ ifit , ip , istate ]
#FitDataChi = [ ifit , ip , istate ]
print ''
print 'Printing Fits to File.'
PrintFitMassSetToFile(np.array(FitData),np.array(FitDataChi),thisMomList,
                      SetList,FitMassList)


# FitData = []
# FitDataAve = []
# sys.stdout.write('ifit = ')
# sys.stdout.flush()
# for thisfitmin,thisfitmax in thisFitList:
#     sys.stdout.write(str((thisfitmin,thisfitmax))+' ')
#     sys.stdout.flush()
#     [FDhold,FDAhold] = FitRFSet(SetData,thisfitmin,thisfitmax)
#     FitData.append(FDhold)
#     FitDataAve.append(FDAhold)
# sys.stdout.write('\n')
# sys.stdout.flush()



# for sFD,sFDA, thisstate in zip(np.rollaxis(np.array(FitData),2),
#                                np.rollaxis(np.array(FitDataAve),2),SetList):
#     print thisstate , sFDA[1][1][0],sFD[1][1].Std, sFDA[1][1][1],abs(sFDA[1][1][0]-sFD[1][1].Avg)



#!/usr/bin/env python
from array import array
from os import walk
import os
import numpy as np
import sys
from Params import *
from RFCalc import CalcRatioFactorCM, CalcRatioFactorSS
from ReadCMCfuns import ReadSet,ReadList
from CMTech import CreateCMCfuns
from OutputData import PrintSSSetToFile,PrintCMSetToFile,PrintLREvecMassToFile
from CreateCombs import CreategiDi, CreateDS
from Fitting import FitRFSet

dir = "/raid/jdragos/scratch/cfun/2ndk12090/"
logfile = '/home/accounts/jdragos/scripts/PythonAnalysis/TryCMPick.log'
to = '18'
dt = '2'

sys.stdout = open(logfile,'a',0)
sys.stderr = sys.stdout
print '----------------------------------------------------------------------------------'
# print "qlow to qhigh is: " , ipTOiq(0), ipTOiq(nmom-1)
# print ''

tsink = 29
ReadMomList = [iqTOip(0),iqTOip(1),qvecTOip([-1,0,0])]
ReadProjGammaList = { 'GMA4' : ['I','g4'] ,'GMA3' :['g1g2','g3g5'] }
ReadProjDerList = {'GMA4' :['g4D4','g1D1','g2D2','g3D3']}
# ReadMomList = [iqTOip(0)]
# ReadProjGammaList = { 'GMA4' : ['g4'] }
# ReadProjDerList = {}

ReadGammaList = []
for iGL,GL in ReadProjGammaList.iteritems():
    ReadGammaList += GL
for iGL,GL in ReadProjDerList.iteritems():
    ReadGammaList += GL


### Read CM Set ###
[data2pt,data3pt,filelist] = ReadSet(ReadMomList,ReadProjGammaList,ReadProjDerList,[tsink],dir,'cm')
###   ###


# ### read CM from List ###
# [data2pt,data3pt,filelist] = ReadList(ReadMomList,ReadProjGammaList,ReadProjDerList,[tsink],conflist,'cm')
# ###   ###

data3pt = data3pt[0]
## data2pt = [ ism , jsm , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
## data3pt = [ DS , ism , jsm , igamma , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)

ncon = np.size(filelist)
print ''
print 'ncon = ' + str(ncon)
print 'nboot = ' + str(nboot)


print ''
print 'Combining Doub and Sing (unspecified remains as doublet)'
data3pt = CreateDS(data3pt[0],data3pt[1],ReadGammaList)

print ''
print 'Creating giDi'
[data3ptgiDi,CombGammaList] = CreategiDi(data3pt,ReadGammaList)
# data3ptgiDi = data3pt
# CombGammaList = ReadGammaList
 

print ['All Opperators: '] + CombGammaList
print ['All Momenta: '] + qvecSet[ReadMomList].tolist()

print ''
print 'Creating Correlation Matrix Correlators'
[CMdata2pt,CMdata3ptgiDi,LREvec,Emass] = CreateCMCfuns(data3ptgiDi,data2pt,to,dt)
## CMdata2pt [ ip , istate , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
## CMdata3ptgiDi  [ ip , igamma , istate , it] = bootstrap1 class (.Avg, .Std, .values, .nboot)

print ''
print 'Calculating Ratios of Diagonally Smeared Correlators'

[RFrff,SqrtFac] = CalcRatioFactorSS(data2pt,data3ptgiDi,str(tsink))

print ''
print 'Calculating Ratios of CM Correlators'

[RFrCM,SqrtFac] = CalcRatioFactorCM(CMdata2pt,CMdata3ptgiDi,str(tsink))

## RFrff = [ igamma , ip , ism , it ] bs1
## RFrCM = [ igamma , ip , istate , it ] bs1

print ''
print 'Writing Data'

outputdir = dir+'/results/'
PrintSSSetToFile(RFrff,,ReadMomList,CombGammaList,tsink)
PrintSSSetToFile([data2pt],outputdir,ReadMomList,['Mass'],0)
PrintCMSetToFile(RFrCM,outputdir,ReadMomList,CombGammaList,tsink,to,dt)
PrintCMMassToFile([CMdata2pt],outputdir,ReadMomList,['Mass'],0,to,dt)
PrintLREvecMassToFile(LREvec,Emass,outputdir,ReadMomList,to,dt)

print ''
print 'finding fits, fit range: ' + str(FitMin) + ' to ' + str(FitMax)

[CMFitBoot, CMFitAvg] = FitRFSet(RFrCM,FitMin,FitMax)
[ffFitBoot, ffFitAvg] = FitRFSet(RFrff,FitMin,FitMax)

#CMFitBoot = [ igamma , ip , istate ] bs1
#CMFitAvg = [  igamma , ip , ism , Avg/Std ]

#ffFitBoot = [ igamma , ip , ism ] bs1
#ffFitAvg = [  igamma , ip , ism , Avg/Std ]

print ''
print 'debug, Fit value for state 1 , ' + CombGammaList[0] , ipTOqvec(ReadMomList[0])
print ffFitAvg[0][0][0][0] , ffFitAvg[0][0][0][1]

#!/usr/bin/env python
from array import array
from os import walk
import os
import numpy as np
from Params import *
from ReadBinaryCfuns import ReadFSCfun, Read2ptCfun
from CfunBoot import BootSet3pt, BootSet2pt,BootSet3ptCM, BootSet2ptCM
from RFCalc import CalcRatioFactor, CalcRatioFactorSS
from ReadCMCfuns import ReadCMSet,ReadCMList
from CMTech import CreateCMCfuns
from OutputData import PrintSSSetToFile
import sys

dir = "/raid/jdragos/scratch/cfun/2ndk12090/"
logfile = '/home/accounts/jdragos/scripts/PythonAnalysis/TrySSRead.log'
sys.stdout = open(logfile, 'a', 0)
sys.sterr = sys.stdout
print "qlow to qhigh is: " , ipTOiq(0), ipTOiq(nmom-1)

### Read CM Set ###
[cmdata2pt,cmdata3pt,filelist] = ReadCMSet('doub','4','',dir,'cm')
###   ###


# ### Read CM from List ###
# [cmdata2pt,cmdata3pt] = ReadCMList('doub','4','',conflist,'cm')
# filelist = conflist
# ###   ###

ncon = np.size(filelist)
print 'ncon = ' + str(ncon)


## cmdata2pt { 'real'/'cmplx' } = [ ism , jsm , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
## cmdata3pt { 'real'/'cmplx' } = [ ism , jsm , itsink , igamma , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)

# [RFi,SqrtFac] = CalcRatioFactorSS(cmdata2pt['real'] ,cmdata3pt['cmplx'])
print 'Creating RatioFactor'
[RFr,SqrtFac] = CalcRatioFactorSS(cmdata2pt['real'] ,cmdata3pt['real'])
##RFr = [ itsink , igamma , ip , ism , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)

outputdir = dir+"results/" 
print "printing results to " + outputdir
mkdir_p(outputdir)

PrintSSSetToFile(RFr,outputdir)


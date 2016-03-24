#!/usr/bin/env python
from array import array
from os import walk
import os
import numpy as np
import sys
from Params import *
from CfunBoot import BootSet3ptCM, BootSet2ptCM
from RFCalc import CalcRatioFactorCM, CalcRatioFactorSS
from ReadCMCfuns import ReadCMSet,ReadCMList
from CMTech import CreateCMCfuns


dir = "/raid/jdragos/scratch/cfun/2ndk12090/"
logfile = '/home/accounts/jdragos/scripts/PythonAnalysis/TryCMRead.log'

sys.stdout = open(logfile,'a',0)
sys.stderr = sys.stdout
print '----------------------------------------------------------------------------------'
print "qlow to qhigh is: " , ipTOiq(0), ipTOiq(nmom-1)


# ### Read CM Set ###
# [cmdata2pt,cmdata3pt,filelist] = ReadCMSet('doub','4','',dir,'cm')
# ###   ###


### Read CM from List ###
[cmdata2pt,cmdata3pt,filelist] = ReadCMList('doub','4','',conflist,'cm')
###   ###

## cmdata2pt = [ r/i , ism , jsm , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
## cmdata3pt = [ r/i , ism , jsm , tsink , igamma , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)

ncon = np.size(filelist)
print 'ncon = ' + str(ncon)


[CMdata2pt,CMdata3ptreal,CMdata3ptcmplx] = CreateCMCfuns(cmdata3pt,cmdata2pt[0],'18','2')
## CMdata2pt [ ip , istate , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
## CMdata3pt(real/cmplx)  [ ip , igamma , istate , it] = bootstrap1 class (.Avg, .Std, .values, .nboot)


# [RFrffcmplx,SqrtFac] = CalcRatioFactorSS(bootdata2ptCM['real'],bootdata3ptCM['cmplx'])
[RFrff,SqrtFac] = CalcRatioFactorSS(bootdata2ptCM['real'],bootdata3ptCM['real'])

# [RFrcmplx,SqrtFac] = CalcRatioFactorCM(CMdata2pt,CMdata3ptcmplx)
[RFr,SqrtFac] = CalcRatioFactorCM(CMdata2pt,CMdata3ptreal,str(tsink))



#!/usr/bin/env python

from Params import *
import numpy as np
from GraphSummary import ReadAndPlotSummary
from SetLists import *
from OppFuns import CreateGammaList
from FFParams import *
import time,datetime
from multiprocessing import Pool
from MultiWrap import *
from InputArgs import *


feedin = InputParams(sys.argv[1:] + ['-noprompt'])

thisGammaList = CreateGammaList(feedin['gamma'],twopt=False)

# MethodList.remove('RF')
# MethodList.remove('OSFCM')
# MethodList.remove('OSFTsink')

print 'SetLists:\n','\n'.join(feedin['set']) + '\n'
# print 'MomList:\n','\n'.join(feedin['mom']) + '\n'
if 'RF' in feedin['method']: feedin['method'].remove('RF')

if DoMulticore and feedin['anaproc'] > 1:
    inputparams = []
    for igamma in thisGammaList:
        if 'doub' not in igamma and 'sing' not in igamma:
            inputparams.append((feedin['method'],['doub'+igamma,'sing'+igamma,'twopt'],feedin['set'],feedin['mom'],feedin['comb']))
        elif igamma.replace('doub','').replace('sing','') not in thisGammaList:
            inputparams.append((feedin['method'],[igamma,'twopt'],feedin['set'],feedin['mom'],feedin['comb']))            
    makeContextFunctions(ReadAndPlotSummary)
    thisPool = Pool(min(len(inputparams),feedin['anaproc']))
    thisPool.map(ReadAndPlotSummary.mapper,inputparams)
    thisPool.close()
    thisPool.join()
else:
    ReadAndPlotSummary(feedin['method'],thisGammaList,feedin['set'],feedin['mom'],feedin['comb'])    
    

# for iin in inputparams: ReadAndPlotSummary(*iin)

print 'Graphing Complete                             '

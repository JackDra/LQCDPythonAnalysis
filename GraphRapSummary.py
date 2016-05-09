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


feedin = InputParams(sys.argv[1:])

thisGammaList = CreateGammaList(feedin['gamma'],twopt=False)

# MethodList.remove('RF')
# MethodList.remove('OSFCM')
# MethodList.remove('OSFTsink')

print 'SetLists:\n','\n'.join(feedin['set']) + '\n'
# print 'MomList:\n','\n'.join(feedin['mom']) + '\n'
if 'RF' in feedin['method']: feedin['method'].remove('RF')

if DoMulticore:
    inputparams = []
    for igamma in thisGammaList:
        inputparams.append((feedin['method'],[igamma,'twopt'],feedin['set'],feedin['mom']))
    makeContextFunctions(ReadAndPlotSummary)
    thisPool = Pool(min(len(inputparams),feedin['anaproc']))
    thisPool.map(ReadAndPlotSummary.mapper,inputparams)
    thisPool.close()
    thisPool.join()
else:
    ReadAndPlotSummary(feedin['method'],thisGammaList,feedin['set'],feedin['mom'])    
    

# for iin in inputparams: ReadAndPlotSummary(*iin)

print 'Graphing Complete                             '

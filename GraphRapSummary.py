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



# thisMomList = ['q = 0 0 0']
thisMomList = ['q = 0 0 0','q = -1 0 0' , 'q = -2 0 -1' , 'q = -2 -2 0' ]
thisGammaList = CreateGammaList(sys.argv[1:],twopt=False)

MethodList.remove('RF')
MethodList.remove('OSFCM')
MethodList.remove('OSFTsink')
print 'SetLists:\n','\n'.join(DefSetList) + '\n'

inputparams = []
for igamma in thisGammaList:
    inputparams.append((MethodList,[igamma,'twopt'],DefSetList,thisMomList))
makeContextFunctions(ReadAndPlotSummary)
thisPool = Pool(min(len(inputparams),AnaProc))
thisPool.map(ReadAndPlotSummary.mapper,inputparams)
thisPool.close()
thisPool.join()


# for iin in inputparams: ReadAndPlotSummary(*iin)

print 'Graphing Complete                             '

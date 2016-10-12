#!/usr/bin/env python

from Params import *
from CreateCombs import *
from InputArgs import *
from OppFuns import *
from ReadDir import *
from copy import deepcopy
import time
import datetime

feedin = InputParams(sys.argv[1:])
thisGammaList = CreateGammaList(feedin['gamma'])

ShowSetLists(feedin['set'])
ShowCombList(feedin['comb'])
ShowMethodList(feedin['method'])

thisCurrList = deepcopy(feedin['current'])
for icurr in feedin['current']:
    if ('doub' not in icurr) and ('sing' not in icurr):
        if 'doub'+icurr not in feedin['current']:
            thisCurrList.append('doub'+icurr)
        if 'sing'+icurr not in feedin['current']:
            thisCurrList.append('sing'+icurr)
        thisCurrList.remove(icurr)


for iDS in DefDSList:
    if iDS in feedin['comb']: feedin['comb'].remove(iDS)
        
starttime = time.time()
ReadAndComb(feedin,[CombFunsDict[iCombType] for iCombType in feedin['comb']],feedin['comb'])
    
if feedin['DoCurr']:
    ReadAndCombFF(GetCurrDict(thisCurrList)[kappa],[CombFunsDict[iCombType] for iCombType in feedin['comb']],feedin['comb'])

print 'Combining Complete, time taken:', str(datetime.timedelta(seconds=time.time()-starttime)) , ' h:m:s '


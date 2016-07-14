#!/usr/bin/env python

from Params import *
from CreateCombs import *
from InputArgs import *
from OppFuns import *
from ReadDir import *
from copy import deepcopy

feedin = InputParams(sys.argv[1:])
thisGammaList = CreateGammaList(feedin['gamma'])

ShowSetLists(feedin['set'])
ShowCombList(feedin['comb'])

thisCurrList = deepcopy(feedin['current'])
for icurr in feedin['current']:
    if ('doub' not in icurr) and ('sing' not in icurr):
        if 'doub'+icurr not in feedin['current']:
            thisCurrList.append('doub'+icurr)
        if 'sing'+icurr not in feedin['current']:
            thisCurrList.append('sing'+icurr)
        thisCurrList.remove(icurr)

for iCombType in feedin['comb']:
    ReadAndComb(feedin,CombFunsDict[iCombType],iCombType)
    ReadAndCombFF(GetCurrDict(thisCurrList),CombFunsDict[iCombType],iCombType)

#!/usr/bin/env python

from Params import *
from CreateCombs import *
from InputArgs import *
from OppFuns import *




feedin = InputParams(sys.argv[1:])
thisGammaList = CreateGammaList(feedin['gamma'])

ShowSetLists(feedin['set'])
ShowCombList(feedin['comb'])

for iCombType in feedin['comb']:
    ReadAndComb(feedin,CombFunsDict[iCombType],iCombType)
    ReadAndCombFF(GetCurrDict(feedin['current']),CombFunsDict[iCombType],iCombType)

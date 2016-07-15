#!/usr/bin/env python

from Params import *
from CreateCombs import *
from InputArgs import *
from OppFuns import *
from ReadDir import *
from MiscFuns import *
from copy import deepcopy

feedin = InputParams(sys.argv[1:])
thisGammaList = CreateGammaList(feedin['gamma'])

ShowSetLists(feedin['set'])
ShowCombList(feedin['comb'])
ShowMethodList(feedin['method'])

thisCurrList = ElongateName(feedin['comb'],feedin['current'])

funlist,comblist = [].[]
for iFFcomb in feedin['FFcomb']:
    print 'Creating ' , iFFcomb
    funlist.append(CombFFFunsDict[iFFcomb])
    comblist.append(iFFcomb)
    
ReadAndCombTheFFs(GetCurrDict(thisCurrList),funlist,comblist)

print 'FF Combining Complete'

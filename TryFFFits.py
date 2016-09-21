#!/usr/bin/env python

import numpy as np
from Params import *
from MiscFuns import *
from ReadTxt import ExtractValues
from ReadDir import GetCurrDict
from ReadTxt import ReadFFDict
from SetLists import *
from FormFactors import CreateFF
from OutputData import PrintFFSet
from OppFuns import PrintOpps
from FFParams import *
import sys
from multiprocessing import Pool
from MultiWrap import *
import time
import datetime
# from guppy import hpy; h=hpy()
import resource
from InputArgs import *



feedin = InputParams(sys.argv[1:])

print 'CurrList:\n' , '\n'.join(feedin['current'])
print ''

if 'RF' in feedin['method']: feedin['method'].remove('RF')
print 'MethodList:\n' , '\n'.join(feedin['method'])
print ''


thisFFcomb = []
for icomb in feedin['comb']:
    if icomb in CombListForFFComb:
        thisFFcomb.append(icomb)
thisCurrcomb = []
for icomb in feedin['current']:
    if icomb in CurrListForFFComb:
        thisCurrcomb.append(icomb)

thisCurrDict = []
for iFFcomb in feedin['FFcomb']:
    if iFFcomb == '':
        for icurr in ElongateName(feedin['comb'],feedin['current']):
            print 'Looking in ', icurr
            thisCurrDict.append([GetCurrDict([icurr]),feedin['ffgraph']])
    else:
        for icurr in ElongateName(thisFFcomb,thisCurrcomb):
            print 'Looking in ', icurr+iFFcomb
            thisCurrDict.append([GetCurrDict([icurr+iFFcomb]),feedin['ffgraph']])

## datadict { FormFactor } { Set } { Mass:Set/Avg/Std/Chi/Boot , FF#:qsqrd:Avg/Std/Boot , Chi:qsqrd}
datadict = ReadFFDict(outputdir,thisCurrDict)

for iFF,FFdata in datadict.iteritems():
    outputdict = OrderedDict()
    for iSet,Setdata in FFdata.iteritems():
        outputdict[iSet] = OrderedDict()
        iFFloop = Setdata.keys()
        iFFloop.pop('Mass',None)
        iFFloop.pop('Chi',None)
        if Debug: print iFFloop
        for nFF in iFFloop:
            outputdict[iSet][nFF] = OrderedDict()
            nFFdata = Setdata[nFFdata]
            ydatain,xdatain = []
            for iQs,Qsdata in nFFdata.iteritems():
                ydatain.append(Qsdata['Boot'])
                xdatain.append(GetQsqrd(float(qsqrd.replace('qsqrd','')),Phys=PhysicalUnits))
            DPfit,DPfitAvg,DPfitChi = FitBoots(ydatain,xdatain,DPfit)
            outputdict[iSet][nFF]['Boot'],outputdict[iSet][nFF]['Avg'],outputdict[iSet][nFF]['Chi'] = DPfit,DPfitAvg,DPfitChi
    PrintDPfit(iFF,outputdict)

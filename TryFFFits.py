#!/usr/bin/env python

import numpy as np
from Params import *
from MiscFuns import *
from ReadTxt import ExtractValues
from ReadDir import GetCurrDict
from ReadTxt import ReadFFDict
from SetLists import *
from FormFactors import CreateFF
from OutputXmlData import PrintDPfit
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

datadict = OrderedDict()
for iFFcomb in feedin['FFcomb']:
    if iFFcomb == '':
        for icurr in ElongateName(feedin['comb'],feedin['current']):
            print 'Looking in ', icurr
            datadict[icurr] = ReadFFDict(outputdir,GetCurrDict([icurr]))
    else:
        for icurr in ElongateName(thisFFcomb,thisCurrcomb):
            print 'Looking in ', icurr+iFFcomb
            datadict[icurr+iFFcomb.replace('/','')] = ReadFFDict(outputdir,GetCurrDict([icurr+iFFcomb]))

if Debug: print datadict[datadict.keys()[0]].keys()

## datadict { FormFactor } { Set } { Mass:Set/Avg/Std/Chi/Boot , FF#:qsqrd:Avg/Std/Boot , Chi:qsqrd}
    
for iFF,FFdata in datadict.iteritems():
    outputdict = OrderedDict()
    for iSet,Setdata in FFdata.iteritems():
        outputdict[iSet] = OrderedDict()
        iFFloop = Setdata.keys()
        if 'Mass' in iFFloop: iFFloop.remove('Mass')
        if 'Chi' in iFFloop: iFFloop.remove('Chi')
        if Debug: print iFFloop
        for nFF in iFFloop:
            outputdict[iSet][nFF] = OrderedDict()
            nFFdata = Setdata[nFF]
            ydatain,xdatain = [],[]
            for iQs,Qsdata in nFFdata.iteritems():
                ydatain.append(Qsdata['Boot'])
                xdatain.append(GetQsqrd(float(qsqrd.replace('qsqrd','')),Phys=PhysicalUnits))
            DPfit,DPfitAvg,DPfitChi = FitBoots(ydatain,xdatain,DPfit)
            outputdict[iSet][nFF]['Boot'],outputdict[iSet][nFF]['Avg'],outputdict[iSet][nFF]['Chi'] = DPfit,DPfitAvg,DPfitChi
    PrintDPfit(iFF,outputdict)

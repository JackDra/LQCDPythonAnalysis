#!/usr/bin/env python

import numpy as np
from Params import *
from MiscFuns import *
from FitFunctions import DPfitfun,DPfitfunDer
from LLSBoot import *
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
            datadict[icurr] = datadict[icurr][datadict[icurr].keys()[0]]
    else:
        for icurr in ElongateName(thisFFcomb,thisCurrcomb):
            print 'Looking in ', icurr+iFFcomb
            datadict[icurr+iFFcomb.replace('/','')] = ReadFFDict(outputdir,GetCurrDict([icurr+iFFcomb]))
            datadict[icurr+iFFcomb.replace('/','')] = datadict[icurr+iFFcomb.replace('/','')][datadict[icurr+iFFcomb.replace('/','')].keys()[0]]
# if Debug: print datadict.keys()
# if Debug: print datadict[datadict.keys()[0]].keys()

## datadict { Set } { Mass:Set/Avg/Std/Chi/Boot , FF#:qsqrd:Avg/Std/Boot , Chi:qsqrd}
    
  
def CurrFFDPfit(iCurr,Currdata):
    outputdict = OrderedDict()
    CurrInfo = False
    if Debug: print 'iCurr' , icurr
    for iSet,Setdata in Currdata.iteritems():
        if Debug: print 'iSet', iSet
        outputdict[iSet] = OrderedDict()
        iFFloop = Setdata.keys()
        if 'Mass' in iFFloop: iFFloop.remove('Mass')
        if 'Chi' in iFFloop: iFFloop.remove('Chi')
        if 'Info' in iFFloop:
            CurrSetInfo = Setdata['Info']
            iFFloop.remove('Info')
        for nFF in iFFloop:
            outputdict[iSet][nFF] = OrderedDict()
            nFFdata = Setdata[nFF]
            ydatain,xdatain = [],[]
            for iQs,Qsdata in nFFdata.iteritems():
                if 'Boot' in Qsdata:
                    ydatain.append(Qsdata['Boot'])
                    xdatain.append([GetQsqrd(float(iQs.replace('qsqrd','')),Phys=PhysicalUnits)])
                else:
                    print 'Warning, Boot not found in', iCurr, iSet, nFF, iQs 
            # if Debug:
            #     print 'Fitting to points:'
            #     for ix,iy in zip(xdatain, ydatain):
            #         print ix, iy.Avg
            if len(ydatain) < 2:
                print "too short ydata, skipping",iCurr, iSet, nFF, iQs 
            else:
                DPfit,DPfitAvg,DPfitChi = FitBoots(ydatain,xdatain,DPfitfun)
                outputdict[iSet][nFF]['Boot'],outputdict[iSet][nFF]['Avg'],outputdict[iSet][nFF]['Chi'] = DPfit,DPfitAvg,DPfitChi
    PrintDPfit(iCurr,outputdict,CurrSetInfo)


inputparams = []
for iCurr,Currdata in datadict.iteritems():
    inputparams.append([iCurr,Currdata])


starttime = time.time()
feedin['anaproc'] = min(feedin['anaproc'],len(inputparams))
if DoMulticore and feedin['anaproc'] > 1:
    makeContextFunctions(CurrFFDPfit)
    thisPool = Pool(feedin['anaproc'])
    thisPool.map(CurrFFDPfit.mapper,inputparams)
    thisPool.close()
    thisPool.join()
else:
    for ip,iparam in enumerate(inputparams):
        print 'Total percent: ' ,GetPercent(ip,len(inputparams))
        print 'Time Left:' , GetTimeLeftStr(ip,len(inputparams),time.time()-starttime)
        CurrFFDPfit(*iparam)
        print 'Form Factor Creation Complete, time taken:', str(datetime.timedelta(seconds=time.time()-starttime)) , ' h:m:s '
                                        

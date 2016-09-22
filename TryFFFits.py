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
from BootTest import BootStrap1


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
    
  
def CurrFFDPfit(iCurr,Currdata,thisSetList,thisMethodList):
    outputdict = OrderedDict()
    CurrInfo = False
    if Debug: print 'iCurr' , icurr
    # if 'tsink29sm32Fitscut5' in Currdata.keys():
    #     iSet = 'tsink29sm32Fitscut5'
    #     Setdata = Currdata[iSet]
    ZeroBoot = BootStrap(nboot,1)
    OneBoot = BootStrap(nboot,1)
    OneBoot.values = np.array([1.0]*nboot)
    OneBoot.Stats()
    TwoBoot = BootStrap(nboot,1)
    TwoBoot.values = np.array([2.0]*nboot)
    TwoBoot.Stats()
    for iSet,Setdata in Currdata.iteritems():    
        if not any([imethod in iSet for imethod in thisMethodList]): continue 
        if 'TSF' in iSet or 'SumMeth' in iSet:
            if not any([RemoveTSink(iset) in iSet for iset in thisSetList]): continue 
        else:
            if not any([iset in iSet for iset in thisSetList]): continue                 
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
                if (iQs == 'qsqrd0') and ('1' not in nFF): continue
                if iQs == 'qsqrd0' and '1' in nFF and ('Ge' in iCurr or ('Vector' in iCurr and 'PsVector' not in iCurr.replace('IsoVector',''))):
                    if 'IsoVector' in iCurr or 'Proton' in iCurr or 'sing' in iCurr:
                        ydatain.append(OneBoot)
                        xdatain.append(0.0)
                    elif 'Neutron' in iCurr:
                        ydatain.append(ZeroBoot)
                        xdatain.append(0.0)                        
                    elif 'doub' in iCurr:
                        ydatain.append(TwoBoot)
                        xdatain.append(0.0)                        
                elif 'Boot' in Qsdata:
                    ydatain.append(Qsdata['Boot'])
                    xdatain.append(GetQsqrd(float(iQs.replace('qsqrd','')),Phys=PhysicalUnits))
                else:
                    print 'Warning, Boot not found in', iCurr, iSet, nFF, iQs 
            if Debug:
                print 'Fitting to points:'
                for ix,iy in zip(xdatain, ydatain):
                    print ix, iy.Avg
            if len(ydatain) < 2:
                print "too short ydata, skipping",iCurr, iSet, nFF, iQs 
            else:
                DPfit,DPfitAvg,DPfitChi = FitBoots(np.array(ydatain),np.array(xdatain),DPfitfun)
                outputdict[iSet][nFF]['Boot'],outputdict[iSet][nFF]['Avg'],outputdict[iSet][nFF]['Chi'] = DPfit,DPfitAvg,DPfitChi
    PrintDPfit(iCurr,outputdict,CurrSetInfo)


inputparams = []
for iCurr,Currdata in datadict.iteritems():
    inputparams.append([iCurr,Currdata,feedin['set'],feedin['method']])


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
                                        

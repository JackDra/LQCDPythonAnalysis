#!/usr/bin/env python

import numpy as np
from Params import *
from MiscFuns import *
from FitFunctions import DPfitfun,DPfitfunDer
from FitFunctions import DPfitfunOnePar,DPfitfunOneParDer
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
    ZeroBoot = BootStrap1(nboot,1)
    OneBoot = BootStrap1(nboot,1)
    OneBoot.values = np.array([1.0]*nboot)
    OneBoot.Stats()
    TwoBoot = BootStrap1(nboot,1)
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
            yZero = False
            for iQs,Qsdata in nFFdata.iteritems():
                if iQs == 'qsqrd0':
                    # continue
                    if '1' not in nFF:
                        yZero = False
                    else:
                        if ('Ge' in iCurr or ('Vector' in iCurr and 'PsVector' not in iCurr.replace('IsoVector',''))):
                            if 'IsoVector' in iCurr or 'Proton' in iCurr or 'sing' in iCurr:
                                yZero = OneBoot
                            elif 'Neutron' in iCurr:
                                yZero = ZeroBoot
                            elif 'doub' in iCurr:
                                yZero = TwoBoot
                            else:
                                yZero = OneBoot
                        elif 'Boot' in Qsdata:
                            ydatain.append(Qsdata['Boot'])
                            xdatain.append(GetQsqrd(float(iQs.replace('qsqrd','')),Phys=PhysicalUnits))
                elif 'Boot' in Qsdata:
                    ydatain.append(Qsdata['Boot'])
                    xdatain.append(GetQsqrd(float(iQs.replace('qsqrd','')),Phys=PhysicalUnits))
                else:
                    print 'Warning, Boot not found in', iCurr, iSet, nFF, iQs 
            if len(ydatain) < 2:
                print "too short ydata, skipping",iCurr, iSet, nFF, iQs 
            else:
                if yZero == False:
                    if Debug:
                        print 'Fitting to points using two parameter fit:'
                        for ix,iy in zip(xdatain, ydatain):
                            print ix, iy.Avg, iy.Std
                    DPfit,DPfitAvg,DPfitChi = FitBoots(np.array(ydatain),np.array(xdatain),DPfitfun)
                    DPfit[0].values = np.abs(DPfit[0].values)
                    DPfit[0].Stats()
                    outputdict[iSet][nFF]['Boot'],outputdict[iSet][nFF]['Avg'],outputdict[iSet][nFF]['Chi'] = [DPfit[0],DPfit[1]],[abs(DPfitAvg[0]),DPfitAvg[1]],DPfitChi
                else:
                    if Debug:
                        print 'Fitting to points using one parameter fit:'
                        for ix,iy in zip(xdatain, ydatain):
                            print ix, iy.Avg, iy.Std
                    DPfit,DPfitAvg,DPfitChi = FitBoots(np.array(ydatain),np.array(xdatain),DPfitfunOnePar)
                    outputdict[iSet][nFF]['Boot'],outputdict[iSet][nFF]['Avg'],outputdict[iSet][nFF]['Chi'] = [yZero,DPfit[0]],[yZero.Avg,DPfitAvg[0]],DPfitChi*2 
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
                                        

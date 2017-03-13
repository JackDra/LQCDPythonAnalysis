#!/usr/bin/env python

import numpy as np
from Params import *
from MiscFuns import *
from ReadTxt import ExtractValues,ReadAlphaList
from SetLists import *
from FormFactors import CreateFF
from OutputData import PrintFFSet
from FFFuns import *
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
from FitFunctions import *
from collections import OrderedDict
from BootTest import BootStrap1
from LLSBoot import *




# ## data { { gamma } { mom } { Fit(Boot/Avg/Std/Chi) } }
# ## dataout { { momsqrd } { Boot/Avg/Chi }
# ##REMEBER deal with cmplx signals
# def CreateFF(data,mass,iCurr,gammaflag='',Rfac=True,alphalist = [1.0]):
#     thisdataout = {}
#     DoTop = 'Top' in iCurr
#     baseCurr = iCurr.replace(gammaflag,'')
#     Opps = CurrOpps[baseCurr]
#     thisdataout = OrderedDict()
#     infodict = OrderedDict()
#     for iqsqrd in MomSqrdSet:
#         momstart = time.time()
#         iqs = 'qsqrd'+str(iqsqrd)
#         thisdataout['qsqrd'+str(iqsqrd)] = {}
#         datavals,FFcoeff = [],[]
#         for ica in xrange(len(alphalist)):
#             FFcoeff.append([])
#             for iff in xrange(NoFFPars[baseCurr]):
#                 FFcoeff[ica].append([])
#         opplist = []
#         for iopp in Opps:                    
#             flagopp = gammaflag+iopp
#             RealVal,CmplxVal = False,False
#             if flagopp in data.keys(): RealVal = True
#             if flagopp+'cmplx' in data.keys(): CmplxVal = True
#             if not RealVal and not CmplxVal: continue
#             start = time.time()
#             for iq in qvecINqsqrd(int(iqsqrd)):
#                 if RealVal: 
#                     if iq not in data[flagopp].keys(): continue 
#                 if CmplxVal: 
#                     if iq not in data[flagopp+'cmplx'].keys(): continue 
#                 for ica,ialpha in enumerate(alphalist):
#                     FFcoeffhold,rcheck,ccheck = CurrFFs[baseCurr](iopp,np.array(qstrTOqvec(iq))*qunit,[0,0,0],mass,Rfac=Rfac,alpha=ialpha)
#                     if CmplxVal and ccheck:
#                         for iFF,iFFcof in enumerate(FFcoeffhold):
#                             FFcoeff[ica][iFF].append(iFFcof.imag)
#                     if RealVal and rcheck:
#                         for iFF,iFFcof in enumerate(FFcoeffhold):
#                             FFcoeff[ica][iFF].append(iFFcof.real)
#                     if ica == 0:
#                         if CmplxVal and ccheck:
#                             datavals.append(data[flagopp+'cmplx'][iq]['Boot'])
#                             infodict[iqs] = data[flagopp+'cmplx'][iq]['Info']
#                             opplist.append(flagopp+'cmplx '+ iq)
#                         if RealVal and rcheck:
#                             datavals.append(data[flagopp][iq]['Boot'])
#                             infodict[iqs] = data[flagopp][iq]['Info']
#                             opplist.append(flagopp +' '+ iq)
#             # print 'PullOutLHS iopp',iopp,', time taken:' , GetTimeStr(time.time()-start)
#         if len(datavals) == 0: continue
#         thisdataout[iqs]['Chi'] = FFChihold[0]
#         # print 'All Done',iqs,', time taken:' , GetTimeStr(time.time()-momstart)
#     return thisdataout,infodict



## WORK IN PROGRESS

## topdata { mom : flow }
## F1F2data { FF# : qsqrd : Boot//Info}
def ExtratF3(topdata,F1F2data,alphalist,thisInfo):    
    for iqsqrd,ff1data in F1F2data['FF1'].iteritems():
        if iqsqrd in F1F2data['FF2'].keys():
            ff2data = F1F2data['FF2'][iqsqrd]
            FFcoeffhold,rcheck,ccheck = CurrFFs[baseCurr](iopp,np.array(qstrTOqvec(iq))*qunit,[0,0,0],mass,Rfac=Rfac,alpha=ialpha) 
           


feedin = InputParams(sys.argv[1:])

if 'RF' in feedin['method']: feedin['method'].remove('RF')
print 'MethodList:\n' , '\n'.join(feedin['method'])
print ''

start = time.time()
## data { set } { gamma } { mom }
data,MassSet = ExtractValues(outputdir[0],['doubP3g4Top','singP3g4Top'],feedin['set'],feedin['method'],thisMomList=RunMomList,TopRead=True)
# FFdata { Set } { Mass:Set/Avg/Std/Chi/Boot , FF#:qsqrd:Avg/Std/Boot/Info , Chi:qsqrd}
FFdata = ReadFFDict(outputdir[0],GetCurrDict(['Vector']))
## alphalist = [iAvg,boot1,boot2,...,bootn]
alphalist,alphainfo = ReadAlphaList(theset)
thisInfo = alphainfo
thisInfo['Mass'] = DefMassVal[massfitr]
thisInfo['Mass']['Set'] = 'Default'
thisInfo['Mass']['fit_range'] = massfitr
doubFFdata = FFdata['doubVector']
singFFdata = FFdata['singVector']
print 'data read, time taken:', str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s '

for iset,setdata in data.iteritems():
    for doubffset,doubffsetdata in doubFFdata.iteritems():
        if iset in doubffset and doubffset.replace('doub','sing') in singFFdata.keys() :
            singffsetdata = singFFdata[doubffset.replace('doub','sing')]
            print 'solving set ',iset,' with ffset ' , doubffset.replace('doub','')
            thisdataF3doub,infodictdoub = ExtractF3(setdata['doubP3g4Top'],doubffsetdata,thisInfo)
            thisdataF3sing,infodictsing = ExtractF3(setdata['singP3g4Top'],singffsetdata,thisInfo)
            


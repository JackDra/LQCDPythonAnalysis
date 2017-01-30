#!/usr/bin/env python

import numpy as np
from Params import *
from MiscFuns import *
from ReadTxt import ExtractValues,ReadAlphaList
from SetLists import *
from FormFactors import CreateFF
from OutputData import PrintFFSet
from FFFuns import RenormFF,CombineVector
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



def PickMassSet(MassSet,theset):
    theMass = DefMassVal[DefMassVal.keys()[0]]
    thesetmass = 'Default'
    # if 'SF' in theset and RemoveTSink(theset) in MassSet.keys():
    #     thesetmass = RemoveTSink(theset)
    #     theMass = MassSet[thesetmass]['Boot'].Avg
    # elif 'SumMeth' in theset and 'sm32OSFTsink' in MassSet.keys():
    #     thesetmass = 'sm32OSFTsink'
    #     theMass = MassSet[thesetmass]['Boot'].Avg
    return theMass,thesetmass

def CreateFFWrap(thisMass,thesetmass,theset,setdict,thisCurr,Rfac):
    # mprint( 'Set:' + theset + ' MassSetPicked:'+thesetmass)
## FF { { momsqrd } { Boot/Avg/Chi } }
    thisstart = time.time()
    thisDS,baseCurr,dump = SplitDSCurr(thisCurr)
    if 'Top' in thisCurr:
        baseCurr = baseCurr+'Top'
        alphalist,alphadata = ReadAlphaList(theset)
    else:
        alphalist = [1.0]
        alphainfo = OrderedDict()
        alphainfo['Avg'] = 1.0
        alphainfo['Std'] = 0.0
        alphainto['Chi'] = 0.0
        alphainfo['fit_range'] = 'fitr0-0'
        alphainfo['File'] = ''
    if thisDS == '':
        thisCurrList = [thisCurr for ids in DefDSList]
        thisgflist = DefDSList
    else:
        thisCurrList = [thisCurr.replace(thisDS,'')]
        thisgflist = [thisDS]
    for iCurr,igf in zip(thisCurrList,thisgflist):
        combCurr = igf+iCurr
        if NoSFRfacScale:
            FF,infodict = CreateFF(setdict,thisMass['Avg'],iCurr,gammaflag=igf,Rfac=Rfac,alphalist=alphalist)
        else:
            FF,infodict = CreateFF(setdict,thisMass['Avg'],iCurr,gammaflag=igf,Rfac=True,alphalist=alphalist)            
        infodict['Mass'] = thisMass
        infodict['Mass']['Set'] = thesetmass
        if 'Vector' in thisCurr and 'Top' not in thisCurr and 'IsoVector' not in thisCurr and 'PsVector' not in thisCurr:
            if ForceVecNorm: FF = RenormFF(FF,FF['qsqrd0']['Boot'][0].Avg,igf)
            PrintFFSet(FF,theset,thisMass,thesetmass,combCurr,infodict)
            NewFF = CombineVector(FF,thisMass)
            PrintFFSet(NewFF,theset,thisMass,thesetmass,combCurr.replace('Vector','GeGm'),infodict)
        else:
            PrintFFSet(FF,theset,thisMass,thesetmass,combCurr,infodict)

    print 'Fit and Print for ' , theset , ' took: ',str(datetime.timedelta(seconds=time.time()-thisstart)) , ' h:m:s'   


#FitMasses later
# def DoFF(thisMethodList,thisCurr,thisSetList,thisGammaList,thisMomList):
def DoFF(inputlist):
    datalist = []
    MassSetlist = []
    for iin in inputlist:
        thisMethodList,thisCurr,thisSetList,thisGammaList,thisMomList = iin
        data,MassSet = ExtractValues(outputdir[0],thisGammaList,thisSetList,thisMethodList,thisMomList=thisMomList,TopRead='Top' in thisCurr)
        print 'data Collected:'
        for iCol in data.keys():
            print thisMethodList[0], thisCurr, iCol
        print ''
        datalist.append(data)
        MassSetlist.append(MassSet)
    if len(data.keys()) == 0: return

    inputparams = []
    for data,MassSet in zip(datalist,MassSetlist):
        for theset,setdict in data.iteritems():
            # print outputdir[0] +'/FormFactors/'+DefDSList[0]+thisCurr+'/'+DefDSList[0]+thisCurr+theset+'.xml'
            if (not all([os.path.isfile(outputdir[0] +'/FormFactors/'+iDS+thisCurr+'/'+iDS+thisCurr+theset+'.xml') for iDS in DefDSList])) or DefWipe:
                inputparams.append(PickMassSet(MassSet,theset)+(theset,setdict,thisCurr,'SF' not in theset))


    starttime = time.time()
    print 'Total Paralell:', len(inputparams)
    feedin['anaproc'] = min(feedin['anaproc'],len(inputparams))
    if DoMulticore and feedin['anaproc'] > 1:
        makeContextFunctions(CreateFFWrap)
        thisPool = Pool(feedin['anaproc'])
        thisPool.map(CreateFFWrap.mapper,inputparams)
        thisPool.close()
        thisPool.join()
    else:
        for ip,iparam in enumerate(inputparams):
            print 'Total percent: ' ,GetPercent(ip,len(inputparams))
            print 'Time Left:' , GetTimeLeftStr(ip,len(inputparams),time.time()-starttime)
            CreateFFWrap(*iparam)

feedin = InputParams(sys.argv[1:])

print 'CurrList:\n' , '\n'.join(feedin['current'])
print ''

if 'RF' in feedin['method']: feedin['method'].remove('RF')
print 'MethodList:\n' , '\n'.join(feedin['method'])
print ''

inputparams = []
thisGammaList = []
if 'GeGm' in feedin['current']: feedin['current'].remove('GeGm')
for iDS in DefDSList:
    if iDS+'GeGm' in feedin['current']: feedin['current'].remove(iDS+'GeGm')
for thisCurr in feedin['current']:
    thisDS,baseCurr,dump = SplitDSCurr(thisCurr)
    if thisDS == '':
        thisGammaList = CurrOpps[baseCurr] + [iC+'cmplx' for iC in CurrOpps[baseCurr]]
        thisGammaList = ['doub'+ig for ig in thisGammaList] + ['sing'+ig for ig in thisGammaList]+ ['twopt']
    else:
        thisGammaList = CurrOpps[baseCurr] + [iC+'cmplx' for iC in CurrOpps[baseCurr]]
        thisGammaList = [thisDS+ig for ig in thisGammaList] + ['twopt']
        
    for imeth in feedin['method']:
        if 'Fits' in imeth or 'OSF' in imeth:
            for iSet in feedin['set']:
                print 'Adding to queue FF: ' , imeth , thisCurr , iSet
                inputparams.append(([imeth],thisCurr,[iSet],thisGammaList,feedin['mom']))
        else:
            for iSet in ReduceTsink(feedin['set']):
                print 'Adding to queue FF: ' , imeth , thisCurr , iSet
                inputparams.append(([imeth],thisCurr,[iSet],thisGammaList,feedin['mom']))
            
starttime = time.time()
DoFF(inputparams)
print 'Form Factor Creation Complete, time taken:', str(datetime.timedelta(seconds=time.time()-starttime)) , ' h:m:s '

#!/usr/bin/env python


from OppFuns import CreateGammaList,WipeSet
from Params import *
import numpy as np
from BootTest import BootStrap1
from ReadTxt import *
from SetLists import *
from Fitting import FitAlpha
from OutputXmlData import PrintAlphaFitFile
from FitParams import *
from multiprocessing import Pool
from MultiWrap import *
from FFParams import *
import itertools as it
import sys
import time
import datetime
from InputArgs import *
from CheckXml import *



def TryAlphaFitsFun(thisSetList,thisMomList):
    # dataRF = { iset , imom , Info/ Boots: itflow , tsink } bs
    dataRF = ReadAlphaSet(thisSetList,thisMomList)
    start = time.time()
    FitData = OrderedDict()
    for iset,setdata in dataRF.iteritems():
        FitData[iset] = OrderedDict()
        for thismom,momdata in setdata.iteritems():
            # print 'Fitting ' , thisgamma , thismom , '       \r',
            start = time.time()
            FitData[iset][thismom] = OrderedDict()            
            FitData[iset][thismom]['Info'] = momdata['Info']
            FitData[iset][thismom]['Avg'] = OrderedDict()
            FitData[iset][thismom]['Boots'] = OrderedDict()
            FitData[iset][thismom]['Chi'] = OrderedDict()
            for iflow,flowdata in momdata['Boots'].iteritems():
                FitData[iset][thismom]['Avg'][iflow] = OrderedDict()
                FitData[iset][thismom]['Boots'][iflow] = OrderedDict()
                FitData[iset][thismom]['Chi'][iflow] = OrderedDict()
                fitdata = CreateBootClass(flowdata,nboot)
                for fitr,ifit in zip(FitAlphaList,FitAlphaArgs):
                    FitDatahold = FitAlpha(fitdata,fitr)
                    FitData[iset][thismom]['Avg'][iflow][ifit] = FitDatahold['Avg'][0]
                    FitData[iset][thismom]['Boots'][iflow][ifit] = FitDatahold['Boots'][0]
                    FitData[iset][thismom]['Chi'][iflow][ifit] = FitDatahold['Chi'][0]
            print ' Took: ',GetTimeStr(time.time()-start) , 'For ',thismom  , '   ', iset
    #FitData = { iset , ip , Boot/Avg/Chi , iflow , ifit   }
    return FitData,thisSetList



feedin = InputParams(sys.argv[1:])

feedin['set'] = ReduceTooMassSet(feedin['set'])

ShowSetLists(feedin['set'])

totstart = time.time()
inputparams = []
# RunGammaList = []
thisMomList = RunAvgMomList
for iSet in feedin['set']:
    for imom in thisMomList:
        inputparams.append(([iSet],[imom]))

if len(inputparams) > 0:
    if DoMulticore and feedin['anaproc'] > 1 and len(inputparams) > 1:
        print 'Running Multicore'
        makeContextFunctions(TryAlphaFitsFun)
        thisPool = Pool(min(len(inputparams),feedin['anaproc']))
        output = thisPool.map(TryAlphaFitsFun.mapper,inputparams)
        thisPool.close()
        thisPool.join()
    else:
        print 'Running Single core'
        output = []
        for icount,iin in enumerate(inputparams):
            output.append(TryAlphaFitsFun(*iin))
            print int((icount*100)/float(len(inputparams))) , '% done' + ' '*50 + '\r',
else:
    print 'nothing to calculate'        
    output = []
    
# WipeSet(outputdir[0],RunGammaList,feedin['set'],filepref='Fits/')
print 'Done Fits ' + ' '*50
#FitData = { iset , ip , iflow , ifit  , Boot/Avg/Chi }
topfitdir = outputdir[0] + 'Top/Alpha/Fits'
mkdir_p(topfitdir)

for iout in output:
    thisset = iout[1][0]
    PrintAlphaFitFile(iout[0][thisset],thisset,topfitdir)

print 'Total fit time took: ' , str(datetime.timedelta(seconds=time.time()-totstart)) , ' h:m:s '

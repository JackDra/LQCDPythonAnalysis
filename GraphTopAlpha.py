#!/usr/bin/env python

from Params import *
import numpy as np
from ReadXml import *
from MiscFuns import *
from FitParams import *
from StringFix import *
from SetLists import *
from OppFuns import CreateGammaList
import time,datetime
from MultiWrap import *
from multiprocessing import Pool
from InputArgs import *
from CreateCombs import CreateDictOldCombs
from CombParams import *
import cPickle as pickle


##datadict = { gamma } { mom } { method } { set }
##thisdatadict = { method } { set }
##thisGammaList format: D/S P4/3 Opp (e.g. doubP4g4)
##thisMomList in qstr format (e.g. 'q = 0 0 0')
##NB q = 0 0 0 MUST BE INCLUDED FOR TSF ETC



def ReadAndPlotDis(thisSetList,thisMomList):
    # thisAllSetList = thisSmearList+thisSetList
    # for isetlist,dump in thisSetPoFLists[:len(thisSetPoFLists)/2]:
    #     thisAllSetList += isetlist
    iterSetList = SortMySet(ReduceTooMassSet(thisSetList))[0]
    for imom in thisMomList:
        datadict = ReadTopFile(outputdir[0],iset,thisMomList=thisMomList)
        setstart = time.time()
        PlotTopSetCharge(datadict,iterSetList,imom,feedin['ForceTitle'])
        print 'Getting and plotting ' , imom, 'Took: ' , str(datetime.timedelta(seconds=(time.time()-setstart))) ,' h:m:s                      '


# if all('-m' not in iin for iin in sys.argv[1:]):
#     feedin = InputParams(sys.argv[1:] + ['-noprompt'] + ['-m=Fits,TSFTsink,TSFtest32,OSF'])
# else:
feedin = InputParams(sys.argv[1:] + ['-noprompt'])
    


feedin['set'] = ReduceTooMassSet(feedin['set'])
ShowSetLists(feedin['set'])
feedin['mom'] = GetAvgMomList(feedin['mom'])

if DoMulticore and len(feedin['set']) > 1  and feedin['anaproc'] > 1:
    # inputparams = [([iset],feedin['mom'],feedin['method']) for iset in feedin['set']]
    inputparams = []
    for imom in feedin['mom']:
        inputparams.append((feedin['set'],[imom]))
    makeContextFunctions(ReadAndPlotDis)
    thisPool = Pool(min(len(inputparams),feedin['anaproc']))
    thisPool.map(ReadAndPlotDis.mapper,inputparams)
    thisPool.close()
    thisPool.join()
else:
    ReadAndPlotDis(feedin['set'],feedin['mom'])


    

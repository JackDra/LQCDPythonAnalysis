#!/usr/bin/env python

from Params import *
import numpy as np
from ReadTxt import *
from MiscFuns import *
from FitParams import *
from StringFix import *
from GraphDataNew import *
from SetLists import *
from OppFuns import CreateGammaList
import time,datetime
from MultiWrap import *
from multiprocessing import Pool
from InputArgs import *
from CreateCombs import CreateDictOldCombs
from CombParams import *

DoDS = True
feedin = InputParams(sys.argv[1:] + ['-noprompt'])
    

##datadict = { gamma } { mom } { method } { set }
##thisdatadict = { method } { set }
##thisGammaList format: D/S P4/3 Opp (e.g. doubP4g4)
##thisMomList in qstr format (e.g. 'q = 0 0 0')
##NB q = 0 0 0 MUST BE INCLUDED FOR TSF ETC
def progprint(numb,starttime,igamma):
    print 'Graphing Operator: ' + igamma , int(numb*100/float(9)) , '% time taken:' , str(datetime.timedelta(seconds=time.time()-starttime)) ,' h:m:s          '


def ReadAndPlotDis(thisSetList,thisMomList,thisMethodList):
    # thisAllSetList = thisSmearList+thisSetList
    # for isetlist,dump in thisSetPoFLists[:len(thisSetPoFLists)/2]:
    #     thisAllSetList += isetlist
    datadictlist = []
    thiskappalist = []
    if 'RF' in thisMethodList: thisMethodList.remove('RF')
    for ioutput,ikappa in zip(outputdir,kappalist):
        hold = ReadSetFitRFDict(ioutput,thisSetList,['twopt'],thisMethodList,thisMomList=thisMomList)
        if CheckDict(hold,'twopt','q = 0 0 0'):
            datadictlist.append(hold)
            # thiskappalist.append(ikappa.replace('k'+str(ikappa),''))
            thiskappalist.append(ikappa)
    datadict = datadictlist[0]
    iterSetList = SortMySet(ReduceTooMassSet(thisSetList))[0]
    for iset in iterSetList:
        setstart = time.time()
        PlotDispersion(datadict['twopt'],iset,feedin['ForceTitle'])
        print 'Plotting ' , iset, 'Took: ' , str(datetime.timedelta(seconds=(time.time()-setstart))) ,' h:m:s                      '
    
                    


# if all('-m' not in iin for iin in sys.argv[1:]):
#     feedin = InputParams(sys.argv[1:] + ['-noprompt'] + ['-m=Fits,TSFTsink,TSFtest32,OSF'])
# else:
feedin = InputParams(sys.argv[1:] + ['-noprompt'])
    
thisGammaList = CreateGammaList(feedin['gamma'],twopt=True)


feedin['set'] = ReduceTooMassSet(feedin['set'])
ShowSetLists(feedin['set'])
ShowMethodList(feedin['method'])
feedin['mom'] = GetAvgMomList(feedin['mom'])

if DoMulticore and len(feedin['set']) > 1  and feedin['anaproc'] > 1:
    inputparams = [([iset],feedin['mom'],feedin['method']) for iset in feedin['set']]
    makeContextFunctions(ReadAndPlotDis)
    thisPool = Pool(min(len(inputparams),feedin['anaproc']))
    output = thisPool.map(ReadAndPlotDis.mapper,inputparams)
    thisPool.close()
    thisPool.join()
else:
    ReadAndPlotDis(feedin['set'],feedin['mom'],feedin['method'])
        
print 'Graphing all complete'
    

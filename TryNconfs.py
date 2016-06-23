#!/usr/bin/env python


from OppFuns import CreateGammaList
from Params import *
import sys
from InputArgs import *
from CheckXml import *
from OutputXmlData import WriteXml
from SetLists import ReduceTooMassSet

minmax = 'max'
if len(sys.argv) > 1:
    if sys.argv[1] == 'min' or sys.argv[1] == 'max':
        minmax = sys.argv[1]
        sys.argv = [sys.argv[0]] + sys.argv[2:]
        
feedin = InputParams(sys.argv[1:] + ['-noprompt'])
thisGammaList = CreateGammaList(feedin['gamma'])

ShowSetLists(feedin['set'])

ShowMethodList(feedin['method'])

def PickSetForMethod(thismethod,thisSetList):
    if 'Tsink' in thismethod:
        outSetList = []
        for itsink in AllTSinkStrList:
            if itsink+SingSmList[0] in thisSetList:
                outSetList.append(itsink+SingSmList[0])
    ##HERE NOT DONE##
    if 'TSF' in thismethod:
        outSetList = ReduceTsink(outSetList,NoPoF=True)
    elif 'SumMeth' in thismethod:
        outSetList = SingSmList
    else:
        outSetList = thisSetList
    

for imethod in feedin['method']:
    PickSetForMethod(thismethod,feedin['set']
    for iset in thisSetList:
        nconf,nconfDict = CheckNconf(thisGammaList,[iset],thisMomList=feedin['mom'],CheckList=[imethod],minmax=minmax)
        if not isinstance(nconf, str):
            print imethod, iset , 'Nconfs =', nconf
        else:
            print imethod, iset , nconf
        mkdir_p(outputdir+'/Nconfs')
        if len(nconfDict.keys()) > 0: WriteXml(outputdir+'/Nconfs/'+imethod+iset,{'Data':nconfDict})

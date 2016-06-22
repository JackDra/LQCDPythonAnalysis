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
RedSetList = ReduceTooMassSet(feedin['set'])

for imethod in feedin['method']:
    if 'TSF' in imethod:
        thisSetList = RedSetList
    if 'SumMeth' in imethod:
        thisSetList = SingSmList
    else:
        thisSetList = feedin['set']
    for iset in thisSetList:
        nconf,nconfDict = CheckNconf(thisGammaList,[iset],thisMomList=feedin['mom'],CheckList=[imethod],minmax=minmax)
        if not isinstance(nconf, str):
            if nconf > 10e9 or nconf < 0:
                print imethod, iset ,'No nconfs in file (depreceated code)'
            else:
                print imethod, iset , 'Nconfs=', nconf
        else:
            print nconf
        mkdir_p(outputdir+'/Nconfs')
        if len(nconfDict.keys()) > 0: WriteXml(outputdir+'/Nconfs/'+imethod+iset,{'Data':nconfDict})

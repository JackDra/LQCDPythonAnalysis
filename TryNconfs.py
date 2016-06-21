#!/usr/bin/env python


from OppFuns import CreateGammaList
from Params import *
import sys
from InputArgs import *
from CheckXml import *
from OutputXmlData import WriteXml



feedin = InputParams(sys.argv[1:] + ['-noprompt'])
thisGammaList = CreateGammaList(feedin['gamma'])

ShowSetLists(feedin['set'])

ShowMethodList(feedin['method'])

for imethod in feedin['method']:
    for iset in feedin['set']:
        nconf,nconfDict = CheckNconf(thisGammaList,[iset],thisMomList=feedin['mom'],CheckList=[imethod])
        if not isinstance(nconf, str):
            if nconf > 10e9:
                print imethod, iset ,'No nconfs in file (depreceated code)'
            else:
                print imethod, iset , 'Nconfs=', nconf
        else:
            print nconf
        mkdir_p(outputdir+'/Nconfs')
        if len(nconfDict.keys()) > 0: WriteXml(outputdir+'/Nconfs/'+imethod+iset,{'Data':nconfDict})

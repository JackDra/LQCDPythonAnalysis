#!/usr/bin/env python


from OppFuns import CreateGammaList
from Params import *
import sys
from InputArgs import *
from CheckXml import *
from OutputXmlData import WriteXml
from SetLists import PickSetForMethod

minmax = 'max'
if len(sys.argv) > 1:
    if sys.argv[1] == 'min' or sys.argv[1] == 'max':
        minmax = sys.argv[1]
        sys.argv = [sys.argv[0]] + sys.argv[2:]
        
feedin = InputParams(sys.argv[1:] + ['-noprompt'])
thisGammaList = CreateGammaList(feedin['gamma'])
ShowSetLists(feedin['set'])

ShowMethodList(feedin['method'])

thisNconfList = []
mkdir_p(outputdir+'/Nconfs')
for imethod in feedin['method']:
    thisSetList = PickSetForMethod(imethod,feedin['set'])
    for iset in thisSetList:
        nconf,nconfDict = CheckNconf(thisGammaList,[iset],thisMomList=feedin['mom'],CheckList=[imethod],minmax=minmax)
        if not isinstance(nconf, str):
            print imethod, iset , 'Nconfs =', nconf
            thisNconfList.append(nconf)
        else:
            print imethod, iset , nconf
        if len(nconfDict.keys()) > 0: WriteXml(outputdir+'/Nconfs/'+imethod+iset,{'Data':nconfDict})

if len(thisNconfList) > 0:
    TotNconf = min(thisNconfList)
else:
    TotNconf = False
with open(scriptdir+'/setup.cfg','w') as f:
    f.write('\nscriptdir:\n')
    f.write(scriptdir+'/\n')
    f.write('\ndatadir:\n')
    f.write(datadir+'\n')
    f.write('\nAnaProc:\n')
    f.write(str(AnaProc)+'\n')
    f.write('\nListOrSet:\n')
    f.write(ListOrSet+'\n')
    f.write('\nkappa:\n')
    f.write(str(kappa)+'\n')
    f.write('\nPoFShifts:\n')
    f.write(str(PoFShifts)+'\n')
    f.write('\nNconfs:\n')
    f.write(str(TotNconf)+'\n')

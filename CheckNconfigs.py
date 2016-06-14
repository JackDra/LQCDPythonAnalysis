#!/usr/bin/env python


from Params import *
from InputArgs import *
from CheckXml import CheckNconf
from OppFuns import CreateGammaList


feedin = InputParams(sys.argv[1:])
DefWipeWarning()
thisGammaList = CreateGammaList(feedin['gamma'])

ShowSetLists(feedin['set'])

Nconf = CheckNconf(thisGammaList,feedin['set'],thisMomList=feedin['mom'],CheckList=feedin['method'],cfuns=False)

print Nconf

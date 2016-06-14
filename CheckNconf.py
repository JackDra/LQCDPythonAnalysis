#!/usr/bin/env python


from InputArgs import *
from ReadXml import CheckNconfFile
from OppFuns import CreateGammaList


feedin = InputParams(sys.argv[1:])
DefWipeWarning()
thisGammaList = CreateGammaList(feedin['gamma'])

ShowSetLists(feedin['set'])

Nconf = CheckNconfFile(thisGammaList,feedin['set'],thisMomList=feedin['mom'],CheckList=feedin['method'],cfuns=False)

print Nconf

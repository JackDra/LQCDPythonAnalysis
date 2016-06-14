#!/usr/bin/env python


from OppFuns import CreateGammaList,WipeSet
from Params import *
import numpy as np
from BootTest import BootStrap1
from ReadTxt import *
from SetLists import *
from Fitting import FitRFSet
from OutputData import PrintFitSetToFile
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



feedin = InputParams(sys.argv[1:])
DefWipeWarning()
thisGammaList = CreateGammaList(feedin['gamma'])

ShowSetLists(feedin['set'])

ShowMethodList(feedin['method'])

for imethod in feedin['method']:
    for iset in feedin['set']:
        nconf = CheckNconf(thisGammaList,[iset],thisMomList=feedin['mom'],CheckList=[imethod])
        if not isinstance(nconf, str):
            if nconf > 10e9:
                print 'No nconfs in file (depreceated code)'
            else:
                print imethod, iset , 'Smallest number of configurations used is', nconf
        else:
            print nconf

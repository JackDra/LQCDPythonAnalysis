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

nconf = CheckNconf(thisGammaList,feedin['set'],thisMomList=feedin['mom'],CheckList=feedin['method'])

print nconf

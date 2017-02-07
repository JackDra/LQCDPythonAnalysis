#!/usr/bin/env python

from Params import *
import numpy as np
from ReadXml import *
from MiscFuns import *
from FitParams import *
from StringFix import *
from SetLists import *
import time,datetime
from MultiWrap import *
from multiprocessing import Pool
from InputArgs import *
from CreateCombs import CreateDictOldCombs
from CombParams import *
from GraphDataNew import *

## cfglistout [ icfg ]
## tflow [ icfg , itflow ]
## topcharge [ icfg , itflow ]
kappatopc,kappatflow = [],[]
for ikappa in kappalist:    
    print 'Reading kappa=', ikappa
    if Debug: print WeinDir.replace('Kud0'+str(kappa),ikappa.replace('k','Kud0'))
    filelist,topcharge,tflow = ReadTopAll(WeinDir.replace('Kud0'+str(kappa),ikappa.replace('k','Kud0')))
    if ikappa == kappa:
        GraphWExp(topcharge,tflow[0])
        GraphWLines(topcharge,tflow[0],np.arange(0,200,33))
        GraphWchit(topcharge,tflow[0])
        GraphW2Hist(topcharge,tflow[0].tolist().index(3.01))


    kappatopc.append(topcharge)
    kappatflow.append(tflow[0])
    

# GraphWchitKappas(kappatopc,kappatflow)
    

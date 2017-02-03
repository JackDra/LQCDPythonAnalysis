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
    filelist,topcharge,tflow = ReadTopAll(TCDir.replace('k'+str(kappa),ikappa))
    if ikappa == kappa:
        GraphQExp(topcharge,tflow[0])
        GraphQLines(topcharge,tflow[0],np.arange(0,200,33))
        Graphchit(topcharge,tflow[0])
        GraphQ2Hist(topcharge,tflow[0].tolist().index(3.01))


    kappatopc.append(topcharge)
    kappatflow.append(tflow[0])
    

GraphchitKappas(kappatopc,kappatflow)
    

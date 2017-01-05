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
filelist,topcharge,tflow = ReadTopAll(TCDir)

GraphQExp(topcharge,tflow[0])
GraphQLines(topcharge,tflow[0],np.arange(0,200,33))
Graphchit(topcharge,tflow[0])
GraphQ2Hist(topcharge,tflow[0].tolist().index(3.01))

    

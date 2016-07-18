#!/usr/bin/env python

import numpy as np
from Params import *
from GraphRapNew import SummaryFeed
import sys

thisSmearList = ['32','64','128']
thisStateList = ['1']
thisTSinkList = ['26','29','32','35','38']
thisTvarList = ['to20dt2']
thisGammaList = DefGammaList

SummaryFeed(thisGammaList,thisSmearList,thisStateList,thisTvarList,thisTSinkList,'Summary')
    



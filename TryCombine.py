#!/usr/bin/env python

from Params import *
from CreateCombs import *
from InputArgs import *
from OppFuns import *




feedin = InputParams(sys.argv[1:])
thisGammaList = CreateGammaList(feedin['gamma'])

ShowSetLists(feedin['set'])

CombType = raw_input("What Combination do you want? (Proton, Neutron, IsoVector) \n")
if CombType not in ['Proton','Neutron','IsoVector']:
    raise IOError("choose Proton, Neutron, IsoVector")

if CombType == 'Proton':
    ReadAndComb(feedin,FFProton,CombType)
if CombType == 'Neutron':
    ReadAndComb(feedin,FFNeutron,CombType)
if CombType == 'IsoVector':
    ReadAndComb(feedin,IsoVector,'')

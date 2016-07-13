#!/usr/bin/env python

from Params import *
from CreateCombs import *





feedin = InputParams(sys.argv[1:])
thisGammaList = CreateGammaList(feedin['gamma'])

ShowSetLists(feedin['set'])

CombType = raw_input("What Combination do you want? (Proton, Neutron, IsoVector)")
if CombType not in ['Proton','Neutron','IsoVector']:
    raise IOError("choose Proton, Neutron, IsoVector")

if CombType == 'Proton':
    ReadAndComb(inputargs,FFProton,CombType)
if CombType == 'Neutron':
    ReadAndComb(inputargs,FFNeutron,CombType)
if CombType == 'IsoVector':
    ReadAndComb(inputargs,IsoVector,'')

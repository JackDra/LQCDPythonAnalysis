#!/usr/bin/env python
from array import array
from os import walk
import os.path
import numpy as np
from Params import *
from ReadBinaryCfuns import ReadFSCfunPick, Read2ptCfunPick
from CfunBoot import BootSet3pt, BootSet2pt,BootSet3ptCM, BootSet2ptCM
from RFCalc import CalcRatioFactor, CalcRatioFactorCM
from ReadCMCfuns import ReadCMSet,ReadCMList


dir = "/raid/jdragos/scratch/cfun/2ndk12090/"
group3pt = ".t16sm32Xsm32GMA4t29p000.doubk12090.ifms.3cf"
group2pt = ".t16sm32k12090.ifmssi32.nucleon.u.2cf"

print "qlow to qhigh is: " , ipTOiq(0), ipTOiq(nmom-1)

### Read specific list ###
data2pt = []
data3pt = []
for iconf in conflist:
    file = filelist.replace('*',iconf)
    # print file
    data2pt.append(Read2ptCfunPick(file.replace('@','twoptsm32si32')+group2pt,[iqTOip(0)]))
    data3pt.append(ReadFSCfunPick(file.replace('@','cmsm32Xsm32GMA4t29p000.doub')+group3pt,[iqTOip(0)],['g4']))

### Read whole set ###
# for (dirname,dirs,files) in walk(dir+'source7/'):
#     for file in files:
#         if group3pt not in file or ".metadata" in file: continue    
#         if not os.path.isfile(dir+'source7/'+file.replace(group3pt,group2pt)): continue
#         print file
#         data3pt.append(ReadFSCfun(dir+'source7/'+file))
#         data2pt.append(Read2ptCfun(dir+'source7/'+file.replace(group3pt,group2pt)))
###   ###


ncon = np.size(data3pt)
print 'ncon = ' + str(ncon)

print data3pt

[RFr,SqrtFac] = CalcRatioFactorCM(bootdata2pt['real'] ,bootdata3pt['real'])


# # RFi = CalcRatioFactor(bootdata2pt['real'] ,bootdata3pt['cmplx'])
# [RFr,SqrtFac] = CalcRatioFactor(bootdata2pt['real'] ,bootdata3pt['real'])

for ictsink,itsink in enumerate(TSinkSet):
    for icsm,ism in enumerate(SmearSet):
        with open('results_t'+itsink+'sm'+ism+'.txt', 'w') as f:
            for it in range(int(itsink)-tsource+1):
                # RFr[4][iqTOip(0)][it].Stats()
                f.write( '{0}{1:20.10f}{2:20.10f}'.format(repr(it+tsource).rjust(3), RFr[ictsink][4][iqTOip(0)][icsm][it].Avg, RFr[ictsink][4][iqTOip(0)][icsm][it].Std) + '\n' )


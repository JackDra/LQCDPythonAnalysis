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
# for iconf in conflist:
#     file = filelist.replace('*',iconf)
#     print file
#     data3pt.append(ReadFSCfun(file))
#     data2pt.append(Read2ptCfun(file.replace(group3pt,group2pt)))

### Read whole set ###
# for (dirname,dirs,files) in walk(dir+'source7/'):
#     for file in files:
#         if group3pt not in file or ".metadata" in file: continue    
#         if not os.path.isfile(dir+'source7/'+file.replace(group3pt,group2pt)): continue
#         print file
#         data3pt.append(ReadFSCfun(dir+'source7/'+file))
#         data2pt.append(Read2ptCfun(dir+'source7/'+file.replace(group3pt,group2pt)))
###   ###

# ### Read CM Set ###
# [cmdata2pt,cmdata3pt,filelist] = ReadCMSet('doub','4','',dir)
# ###   ###


### Read CM Set ###
[cmdata2pt,cmdata3pt,filelist] = ReadCMList('doub','4','',conflist)
###   ###

ncon = np.size(data3pt)
print 'ncon = ' + str(ncon)

# bootdata3pt = BootSet3pt(data3pt,nboot)
# bootdata2pt = BootSet2pt(data2pt,nboot)
##bootdata3pt { 'real'/'cmplx' } [ gamma , ip , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
##bootdata2pt { 'real'/'cmplx' } [ ip , it ]  = bootstrap1 class (.Avg, .Std, .values, .nboot)


bootdata3ptCM = BootSet3ptCM(data3pt,nboot)
bootdata2ptCM = BootSet2ptCM(data2pt,nboot)
##bootdata3ptCM { 'real'/'cmplx' } [ itsink , gamma , ip , ism , jsm, it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
##bootdata2ptCM { 'real'/'cmplx' } [ ip , ism , jsm , it ]  = bootstrap1 class (.Avg, .Std, .values, .nboot)

# RFi = CalcRatioFactor(bootdata2pt['real'] ,bootdata3pt['cmplx'])
[RFr,SqrtFac] = CalcRatioFactorCM(bootdata2pt['real'] ,bootdata3pt['real'])


# # RFi = CalcRatioFactor(bootdata2pt['real'] ,bootdata3pt['cmplx'])
# [RFr,SqrtFac] = CalcRatioFactor(bootdata2pt['real'] ,bootdata3pt['real'])

for ictsink,itsink in enumerate(TSinkSet):
    for icsm,ism in enumerate(SmearSet):
        with open('results_t'+itsink+'sm'+ism+'.txt', 'w') as f:
            for it in range(int(itsink)-tsource+1):
                # RFr[4][iqTOip(0)][it].Stats()
                f.write( '{0}{1:20.10f}{2:20.10f}'.format(repr(it+tsource).rjust(3), RFr[ictsink][4][iqTOip(0)][icsm][it].Avg, RFr[ictsink][4][iqTOip(0)][icsm][it].Std) + '\n' )


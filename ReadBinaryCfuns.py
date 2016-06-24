#!/usr/bin/env python
from array import array
import numpy as np
from Params import *

SeekIncSize = 8

## Der[igamma , ip , it]
class ReadFSDerCfunPick:
    def __init__(self,file,thisMomList,thisDGList):
        self.data = []
        f = open(file,'rb')
        for igd,thisgamma in enumerate(thisDGList):
            cmplxflag = False
            if 'cmplx' in thisgamma: 
                cmplxflag = True
                thisgamma = thisgamma.replace('cmplx','')
            igdloc = DGSet.index(thisgamma)
            self.data.append([])
            for ip,iploc in enumerate(thisMomList):         
                self.data[igd].append([])
                for it in range(nt):
                    loc = 2*SeekIncSize*(igdloc + iploc*DGBSize + it*DGBSize*nmom)
                    if cmplxflag: loc += SeekIncSize
                    f.seek(loc)
                    tmpdata = array('d')
                    tmpdata.read(f,1)
                    tmpdata.byteswap()
                    self.data[igd][ip].append(tmpdata[0])
                    if np.isnan(self.data[igd][ip][it]):
                        raise NaNCfunError('NaN Values: '+thisgamma+' ' +qvecSet[iploc] + ' it='+str(it+1) )
        f.close()

## noDer[igamma , ip , it]
class ReadFSCfunPick:
    def __init__(self,file,thisMomList,thisGammaList):
        self.data = []
        f = open(file,'rb')
        for igamma,thisgamma in enumerate(thisGammaList):
            cmplxflag = False
            if 'cmplx' in thisgamma: 
                cmplxflag = True
                thisgamma = thisgamma.replace('cmplx','')
            igammaloc = GammaSet.index(thisgamma)
            self.data.append([])
            for ip,iploc in enumerate(thisMomList):      
                self.data[igamma].append([])
                for it in range(nt):
                    loc = 2*SeekIncSize*(igammaloc + iploc*GBSize + it*GBSize*nmom)
                    if cmplxflag: loc += SeekIncSize
                    f.seek(loc)
                    tmpdata = array('d')
                    tmpdata.read(f,1)
                    tmpdata.byteswap()
                    self.data[igamma][ip].append(tmpdata[0])
                    if np.isnan(self.data[igamma][ip][it]):
                        if Debug: self.data[igamma][ip][it] = 0.0
                        raise NaNCfunError('NaN Values: '+thisgamma+ ' ' +qvecSet[iploc] + ' it='+str(it+1) )
        f.close()



class Read2ptCfunPick:
    def __init__(self,file,thisMomList):
        self.data = []
        f = open(file,'rb')
        for ip,iploc in enumerate(thisMomList):
            self.data.append([])
            for it in range(nt):
        ##11 + 22 projected spin component ns = 4 works
                loc = SeekIncSize*2*(ns**2)*(it + iploc*nt)
                f.seek(loc)
                tmpdata = array('d')
                tmpdata.read(f,1)
                tmpdata.byteswap()
                val = tmpdata[0]
                f.seek(loc+10*SeekIncSize)
                tmpdata = array('d')
                tmpdata.read(f,1)
                tmpdata.byteswap()
                val += tmpdata[0]
                self.data[ip].append(val)
                if np.isnan(self.data[ip][it]):
                    raise NaNCfunError('NaN Values: ' +qvecSet[iploc] + ' it='+str(it+1) )
                # if 'cmplx' in RI:
                #     self.datai[ipread].append(tmpdata[loc+1]+tmpdata[loc+11])
            if self.data[ip][tsource-1] < 0.0:
                self.data[ip] = np.negative(self.data[ip])
        f.close()
 
class NaNCfunError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)












# ## Der[igamma , ip , it]
# class ReadFSDerCfunPick:
#     def __init__(self,file,thisMomList,thisDGList):
#         self.data = []
#         f = open(file,'rb')
#         tmpdata = array('d')
#         tmpdata.read(f,DGBSize*nmom*nt*2)
#         tmpdata.byteswap()
#         for igd,thisgamma in enumerate(thisDGList):
#             cmplxflag = False
#             if 'cmplx' in thisgamma: 
#                 cmplxflag = True
#                 thisgamma = thisgamma.replace('cmplx','')
#             igdloc = DGSet.index(thisgamma)
#             self.data.append([])
#             for ip,iploc in enumerate(thisMomList):         
#                 self.data[igd].append([])
#                 for it in range(nt):
#                     loc = igdloc*2 + iploc*DGBSize*2 + it*DGBSize*nmom*2
#                     if cmplxflag: loc += 1
#                     self.data[igd][ip].append(tmpdata[loc])
#                     if np.isnan(self.data[igd][ip][it]):
#                         raise NaNCfunError('NaN Values: '+thisgamma+' ' +qvecSet[iploc] + ' it='+str(it+1) )

#         f.close()

# ## noDer[igamma , ip , it]
# class ReadFSCfunPick:
#     def __init__(self,file,thisMomList,thisGammaList):
#         self.data = []
#         f = open(file,'rb')
#         tmpdata = array('d')
#         try:
#             tmpdata.read(f,GBSize*nmom*nt*2)
#         except:
#             print file
#         tmpdata.byteswap()
#         for igamma,thisgamma in enumerate(thisGammaList):
#             cmplxflag = False
#             if 'cmplx' in thisgamma: 
#                 cmplxflag = True
#                 thisgamma = thisgamma.replace('cmplx','')
#             igammaloc = GammaSet.index(thisgamma)
#             self.data.append([])
#             for ip,iploc in enumerate(thisMomList):      
#                 self.data[igamma].append([])
#                 for it in range(nt):
#                     loc = igammaloc*2 + iploc*GBSize*2 + it*GBSize*nmom*2
#                     if cmplxflag: loc += 1
#                     self.data[igamma][ip].append(tmpdata[loc])
#                     if np.isnan(self.data[igamma][ip][it]):
#                         raise NaNCfunError('NaN Values: '+thisgamma+ ' ' +qvecSet[iploc] + ' it='+str(it+1) )

#         f.close()

# ##[ip,it]
# class Read2ptCfunPickOld:
#     def __init__(self,file,thisMomList):
#         self.data = []
#         f = open(file,'rb')
#         tmpdata = array('d')
#         tmpdata.read(f,nmom*ns*ns*nt*2)
#         tmpdata.byteswap()
#         for ip,iploc in enumerate(thisMomList):
#             self.data.append([])
#             for it in range(nt):
#         ##11 + 22 projected spin component ns = 4 works
#                 loc = it*2*ns**2 + iploc*nt*2*ns**2
#                 self.data[ip].append(tmpdata[loc]+tmpdata[loc+10])
#                 if np.isnan(self.data[ip][it]):
#                     raise NaNCfunError('NaN Values: ' +qvecSet[iploc] + ' it='+str(it+1) )
#                 # if 'cmplx' in RI:
#                 #     self.datai[ipread].append(tmpdata[loc+1]+tmpdata[loc+11])
#             if self.data[ip][tsource-1] < 0.0:
#                 self.data[ip] = np.negative(self.data[ip])
#         f.close()

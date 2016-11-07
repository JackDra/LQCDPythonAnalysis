#!/usr/bin/env python

from Params import *
import glob,os
from array import array

SeekIncSize = 8

# data [it]

def ReadCfun(thisfile):
    data = []
    with open(thisfile,'rb') as f:
        ip,iploc = 0,qstrTOip('q = 0 0 0')
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
            data.append(val)
            if np.isnan(data[it]):
                raise NaNCfunError('NaN Values: ' +qvecSet[iploc] + ' it='+str(it+1) )
            # if 'cmplx' in RI:
            #     datai[ipread].append(tmpdata[loc+1]+tmpdata[loc+11])
        if data[tsource-1] < 0.0:
            data = np.negative(data)
    return data

def WriteCfun(thisfile,data):
    with open(thisfile,'wb') as f:
        for tdata in data:
            f.write(bin(tdata))

srclist = ['1','2','3','7','8','9']

cfundir = '/raid/jdragos/data/cfuns/k12090/source@/twoptsm32si32/'

cfundirlist = [cfundir.replace('@',isrc) for isrc in srclist]
cfundirout = '/raid/jdragos/data/cfuns/Yousifah/'

totfilelist = []
outfilelist = []
for idir in cfundirlist:
    print 'Reading: ', idir
    totfilelist += glob.glob(idir+'*.665.*t16*.2cf')
    outfilelist += [ifile.replace(idir,cfundirout) for ifile in totfilelist]

print 'files found:'
for ifile in totfilelist:
    print ifile

for outfile,ifile in zip(totfilelist,outfilelist):
    if not os.path.isfile(ifile): continue
    idata = ReadCfun(ifile)
    WriteCfun(outfile,idata)
    
        

#!/usr/bin/env python

from Params import *
import glob,os,re
from array import array
import numpy as np
from ReadBinaryCfuns import *
from collections import OrderedDict

def R2CCOut(data,thisfile):
    with open(thisfile,'w') as f:
        for it,idata in enumerate(data):
            '{:5d} {:10.20e}  \n'.format(it,idata)


cfundir = '/mnt/research/lqcd/CfunAnalysis/cfunPChroma/Kud01370000Ks01364000/twoptRandT/twoptsm64si64/'
cfunOutdir = '/mnt/research/lqcd/CfunAnalysis/cfunForAhmed/Kud01370000Ks01364000/twoptRandT/twoptsm64si64/'
cfunOutdirg5 = '/mnt/research/lqcd/CfunAnalysis/cfunForAhmed/Kud01370000Ks01364000/twoptRandT/twoptsm64si64g5/'

mkdir_p(cfunOutdir)
mkdir_p(cfunOutdirg5)

Totfilelist = glob.glob(cfundir+'*')

print 'Found ' , len(Totfilelist) , ' configurations'

FileListIn = []

def SplitCfgSrc(ifile):
    isrc = re.search('xsrc.*_k',ifile).group().replace('_k','')
    icfg = re.search('-.-.*_xsrc',ifile).group().replace('_xsrc','')
    return isrc,icfg

FileListIn,FileListOut,FileListOut = OrderedDict(),[],[]
for ifile in Totfilelist:
    isrc,icfg = SplitCfgSrc(ifile)
    if icfg not in FileListIn.keys():        
        FileListIn[icfg] = []
        FileListOut.append(ifile.replace(cfundir,cfunOutdir).replace(isrc,''))
        FileListOutg5.append(ifile.replace(cfundir,cfunOutdirg5).replace(isrc,''))
    FileListIn[icfg].append(ifile)
    


for cfgFileOut,cfgFileOutg5,(icfg,cfgFileListIn) in zip(FileListOut,FileListOutg5,FileListIn.iteritems()):
    print 'Reading and Writing Config' , icfg
    data = R2CChromaXMLFileList(cfgFileListIn,[0],Dog5=True)
    R2CCOut(data.data[0],cfgFileOut)
    R2CCOut(data.datag5[0],cfgFileOutg5)


    
### SeekIncSize = 8

### # data [it]

### def ReadCfun(thisfile):
###     data = []
###     f = open(thisfile,'rb')
###     # with open(thisfile,'rb') as f:
###     ip,iploc = 0,qstrTOip('q = 0 0 0')
###     for it in range(nt):
###         ##11 + 22 projected spin component ns = 4 works
###         try:
###             loc = SeekIncSize*2*(ns**2)*(it + iploc*nt)
###             f.seek(loc)
###             tmpdata = array('d')
###             tmpdata.read(f,1)
###             tmpdata.byteswap()
###             val = tmpdata[0]
###             f.seek(loc+10*SeekIncSize)
###             tmpdata = array('d')
###             tmpdata.read(f,1)
###             tmpdata.byteswap()
###             val += tmpdata[0]
###             data.append(val)
###         except:
###             return False
###         if np.isnan(data[it]):
###             raise NaNCfunError('NaN Values: ' +qvecSet[iploc] + ' it='+str(it+1) )
###         # if 'cmplx' in RI:
###         #     datai[ipread].append(tmpdata[loc+1]+tmpdata[loc+11])
###     if data[tsource-1] < 0.0:
###         data = np.negative(data)
###     f.close()
###     return data

### def WriteCfun(thisfile,data):
###     np.array(data).tofile(thisfile)
###     # with open(thisfile,'wb') as f:
###     #     for tdata in data:
###     #         f.write(tdata)

### srclist = ['1','2','3','6','7','8']
### # srclist = ['1','3','7','8','9']


### cfundir = '/raid/jdragos/data/cfuns/k12090/source@/twoptsm32si32/'

### cfundirlist = [cfundir.replace('@',isrc) for isrc in srclist]
### cfundirout = '/raid/jdragos/data/cfuns/Yousifah/'

### totfilelist = []
### outfilelist = []
### for idir in cfundirlist:
###     print 'Reading: ', idir
###     thisfilelist = glob.glob(idir+'*.665.*t15*.2cf')
###     # thisfilelist = glob.glob(idir+'*kp120900.0*t16*.2cf')
###     totfilelist += thisfilelist
###     outfilelist += [ifile.replace(idir,cfundirout) for ifile in thisfilelist]

### print 'files IO:'

### for ifile,outfile in zip(totfilelist,outfilelist):
###     print ifile
###     if not os.path.isfile(ifile): continue
###     idata = ReadCfun(ifile)
###     if idata != False:
###         print outfile
###         print ''
###         WriteCfun(outfile,idata)
    
        

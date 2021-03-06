#!/usr/bin/env python
from array import array
import numpy as np
import BootTest as bt
from Params import *
from MiscFuns import *
from ReadBinaryCfuns import *

##NB CreateBoot take array [ iconf , it ]

#data = [ iconf , igamma , ip , it ]
#dataout = [igamma , ip , it ] bs
def BootSet3pt(data,thisMomList,thisGammaList,nboot,printstr='',randlist=[]):
    dataout = []
    for ig,igamma in enumerate(thisGammaList):
        dataout.append([])
        for ip,imom in enumerate(thisMomList):
            # print 'Booting '+printstr+igamma+' {}%               \r'.format(int((ip*100)/float(len(thisMomList)))),
            bootdata,randlist = bt.CreateBoot(np.array(data)[:,ig,ip,:],nboot,0,randlist=randlist)
            dataout[ig].append(bootdata)
    # print '                             \r',
    return dataout


#data = [ iconf , ip , it]
#dataout = [ ip , it ]. bs
def BootSet2pt(data,thisMomList,nboot,randlist=[]):
    dataout = []
    randlist = []
    for ip,imom in enumerate(thisMomList):
        # print 'Booting {}%  \r'.format(int((ip*100)/float(len(thisMomList)))),
        bootdata,randlist = bt.CreateBoot(np.array(data)[:,ip,:],nboot,0,randlist=randlist)
        dataout.append(bootdata)
    # print '                              \r',
    return dataout,randlist

#dataout = [ ip , it ]. bs
def ReadAndBoot2pt(readfilelist,thisMomList,thisnboot,randlist=[]):
    tempdata = []
    for iconf,ifile in enumerate(readfilelist):
        # print 'Reading {}%  \r'.format(int((iconf*100)/float(len(readfilelist)))),
        try:
            tempdata.append(Read2ptCfunPick(ifile,thisMomList).data)
        except NaNCfunError as e:
            print 'Deleting file ' + ifile
            print e
            os.remove(ifile)
    return BootSet2pt(tempdata,thisMomList,thisnboot,randlist=randlist)


#dataout [ igamma , ip , it ]. bs
def ReadAndBoot3pt(readfilelist,thisMomList,thisGammaList,thisDerGammaList,thisnboot,printstr='',randlist=[]):
    tempdata = []
    for iconf,ifile in enumerate(readfilelist):
        # print 'Reading '+printstr+' {}%            \r'.format(int((iconf*100)/float(len(readfilelist)))),
        try:
            if len(thisGammaList) > 0:
                tempdata.append(ReadFSCfunPick(ifile,thisMomList,thisGammaList).data)
            if len(thisDerGammaList) > 0:
                tempdata.append(ReadFSDerCfunPick(ifile,thisMomList,thisDerGammaList).data)
        except NaNCfunError as e:
            print 'Deleting file ' + ifile
            print e
            readfilelist.remove(ifile)
    if len(thisGammaList) > 0:
        return BootSet3pt(tempdata,thisMomList,thisGammaList,thisnboot,printstr='',randlist=randlist)
    elif len(thisDerGammaList) > 0:
        return BootSet3pt(tempdata,thisMomList,thisDerGammaList,thisnboot,printstr='',randlist=randlist)

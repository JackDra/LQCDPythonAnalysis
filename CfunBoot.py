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
            # for it in range(16):
            #     print 'pie3pt'
            #     for idata in np.array(data)[:,ig,ip,:]:
            #         print it,idata[it]
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
        # for icfg,cfgdata in enumerate(np.array(data)[:,ip,:]):
        #     print ''
        #     print 'icfg=',icfg
        #     for it,tdata in enumerate(cfgdata):
        #         print tdata
        # print 'pie2pt'
        # for iboot,idata in enumerate(np.array(data)[:,ip,12]):
        #     print iboot,idata
        bootdata,randlist = bt.CreateBoot(data[:,ip,:],nboot,0,randlist=randlist)
        dataout.append(bootdata)
    # print '                              \r',
    return dataout,randlist

#data = [ iconf , iflow , ip , it]
#dataout = [ iflow, ip , it ]. bs
def BootSet2ptTC(data,thisMomList,nboot,tflowlist,randlist=[]):
    dataout = []
    randlist = []
    for icf,iflow in enumerate(tflowlist):
        dataout.append([])
        for ip,imom in enumerate(thisMomList):
        # print 'Booting {}%  \r'.format(int((ip*100)/float(len(thisMomList)))),
        # for icfg,cfgdata in enumerate(np.array(data)[:,ip,:]):
        #     print ''
        #     print 'icfg=',icfg
        #     for it,tdata in enumerate(cfgdata):
        #         print tdata
        # print 'pie2pt'
        # for iboot,idata in enumerate(np.array(data)[:,ip,12]):
        #     print iboot,idata
            bootdata,randlist = bt.CreateBoot(data[:,icf,ip,:],nboot,0,randlist=randlist)
            dataout[-1].append(bootdata)
    # print '                              \r',
    return dataout

#dataout = [ ip , it ]. bs
def ReadAndBoot2pt(readfilelist,thisMomList,thisnboot,randlist=[]):
    tempdata = []
    shiftlist = []
    for iconf,ifile in enumerate(readfilelist):
        # print 'Reading {}%  \r'.format(int((iconf*100)/float(len(readfilelist)))),
        try:
            if CHROMA:
                if xsrcList[0] in ifile or not XAvg:
                    data = Read2ptCfunChromaXML(ifile,thisMomList)
                    tempdata.append(data.data)
                    shiftlist.append(data.tshiftlist)
            else:
                tempdata.append(Read2ptCfunPick(ifile,thisMomList).data)
        except NaNCfunError as e:
            print 
            print 'Deleting file ' + ifile
            print e
            print 'MUST RE-RUN AFTER THIS TO EXCLUDE BAD CONFIGS'
            print
            # os.remove(ifile)
    return BootSet2pt(np.array(tempdata),thisMomList,thisnboot,randlist=randlist),shiftlist


#dataout [ igamma , ip , it ]. bs
def ReadAndBoot3pt(readfilelist,thisMomList,thisGammaList,thisDerGammaList,thisnboot,shiftlist,printstr='',randlist=[]):
    tempdata = []
    counter = -1
    for iconf,ifile in enumerate(readfilelist):
        # print 'Reading '+printstr+' {}%            \r'.format(int((iconf*100)/float(len(readfilelist)))),
        try:
            if CHROMA:
                if xsrcList[0] in ifile or not XAvg:
                    counter += 1
                    if len(thisGammaList) > 0:
                        tempdata.append(ReadFSCfunPickCHROMA(ifile,thisMomList,thisGammaList,srcshift=shiftlist[counter]).data)
                    if len(thisDerGammaList) > 0:
                        raise IOError('Chroma version does not do derivatives, make DerList in Params.py be empty')
            else:
                if len(thisGammaList) > 0:
                    tempdata.append(ReadFSCfunPick(ifile,thisMomList,thisGammaList).data)
                if len(thisDerGammaList) > 0:
                    tempdata.append(ReadFSDerCfunPick(ifile,thisMomList,thisDerGammaList).data)
        except NaNCfunError as e:
            print 
            print 'Deleting file ' + ifile
            print e
            print 'MUST RE-RUN AFTER THIS TO EXCLUDE BAD CONFIGS'
            print
            # os.remove(ifile)
    if len(thisGammaList) > 0:
        return BootSet3pt(tempdata,thisMomList,thisGammaList,thisnboot,printstr='',randlist=randlist)
    elif len(thisDerGammaList) > 0:
        return BootSet3pt(tempdata,thisMomList,thisDerGammaList,thisnboot,printstr='',randlist=randlist)

    
def ReadAndBoot2ptTop(readfilelist,thisMomList,thisnboot,chargedata,chargecfglist,flowlist,randlist=[]):
    tempdata = []
    tempdataTop = []
    shiftlist = []
    for ifilepref,ifileList in readfilelist.iteritems():
        # print 'Reading {}%  \r'.format(int((iconf*100)/float(len(readfilelist)))),
        try:
            if CHROMA:
                chargeindex = FileToChargeCfg(ifilepref,chargecfglist)
                data = R2CChromaXMLFileList(ifileList,thisMomList,Dog5=True)
                tempdataTop.append([])
                tempdata.append(data.data)
                for iflowdata in chargedata[chargeindex]:                      
                    tempdataTop[-1].append(np.array(data.datag5)*iflowdata)
                shiftlist.append(data.tshiftlist)
            else:
                raise IOError('Top Charge not implemented for non chroma results')
        except NaNCfunError as e:
            print 
            print 'Deleting file ' + ifile
            print e
            print 'MUST RE-RUN AFTER THIS TO EXCLUDE BAD CONFIGS'
            print
            # os.remove(ifile)
    return (BootSet2ptTC(np.array(tempdataTop),thisMomList,thisnboot,flowlist,randlist=randlist),
            BootSet2pt(np.array(tempdata),thisMomList,thisnboot,randlist=randlist),shiftlist)


def FileToChargeCfg(ifile,chargecfglist):
    for ic,icharge in enumerate(chargecfglist):
        if icharge in ifile:
           return ic
    print 'charge cfg not found'
    return None

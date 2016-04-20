#!/usr/bin/env python
import os
from os import walk
from Params import *
from ReadBinaryCfuns import ReadFSCfunPick, ReadFSDerCfunPick, Read2ptCfunPick, NaNCfunError
from CfunBoot import *
from MiscFuns import *
import time,datetime


def TestMomList(thisMomList):
    try:
        PrintZMomI(np.where(np.array(thisMomList)==iqTOip(0))[0][0])
    except IndexError:
        print "Must Read Zero Momentum value " + iqTOip(0)


def ReadList(thisSmearList,thisMomList,thisProjGammaList,thisProjDerList,thisDSList,
             thisTSinkList,thisconflist,Flag,Interps=['nucleon']):
    thisfilelist = []
    f = open('./cfglistlist.txt','w')
    for iconf in thisconflist:
        file = filelist.replace('*',iconf)
        f.write(file+'\n')
        thisfilelist.append(file)
    f.close()
    data2pt,randlist = Read2ptSet(thisfilelist,thisSmearList,thisMomList,Interps)
    # print ''
    data3pt = []
    for iFlag,itsink in zip(Flag,thisTSinkList):
        data3pt.append(Read3ptSet(thisfilelist,thisSmearList,thisMomList,thisProjGammaList,
                                  thisProjDerList,thisDSList,itsink,iFlag,randlist=randlist))
    # print ''
    return [data2pt,data3pt,thisfilelist]

def ReadSet(thisSmearList,thisMomList,thisProjGammaList,thisProjDerList, thisDSList,
            thisTSinkList,directory,Flag,Interps=['nucleon']):
    # if len(thisTSinkList) > 0:
    #     TestMomList(thisMomList)
    thisfilelist = []
    f = open('./cfglistset.txt','w')
    for isource in SourceList:
        for (dirname,dirs,files) in walk(directory+'/'+isource):
            for file in files:
                fileend2pt = CreateEnd2pt(thisSmearList[0],thisSmearList[0],Interps[0],Interps[0])
                if fileend2pt not in file or ".metadata" in file: continue
                fileprefix = file.replace(fileend2pt,'')
                if CheckSet(fileprefix,dirname+'/',thisSmearList,thisProjGammaList,
                            thisProjDerList,thisDSList,thisTSinkList,Flag,Interps):
                    f.write(directory+'/'+isource+'/@/'+fileprefix+'\n')
                    thisfilelist.append(directory+'/'+isource+'/@/'+fileprefix)
    f.close()
    if ShowConfNum:
        print 'number of configs = ' , len(thisfilelist)
        print ''
        for ifile in thisfilelist:
            print ifile
        print ''
        
    data2pt,randlist = Read2ptSet(thisfilelist,thisSmearList,thisMomList,Interps)
    # print ''
    data3pt = []
    for iFlag,itsink in zip(Flag,thisTSinkList):
        data3pt.append(Read3ptSet(thisfilelist,thisSmearList,thisMomList,thisProjGammaList
                                  ,thisProjDerList,thisDSList,itsink,iFlag,randlist=randlist))
        # print ''
    return [data2pt,data3pt,thisfilelist]

def CreateEnd3pt(ism,jsm,Tsink,DS,Proj,Der):
    return '.t'+str(tsource)+'sm'+ism+jsm+Proj+'t'+str(Tsink)+'p000.'+DS+Der+'k'+str(kappa)+'.ifms.3cf'

def CreateEnd2pt(ism,jsm,iterp,jterp):
    if iterp == jterp: return '.t'+str(tsource)+'sm'+ism+'k'+str(kappa)+'.ifmssi'+jsm+'.'+iterp+'.u.2cf'
    else: return '.t'+str(tsource)+'sm'+ism+'k'+str(kappa)+'.ifmssi'+jsm+'.'+iterp+'_'+jterp+'.u.2cf'

def CreateDir3pt(ism,jsm,Tsink,DS,Proj,Flag):
    return Flag+'sm'+ism+jsm+Proj+'t'+str(Tsink)+'p000.'+DS

def CreateDir2pt(ism,jsm):
    return 'twoptsm'+ism+'si'+jsm

def CheckSet(FilePrefix,directory,thisSmearList,thisProjGammaList,thisProjDerList,thisDSList,
             thisTSinkList,Flag,Interps):
    CompSet = True
    for iterp,ism in Elongate(Interps,thisSmearList):
        for jcsm,(jterp,jsm) in enumerate(Elongate(Interps,thisSmearList)):
            testfile2pt = (directory.replace(CreateDir2pt(thisSmearList[0],thisSmearList[0]),
                                             CreateDir2pt(ism,jsm))
                           +FilePrefix+CreateEnd2pt(ism,jsm,iterp,jterp))
            if not os.path.isfile(testfile2pt): CompSet = False
        for iFlag,itsink in zip(Flag,thisTSinkList):
            thisFlag = iFlag
            C2C3Dis = ''
            if iFlag == 'REvec':
                Jsmlist = ['REvec']
            elif 'PoF' in iFlag:
                Jsmlist = ['RE'+PoFReadTvarList[0]]
                thisFlag = 'RE'+PoFDirTvarList[0]
                C2C3Dis = PoFC2C3Dis
            else:
                Jsmlist = ['Xsm'+jsm for jsm in thisSmearList]
            for jcsm,(jsm3pt,jsm) in enumerate(zip(Jsmlist,thisSmearList)):
                for iDS in thisDSList:
                    for iProj in thisProjGammaList:
                        testfile3pt = (directory.replace(CreateDir2pt(thisSmearList[0],thisSmearList[0]),
                                                         CreateDir3pt(ism,jsm3pt,itsink,iDS,iProj,thisFlag))
                                       +FilePrefix+CreateEnd3pt(ism,jsm3pt,itsink,iDS,iProj,''))
                        if not os.path.isfile(testfile3pt.replace(FileStruct,FileStruct+C2C3Dis)): CompSet = False
                    for iProj in thisProjDerList:
                        testfile3pt = (directory.replace(CreateDir2pt(thisSmearList[0],thisSmearList[0]),
                                                         CreateDir3pt(ism,jsm3pt,itsink,iDS,iProj,thisFlag))
                                       +FilePrefix+CreateEnd3pt(ism,jsm3pt,itsink,iDS,iProj,'D'))
                        if not os.path.isfile(testfile3pt.replace(FileStruct,FileStruct+C2C3Dis)): CompSet = False
    return CompSet

##tempdata [ iconf ] .data [ ip , it ]
##thisdata2pt [ ism , jsm , ip , it ] bootdataclas
def Read2ptSet(readfilelist,thisSmearList,thisMomList,Interps):
    thisdata2pt = []
    start = time.time()
    icount = -1
    randlist = []
    for icsm,(iterp,ism) in enumerate(Elongate(Interps,thisSmearList)):
        thisdata2pt.append([])
        for jcsm,(jterp,jsm) in enumerate(Elongate(Interps,thisSmearList)):
            icount += 1
            timeleft = GetTimeLeft(icount,len(Elongate(Interps,thisSmearList))**2,time.time()-start)
            print 'Read 2pt: sm' + ism + 'Xsm'+jsm+' Time Left: ' +str(datetime.timedelta(seconds=timeleft)) , ' h:m:s \r',
            thisstart = time.time()
            thisfilelist = [ifile.replace('@',CreateDir2pt(ism,jsm))+CreateEnd2pt(ism,jsm,iterp,jterp) for ifile in readfilelist]
            dataout,randlist = ReadAndBoot2pt(thisfilelist,thisMomList,nboot) 
            thisdata2pt[icsm].append(dataout)
            print 'Read 2pt: sm' + ism + 'Xsm'+jsm + ' took: ' +str(datetime.timedelta(seconds=time.time()-thisstart)) , ' h:m:s        '
    print 'Read 2pt total time: ' + str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s'
    return thisdata2pt,randlist


##thisdata3pt [ ism , jsm , igamma , ip ] bootdataclas
def Read3ptSet(readfilelist,thisSmearList,thisMomList,thisProjGammaList,thisProjDerList,thisDSList,tsink,iFlag,randlist=[]):
    thisdata3pt = []
    start = time.time()
    icount = -1
    for icsm,ism in enumerate(thisSmearList):
        thisdata3pt.append([])
        thisFlag = iFlag
        C2C3Dis,rownumb = '',''
        if iFlag == 'REvec':
            Jsmlist = ['REvec']
        elif iFlag == 'PoF':
            Jsmlist = ['RE'+PoFReadTvarList[0]]
            thisFlag = 'RE'+PoFDirTvarList[0]
            C2C3Dis = PoFC2C3Dis
        else:
            Jsmlist = ['Xsm'+jsm for jsm in thisSmearList]
        for jcsm,jsm in enumerate(Jsmlist):
            thisdata3pt[icsm].append([])
            thisstart = time.time()
            icount += 1
            timeleft = GetTimeLeft(icount,len(thisSmearList)*len(Jsmlist),time.time()-start)
            print 'Read 3pt: sm' + ism + jsm + ' ts'+str(tsink)+' Time Remaining: '+ str(datetime.timedelta(seconds=timeleft)) , ' h:m:s \r',
            for iDS in thisDSList:
                for iProj,thisGammaList in thisProjGammaList.iteritems():
                    thisfilelist = [ifile.replace('@',CreateDir3pt(ism,jsm,tsink,iDS,iProj,thisFlag))
                                    .replace(FileStruct,FileStruct+C2C3Dis)
                                    +CreateEnd3pt(ism,jsm,tsink,iDS,iProj,'') for ifile in readfilelist]
                    holdres = ReadAndBoot3pt(thisfilelist,thisMomList,thisGammaList,[],nboot,printstr=iDS+iProj,randlist=randlist)
                    for gammares in holdres:
                        thisdata3pt[icsm][jcsm].append(gammares)
                for iProj,thisDerList in thisProjDerList.iteritems():
                    thisfilelist = [ifile.replace('@',CreateDir3pt(ism,jsm,tsink,iDS,iProj,thisFlag))
                                    .replace(FileStruct,FileStruct+C2C3Dis)
                                    +CreateEnd3pt(ism,jsm,tsink,iDS,iProj,'D') for ifile in readfilelist]
                    holdres = ReadAndBoot3pt(thisfilelist,thisMomList,[],thisDerList,nboot,printstr=iDS+iProj,randlist=randlist)
                    for gammares in holdres:
                        thisdata3pt[icsm][jcsm].append(gammares)
            print 'Read 3pt: sm' + ism + jsm + ' ts'+str(tsink)+ ' took: ' + str(datetime.timedelta(seconds=time.time()-thisstart)) , ' h:m:s        '
    print 'Read 3pt total time: ' + str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s'
    return thisdata3pt


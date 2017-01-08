#!/usr/bin/env python
import os
from os import walk
from Params import *
from ReadBinaryCfuns import ReadFSCfunPick, ReadFSDerCfunPick, Read2ptCfunPick, NaNCfunError
from CfunBoot import *
from ReadTxt import *
from MiscFuns import *
import time,datetime


def TestMomList(thisMomList):
    try:
        PrintZMomI(np.where(np.array(thisMomList)==iqTOip(0))[0][0])
    except IndexError:
        print "Must Read Zero Momentum value " + iqTOip(0)


def ReadListTopCharge(thisSmearList,thisMomList,thisconflist,Interps=['nucleon'],thistsourceList=[tsource]):
    thisfilelist = OrderedDict()
    f = open('./cfglistlist.txt','w')
    for iconf in thisconflist:
        ifile = filelist.replace('*',iconf)
        prefnosrc = re.sub('_xsrc.*','',ifile)
        if prefnosrc not in thisfilelist.keys():
            thisfilelist[prefnosrc] = [ifile]
        else:
            thisfilelist[prefnosrc].append(ifile)
        f.write(ifile+'\n')
    f.close()

    
    print thisfilelist.keys()
    cfglistout,topcharge,tflow = ReadTopList(TCDir,StripSrc(thisfilelist.keys()))
    if not np.all([x==tflow[0] for x in tflow]):
        print 'warning, files had different flow times'
    thisTflowList = tflow[0]
    # if len(cfglistout) != len(thiscfglist):
    #     print 'Not All Top Charge Found'
        
    dataTop,data2pt,randlist,shiftlist = Read2ptSetTop(topcharge,thisfilelist,thisSmearList,GetAvgMomListip(thisMomList),
                                                       Interps,cfglistout,thisTflowList,tsourceList=thistsourceList)
    return [data2pt,dataTop,thisTflowList,thisfilelist]



def ReadSetTopCharge(thisSmearList,thisMomList,directory,Interps=['nucleon'],thistsourceList=[tsource]):
    # if len(thisTSinkList) > 0:
    #     TestMomList(thisMomList)
    thisfilelist = OrderedDict()
    f = open('./cfglistset.txt','w')
    
    fileend2pt = CreateEnd2pt(DefSmearList[0],DefSmearList[0],thistsourceList[0],Interps[0],Interps[0])
    if Debug: print ' comparing to: ',fileend2pt
    for isource in SourceList:
        print 'reading directory: ' , directory+'/'+isource
        for (dirname,dirs,files) in walk(directory+'/'+isource):
            for ifile in files:
                if fileend2pt not in ifile or ".metadata" in ifile: continue
                # ## FOR DEBUGGING
                # if xsrcList[0] not in ifile: continue
                if 'xsrc' in ListOrSet:
                    if ListOrSet.replace('ReadSet','').replace('ReadList','')+'_' not in ifile: continue
                fileprefix = ifile.replace(fileend2pt,'')
                if CheckSet(fileprefix,dirname+'/',thisSmearList,{},{},[],[],'',Interps,tsourceList=thistsourceList):
                    prefnosrc = re.sub('_xsrc.*','',fileprefix)
                    if prefnosrc not in thisfilelist.keys():
                        thisfilelist[prefnosrc] = [directory+'/'+isource+'/@/'+fileprefix]
                    else:
                        if XSrcLen <= len(thisfilelist[prefnosrc]) and ForceXSrcLen: continue
                        thisfilelist[prefnosrc].append(directory+'/'+isource+'/@/'+fileprefix)
                    f.write(directory+'/'+isource+'/@/'+fileprefix+'\n')
    f.close()
    if ExactXSrcNumber:
        maxlen = np.max([len(ifilelist) for ifilelist in thisfilelist.itervalues()])
        for ikey in thisfilelist.iterkeys():
            if len(thisfilelist[ikey]) != maxlen:
                del thisfilelist[ikey]
    print 'number of configs = ' , len(thisfilelist.keys())
    print 'average number of sources per cfg = ' ,np.mean([len(ifilelist) for ifilelist in thisfilelist.itervalues()])
    print 'total number of measurements = ' , np.sum([len(ifilelist) for ifilelist in thisfilelist.itervalues()])
    print ''
    if len(thisfilelist) == 0: raise IOError('No Configurations Found')
    if ShowConfNum:
        for ifile in thisfilelist:
            print ifile
        print ''
        
    cfglistout,topcharge,tflow = ReadTopList(TCDir,StripSrc(thisfilelist.keys()))
    if not np.all([x==tflow[0] for x in tflow]):
        print 'warning, files had different flow times'
    thisTflowList = tflow[0]
    # if len(cfglistout) != len(thiscfglist):
    #     print 'Not All Top Charge Found'

    dataTop,data2pt,randlist,shiftlist = Read2ptSetTop(topcharge,thisfilelist,thisSmearList,
                                                       GetAvgMomListip(thisMomList),Interps,cfglistout,thisTflowList,tsourceList=thistsourceList)
    
    return [data2pt,dataTop,thisTflowList,thisfilelist]

        
def ReadList(thisSmearList,thisMomList,thisProjGammaList,thisProjDerList,thisDSList,
             thisTSinkList,thisconflist,Flag,Interps=['nucleon'],thistsourceList=[tsource]):
    thisfilelist = []
    f = open('./cfglistlist.txt','w')
    for iconf in thisconflist:
        file = filelist.replace('*',iconf)
        f.write(file+'\n')
        thisfilelist.append(file)
    f.close()
    data2pt,randlist,shiftlist = Read2ptSet(thisfilelist,thisSmearList,GetAvgMomListip(thisMomList),Interps,tsourceList=thistsourceList)
    # print ''
    data3pt = []
    for iFlag,itsink in zip(Flag,thisTSinkList):
        data3pt.append(Read3ptSet(thisfilelist,thisSmearList,thisMomList,thisProjGammaList,
                                  thisProjDerList,thisDSList,itsink,iFlag,shiftlist,randlist=randlist,tsourceList=thistsourceList))
    # print ''
    return [data2pt,data3pt,thisfilelist]

def ReadSet(thisSmearList,thisMomList,thisProjGammaList,thisProjDerList, thisDSList,
            thisTSinkList,directory,Flag,Interps=['nucleon'],thistsourceList=[tsource]):
    # if len(thisTSinkList) > 0:
    #     TestMomList(thisMomList)
    thisfilelist = []
    f = open('./cfglistset.txt','w')
    
    fileend2pt = CreateEnd2pt(DefSmearList[0],DefSmearList[0],thistsourceList[0],Interps[0],Interps[0])
    if Debug: print ' comparing to: ',fileend2pt
    for isource in SourceList:
        print 'reading directory: ' , directory+'/'+isource
        for (dirname,dirs,files) in walk(directory+'/'+isource):
            for ifile in files:
                if CfunConfigCheck:
                    if xsrcList[0] not in ifile and XAvg : continue
                    # ## FOR DEBUGGING
                    # if xsrcList[0] not in ifile: continue
                    if fileend2pt not in ifile or ".metadata" in ifile: continue
                    if 'xsrc' in ListOrSet:
                        if ListOrSet.replace('ReadSet','').replace('ReadList','')+'_' not in ifile: continue
                    fileprefix = ifile.replace(fileend2pt,'')
                    if CheckAllSet(fileprefix,dirname+'/',Interps,tsourceList=thistsourceList):
                        f.write(directory+'/'+isource+'/@/'+fileprefix+'\n')
                        thisfilelist.append(directory+'/'+isource+'/@/'+fileprefix)                    
                else:
                    if fileend2pt not in ifile or ".metadata" in ifile: continue
                    # ## FOR DEBUGGING
                    # if xsrcList[0] not in ifile: continue
                    if xsrcList[0] not in ifile and XAvg: continue
                    if 'xsrc' in ListOrSet:
                        if ListOrSet.replace('ReadSet','').replace('ReadList','')+'_' not in ifile: continue
                    fileprefix = ifile.replace(fileend2pt,'')
                    if CheckSet(fileprefix,dirname+'/',thisSmearList,thisProjGammaList,
                                thisProjDerList,thisDSList,thisTSinkList,Flag,Interps,tsourceList=thistsourceList):
                        f.write(directory+'/'+isource+'/@/'+fileprefix+'\n')
                        thisfilelist.append(directory+'/'+isource+'/@/'+fileprefix)
    f.close()
    print 'number of configs = ' , len(thisfilelist)
    print ''
    if len(thisfilelist) == 0: raise IOError('No Configurations Found')
    if ShowConfNum:
        for ifile in thisfilelist:
            print ifile
        print ''
        
    data2pt,randlist,shiftlist = Read2ptSet(thisfilelist,thisSmearList,GetAvgMomListip(thisMomList),Interps,tsourceList=thistsourceList)
    # print ''
    data3pt = []
    for iFlag,itsink in zip(Flag,thisTSinkList):
        data3pt.append(Read3ptSet(thisfilelist,thisSmearList,thisMomList,thisProjGammaList
                                  ,thisProjDerList,thisDSList,itsink,iFlag,shiftlist,randlist=randlist,tsourceList=thistsourceList))
        # print ''
    return [data2pt,data3pt,thisfilelist]

def CreateEnd3pt(ism,jsm,thists,Tsink,DS,Proj,Der):
    if CHROMA:
        return '_k'+str(kappa)+'_tsrc'+str(thists)+'sm'+ism+'_'+jsm+Proj+'tsink'+str(Tsink)+'p000'+DS+'NDer0.3cf'
    else:
        return '.t'+str(thists)+'sm'+ism+jsm+Proj+'t'+str(Tsink)+'p000.'+DS+Der+'k'+str(kappa)+'.ifms.3cf'

def CreateEnd2pt(ism,jsm,thists,iterp,jterp):
    # if iterp == jterp: return '.t'+str(thists)+'sm'+ism+'k'+str(kappa)+'.ifmssi'+jsm+'.'+iterp+'.u.2cf'
    # else: return '.t'+str(thists)+'sm'+ism+'k'+str(kappa)+'.ifmssi'+jsm+'.'+iterp+'_'+jterp+'.u.2cf'
    return '_k'+str(kappa)+'_tsrc'+str(thists)+'sm'+ism+'si'+jsm+'_'+iterp+'.2cf.xml'

def CreateDir3pt(ism,jsm,Tsink,DS,Proj,Flag):    
    if CHROMA:
        return 'threept/sm'+ism+jsm+Proj+'tsink'+str(Tsink)+'p000'+DS
    else:
        return Flag+'sm'+ism+jsm+Proj+'t'+str(Tsink)+'p000.'+DS

def CreateDir2pt(ism,jsm):
    return 'twoptRandT/twoptsm'+ism+'si'+jsm

def CheckAllSet(FilePrefix,directory,Interps,tsourceList=[tsource]):
    if XAvg: thisxsrcList = xsrcList
    else: thisxsrcList = [xsrcList[0]]
    for xsrc in thisxsrcList:
        for iterp,ism in Elongate(Interps,DefSmearList):
            for jcsm,(jterp,jsm) in enumerate(Elongate(Interps,DefSmearList)):
                for thists in tsourceList:
                    testfile2pt = (directory.replace(CreateDir2pt(DefSmearList[0],DefSmearList[0]),
                                                     CreateDir2pt(ism,jsm))
                                   +FilePrefix+CreateEnd2pt(ism,jsm,thists,iterp,jterp))
                    testfile2pt = testfile2pt.replace(xsrcList[0],xsrc)
                    if Debug: print 'Checking!! ' ,testfile2pt
                    if not os.path.isfile(testfile2pt): return False
        for iFlag in CfunCheckList:
            if iFlag not in TSinkDictList: continue
            for thists in tsourceList:
                for iterp,ism in Elongate(Interps,SmearDictList[iFlag]):
                    for itsink in TSinkDictList[iFlag]:
                        thisFlag = iFlag.replace('Read','')
                        C2C3Dis = ''
                        if 'REvec' in iFlag:
                            Jsmlist = ['REvec']
                        elif 'PoF' in iFlag:
                            if CHROMA:
                                Jsmlist = ['PoF'+str(PoFShifts)+'D'+str(PoFDelta)]
                                thisFlag = ''
                                C2C3Dis = ''
                            else:
                                Jsmlist = ['RE'+PoFReadTvarList[0]]
                                thisFlag = 'RE'+PoFDirTvarList[0]
                                C2C3Dis = PoFC2C3Dis
                        else:
                            Jsmlist = ['Xsm'+jsm for jsm in SmearDictList[iFlag]]
                        for jcsm,jsm3pt in enumerate(Jsmlist):
                            for iDS in DefDSList:
                                for iProj in DefProjGammaList:
                                    testfile3pt = (directory.replace(CreateDir2pt(DefSmearList[0],DefSmearList[0]),
                                                                     CreateDir3pt(ism,jsm3pt,itsink,iDS,iProj,thisFlag))
                                                   +FilePrefix+CreateEnd3pt(ism,jsm3pt,thists,itsink,iDS,iProj,''))
                                    if Debug: print 'Checking!!! ' ,  testfile3pt.replace(FileStruct,FileStruct+C2C3Dis)
                                    testfile3pt = testfile3pt.replace(xsrcList[0],xsrc)
                                    if not os.path.isfile(testfile3pt.replace(FileStruct,FileStruct+C2C3Dis)): return False
                                for iProj in DefProjDerList:
                                    testfile3pt = (directory.replace(CreateDir2pt(DefSmearList[0],DefSmearList[0]),
                                                                     CreateDir3pt(ism,jsm3pt,itsink,iDS,iProj,thisFlag))
                                                   +FilePrefix+CreateEnd3pt(ism,jsm3pt,thists,itsink,iDS,iProj,'D'))
                                    if Debug: print 'Checking!!! ' ,testfile3pt.replace(FileStruct,FileStruct+C2C3Dis)
                                    testfile3pt = testfile3pt.replace(xsrcList[0],xsrc)
                                    if not os.path.isfile(testfile3pt.replace(FileStruct,FileStruct+C2C3Dis)): return False
    return True
                


def CheckSet(FilePrefix,directory,thisSmearList,thisProjGammaList,thisProjDerList,thisDSList,
             thisTSinkList,Flag,Interps,tsourceList=[tsource]):
    for iterp,ism in Elongate(Interps,thisSmearList):
        for thists in tsourceList:
            for jcsm,(jterp,jsm) in enumerate(Elongate(Interps,thisSmearList)):
                testfile2pt = (directory.replace(CreateDir2pt(thisSmearList[0],thisSmearList[0]),
                                                 CreateDir2pt(ism,jsm))
                               +FilePrefix+CreateEnd2pt(ism,jsm,thists,iterp,jterp))
                # testfile2pt = testfile2pt.replace(xsrcList[0],xsrc)
                if Debug: print 'Checking~!! ' ,testfile2pt
                if not os.path.isfile(testfile2pt): return False
            for iFlag,itsink in zip(Flag,thisTSinkList):
                thisFlag = iFlag
                C2C3Dis = ''
                if iFlag == 'REvec':
                    Jsmlist = ['REvec']
                elif 'PoF' in iFlag:
                    if CHROMA:
                        Jsmlist = ['PoF'+str(PoFShifts)+'D'+str(PoFDelta)]
                        thisFlag = ''
                        C2C3Dis = ''
                    else:
                        Jsmlist = ['RE'+PoFReadTvarList[0]]
                        thisFlag = 'RE'+PoFDirTvarList[0]
                        C2C3Dis = PoFC2C3Dis
                else:
                    Jsmlist = ['Xsm'+jsm for jsm in thisSmearList]
                for jcsm,jsm3pt in enumerate(Jsmlist):
                    for iDS in thisDSList:
                        for iProj in thisProjGammaList:
                            testfile3pt = (directory.replace(CreateDir2pt(thisSmearList[0],thisSmearList[0]),
                                                             CreateDir3pt(ism,jsm3pt,itsink,iDS,iProj,thisFlag))
                                           +FilePrefix+CreateEnd3pt(ism,jsm3pt,thists,itsink,iDS,iProj,''))
                            # testfile3pt = testfile3pt.replace(xsrcList[0],xsrc)
                            if Debug: print 'Checking~!!! ' ,  testfile3pt.replace(FileStruct,FileStruct+C2C3Dis)
                            if not os.path.isfile(testfile3pt.replace(FileStruct,FileStruct+C2C3Dis)): return False
                        for iProj in thisProjDerList:
                            testfile3pt = (directory.replace(CreateDir2pt(thisSmearList[0],thisSmearList[0]),
                                                             CreateDir3pt(ism,jsm3pt,itsink,iDS,iProj,thisFlag))
                                           +FilePrefix+CreateEnd3pt(ism,jsm3pt,thists,itsink,iDS,iProj,'D'))
                            # testfile3pt = testfile3pt.replace(xsrcList[0],xsrc)
                            if Debug: print 'Checking~!!! ' ,  testfile3pt.replace(FileStruct,FileStruct+C2C3Dis)
                            if not os.path.isfile(testfile3pt.replace(FileStruct,FileStruct+C2C3Dis)): return False
    return True

##tempdata [ iconf ] .data [ ip , it ]
##thisdata2pt [ tsource, ism , jsm , ip , it ] bootdataclas
##shiftlist = [ tsource , ism, jsm , icfg, isrc_number ]
def Read2ptSet(readfilelist,thisSmearList,thisMomList,Interps,tsourceList=[tsource]):
    thisdata2pt = []
    start = time.time()
    icount = -1
    randlist = []
    shiftlist = []
    for thists in tsourceList:
        thisdata2pt.append([])
        shiftlist.append([])
        for icsm,(iterp,ism) in enumerate(Elongate(Interps,thisSmearList)):
            thisdata2pt[-1].append([])
            shiftlist[-1].append([])
            for jcsm,(jterp,jsm) in enumerate(Elongate(Interps,thisSmearList)):
                icount += 1
                timeleft = GetTimeLeft(icount,len(tsourceList)*len(Elongate(Interps,thisSmearList))**2,time.time()-start)
                print 'Read 2pt: sm' + ism + 'Xsm'+jsm+' Time Left: ' +str(datetime.timedelta(seconds=timeleft)) , ' h:m:s \r',
                thisstart = time.time()
                thisfilelist = [ifile.replace('@',CreateDir2pt(ism,jsm))+CreateEnd2pt(ism,jsm,thists,iterp,jterp) for ifile in readfilelist]
                (dataout,randlist),thisshiftlist = ReadAndBoot2pt(thisfilelist,thisMomList,nboot) 
                shiftlist[-1][icsm].append(thisshiftlist)
                thisdata2pt[-1][icsm].append(dataout)
                print 'Read 2pt: sm' + ism + 'Xsm'+jsm + ' t_src'+str(thists)+' took: ' +str(datetime.timedelta(seconds=time.time()-thisstart)) , ' h:m:s        '
    print 'Read 2pt total time: ' + str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s'
    return thisdata2pt,randlist,shiftlist


##thisdata3pt [ thistsource, ism , jsm , igamma , ip ] bootdataclas
def Read3ptSet(readfilelist,thisSmearList,thisMomList,thisProjGammaList,thisProjDerList,thisDSList,tsink,iFlag,shiftlist,randlist=[],tsourceList=[tsource]):
    thisdata3pt = []
    start = time.time()
    icount = -1
    for itsource,thists in enumerate(tsourceList):
        thisdata3pt.append([])
        for icsm,ism in enumerate(thisSmearList):
            thisdata3pt[-1].append([])
            thisFlag = iFlag
            C2C3Dis,rownumb = '',''
            if iFlag == 'REvec':
                Jsmlist = ['REvec']
            elif iFlag == 'PoF':
                if CHROMA:
                    Jsmlist = ['PoF'+str(PoFShifts)+'D'+str(PoFDelta)]
                    thisFlag = ''
                    C2C3Dis = ''
                else:
                    Jsmlist = ['RE'+PoFReadTvarList[0]]
                    thisFlag = 'RE'+PoFDirTvarList[0]
                    C2C3Dis = PoFC2C3Dis
            else:
                Jsmlist = ['Xsm'+jsm for jsm in thisSmearList]
            for jcsm,jsm in enumerate(Jsmlist):
                thisdata3pt[-1][icsm].append([])
                thisstart = time.time()
                icount += 1
                timeleft = GetTimeLeft(icount,len(tsourceList)*len(thisSmearList)*len(Jsmlist),time.time()-start)
                print 'Read 3pt: sm' + ism + jsm + ' ts'+str(tsink)+' Time Remaining: '+ str(datetime.timedelta(seconds=timeleft)) , ' h:m:s \r',
                thisshiftlist = shiftlist[itsource][icsm][jcsm]
                for iDS in thisDSList:
                    for iProj,thisGammaList in thisProjGammaList.iteritems():
                        thisfilelist = [ifile.replace('@',CreateDir3pt(ism,jsm,tsink,iDS,iProj,thisFlag))
                                        .replace(FileStruct,FileStruct+C2C3Dis)
                                        +CreateEnd3pt(ism,jsm,thists,tsink,iDS,iProj,'') for ifile in readfilelist]

                        holdres = ReadAndBoot3pt(thisfilelist,thisMomList,thisGammaList,[],nboot,thisshiftlist,printstr=iDS+iProj,randlist=randlist)
                        for gammares in holdres:
                            thisdata3pt[-1][icsm][jcsm].append(gammares)
                    for iProj,thisDerList in thisProjDerList.iteritems():
                        thisfilelist = [ifile.replace('@',CreateDir3pt(ism,jsm,tsink,iDS,iProj,thisFlag))
                                        .replace(FileStruct,FileStruct+C2C3Dis)
                                        +CreateEnd3pt(ism,jsm,thists,tsink,iDS,iProj,'D') for ifile in readfilelist]
                        holdres = ReadAndBoot3pt(thisfilelist,thisMomList,[],thisDerList,nboot,thisshiftlist,printstr=iDS+iProj,randlist=randlist)
                        for gammares in holdres:
                            thisdata3pt[-1][icsm][jcsm].append(gammares)
                print 'Read 3pt: sm' + ism + jsm + 't_src'+str(thists)+' t_sink'+str(tsink)+ ' took: ' + str(datetime.timedelta(seconds=time.time()-thisstart)) , ' h:m:s        '
    print 'Read 3pt total time: ' + str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s'
    return thisdata3pt





##tempdata [ iconf ] .data [ ip , it ]
##thisdata2pt [ tflow, tsource, ism , jsm , ip , it ] bootdataclas
##shiftlist = [ tsource , ism, jsm , icfg, isrc_number ]
##TopCharge = [ iconf , tflow ] 
def Read2ptSetTop(TopCharge,readfilelist,thisSmearList,thisMomList,Interps,Topcfglist,thisTflowList,tsourceList=[tsource]):
    thisdata2pt = []
    thisdataTop = []
    start = time.time()
    icount = -1
    randlist = []
    shiftlist = []
    for thists in tsourceList:
        thisdata2pt.append([])
        thisdataTop.append([])
        shiftlist.append([])
        for icsm,(iterp,ism) in enumerate(Elongate(Interps,thisSmearList)):
            thisdata2pt[-1].append([])
            thisdataTop[-1].append([])
            shiftlist[-1].append([])
            for jcsm,(jterp,jsm) in enumerate(Elongate(Interps,thisSmearList)):
                icount += 1
                timeleft = GetTimeLeft(icount,len(tsourceList)*len(Elongate(Interps,thisSmearList))**2,time.time()-start)
                print 'Read 2pt: sm' + ism + 'Xsm'+jsm+' Time Left: ' +str(datetime.timedelta(seconds=timeleft)) , ' h:m:s \r',
                thisstart = time.time()
                thisfilelist = OrderedDict()
                for ireadkey,ireadfl in readfilelist.iteritems():
                    thisfilelist[ireadkey] = [ifile.replace('@',CreateDir2pt(ism,jsm))+CreateEnd2pt(ism,jsm,thists,iterp,jterp) for ifile in ireadfl]
                dataoutTop,(dataout,randlist),thisshiftlist = ReadAndBoot2ptTop(thisfilelist,thisMomList,nboot,TopCharge,Topcfglist,thisTflowList) 
                shiftlist[-1][icsm].append(thisshiftlist)
                thisdata2pt[-1][icsm].append(dataout)
                thisdataTop[-1][icsm].append(dataoutTop)
                print 'Read 2pt: sm' + ism + 'Xsm'+jsm + ' t_src'+str(thists)+' took: ' +str(datetime.timedelta(seconds=time.time()-thisstart)) , ' h:m:s        '
    print 'Read 2pt total time: ' + str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s'
    return np.rollaxis(np.array(thisdataTop),3),thisdata2pt,randlist,shiftlist
##thisdata2pt [ tflow, tsource, ism , jsm , ip , it ] bootdataclas

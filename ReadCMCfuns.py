#!/usr/bin/env python
import os
from os import walk
from Params import *
from ReadBinaryCfuns import ReadFSCfunPick, ReadFSDerCfunPick, Read2ptCfunPick, NaNCfunError
from CfunBoot import *
from ReadTxt import *
from MiscFuns import *
import time,datetime

def RemoveCfgs(thisfilelist,newcfglist):
    filelistout = OrderedDict()
    for ioldcfg in thisfilelist.iterkeys():
        print ioldcfg
    for icfg in newcfglist:
        print icfg
        for ioldcfg,icfgvalues in thisfilelist.iteritems():
            if icfg in ioldcfg:
                filelistout[ioldcfg] = icfgvalues
    return filelistout


def ReadAndCheckTop(thisWein,thisfilelist):
    thiscfglist = StripSrc(thisfilelist.keys())
    if thisWein:
        cfglistout,topcharge,tflow = ReadTopList(WeinDir,thiscfglist)
        thisfilelist = RemoveCfgs(thisfilelist,cfglistout)
        if QoppConfigCheck:
            cfglistout,dump,dump2 = ReadTopList(TCDir,thiscfglist,OnlyCheck=True)
            thisfilelist = RemoveCfgs(thisfilelist,cfglistout)                
    else:
        cfglistout,topcharge,tflow = ReadTopList(TCDir,thiscfglist)
        thisfilelist = RemoveCfgs(thisfilelist,cfglistout)
        if WoppConfigCheck:
            cfglistout,dump,dump2 = ReadTopList(WeinDir,thiscfglist,OnlyCheck=True)
            thisfilelist = RemoveCfgs(thisfilelist,cfglistout)
    return thisfilelist,topcharge,tflow


def CheckTop(thisfilelist):
    thiscfglist = StripSrc(thisfilelist.keys())
    if WoppConfigCheck:
        cfglistout,topcharge,tflow = ReadTopList(WeinDir,thiscfglist,OnlyCheck=True)
        thisfilelist = RemoveCfgs(thisfilelist,cfglistout)
    if QoppConfigCheck:
        cfglistout,topcharge,tflow = ReadTopList(TCDir,thiscfglist,OnlyCheck=True)
        thisfilelist = RemoveCfgs(thisfilelist,cfglistout)                
    return thisfilelist


def DoForceChecks(thisfilelist):
    if FlipConfs: thisfilelist = OrderedDict(reversed(list(thisfilelist.items())))
    if ExactXSrcNumber:
        maxlen = np.max([len(ifilelist) for ifilelist in thisfilelist.itervalues()])
        for ikey in thisfilelist.iterkeys():
            if len(thisfilelist[ikey]) != maxlen:
                del thisfilelist[ikey]
    if ForceMinXSrcLen:
        for ikey in thisfilelist.iterkeys():
            if len(thisfilelist[ikey]) < MinXSrcLen:
                del thisfilelist[ikey]
    if ForceNConf:
        for inum,ikey in enumerate(thisfilelist.iterkeys()):
            if inum >= NConfLen:
                del thisfilelist[ikey]
    return thisfilelist

def ReadSetTopCharge(thisiSmearList,thisjSmearList,thisMomList,thisProjGammaList,thisProjDerList,thisDSList,
                     thisTSinkList,directory,Flag,Interps=['nucleon'],thistsourceList=[tsource],Wein=False):
    # if len(thisTSinkList) > 0:
    #     TestMomList(thisMomList)
        
    thisfilelist = OrderedDict()
    f = open('./cfglistset.txt','w')
    
    fileend2pt = CreateEnd2pt(DefiSmearList[0],DefjSmearList[0],thistsourceList[0],Interps[0],Interps[0])
    if Debug: print ' comparing to: ',fileend2pt

    for isource in SourceList:
        print 'reading directory: ' , directory+'/'+isource
        for (dirname,dirs,files) in walk(directory+'/'+isource):
            if 'Temp' in dirname: continue
            for ifile in filesortCfuns(files):
                if CfunConfigCheck:
                    # ## FOR DEBUGGING
                    # if xsrcList[0] not in ifile: continue
                    if fileend2pt not in ifile or ".metadata" in ifile: continue
                    if 'xsrc' in ListOrSet:
                        if ListOrSet.replace('ReadSet','').replace('ReadList','')+'_' not in ifile: continue
                    fileprefix = ifile.replace(fileend2pt,'')
                    if CheckAllSet(fileprefix,dirname+'/',Interps,tsourceList=thistsourceList):
                        prefnosrc = re.sub('_xsrc.*','',fileprefix)
                        if prefnosrc not in thisfilelist.keys():
                            thisfilelist[prefnosrc] = [directory+'/'+isource+'/@/'+fileprefix]
                        else:
                            if XSrcLen <= len(thisfilelist[prefnosrc]) and ForceXSrcLen: continue
                            thisfilelist[prefnosrc].append(directory+'/'+isource+'/@/'+fileprefix)
                        f.write(directory+'/'+isource+'/@/'+fileprefix+'\n')
                else:
                    if fileend2pt not in ifile or ".metadata" in ifile: continue
                    # ## FOR DEBUGGING
                    # if xsrcList[0] not in ifile: continue
                    if 'xsrc' in ListOrSet:
                        if ListOrSet.replace('ReadSet','').replace('ReadList','')+'_' not in ifile: continue
                    fileprefix = ifile.replace(fileend2pt,'')
                    if CheckSet(fileprefix,dirname+'/',thisiSmearList,thisjSmearList,thisProjGammaList,
                                thisProjDerList,thisDSList,thisTSinkList,Flag,Interps,tsourceList=thistsourceList):
                        prefnosrc = re.sub('_xsrc.*','',fileprefix)
                        if prefnosrc not in thisfilelist.keys():
                            thisfilelist[prefnosrc] = [directory+'/'+isource+'/@/'+fileprefix]
                        else:
                            if XSrcLen <= len(thisfilelist[prefnosrc]) and ForceXSrcLen: continue
                            thisfilelist[prefnosrc].append(directory+'/'+isource+'/@/'+fileprefix)
                        f.write(directory+'/'+isource+'/@/'+fileprefix+'\n')
    f.close()
    thisfilelist = DoForceChecks(thisfilelist)
    thisfilelist,topcharge,tflow = ReadAndCheckTop(Wein,thisfilelist)
    print 'number of configs = ' , len(thisfilelist.keys())
    print 'average number of sources per cfg = ' ,np.mean([len(ifilelist) for ifilelist in thisfilelist.itervalues()])
    print 'total number of measurements = ' , np.sum([len(ifilelist) for ifilelist in thisfilelist.itervalues()])
    print ''
    if len(thisfilelist.keys()) == 0: raise IOError('No Configurations Found')
    if ShowConfNum:
        for ifile in thisfilelist.itervalues():
            for iifile in ifile:
                print iifile
        print ''
    # print 'number of configs = ' , len(thisfilelist)
    # print ''
    # if len(thisfilelist) == 0: raise IOError('No Configurations Found')
    # if ShowConfNum:
    #     for ifile in thisfilelist:
    #         print ifile
    #     print ''
    data2pt,randlist,shiftlist = Read2ptSet(thisfilelist,thisiSmearList,thisjSmearList,GetAvgMomListip(thisMomList),Interps,tsourceList=thistsourceList)
    # print ''
    #thisfilelist = { icfg: [isrc] }
    #newcfglist = [icfg]
    
    if not np.all([x==tflow[0] for x in tflow]):
        print 'warning, files had different flow times'
    thisTflowList = tflow[0]

    data3pt,dataTop = [],[]
    for iFlag,itsink in zip(Flag,thisTSinkList):
        holdTop,hold3pt = Read3ptSetTop(topcharge,thisfilelist,thisiSmearList,thisjSmearList,thisMomList,thisProjGammaList,
                                        thisProjDerList,thisDSList,itsink,iFlag,shiftlist,cfglistout,thisTflowList,
                                        randlist=randlist,tsourceList=thistsourceList)
        data3pt.append(hold3pt)
        dataTop.append(holdTop)
    # print ''
    return [data2pt,data3pt,np.rollaxis(np.array(dataTop),1),thisTflowList,thisfilelist]


def ReadListTopCharge(thisiSmearList,thisjSmearList,thisMomList,thisProjGammaList,thisProjDerList,thisDSList,
                      thisTSinkList,thisconflist,Flag,Interps=['nucleon'],thistsourceList=[tsource],Wein=False):
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
    thisfilelist,topcharge,tflow = ReadAndCheckTop(Wein,thisfilelist)
    data2pt,randlist,shiftlist = Read2ptSet(thisfilelist,thisiSmearList,thisjSmearList,GetAvgMomListip(thisMomList),Interps,tsourceList=thistsourceList)
    # print ''
    if not np.all([x==tflow[0] for x in tflow]):
        print 'warning, files had different flow times'
    thisTflowList = tflow[0]

    data3pt,dataTop = [],[]
    for iFlag,itsink in zip(Flag,thisTSinkList):
        holdTop,hold3pt = Read3ptSetTop(topcharge,thisfilelist,thisiSmearList,thisjSmearList,thisMomList,thisProjGammaList,
                                        thisProjDerList,thisDSList,itsink,iFlag,shiftlist,cfglistout,thisTflowList,
                                        randlist=randlist,tsourceList=thistsourceList)
        data3pt.append(hold3pt)
        dataTop.append(holdTop)
    # print ''
    return [data2pt,data3pt,np.rollaxis(np.array(dataTop),1),thisTflowList,thisfilelist]


def TestMomList(thisMomList):
    try:
        PrintZMomI(np.where(np.array(thisMomList)==iqTOip(0))[0][0])
    except IndexError:
        print "Must Read Zero Momentum value " + iqTOip(0)

        
def ReadListAlpha(thisiSmearList,thisjSmearList,thisMomList,thisconflist,Interps=['nucleon'],thistsourceList=[tsource],Wein=False):
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

    
    thisfilelist,topcharge,tflow = ReadAndCheckTop(Wein,thisfilelist)
    print thisfilelist.keys()
    if not np.all([x==tflow[0] for x in tflow]):
        print 'warning, files had different flow times'
    thisTflowList = tflow[0]
    # if len(cfglistout) != len(thiscfglist):
    #     print 'Not All Top Charge Found'
        
    dataTop,data2pt,randlist,shiftlist = Read2ptSetTop(topcharge,thisfilelist,thisiSmearList,thisjSmearList,GetAvgMomListip(thisMomList),
                                                       Interps,cfglistout,thisTflowList,tsourceList=thistsourceList)
    return [data2pt,dataTop,thisTflowList,thisfilelist]



def ReadSetAlpha(thisiSmearList,thisjSmearList,thisMomList,directory,Interps=['nucleon'],thistsourceList=[tsource],Wein=False):
    # if len(thisTSinkList) > 0:
    #     TestMomList(thisMomList)
    thisfilelist = OrderedDict()
    f = open('./cfglistset.txt','w')
    
    fileend2pt = CreateEnd2pt(DefiSmearList[0],DefjSmearList[0],thistsourceList[0],Interps[0],Interps[0])
    if Debug: print ' comparing to: ',fileend2pt
    for isource in SourceList:
        print 'reading directory: ' , directory+'/'+isource
        for (dirname,dirs,files) in walk(directory+'/'+isource):
            if 'Temp' in dirname: continue
            for ifile in filesortCfuns(files):
                if fileend2pt not in ifile or ".metadata" in ifile: continue
                # ## FOR DEBUGGING
                # if xsrcList[0] not in ifile: continue
                if 'xsrc' in ListOrSet:
                    if ListOrSet.replace('ReadSet','').replace('ReadList','')+'_' not in ifile: continue
                fileprefix = ifile.replace(fileend2pt,'')
                if CfunConfigCheck:
                    if CheckAllSet(fileprefix,dirname+'/',Interps,tsourceList=thistsourceList):
                        prefnosrc = re.sub('_xsrc.*','',fileprefix)
                        if prefnosrc not in thisfilelist.keys():
                            thisfilelist[prefnosrc] = [directory+'/'+isource+'/@/'+fileprefix]
                        else:
                            if XSrcLen <= len(thisfilelist[prefnosrc]) and ForceXSrcLen: continue
                            thisfilelist[prefnosrc].append(directory+'/'+isource+'/@/'+fileprefix)
                        f.write(directory+'/'+isource+'/@/'+fileprefix+'\n')
                else:
                    if CheckSet(fileprefix,dirname+'/',thisiSmearList,thisjSmearList,{},{},[],[],'',Interps,tsourceList=thistsourceList):
                        prefnosrc = re.sub('_xsrc.*','',fileprefix)
                        if prefnosrc not in thisfilelist.keys():
                            thisfilelist[prefnosrc] = [directory+'/'+isource+'/@/'+fileprefix]
                        else:
                            if XSrcLen <= len(thisfilelist[prefnosrc]) and ForceXSrcLen: continue
                            thisfilelist[prefnosrc].append(directory+'/'+isource+'/@/'+fileprefix)
                        f.write(directory+'/'+isource+'/@/'+fileprefix+'\n')
    f.close()
    thisfilelist = DoForceChecks(thisfilelist)
    thisfilelist,topcharge,tflow = ReadAndCheckTop(Wein,thisfilelist)
    print 'number of configs = ' , len(thisfilelist.keys())
    print 'average number of sources per cfg = ' ,np.mean([len(ifilelist) for ifilelist in thisfilelist.itervalues()])
    print 'total number of measurements = ' , np.sum([len(ifilelist) for ifilelist in thisfilelist.itervalues()])
    print ''
    if len(thisfilelist) == 0: raise IOError('No Configurations Found')
    if ShowConfNum:
        for ifile in thisfilelist.itervalues():
            for iifile in ifile:
                print iifile
        print ''
        
    if len(tflow) == 0:
        if Wein:
            raise IOError('Weinerg operator data not found in ' + TopReadDir )
        else:
            raise IOError('Topological charge data not found in ' + TopReadDir )            
    if not np.all([x==tflow[0] for x in tflow]):
        print 'warning, files had different flow times'
    thisTflowList = tflow[0]
    # if len(cfglistout) != len(thiscfglist):
    #     print 'Not All Top Charge Found'

    dataTop,data2pt,randlist,shiftlist = Read2ptSetTop(topcharge,thisfilelist,thisiSmearList,thisjSmearList,
                                                       GetAvgMomListip(thisMomList),Interps,cfglistout,thisTflowList,tsourceList=thistsourceList)
    
    return [data2pt,dataTop,thisTflowList,thisfilelist]

        
def ReadList(thisiSmearList,thisjSmearList,thisMomList,thisProjGammaList,thisProjDerList,thisDSList,
             thisTSinkList,thisconflist,Flag,Interps=['nucleon'],thistsourceList=[tsource]):
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
    thisfilelist = CheckTop(thisfilelist)
    data2pt,randlist,shiftlist = Read2ptSet(thisfilelist,thisiSmearList,thisjSmearList,GetAvgMomListip(thisMomList),Interps,tsourceList=thistsourceList)
    # print ''
    data3pt = []
    for iFlag,itsink in zip(Flag,thisTSinkList):
        data3pt.append(Read3ptSet(thisfilelist,thisiSmearList,thisjSmearList,thisMomList,thisProjGammaList,
                                  thisProjDerList,thisDSList,itsink,iFlag,shiftlist,randlist=randlist,tsourceList=thistsourceList))
    # print ''
    return [data2pt,data3pt,thisfilelist]

def ReadSet(thisiSmearList,thisjSmearList,thisMomList,thisProjGammaList,thisProjDerList, thisDSList,
            thisTSinkList,directory,Flag,Interps=['nucleon'],thistsourceList=[tsource]):
    # if len(thisTSinkList) > 0:
    #     TestMomList(thisMomList)
    thisfilelist = OrderedDict()
    f = open('./cfglistset.txt','w')
    
    fileend2pt = CreateEnd2pt(DefiSmearList[0],DefjSmearList[0],thistsourceList[0],Interps[0],Interps[0])
    if Debug: print ' comparing to: ',fileend2pt

    for isource in SourceList:
        print 'reading directory: ' , directory+'/'+isource
        for (dirname,dirs,files) in walk(directory+'/'+isource):
            if 'Temp' in dirname: continue
            for ifile in filesortCfuns(files):
                if fileend2pt not in ifile or ".metadata" in ifile: continue
                if 'xsrc' in ListOrSet:
                    if ListOrSet.replace('ReadSet','').replace('ReadList','')+'_' not in ifile: continue
                fileprefix = ifile.replace(fileend2pt,'')
                if CfunConfigCheck:
                    # ## FOR DEBUGGING
                    # if xsrcList[0] not in ifile: continue
                    if CheckAllSet(fileprefix,dirname+'/',Interps,tsourceList=thistsourceList):
                        prefnosrc = re.sub('_xsrc.*','',fileprefix)
                        if prefnosrc not in thisfilelist.keys():
                            thisfilelist[prefnosrc] = [directory+'/'+isource+'/@/'+fileprefix]
                        else:
                            if XSrcLen <= len(thisfilelist[prefnosrc]) and ForceXSrcLen: continue
                            thisfilelist[prefnosrc].append(directory+'/'+isource+'/@/'+fileprefix)
                        f.write(directory+'/'+isource+'/@/'+fileprefix+'\n')
                else:
                    if CheckSet(fileprefix,dirname+'/',thisiSmearList,thisjSmearList,thisProjGammaList,
                                thisProjDerList,thisDSList,thisTSinkList,Flag,Interps,tsourceList=thistsourceList):
                        prefnosrc = re.sub('_xsrc.*','',fileprefix)
                        if prefnosrc not in thisfilelist.keys():
                            thisfilelist[prefnosrc] = [directory+'/'+isource+'/@/'+fileprefix]
                        else:
                            if XSrcLen <= len(thisfilelist[prefnosrc]) and ForceXSrcLen: continue
                            thisfilelist[prefnosrc].append(directory+'/'+isource+'/@/'+fileprefix)
                        f.write(directory+'/'+isource+'/@/'+fileprefix+'\n')
    f.close()
    thisfilelist = DoForceChecks(thisfilelist)
    thisfilelist = CheckTop(thisfilelist)
    print 'number of configs = ' , len(thisfilelist.keys())
    print 'average number of sources per cfg = ' ,np.mean([len(ifilelist) for ifilelist in thisfilelist.itervalues()])
    print 'total number of measurements = ' , np.sum([len(ifilelist) for ifilelist in thisfilelist.itervalues()])
    print ''
    if len(thisfilelist.keys()) == 0: raise IOError('No Configurations Found')
    if ShowConfNum:
        for ifile in thisfilelist.itervalues():
            for iifile in ifile:
                print iifile
        print ''
    # print 'number of configs = ' , len(thisfilelist)
    # print ''
    # if len(thisfilelist) == 0: raise IOError('No Configurations Found')
    # if ShowConfNum:
    #     for ifile in thisfilelist:
    #         print ifile
    #     print ''
        
    data2pt,randlist,shiftlist = Read2ptSet(thisfilelist,thisiSmearList,thisjSmearList,GetAvgMomListip(thisMomList),Interps,tsourceList=thistsourceList)
    # print ''
    data3pt = []
    for iFlag,itsink in zip(Flag,thisTSinkList):
        data3pt.append(Read3ptSet(thisfilelist,thisiSmearList,thisjSmearList,thisMomList,thisProjGammaList
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
    for iterp,ism in Elongate(Interps,DefiSmearList):
        for jcsm,(jterp,jsm) in enumerate(Elongate(Interps,DefjSmearList)):
            for thists in tsourceList:
                testfile2pt = (directory.replace(CreateDir2pt(DefiSmearList[0],DefjSmearList[0]),
                                                 CreateDir2pt(ism,jsm))
                               +FilePrefix+CreateEnd2pt(ism,jsm,thists,iterp,jterp))
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
                                testfile3pt = (directory.replace(CreateDir2pt(DefiSmearList[0],DefjSmearList[0]),
                                                                 CreateDir3pt(ism,jsm3pt,itsink,iDS,iProj,thisFlag))
                                               +FilePrefix+CreateEnd3pt(ism,jsm3pt,thists,itsink,iDS,iProj,''))
                                if Debug: print 'Checking!!! ' ,  testfile3pt.replace(FileStruct,FileStruct+C2C3Dis)
                                if not os.path.isfile(testfile3pt.replace(FileStruct,FileStruct+C2C3Dis)): return False
                            for iProj in DefProjDerList:
                                testfile3pt = (directory.replace(CreateDir2pt(DefiSmearList[0],DefjSmearList[0]),
                                                                 CreateDir3pt(ism,jsm3pt,itsink,iDS,iProj,thisFlag))
                                               +FilePrefix+CreateEnd3pt(ism,jsm3pt,thists,itsink,iDS,iProj,'D'))
                                if Debug: print 'Checking!!! ' ,testfile3pt.replace(FileStruct,FileStruct+C2C3Dis)
                                if not os.path.isfile(testfile3pt.replace(FileStruct,FileStruct+C2C3Dis)): return False
    return True
                


def CheckSet(FilePrefix,directory,thisiSmearList,thisjSmearList,thisProjGammaList,thisProjDerList,thisDSList,
             thisTSinkList,Flag,Interps,tsourceList=[tsource]):
    for iterp,ism in Elongate(Interps,thisiSmearList):
        for thists in tsourceList:
            for jcsm,(jterp,jsm) in enumerate(Elongate(Interps,thisjSmearList)):
                testfile2pt = (directory.replace(CreateDir2pt(thisiSmearList[0],thisjSmearList[0]),
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
                    Jsmlist = ['Xsm'+jsm for jsm in thisjSmearList]
                for jcsm,jsm3pt in enumerate(Jsmlist):
                    for iDS in thisDSList:
                        for iProj in thisProjGammaList:
                            testfile3pt = (directory.replace(CreateDir2pt(thisiSmearList[0],thisjSmearList[0]),
                                                             CreateDir3pt(ism,jsm3pt,itsink,iDS,iProj,thisFlag))
                                           +FilePrefix+CreateEnd3pt(ism,jsm3pt,thists,itsink,iDS,iProj,''))
                            # testfile3pt = testfile3pt.replace(xsrcList[0],xsrc)
                            thisfile3pt = testfile3pt.replace(FileStruct,FileStruct+C2C3Dis)
                            if Debug: print 'Checking~!!! ' ,  thisfile3pt
                            if not os.path.isfile(thisfile3pt): return False
                        for iProj in thisProjDerList:
                            testfile3pt = (directory.replace(CreateDir2pt(thisiSmearList[0],thisjSmearList[0]),
                                                             CreateDir3pt(ism,jsm3pt,itsink,iDS,iProj,thisFlag))
                                           +FilePrefix+CreateEnd3pt(ism,jsm3pt,thists,itsink,iDS,iProj,'D'))
                            # testfile3pt = testfile3pt.replace(xsrcList[0],xsrc)
                            thisfile3pt = testfile3pt.replace(FileStruct,FileStruct+C2C3Dis)
                            if Debug: print 'Checking~!!! ' ,  thisfile3pt
                            if not os.path.isfile(testfile3pt): return False
    return True

##tempdata [ iconf ] .data [ ip , it ]
##thisdata2pt [ tsource, ism , jsm , ip , it ] bootdataclas
##shiftlist = [ tsource , ism, jsm , icfg, isrc_number ]
def Read2ptSet(readfilelist,thisiSmearList,thisjSmearList,thisMomList,Interps,tsourceList=[tsource]):
    thisdata2pt = []
    start = time.time()
    icount = -1
    randlist = []
    shiftlist = []
    for thists in tsourceList:
        thisdata2pt.append([])
        shiftlist.append([])
        for icsm,(iterp,ism) in enumerate(Elongate(Interps,thisiSmearList)):
            thisdata2pt[-1].append([])
            shiftlist[-1].append([])
            for jcsm,(jterp,jsm) in enumerate(Elongate(Interps,thisjSmearList)):
                icount += 1
                timeleft = GetTimeLeft(icount,len(tsourceList)*len(thisiSmearList)*len(thisjSmearList)*len(Interps)**2,time.time()-start)
                print 'Read 2pt: sm' + ism + 'Xsm'+jsm+' Time Left: ' +str(datetime.timedelta(seconds=timeleft)) , ' h:m:s \r',
                thisstart = time.time()
                thisfilelist = OrderedDict()
                for ireadkey,ireadfl in readfilelist.iteritems():
                    thisfilelist[ireadkey] = [ifile.replace('@',CreateDir2pt(ism,jsm))+CreateEnd2pt(ism,jsm,thists,iterp,jterp) for ifile in ireadfl]
                # thisfilelist = [ifile.replace('@',CreateDir2pt(ism,jsm))+CreateEnd2pt(ism,jsm,thists,iterp,jterp) for ifile in readfilelist]
                (dataout,randlist),thisshiftlist = ReadAndBoot2pt(thisfilelist,thisMomList,nboot) 
                shiftlist[-1][icsm].append(thisshiftlist)
                thisdata2pt[-1][icsm].append(dataout)
                print 'Read 2pt: sm' + ism + 'Xsm'+jsm + ' t_src'+str(thists)+' took: ' +str(datetime.timedelta(seconds=time.time()-thisstart)) , ' h:m:s        '
    print 'Read 2pt total time: ' + str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s'
    return thisdata2pt,randlist,shiftlist


##thisdata3pt [ thistsource, ism , jsm , igamma , ip ] bootdataclas
def Read3ptSet(readfilelist,thisiSmearList,thisjSmearList,thisMomList,thisProjGammaList,thisProjDerList,thisDSList,tsink,iFlag,shiftlist,randlist=[],tsourceList=[tsource]):
    thisdata3pt = []
    start = time.time()
    icount = -1
    for itsource,thists in enumerate(tsourceList):
        thisdata3pt.append([])
        for icsm,ism in enumerate(thisiSmearList):
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
                Jsmlist = ['Xsm'+jsm for jsm in thisjSmearList]
            for jcsm,jsm in enumerate(Jsmlist):
                thisdata3pt[-1][icsm].append([])
                thisstart = time.time()
                icount += 1
                timeleft = GetTimeLeft(icount,len(tsourceList)*len(thisiSmearList)*len(Jsmlist),time.time()-start)
                print 'Read 3pt: sm' + ism + jsm + ' ts'+str(tsink)+' Time Remaining: '+ str(datetime.timedelta(seconds=timeleft)) , ' h:m:s \r',
                thisshiftlist = shiftlist[itsource][icsm][jcsm]
                for iDS in thisDSList:
                    for iProj,thisGammaList in thisProjGammaList.iteritems():
                        thisfilelist = OrderedDict()
                        for ireadkey,ireadfl in readfilelist.iteritems():
                            thisfilelist[ireadkey] = [ifile.replace('@',CreateDir3pt(ism,jsm,tsink,iDS,iProj,thisFlag))
                                                      .replace(FileStruct,FileStruct+C2C3Dis)
                                                      +CreateEnd3pt(ism,jsm,thists,tsink,iDS,iProj,'') for ifile in ireadfl]

                        holdres = ReadAndBoot3pt(thisfilelist,thisMomList,thisGammaList,[],nboot,printstr=iDS+iProj,randlist=randlist)
                        for gammares in holdres:
                            thisdata3pt[-1][icsm][jcsm].append(gammares)
                    for iProj,thisDerList in thisProjDerList.iteritems():
                        thisfilelist = OrderedDict()
                        for ireadkey,ireadfl in readfilelist.iteritems():
                            thisfilelist = [ifile.replace('@',CreateDir3pt(ism,jsm,tsink,iDS,iProj,thisFlag))
                                            .replace(FileStruct,FileStruct+C2C3Dis)
                                            +CreateEnd3pt(ism,jsm,thists,tsink,iDS,iProj,'D') for ifile in ireadfl]
                        holdres = ReadAndBoot3pt(thisfilelist,thisMomList,[],thisDerList,nboot,printstr=iDS+iProj,randlist=randlist)
                        for gammares in holdres:
                            thisdata3pt[-1][icsm][jcsm].append(gammares)
                print 'Read 3pt: sm' + ism + jsm + 't_src'+str(thists)+' t_sink'+str(tsink)+ ' took: ' + str(datetime.timedelta(seconds=time.time()-thisstart)) , ' h:m:s        '
    print 'Read 3pt total time: ' + str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s'
    return thisdata3pt





##tempdata [ iconf ] .data [ ip , it ]
##thisdata2pt [ tflow, tsource, ism , jsm , ip , it ] bootdataclas
##shiftlist = [ tsource , ism, jsm , icfg, isrc_number ]
##TopCharge = [ iconf , tflow ] 
def Read2ptSetTop(TopCharge,readfilelist,thisiSmearList,thisjSmearList,thisMomList,Interps,Topcfglist,thisTflowList,tsourceList=[tsource]):
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
        for icsm,(iterp,ism) in enumerate(Elongate(Interps,thisiSmearList)):
            thisdata2pt[-1].append([])
            thisdataTop[-1].append([])
            shiftlist[-1].append([])
            for jcsm,(jterp,jsm) in enumerate(Elongate(Interps,thisjSmearList)):
                icount += 1
                timeleft = GetTimeLeft(icount,len(tsourceList)*len(thisiSmearList)*len(thisjSmearList)*len(Interps)**2,time.time()-start)
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



##thisdata3pt [ thistsource, ism , jsm , igamma , ip ] bootdataclas
##thisdataTop [ tflow, thistsource, ism , jsm , igamma , ip ] bootdataclas
def Read3ptSetTop(TopCharge,readfilelist,thisiSmearList,thisjSmearList,thisMomList,thisProjGammaList,thisProjDerList,thisDSList,tsink,iFlag,shiftlist,
                  Topcfglist,thisTflowList,randlist=[],tsourceList=[tsource]):
    thisdata3pt = []
    thisdataTop = []
    start = time.time()
    icount = -1
    for itsource,thists in enumerate(tsourceList):
        thisdata3pt.append([])
        thisdataTop.append([])
        for icsm,ism in enumerate(thisiSmearList):
            thisdata3pt[itsource].append([])
            thisdataTop[itsource].append([])
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
                Jsmlist = ['Xsm'+jsm for jsm in thisjSmearList]
            for jcsm,jsm in enumerate(Jsmlist):
                thisdata3pt[itsource][icsm].append([])
                thisdataTop[itsource][icsm].append([])
                thisstart = time.time()
                icount += 1
                timeleft = GetTimeLeft(icount,len(tsourceList)*len(thisiSmearList)*len(Jsmlist),time.time()-start)
                print 'Read 3pt: sm' + ism + jsm + ' ts'+str(tsink)+' Time Remaining: '+ str(datetime.timedelta(seconds=timeleft)) , ' h:m:s \r',
                thisshiftlist = shiftlist[itsource][icsm][jcsm]
                for iDS in thisDSList:
                    for iProj,thisGammaList in thisProjGammaList.iteritems():
                        thisfilelist = OrderedDict()
                        for ireadkey,ireadfl in readfilelist.iteritems():
                            thisfilelist[ireadkey] = [ifile.replace('@',CreateDir3pt(ism,jsm,tsink,iDS,iProj,thisFlag))
                                                      .replace(FileStruct,FileStruct+C2C3Dis)
                                                      +CreateEnd3pt(ism,jsm,thists,tsink,iDS,iProj,'') for ifile in ireadfl]

                        holdresTop,holdres = ReadAndBoot3ptTop(thisfilelist,thisMomList,thisGammaList,[],nboot,TopCharge,Topcfglist,thisTflowList,
                                                               printstr=iDS+iProj,randlist=randlist)
                        for ig,(gammares,gammaTopres) in enumerate(zip(holdres,holdresTop)):
                            thisdata3pt[itsource][icsm][jcsm].append(gammares)
                            thisdataTop[itsource][icsm][jcsm].append(gammaTopres)

                    # for iProj,thisDerList in thisProjDerList.iteritems():
                    #     thisfilelist = OrderedDict()
                    #     for ireadkey,ireadfl in readfilelist.iteritems():
                    #         thisfilelist = [ifile.replace('@',CreateDir3pt(ism,jsm,tsink,iDS,iProj,thisFlag))
                    #                         .replace(FileStruct,FileStruct+C2C3Dis)
                    #                         +CreateEnd3pt(ism,jsm,thists,tsink,iDS,iProj,'D') for ifile in ireadfl]
                    #     holdresTop,holdres = ReadAndBoot3ptTop(thisfilelist,thisMomList,[],thisDerList,nboot,TopCharge,Topcfglist,thisTflowList,
                    #                                         printstr=iDS+iProj,randlist=randlist)
                    #     for ig,gammares in enumerate(holdres):
                    #         thisdata3pt[itsource][icsm][jcsm].append(gammares)
                    #         thisdataTop[itsource][icsm][jcsm].append([])
                    #         for iflow,flowres in enumerate(holdresTop):
                    #             thisdataTop[itsource][icsm][jcsm][ig].append(flowres[ig])
                print 'Read 3pt: sm' + ism + jsm + 't_src'+str(thists)+' t_sink'+str(tsink)+' took: ' + str(datetime.timedelta(seconds=time.time()-thisstart)) , ' h:m:s        '
    print 'Read 3pt total time: ' + str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s'
    return np.rollaxis(np.array(thisdataTop),4),thisdata3pt



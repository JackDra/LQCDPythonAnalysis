#!/usr/bin/env python
import os
from os import walk
from Params import *
from ReadBinaryCfuns import ReadFSCfun, Read2ptCfun
from CfunBoot.py import BootSet2pt, BootSet3pt

## data2pt [icon][ism][jsm]
## data3pt [icon][ism][jsm][itsink]

def ReadCMList(iDS,iProj,Der,thisconflist,Flag):
    data2pt = []
    data3pt = []
    for iconf in thisconflist:
        file = filelist.replace('*',iconf)
        print file
        data3pt.append(Read3ptSet(file,iDS,iProj,Der,Flag))
        data2pt.append(Read2ptSet(file))
    return [data2pt,data3pt]

def ReadCMSet(iDS,iProj,Der, directory,Flag):
    thisfilelist = []
    data2pt = []
    data3pt = []
    for isource in SourceList:
        for (dirname,dirs,files) in walk(directory+'/'+isource):
            for file in files:
                fileend3pt = CreateEnd3pt(SmearSet[0],SmearSet[0],TSinkSet[0],iDS,iProj,Der)
                fileend2pt = CreateEnd2pt(SmearSet[0],SmearSet[0])
                if fileend2pt not in file or ".metadata" in file: continue
                fileprefix = file.replace(fileend2pt,'')
                # print dirname + '/'+fileprefix
                
                if CheckSet(fileprefix,dirname+'/',iDS,iProj,Der,Flag):
                    print fileprefix
                    thisfilelist.append(directory+'/'+isource+'/@/'+fileprefix)
                    data3pt.append(Read3ptSet(directory+'/'+isource+'/@/'+fileprefix,iDS,iProj,Der,Flag))
                    data2pt.append(Read2ptSet(directory+'/'+isource+'/@/'+fileprefix))
    return [data2pt,data3pt,thisfilelist]
                
def CreateEnd3pt(ism,jsm,Tsink,DS,Proj,Der):
    return '.t'+str(tsource)+'sm'+ism+'Xsm'+jsm+'GMA'+Proj+'t'+Tsink+'p000.'+DS+Der+'k'+str(kappa)+'.ifms.3cf'

def CreateEnd2pt(ism,jsm):
    return '.t'+str(tsource)+'sm'+ism+'k'+str(kappa)+'.ifmssi'+jsm+'.nucleon.u.2cf'

def CreateDir3pt(ism,jsm,Tsink,DS,Proj,Flag):
    return Flag+'sm'+ism+'Xsm'+jsm+'GMA'+Proj+'t'+Tsink+'p000.'+DS

def CreateDir2pt(ism,jsm):
    return 'twoptsm'+ism+'si'+jsm

def CheckSet(FilePrefix,directory,iDS,iProj,Der,Flag):
    CompSet = True
    for ism in SmearSet:
        for jsm in SmearSet:
            testfile2pt = (directory.replace(CreateDir2pt(SmearSet[0],SmearSet[0]),
                                             CreateDir2pt(ism,jsm))
                           +FilePrefix+CreateEnd2pt(ism,jsm))
            if not os.path.isfile(testfile2pt): CompSet = False
            for itsink in TSinkSet:
                testfile3pt = (directory.replace(CreateDir2pt(SmearSet[0],SmearSet[0]),
                                                 CreateDir3pt(ism,jsm,itsink,iDS,iProj,Flag))
                               +FilePrefix+CreateEnd3pt(ism,jsm,itsink,iDS,iProj,Der))
                if not os.path.isfile(testfile3pt): CompSet = False
    return CompSet

def Read2ptSet(filename):
    thisdata2pt = []
    for icsm,ism in enumerate(SmearSet):
        thisdata2pt.append([])
        for jsm in SmearSet:            
            thisdata2pt[icsm].append(Read2ptCfun(filename.replace('@',CreateDir2pt(ism,jsm))+CreateEnd2pt(ism,jsm)))
    return thisdata2pt
        
def Read3ptSet(filename,iDS,iProj,Der,Flag):
    thisdata3pt = []
    for icsm,ism in enumerate(SmearSet):
        thisdata3pt.append([])
        for jcsm,jsm in enumerate(SmearSet):            
            thisdata3pt[icsm].append([])
            for itsink in TSinkSet:
                thisdata3pt[icsm][jcsm].append(ReadFSCfun(filename.replace('@',CreateDir3pt(ism,jsm,itsink,iDS,iProj,Flag))+CreateEnd3pt(ism,jsm,itsink,iDS,iProj,Der)))
    return thisdata3pt

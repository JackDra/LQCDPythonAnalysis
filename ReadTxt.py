#!/usr/bin/env python

import numpy as np
from BootTest import BootStrap1
from RFCalc import *
from collections import OrderedDict
from SetLists import CutDupSet
from FitParams import *
from SetLists import *
from Params import *
import FitFunctions as ff
from OppFuns import CreateOppDir
from CreateCombs import ops
import os.path
from copy import deepcopy, copy
import re
import time
import datetime
from ReadXml import *


# R/L Evecs [ ip , istate , ival ]
# Emass [ ip , istate ]
def ReadLREM(todtval,thisMomList,filepref):
    filename = REvecDir+filepref+'to'+str(todtval[0])+'dt'+str(todtval[1])+'LREM.txt'
    # filename = REvecDir+VarPref+'TestLREM.txt'

    readmom = False
    LEvec,REvec,Emass = [],[],[]
    MomListOut = []
    if os.path.isfile(filename):
        with open(filename,'r') as f:
            for line in f:
                rdata = line.strip()
                if rdata[0] == 'q':
                    if qstrTOip(rdata) in thisMomList:
                        LEvec.append([])
                        REvec.append([])
                        Emass.append([])
                        MomListOut.append(rdata)
                        readmom = True
                    else:
                        readmom = False
                elif readmom and rdata[0] != 's':
                    rdata = rdata.split()
                    if 'L' in rdata[0] and 'R' not in rdata[0]:
                        LEvec[-1].append(map(complex,rdata[1:-1]))
                    elif 'R' in rdata[0] and 'L' not in rdata[0]:
                        REvec[-1].append(map(complex,rdata[1:-1]))
                        Emass[-1].append(complex(rdata[-1]))
                    elif 'R' in rdata[0] and 'L' in rdata[0]:
                        LEvec[-1].append(map(complex,rdata[1:-1]))
                        REvec[-1].append(map(complex,rdata[1:-1]))
                        Emass[-1].append(complex(rdata[-1]))                        
        return MomOrderLists(MomListOut,[ipTOqstr(ip) for ip in thisMomList],LEvec,REvec,Emass)
    else:
        print 'Warning',filename, 'not found'
        return None,None,None

## readdata { gamma } { mom } { method } { set }
## datadictout { collection } { gamma } { mom }
## datamassout { collection }
def ExtractValues(thisindir,thisGammaList,thisSetList,thisMethodList,thisMomList=qvecSet,thisPrintRead=PrintRead):
    def SetupDict(dictin,thisgamma,thiscol):
        if thiscol not in dictin.keys():
            dictin[thiscol] = {}
        if thisgamma not in dictin[thiscol].keys():
            dictin[thiscol][thisgamma] = {}
        return dictin
    readdata = ReadSetDict(thisindir,thisSetList,thisGammaList,thisMethodList,thisMomList=thisMomList,thisPrintRead=thisPrintRead)
    start = time.time()
    datadictout,datamassout = OrderedDict(),OrderedDict()
    for ig,igamma in enumerate(thisGammaList):
        if thisPrintRead: print ' Extracting data from read: ' , int((ig*100)/float(len(thisGammaList))) ,'%  \r',
        if igamma == 'twopt' or igamma == 'Mass':
            zmom = 'q = 0 0 0' 
            if not CheckDict(readdata,igamma,zmom): continue
            for iMeth,Methdata in readdata[igamma][zmom].iteritems():
                for ism,thismassdict in Methdata.iteritems():
                    if DefTvarPicked in ism: ism = PickedStateStr+DefTvarPicked
                    if 'TSF' in iMeth:
                        if CheckDict(thismassdict,'m0',TSFfitr):
                            datamassout[ism+iMeth] = thismassdict['m0'][TSFfitr]
                    elif 'OSF' in iMeth:
                        if CheckDict(thismassdict,'m0',OSFfitr[RemoveToDt(ism)]):
                            datamassout[ism+iMeth] = thismassdict['m0'][OSFfitr[RemoveToDt(ism)]]
        else:
            for imom,momdata in readdata[igamma].iteritems():
                for iMeth,Methdata in momdata.iteritems():
                    for iSet,thisdict in Methdata.iteritems():
                        itsink,ism = SplitTSinkString(iSet)
                        fitsm = RemoveToDt(ism)
                        if DefTvarPicked in ism: ism = PickedStateStr+DefTvarPicked
                        if 'TSF' in iMeth:
                            for icut in TSFCutList:
                                if CheckDict(thisdict,'B00',TSFfitr,icut):
                                    datadictout = SetupDict(datadictout,igamma,iSet+iMeth+icut)
                                    datadictout[iSet+iMeth+icut][igamma][imom] = thisdict['B00'][TSFfitr][icut]
                        elif 'OSF' in iMeth:
                            for icut in OSFCutList:
                                if CheckDict(thisdict,'B00',OSFfitr[fitsm],icut):
                                    datadictout = SetupDict(datadictout,igamma,iSet+iMeth+icut)
                                    datadictout[iSet+iMeth+icut][igamma][imom] = thisdict['B00'][OSFfitr[fitsm]][icut]
                        elif 'SumMeth' in iMeth:
                            for ifit in SumFitRList:
                                for icut in SumCutList:
                                    if CheckDict(thisdict,icut,ifit):
                                        fitdict = ifit.replace('fit sl ','fitr')
                                        datadictout = SetupDict(datadictout,igamma,iSet+iMeth+icut+fitdict)
                                        datadictout[iSet+iMeth+icut+fitdict][igamma][imom] = thisdict[icut][ifit]
                        elif 'Fits' in iMeth:
                            for icut in FitCutArgs:
                                if icut in thisdict.keys():
                                    datadictout = SetupDict(datadictout,igamma,iSet+iMeth+icut)
                                    datadictout[iSet+iMeth+icut][igamma][imom] = thisdict[icut]
    if thisPrintRead: print 'Extracting data took: ', str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s                  '
    return datadictout,datamassout


def Get2ptSetMoms(outputdir,MomListIn,statelist=[],todtlist=[],smlist=[]):
    momlist = set([])
    xmlMomList = map(qstrTOqcond,MomListIn)
    for iflag in ['cfuns/twopt','Mass']:
        for ip in xmlMomList:
            thisdir = outputdir+iflag+MakeMomDir(ip)
            for itodt in todtlist:
                if iflag == 'Mass':
                    ifile = thisdir+itodt+'LREM'+ip+'.xml'
                    if not CheckMomFile(ifile): momlist.add(qcondTOqstr(ip))
                for istate in statelist:
                    ifile = thisdir+MakeMomDir(ip)+'state'+istate+itodt+iflag.replace('cfuns/','')+ip+'.xml'
                    if not CheckMomFile(ifile): momlist.add(qcondTOqstr(ip))
            for ism in smlist:
                ifile = thisdir+MakeMomDir(ip)+ism+iflag.replace('cfuns/','')+ip+'.xml'
                if not CheckMomFile(ifile): momlist.add(qcondTOqstr(ip))
    return OrderMomList(momlist)


def Get3SM(outputdir,thisGammaList,MomListIn,setlist):
    momlist = set([])
    xmlMomList = map(qstrTOqcond,MomListIn)
    for igamma in thisGammaList:
        thisdir = outputdir+CreateOppDir(igamma)
        for iset in setlist:
            for ip in xmlMomList:
                ifile = thisdir+MakeMomDir(ip)+iset+igamma+ip+'.xml'
                if not CheckMomFile(ifile): momlist.add(qcondTOqstr(ip))
    return OrderMomList(momlist)
            

def Get3ptSetMoms(outputdir,thisGammaList,MomListIn,setlist):
    return OrderMomList(set(Get3SM(outputdir,thisGammaList,MomListIn,setlist)) |
                        set(Get3SM(outputdir+'cfuns/',thisGammaList,MomListIn,setlist)))


##NEW FUNCTIONS##

#data2pt = [ mom , set , it ] bs
#data3pt = [ gamma , mom , set , it ] bs
def ReadCfunsnp(thisGammaList,thisSetList,thisMomList=RunMomList,thisPrintRead=PrintRead):
    datadict = ReadCfunsDict(outputdir,thisSetList,thisGammaList,thisMomList=thisMomList,thisPrintRead=thisPrintRead)
    return SetRFDictToList(datadict,thisPrintRead=thisPrintRead)

def ReadRFnp(thisGammaList,thisSetList,thisMomList=RunMomList,thisPrintRead=PrintRead):
    datadict = ReadSetDict(outputdir,thisSetList,thisGammaList,['RF'],thisMomList=thisMomList,thisPrintRead=thisPrintRead)
    return SetRFDictToList(datadict,thisPrintRead=thisPrintRead)

def ReadSetAndCfunsDict(thisindir,thisSetList,thisGammaList,thisMethodList,thisMomList=RunMomList,thisPrintRead=PrintRead):
    datadict = ReadSetDict(thisindir,thisSetList,thisGammaList,thisMethodList,thisMomList=thisMomList,thisPrintRead=thisPrintRead)
    data3ptdict = ReadCfunsDict(thisindir,thisSetList,thisGammaList,thisMomList=thisMomList)
    for igamma in datadict.keys():
        for imom in datadict[igamma].keys():
            datadict[igamma][imom]['C3pt'] = data3ptdict[igamma][imom]['RF']
    return datadict

def RewriteRF(RFdict,threeptdict,thisopp,thismom):
    if CheckDict(threeptdict,thisopp,thismom,'RF'):
        for igamma in RFdict.keys():
            if 'twopt' in igamma or 'P4g4' in igamma: continue
            if CheckDict(RFdict,igamma,thismom,'RF'):
                for iset in RFdict[igamma][thismom]['RF'].keys():
                    for it,(i3pt,i3ptopp) in enumerate(zip(threeptdict[igamma][thismom]['RF'][iset]['Boot'][tsource-1:GetintTSink(iset)],
                                                           threeptdict[thisopp][thismom]['RF'][iset]['Boot'][tsource-1:GetintTSink(iset)])):
                        RFdict[igamma][thismom]['RF'][iset]['Boot'][it].values = i3pt.values/i3ptopp.values
                        GetBootStats(RFdict[igamma][thismom]['RF'][iset]['Boot'][it])
                        RFdict[igamma][thismom]['RF'][iset]['Vals'][it] = RFdict[igamma][thismom]['RF'][iset]['Boot'][it].Avg
                        RFdict[igamma][thismom]['RF'][iset]['Valserr'][it] = RFdict[igamma][thismom]['RF'][iset]['Boot'][it].Std
    return RFdict
                
def ReadSetFitRFDict(thisindir,thisSetList,thisGammaList,thisMethodList,thisMomList=RunMomList,thisPrintRead=PrintRead):
    if 'twopt' not in thisGammaList: raise IOError('twopt needed in set list')
    data3ptdict = ReadCfunsDict(thisindir,thisSetList,thisGammaList,thisMomList=thisMomList)
    datadict = ReadSetDict(thisindir,thisSetList,thisGammaList,thisMethodList,thisMomList=thisMomList,thisPrintRead=thisPrintRead)
    zmomstr = 'q = 0 0 0'
    start = time.time()
    if ScaleByP4g4: datadict = RewriteRF(datadict,data3ptdict,'P4g4',zmomstr)
    for igamma in datadict.keys():
        if igamma == 'twopt': continue
        if thisPrintRead: print 'Constructing Fitted RF Values: ' , igamma , '     \r',
        if zmomstr not in datadict[igamma].keys(): continue
        for iset in datadict[igamma][zmomstr]['RF'].keys():
            if 'RF' not in datadict['twopt'][zmomstr].keys(): continue
            if RemoveTSink(iset) not in datadict['twopt'][zmomstr]['RF'].keys(): continue
                # if thisPrintRead: print RemoveTSink(iset)+' not in two point set list, not constructing RF'
            data3pt = data3ptdict[igamma][zmomstr]['RF'][iset]['Boot']
            for iSF in ['OSF'+iOSF for iOSF in OSFFileFlags]+['TSF'+iTSF for iTSF in TSFFileFlags]:
                if iSF in datadict['twopt'][zmomstr].keys():
                    if RemoveTSink(iset) in datadict['twopt'][zmomstr][iSF].keys():
                        pars2pt = []
                        if 'OSF' in iSF:
                            if not all([iState in datadict['twopt'][zmomstr][iSF][RemoveTSink(iset)].keys() for iState in StateParList['One']['C2']]): continue
                            for ipar in StateParList['One']['C2']:
                                fitrkey = RemoveTSink(iset)
                                if DefTvarPicked in fitrkey: fitrkey = PickedStateStr+DefTvarPicked
                                pars2pt.append(datadict['twopt'][zmomstr][iSF][RemoveTSink(iset)][ipar][OSFfitr[RemoveToDt(fitrkey)]]['Boot'])
                            data2ptZ =  pars2pt[0]*(pars2pt[1]*(-GetintTSink(iset)+tsource)).exp(1)
                        elif 'TSF' in iSF:
                            if not all([iState in datadict['twopt'][zmomstr][iSF][RemoveTSink(iset)].keys() for iState in StateParList['Two']['C2']]): continue
                            for ipar in StateParList['Two']['C2']:
                                pars2pt.append(datadict['twopt'][zmomstr][iSF][RemoveTSink(iset)][ipar][TSFfitr]['Boot'])
                            data2ptZ =  ff.C2TSFLineFun(GetintTSink(iset)-tsource,pars2pt)
                        if 'RF'+iSF not in datadict[igamma][zmomstr].keys(): datadict[igamma][zmomstr]['RF'+iSF] = OrderedDict()
                        if iset not in datadict[igamma][zmomstr]['RF'+iSF].keys(): datadict[igamma][zmomstr]['RF'+iSF][iset] = OrderedDict()
                        datadict[igamma][zmomstr]['RF'+iSF][iset]['Boot'] = [idata3pt/data2ptZ for idata3pt in data3pt[tsource-1:GetintTSink(iset)]]
                        datadict[igamma][zmomstr]['RF'+iSF][iset]['tVals'] = datadict[igamma][zmomstr]['RF'][iset]['tVals']
                        GetBootStats(datadict[igamma][zmomstr]['RF'+iSF][iset]['Boot'])
                        datadict[igamma][zmomstr]['RF'+iSF][iset]['Vals'] = Pullflag(datadict[igamma][zmomstr]['RF'+iSF][iset]['Boot'],'Avg')
                        datadict[igamma][zmomstr]['RF'+iSF][iset]['Valserr'] = Pullflag(datadict[igamma][zmomstr]['RF'+iSF][iset]['Boot'],'Std')
    if thisPrintRead: print 'Constructing Fitted RF Values took: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s '
    return datadict

## DataDict { FormFactor } { Set } { Mass:Set/Avg/Std/Chi/Boot , FF#:qsqrd:Avg/Std/Boot , Chi:qsqrd}
## N.B. Set under Mass is different to Set, it is the set sellected for the mass extraction.
## thisFFDict = { keys=(scalar,vector etc):values=(setlist for current) / Mass:Set/Avg/Std/Chi/Boot / FF#:qsqrd:Avg/Std/Boot / Chi:qsqrd}

def ReadFFDict(thisindir,thisFFDict,thisPrintRead=PrintRead):
    DataDict = OrderedDict()
    if thisPrintRead: print ''
    start = time.time()
    for iFF,(thisFF,FFSetList) in enumerate(thisFFDict.iteritems()):
        startff = time.time()
        DataDict[thisFF] = OrderedDict()
        for iset,thisset in enumerate(FFSetList):
            if thisPrintRead: print 'Reading ',thisFF,'at : ' ,int((iset*100)/float(len(FFSetList))),'%     \r',
            thisdir = thisindir+'FormFactors/'+thisFF+'/'
            thisfile = thisdir+thisFF+thisset+'.xml'
            DataDict[thisFF][thisset] = ReadFFFile(thisfile)            
        if thisPrintRead: print 'Reading ',thisFF,'took : ' , str(datetime.timedelta(seconds=time.time()-startff)) , ' h:m:s                     '
    if thisPrintRead: print 'Reading all FFs took : ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s                     '
    return DataDict

## DataDict { gamma } { mom } { method } { set }
def ReadSetDict(thisindir,thisSetList,thisGammaList,thisMethodList,thisMomList=RunMomList,thisPrintRead=PrintRead):
    DataDict = OrderedDict()
    if thisPrintRead: print ''
    start = time.time()
    massMethodList = copy(thisMethodList)
    if 'SumMeth' in massMethodList:massMethodList.remove('SumMeth')
    if 'Fits' in massMethodList: massMethodList.remove('Fits')
    for ig,thisgamma in enumerate(thisGammaList):
        if thisPrintRead: print 'Reading at: ' ,int((ig*100)/float(len(thisGammaList))),'%     \r',
        if thisgamma == 'twopt':
            gammadirin = thisindir+'cfuns/twopt/'
            DataDict[thisgamma] = MakeMethodsDict(gammadirin,'twopt.xml',
                                                  massMethodList,thisSetList,thisMomList=thisMomList,
                                                  thisPrintRead=thisPrintRead)
        else:
            gammadirin = thisindir+CreateOppDir(thisgamma)+'/'
            DataDict[thisgamma] = MakeMethodsDict(gammadirin,thisgamma+'.xml',
                                                  thisMethodList,thisSetList,thisMomList=thisMomList,
                                                  thisPrintRead=thisPrintRead)
    if 'RF' in thisMethodList: DataDict = CombSetBoot(DataDict,'-',thisPrintRead=thisPrintRead)
    if thisPrintRead: print 'Reading took: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s                     '
    return DataDict

def ReadCfunsDict(thisindir,thisSetList,thisGammaList,thisMomList=RunMomList,thisPrintRead=PrintRead):
    DataDict = OrderedDict()
    if thisPrintRead: print ''
    start = time.time()
    for ig,thisgamma in enumerate(thisGammaList):
        if thisPrintRead: print 'Reading Cfuns at: ' ,int((ig*100)/float(len(thisGammaList))),'%     \r',
        if thisgamma == 'twopt':
            gammadirin = thisindir+'cfuns/twopt/'
            DataDict[thisgamma] = MakeMethodsDict(gammadirin,'twopt.xml',
                                                        ['RF'],thisSetList,thisMomList=thisMomList,
                                                  thisPrintRead=thisPrintRead)
        else:
            gammadirin = thisindir+'cfuns/'+CreateOppDir(thisgamma)+'/'
            DataDict[thisgamma] = MakeMethodsDict(gammadirin,thisgamma+'.xml',
                                                  ['RF'],thisSetList,thisMomList=thisMomList,
                                                  thisPrintRead=thisPrintRead)
    if thisPrintRead: print 'Reading Cfuns took: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s                     '
    DataDict = CombSetBoot(DataDict,'-',thisPrintRead=thisPrintRead)
    return DataDict


def SetRFDictToList(DictData,thisPrintRead=PrintRead):
    dataoutRF,dataout2pt = [],[]
    infolistRF,infolist2pt = [],[]
    gammalistout = OrderedDict()
    BorA = ''
    start = time.time()
    for igamma,(thisgamma,gammadata) in enumerate(DictData.iteritems()):
        if thisPrintRead: print 'Converting DictToList ' + thisgamma+ '          \r',
        gammalistout[thisgamma] = []
        if thisgamma == 'twopt':
            for imom,(thismom,momdata) in enumerate(gammadata.iteritems()):
                gammalistout[thisgamma].append(thismom)
                dataout2pt.append([])
                infolist2pt.append([])
                for iset,setdata in enumerate(momdata['RF'].itervalues()):
                    print setdata.keys()
                    if 'Info' in setdata.keys(): infolist2pt[imom].append(setdata['Info'])
                    if 'Boot' in setdata.keys():
                        dataout2pt[imom].append(setdata['Boot'])
                        if BorA == 'Avg': BorA = 'Mixed'
                        else: BorA = 'Boot'
                    else:
                        dataout2pt[imom].append([setdata['Avg'],setdata['Std']])
                        if BorA == 'Boot': BorA = 'Mixed'
                        else:BorA = 'Avg'
        else:
            dataoutRF.append([])
            infolistRF.append([])
            for imom,(thismom,momdata) in enumerate(gammadata.iteritems()):
                gammalistout[thisgamma].append(thismom)
                dataoutRF[-1].append([])
                infolistRF[-1].append([])
                for iset,(thisset,setdata) in enumerate(momdata['RF'].iteritems()):
                    if 'Info' in setdata.keys(): infolistRF[-1].append(setdata['Info'])
                    if 'Boot' in setdata.keys():
                        dataoutRF[-1][imom].append(setdata['Boot'])
                        if BorA == 'Avg': BorA = 'Mixed'
                        else: BorA = 'Boot'
                    else:
                        dataoutRF[-1][imom].append([setdata['Avg'],setdata['Std']])
                        if BorA == 'Boot': BorA = 'Mixed'
                        else:BorA = 'Avg'
    if 'twopt' in gammalistout.keys():
        if 'q = 0 0 0' in gammalistout['twopt']:
            dataout2pt.insert(0, dataout2pt.pop(gammalistout['twopt'].index('q = 0 0 0')))
            infolist2pt.insert(0, infolist2pt.pop(gammalistout['twopt'].index('q = 0 0 0')))
            gammalistout['twopt'].insert(0, gammalistout['twopt'].pop(gammalistout['twopt'].index('q = 0 0 0')))
    if thisPrintRead: print 'Converting DictToList took: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s                     '
    return dataoutRF,dataout2pt,gammalistout,BorA,infolistRF,infolist2pt
##datadict.keys() is gammalist



def CombBoot(udata,ddata,opp):
    dataout = {'Boot':ops[opp](np.array(udata['Boot']),np.array(ddata['Boot']))}
    GetBootStats(dataout['Boot'])
    dataout['tVals'] = udata['tVals']
    dataout['Vals'] = Pullflag(dataout['Boot'],'Avg')
    dataout['Valserr'] = Pullflag(dataout['Boot'],'Std')
    return dataout

def CombMethBoot(udata,ddata,opp,imom,igamma):
    dataout = OrderedDict()
    for iSet in udata.keys():
        if 'Boot' in udata[iSet].keys() and 'Boot' in ddata[iSet].keys():
            dataout[iSet] = CombBoot(udata[iSet],ddata[iSet],opp)
        else:
            raise IOError(iSet +' does not contain Boot for ' + igamma + ' ' + imom)
    return dataout


def CombSetBoot(data,opp,thisPrintRead=PrintRead):
    dataout = deepcopy(data)
    start = time.time()
    for dgamma,dgammadata in data.iteritems():
        if 'doub' in dgamma:
            combgamma = dgamma.replace('doub','')
            for sgamma,sgammadata in data.iteritems():
                if 'sing'+combgamma == sgamma:
                    if thisPrintRead: print 'Combining ' , dgamma , ' ' , sgamma ,'             \r',
                    if combgamma not in dataout.keys():
                        dataout[combgamma] = OrderedDict()
                    for imom in dgammadata.iterkeys():
                        if imom in sgammadata.keys():
                            if 'RF' in sgammadata[imom].keys():
                                if imom not in dataout[combgamma].keys():
                                    dataout[combgamma][imom] = OrderedDict()
                                dataout[combgamma][imom]['RF'] = CombMethBoot(data[sgamma][imom]['RF'],
                                                                              data[dgamma][imom]['RF'],opp,imom,combgamma)
                            else:
                                print ''
                                print 'Warning: ', imom , ' in ' , dgamma , ' but not in ' , sgamma 
                        else:
                            print ''
                            print 'Warning: ', imom , ' in ' , dgamma , ' but not in ' , sgamma 
    if thisPrintRead: print 'Combining took: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s                     '
    return dataout


def MakeMethodsDict(readdir,readfile,thisMethodList,thisSetList,thisMomList=RunMomList,thisPrintRead=PrintRead):
    MethDict = OrderedDict()
    loopSetList = thisSetList
    if 'twopt' in readdir: loopSetList = ReduceTooMassSet(thisSetList)
    for iMeth in thisMethodList:
        thisDict = OrderedDict()
        if iMeth == 'RF':
            for iSet in loopSetList:
                filename = readdir+iSet+readfile
                thisDict[iSet] = ReadRFFile(readdir,iSet+readfile,thisMomList=thisMomList)
        elif 'TSF' in iMeth:
            for iSet in ReduceTooMassSet(loopSetList):
                thisDict[iSet] = ReadTSFFile(readdir+iMeth+'/',iSet+readfile.replace('.xml','##.xml'),thisMomList=thisMomList)
        elif 'OSF' in iMeth:
            for iSet in loopSetList:
                filename = iSet+readfile.replace('.xml','##.xml')
                thisDict[iSet] = ReadOSFFile(readdir+iMeth+'/',filename,thisMomList=thisMomList)
        elif 'SumMeth' in iMeth:
            for iSet in ReduceTsink(loopSetList,NoCM=True):
                filename = iSet+readfile
                thisDict[iSet] = ReadSumFile(readdir+iMeth+'/',filename,thisMomList=thisMomList)
        elif 'Fits' in iMeth:
            for iSet in loopSetList:
                filename = iSet+readfile
                thisDict[iSet] = ReadFitFile(readdir+iMeth+'/',filename,thisMomList=thisMomList)
        else:
            if thisPrintRead: print iMeth , ' not a known method'
        for iSet in thisDict.keys():
            for imom in thisDict[iSet].keys():
                if imom not in MethDict.keys(): MethDict[imom] = OrderedDict()
                if iMeth not in MethDict[imom].keys(): MethDict[imom][iMeth] = OrderedDict()
                MethDict[imom][iMeth][iSet] = thisDict[iSet][imom]
    return MethDict


def ReadTSFFile(filedir,filename,thisMomList=RunMomList):
    return ReadSFFile(filedir,filename,OneOrTwo='Two',thisMomList=thisMomList)

def ReadOSFFile(filedir,filename,thisMomList=RunMomList):
    return ReadSFFile(filedir,filename,OneOrTwo='One',thisMomList=thisMomList)



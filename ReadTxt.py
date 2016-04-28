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
def ExtractValues(thisindir,thisGammaList,thisSetList,thisMethodList,thisMomList=qvecSet):
    def SetupDict(dictin,thisgamma,thiscol):
        if thiscol not in dictin.keys():
            dictin[thiscol] = {}
        if thisgamma not in dictin[thiscol].keys():
            dictin[thiscol][thisgamma] = {}
        return dictin
    readdata = ReadSetDict(thisindir,thisSetList,thisGammaList,thisMethodList,thisMomList=thisMomList)
    start = time.time()
    datadictout,datamassout = OrderedDict(),OrderedDict()
    for ig,igamma in enumerate(thisGammaList):
        if PrintRead: print ' Extracting data from read: ' , int((ig*100)/float(len(thisGammaList))) ,'%  \r',
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
                        if CheckDict(thismassdict,'m0',OSFfitr[ism]):
                            datamassout[ism+iMeth] = thismassdict['m0'][OSFfitr[ism]]
        else:
            for imom,momdata in readdata[igamma].iteritems():
                for iMeth,Methdata in momdata.iteritems():
                    for iSet,thisdict in Methdata.iteritems():
                        itsink,ism = SplitTSinkString(iSet)
                        if DefTvarPicked in ism: ism = PickedStateStr+DefTvarPicked
                        if 'TSF' in iMeth:
                            for icut in TSFCutList:
                                if CheckDict(thisdict,'B00',TSFfitr,icut):
                                    datadictout = SetupDict(datadictout,igamma,iSet+iMeth+icut)
                                    datadictout[iSet+iMeth+icut][igamma][imom] = thisdict['B00'][TSFfitr][icut]
                        elif 'OSF' in iMeth:
                            for icut in OSFCutList:
                                if CheckDict(thisdict,'B00',OSFfitr[ism],icut):
                                    datadictout = SetupDict(datadictout,igamma,iSet+iMeth+icut)
                                    datadictout[iSet+iMeth+icut][igamma][imom] = thisdict['B00'][OSFfitr[ism]][icut]
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
    if PrintRead: print 'Extracting data took: ', str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s                  '
    return datadictout,datamassout

def GetCompletedMom(thisfile):
    thisqlist = set([])
    with open(thisfile,'r') as f:
        for line in f:
            if len(line.strip()) > 0:
                if line.strip()[0] == 'q':
                    thisqlist = thisqlist | set([line.strip()])
    return thisqlist

def GetTxtAndBootComp(thisfile,thisbootfile):
    txtqlist = GetCompletedMom(thisfile)
    bootqlist = GetCompletedMom(thisbootfile)
    if len(txtqlist ^ bootqlist) != 0:
            raise IOError('momenta missmatch '+thisbootfile+' '+thisfile)
    return txtqlist & bootqlist

def Get2ptSetMoms(outputdir,MomListIn,statelist=[],todtlist=[],smlist=[]):
    momlist = set(MomListIn)
    for iflag in ['cfuns/twopt','Mass']:
        thisdir = outputdir+iflag+'/'
        for itodt in todtlist:
            if iflag == 'Mass':
                ifile = thisdir+itodt+'LREM.txt'
                if os.path.isfile(ifile):
                    tempmomlist = GetCompletedMom(ifile)
                    momlist = momlist & tempmomlist
                else:
                    return MomListIn
            for istate in statelist:
                ifile = thisdir+'state'+istate+itodt+iflag.replace('cfuns/','')+'.txt'
                # ifb = thisdir+'boots/state'+istate+itodt+iflag.replace('cfuns/','')+'.boot.txt'
                # if os.path.isfile(ifile) and os.path.isfile(ifb):
                #     tempmomlist = GetTxtAndBootComp(ifile,ifb)
                if os.path.isfile(ifile):
                    tempmomlist = GetCompletedMom(ifile)
                    momlist = momlist & tempmomlist
                else:
                    return MomListIn
        for ism in smlist:
            ifile = thisdir+'sm'+ism+iflag.replace('cfuns/','')+'.txt'
            # ifb = thisdir+'boots/sm'+ism+iflag.replace('cfuns/','')+'.boot.txt'
            # if os.path.isfile(ifile) and os.path.isfile(ifb):
            #     tempmomlist = GetTxtAndBootComp(ifile,ifb)
            if os.path.isfile(ifile):
                tempmomlist = GetCompletedMom(ifile)
                momlist = momlist & tempmomlist
            else:
                return MomListIn
    return OrderMomList(set(MomListIn) - momlist)


def Get3SM(outputdir,thisGammaList,MomListIn,setlist):
    momlist = set(MomListIn)
    for igamma in thisGammaList:
        thisdir = outputdir+CreateOppDir(igamma)
        for iset in setlist:
            ifile = thisdir+iset+igamma+'.txt'
            if PrintRead: print ifile
            if os.path.isfile(ifile):
                tempmomlist = GetCompletedMom(ifile)
                momlist = momlist & tempmomlist
            else:
                return MomListIn
    return OrderMomList(set(MomListIn) - momlist)
            

def Get3ptSetMoms(outputdir,thisGammaList,MomListIn,setlist):
    return OrderMomList(set(Get3SM(outputdir,thisGammaList,MomListIn,setlist)) |
                        set(Get3SM(outputdir+'cfuns/',thisGammaList,MomListIn,setlist)))


def Check2pt(outputdir,momstr,statelist=[],todtlist=[],smlist=[]):
    DoMom = False
    for iflag in ['cfuns/twopt','Mass']:
        thisdir = outputdir+iflag+'/'
        for itodt in todtlist:
            if iflag == 'Mass':
                ifile = thisdir+itodt+'LREM.txt'
                if os.path.isfile(ifile):
                    tempmomlist = GetTxtAndBootComp(ifile,ifb)
                    momlist = momlist & tempmomlist
                else:
                    return []
            for istate in statelist:
                ifile = thisdir+'state'+istate+itodt+iflag.replace('cfuns/','')+'.txt'
                ifb = thisdir+'boots/state'+istate+itodt+iflag.replace('cfuns/','')+'.boot.txt'
                if os.path.isfile(ifile) and os.path.isfile(ifb):
                    if momstr not in GetTxtAndBootComp(ifile,ifb):
                        DoMom = True
                else:
                    DoMom = True
        for ism in smlist:
            ifile = thisdir+'sm'+ism+iflag.replace('cfuns/','')+'.txt'
            ifb = thisdir+'boots/sm'+ism+iflag.replace('cfuns/','')+'.boot.txt'
            if os.path.isfile(ifile) and os.path.isfile(ifb):
                if momstr not in GetTxtAndBootComp(ifile,ifb):
                    DoMom = True
            else:
                DoMom = True
    return DoMom


def CheckSet(outputdir,momstr,thisGammaList,tlist=[],statelist=[],todtlist=[],smlist=[]):
    DoMom = False
    for igamma in thisGammaList:
        thisdir = outputdir+CreateOppDir(igamma)
        if len(tlist) > 0:
            for it in tlist:
                for istate in statelist:
                    for itodt in todtlist:
                        ifile = thisdir+'tsink'+it+'state'+istate+itodt+igamma+'.txt'
                        # if PrintRead: print ifile
                        ifb = thisdir+'boots/tsink'+it+'state'+istate+itodt+igamma+'.boot.txt'
                        if os.path.isfile(ifile) and os.path.isfile(ifb):
                            if momstr not in GetTxtAndBootComp(ifile,ifb): DoMom = True
                        else:
                            DoMom = True
                for ism in smlist:
                    ifile = thisdir+'tsink'+it+'sm'+ism+igamma+'.txt'
                    # if PrintRead: print ifile
                    ifb = thisdir+'boots/tsink'+it+'sm'+ism+igamma+'.boot.txt'
                    if os.path.isfile(ifile) and os.path.isfile(ifb):
                        if momstr not in GetTxtAndBootComp(ifile,ifb): DoMom = True
                    else:
                        DoMom = True
        else:
            for istate in statelist:
                for itodt in todtlist:
                    ifile = thisdir+'state'+istate+itodt+igamma+'.txt'
                    # if PrintRead: print ifile
                    ifb = thisdir+'boots/state'+istate+itodt+igamma+'.boot.txt'
                    if os.path.isfile(ifile) and os.path.isfile(ifb):
                        if momstr not in GetTxtAndBootComp(ifile,ifb): DoMom = True
                    else:
                        DoMom = True
            for ism in smlist:
                ifile = thisdir+'sm'+ism+igamma+'.txt'
                # if PrintRead: print ifile
                ifb = thisdir+'boots/sm'+ism+igamma+'.boot.txt'
                if os.path.isfile(ifile) and os.path.isfile(ifb):
                    if momstr not in GetTxtAndBootComp(ifile,ifb): DoMom = True
                else:
                    DoMom = True
    return DoMom

##NEW FUNCTIONS##

#data2pt = [ mom , set , it ] bs
#data3pt = [ gamma , mom , set , it ] bs
def ReadCfunsnp(thisGammaList,thisSetList,thisMomList=[]):
    datadict = ReadCfunsDict(outputdir,thisSetList,thisGammaList,thisMomList=thisMomList)
    return SetRFDictToList(datadict)

def ReadRFnp(thisGammaList,thisSetList,thisMomList=[]):
    datadict = ReadSetDict(outputdir,thisSetList,thisGammaList,['RF'],thisMomList=thisMomList)
    return SetRFDictToList(datadict)

def ReadSetAndCfunsDict(thisindir,thisSetList,thisGammaList,thisMethodList,thisMomList=[]):
    datadict = ReadSetDict(thisindir,thisSetList,thisGammaList,thisMethodList,thisMomList=thisMomList)
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
                
def ReadSetFitRFDict(thisindir,thisSetList,thisGammaList,thisMethodList,thisMomList=[]):
    if 'twopt' not in thisGammaList: raise IOError('twopt needed in set list')
    data3ptdict = ReadCfunsDict(thisindir,thisSetList,thisGammaList,thisMomList=thisMomList)
    datadict = ReadSetDict(thisindir,thisSetList,thisGammaList,thisMethodList,thisMomList=thisMomList)
    zmomstr = 'q = 0 0 0'
    start = time.time()
    if ScaleByP4g4: datadict = RewriteRF(datadict,data3ptdict,'P4g4',zmomstr)
    for igamma in datadict.keys():
        if igamma == 'twopt': continue
        if PrintRead: print 'Constructing Fitted RF Values: ' , igamma , '     \r',
        if zmomstr not in datadict[igamma].keys(): continue
        for iset in datadict[igamma][zmomstr]['RF'].keys():
            if 'RF' not in datadict['twopt'][zmomstr].keys(): continue
            if RemoveTSink(iset) not in datadict['twopt'][zmomstr]['RF'].keys(): continue
                # if PrintRead: print RemoveTSink(iset)+' not in two point set list, not constructing RF'
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
                                pars2pt.append(datadict['twopt'][zmomstr][iSF][RemoveTSink(iset)][ipar][OSFfitr[fitrkey]]['Boot'])
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
    if PrintRead: print 'Constructing Fitted RF Values took: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s '
    return datadict

## DataDict { FormFactor } { Set } { Mass:Set/Avg/Std/Chi/Boot , FF#:qsqrd:Avg/Std/Boot , Chi:qsqrd}
## N.B. Set under Mass is different to Set, it is the set sellected for the mass extraction.
## thisFFDict = { keys=(scalar,vector etc):values=(setlist for current) / Mass:Set/Avg/Std/Chi/Boot / FF#:qsqrd:Avg/Std/Boot / Chi:qsqrd}

def ReadFFDict(thisindir,thisFFDict):
    DataDict = OrderedDict()
    if PrintRead: print ''
    start = time.time()
    for iFF,(thisFF,FFSetList) in enumerate(thisFFDict.iteritems()):
        startff = time.time()
        DataDict[thisFF] = OrderedDict()
        for iset,thisset in enumerate(FFSetList):
            if PrintRead: print 'Reading ',thisFF,'at : ' ,int((iset*100)/float(len(FFSetList))),'%     \r',
            thisdir = thisindir+'FormFactors/'+thisFF+'/'
            thisfile = thisdir+thisFF+thisset+'.txt'
            thisbootfile = thisdir + 'boots/'+thisFF+thisset+'.boot.txt'
            DataDict[thisFF][thisset] = ReadFFFile(thisfile,bootfn=thisbootfile)            
        if PrintRead: print 'Reading ',thisFF,'took : ' , str(datetime.timedelta(seconds=time.time()-startff)) , ' h:m:s                     '
    if PrintRead: print 'Reading all FFs took : ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s                     '
    return DataDict

## DataDict { gamma } { mom } { method } { set }
def ReadSetDict(thisindir,thisSetList,thisGammaList,thisMethodList,thisMomList=[]):
    DataDict = OrderedDict()
    if PrintRead: print ''
    start = time.time()
    for ig,thisgamma in enumerate(thisGammaList):
        if PrintRead: print 'Reading at: ' ,int((ig*100)/float(len(thisGammaList))),'%     \r',
        if thisgamma == 'twopt':
            gammadirin = thisindir+'cfuns/twopt/'
            massMethodList = copy(thisMethodList)
            if 'SumMeth' in massMethodList:massMethodList.remove('SumMeth')
            if 'Fits' in massMethodList: massMethodList.remove('Fits')
            DataDict[thisgamma] = MakeMethodsDict(gammadirin,'twopt.txt',
                                                  massMethodList,thisSetList,thisMomList=thisMomList)
        else:
            gammadirin = thisindir+CreateOppDir(thisgamma)+'/'
            DataDict[thisgamma] = MakeMethodsDict(gammadirin,thisgamma+'.txt',
                                                  thisMethodList,thisSetList,thisMomList=thisMomList)
    if PrintRead: print 'Reading took: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s                     '
    if 'RF' in thisMethodList: DataDict = CombSetBoot(DataDict,'-')
    return DataDict

def ReadCfunsDict(thisindir,thisSetList,thisGammaList,thisMomList=[]):
    DataDict = OrderedDict()
    if PrintRead: print ''
    start = time.time()
    for ig,thisgamma in enumerate(thisGammaList):
        if PrintRead: print 'Reading Cfuns at: ' ,int((ig*100)/float(len(thisGammaList))),'%     \r',
        if thisgamma == 'twopt':
            gammadirin = thisindir+'cfuns/twopt/'
            DataDict[thisgamma] = MakeMethodsDict(gammadirin,'twopt.txt',
                                                        ['RF'],thisSetList,thisMomList=thisMomList)
        else:
            gammadirin = thisindir+'cfuns/'+CreateOppDir(thisgamma)+'/'
            DataDict[thisgamma] = MakeMethodsDict(gammadirin,thisgamma+'.txt',
                                                        ['RF'],thisSetList,thisMomList=thisMomList)
    if PrintRead: print 'Reading Cfuns took: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s                     '
    DataDict = CombSetBoot(DataDict,'-')
    return DataDict


def SetRFDictToList(DictData):
    dataoutRF,dataout2pt = [],[]
    gammalistout = OrderedDict()
    BorA = ''
    for igamma,(thisgamma,gammadata) in enumerate(DictData.iteritems()):
        gammalistout[thisgamma] = []
        if thisgamma == 'twopt':
            for imom,(thismom,momdata) in enumerate(gammadata.iteritems()):
                gammalistout[thisgamma].append(thismom)
                dataout2pt.append([])
                for iset,setdata in enumerate(momdata['RF'].itervalues()):
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
            for imom,(thismom,momdata) in enumerate(gammadata.iteritems()):
                gammalistout[thisgamma].append(thismom)
                dataoutRF[-1].append([])
                for iset,(thisset,setdata) in enumerate(momdata['RF'].iteritems()):
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
            gammalistout['twopt'].insert(0, gammalistout['twopt'].pop(gammalistout['twopt'].index('q = 0 0 0')))
    return dataoutRF,dataout2pt,gammalistout,BorA
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
        # elif iMeth == 'SumMeth':
        #     for idict1 in udata[iMeth].keys():
        #         dataout[iMeth][dict1] = OrderedDict()
        #         for idict2 in udata[iMeth][idict1].keys():
        #             dataout[iMeth][idict1][idict2] = CombBoot(udata[iMeth][idict1][idict2],
        #                                                       ddata[iMeth][idict1][idict2],opp)
        # elif 'SF' in iMeth:
        #     for idict1 in udata[iMeth].keys():
        #         for idict2 in udata[iMeth][idict1].keys():
        #             for idict3 in udata[iMeth][idict1][idict2].keys():
        #                 dataout[iMeth][idict1][idict2][idict3] = CombBoot(udata[iMeth][idict1][idict2][idict3],
        #                                                                   ddata[iMeth][idict1][idict2][idict3],
        #                                                                   opp)
    return dataout


def CombSetBoot(data,opp):
    dataout = deepcopy(data)
    start = time.time()
    for dgamma,dgammadata in data.iteritems():
        if 'doub' in dgamma:
            combgamma = dgamma.replace('doub','')
            for sgamma,sgammadata in data.iteritems():
                if 'sing'+combgamma == sgamma:
                    if PrintRead: print 'Combining ' , dgamma , ' ' , sgamma ,'             \r',
                    if combgamma not in dataout.keys():
                        dataout[combgamma] = OrderedDict()
                    for imom in dgammadata.iterkeys():
                        if imom in sgammadata.keys():
                            if imom not in dataout[combgamma].keys():
                                dataout[combgamma][imom] = OrderedDict()
                            dataout[combgamma][imom]['RF'] = CombMethBoot(data[sgamma][imom]['RF'],
                                                                          data[dgamma][imom]['RF'],opp,imom,combgamma)
                        else:
                            print ''
                            print 'Warning: ', imom , ' in ' , dgamma , ' but not in ' , sgamma 
    if PrintRead: print 'Combining took: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s                     '
    return dataout


def MakeMethodsDict(readdir,readfile,thisMethodList,thisSetList,thisMomList=[]):
    MethDict = OrderedDict()
    for iMeth in thisMethodList:
        if iMeth == 'RF':
            thisDict = OrderedDict()
            RFSetList = thisSetList
            if 'twopt' in readfile: RFSetList = ReduceTsink(thisSetList)
            for iSet in RFSetList:
                filename = readdir+iSet+readfile
                bootfn = readdir+'boots/'+iSet+readfile.replace('.txt','.boot.txt')
                thisDict[iSet] = ReadRFFile(filename,bootfn=bootfn,thisMomList=thisMomList)
        elif 'SF' in iMeth:
            thisDict=OrderedDict()
            if 'TSF' in iMeth:
                for iSet in ReduceTsink(thisSetList):
                    filename = readdir+iMeth+'/'+iSet+readfile.replace('.txt','##.txt')
                    bootfn = readdir+iMeth+'/boots/'+iSet+readfile.replace('.txt','##.boot.txt')
                    thisDict[iSet] = ReadTSFFile(filename,bootfn=bootfn,thisMomList=thisMomList)
            elif 'OSF' in iMeth:
                OSFSetList = thisSetList
                if 'twopt' in readfile: OSFSetList = ReduceTsink(thisSetList)
                for iSet in OSFSetList:
                    filename = readdir+iMeth+'/'+iSet+readfile.replace('.txt','##.txt')
                    bootfn = readdir+iMeth+'/boots/'+iSet+readfile.replace('.txt','##.boot.txt')
                    thisDict[iSet] = ReadOSFFile(filename,bootfn=bootfn,thisMomList=thisMomList)
        elif 'SumMeth' in iMeth:
            thisDict=OrderedDict()
            for iSet in ReduceTsink(thisSetList,NoCM=True):
                filename = readdir+iMeth+'/'+iSet+readfile
                bootfn = readdir+iMeth+'/boots/'+iSet+readfile.replace('.txt','.boot.txt')
                thisDict[iSet] = ReadSumFile(filename,bootfn=bootfn,thisMomList=thisMomList)
        elif 'Fits' in iMeth:
            thisDict = OrderedDict()
            for iSet in thisSetList:
                filename = readdir+iMeth+'/'+iSet+readfile
                bootfn = readdir+iMeth+'/boots/'+iSet+readfile.replace('.txt','.boot.txt')
                thisDict[iSet] = ReadFitFile(filename,bootfn=bootfn,thisMomList=thisMomList)
        else:
            if PrintRead: print iMeth , ' not a known method'
        for iSet in thisDict.keys():
            for imom in thisDict[iSet].keys():
                if imom not in MethDict.keys():
                    MethDict[imom] = OrderedDict()
                if iMeth not in MethDict[imom].keys():
                    MethDict[imom][iMeth] = OrderedDict()
                MethDict[imom][iMeth][iSet] = thisDict[iSet][imom]
    return MethDict

## Replaced with ReadXml Version
# ##Also works for cfuns##
# ##outputdict = { thismom , [tVals] / [Vals] / [Valserr] / [Boot] bs }
# def ReadRFFile(filename,bootfn='',thisMomList=[]):
#     ReadMom = False
#     renorm = GetRenorm(filename)
#     outputdict = OrderedDict()
#     # print filename
#     if os.path.isfile(filename):
#         def thisReadFile(thisoutputdict):
#             with open(filename,'r') as f:
#                 for line in f:
#                     rdata = line.strip()
#                     if rdata[0] == 'q':
#                         if len(thisMomList) > 0 and all([imom in thisoutputdict.keys() for imom in thisMomList]): return thisoutputdict
#                         if rdata in thisMomList or len(thisMomList) ==0:
#                             thismom = rdata
#                             ReadMom = True
#                             try:
#                                 thisoutputdict[thismom] = OrderedDict()
#                             except:
#                                 raise IOError("double momenta in file"+filename)
#                             thisoutputdict[thismom]['tVals'] = []
#                             thisoutputdict[thismom]['Vals'] = []
#                             thisoutputdict[thismom]['Valserr'] = []
#                         else:
#                             ReadMom = False
#                     elif ReadMom:
#                         rdata = rdata.split()
#                         thisoutputdict[thismom]['tVals'].append(int(rdata[0]))
#                         thisoutputdict[thismom]['Vals'].append(float(rdata[1])*renorm)
#                         thisoutputdict[thismom]['Valserr'].append(float(rdata[2])*renorm)
#             return thisoutputdict
#         outputdict = thisReadFile(outputdict)

#         if os.path.isfile(bootfn):
#             def thisReadBoot(thisoutputdict):
#                 currMomList = []
#                 ReadMom = False
#                 with open(bootfn,'r') as f:
#                     for line in f:
#                         rdata = line.strip()
#                         if rdata[0] == 'q':
#                             if len(thisMomList) > 0 and all([imom in currMomList for imom in thisMomList]): 
#                                 BootNdimDict(thisoutputdict)
#                                 return thisoutputdict
#                             if rdata in thisMomList or len(thisMomList) ==0:
#                                 thismom = rdata
#                                 currMomList.append(thismom)
#                                 ReadMom = True
#                                 if rdata not in thisoutputdict.keys():
#                                     print 'WARNING: '+rdata+' found in:'
#                                     print bootfn + ' but not in:'
#                                     print filename
#                                     thisoutputdict[thismom] = OrderedDict()
#                                 thisoutputdict[thismom]['Boot'] = []
#                             else:
#                                 ReadMom = False
#                         elif rdata[0] == 't' and ReadMom:
#                             rdata = rdata.split()
#                             if int(rdata[1]) not in thisoutputdict[thismom]['tVals']:
#                                 raise IOError("T value in file but not in Boot file")
#                             thisoutputdict[thismom]['Boot'].append(BootStrap1(nboot,0.9))
#                         elif 'nboot' in rdata:
#                             if nboot != int(rdata.split()[1]): raise IOError("nboot missmatch")
#                         elif ReadMom:
#                             rdata = rdata.split()
#                             iboot = int(rdata[0])
#                             thisoutputdict[thismom]['Boot'][-1].values[iboot] = float(rdata[1])*renorm
#                 BootNdimDict(thisoutputdict)
#                 return thisoutputdict
#             outputdict = thisReadBoot(outputdict)
#     return outputdict

def ReadFitFile(filename,bootfn='',thisMomList=[]):
    ReadMom = False
    renorm = 1.
    outputdict = OrderedDict()
    # print filename
    if os.path.isfile(filename):
        def thisReadFile(thisoutputdict):
            with open(filename,'r') as f:
                for line in f:
                    rdata = line.strip()
                    if len(rdata) > 0:
                        if rdata[0] == 'q':
                            if len(thisMomList) > 0 and all([imom in thisoutputdict.keys() for imom in thisMomList]): return thisoutputdict
                            if rdata in thisMomList or len(thisMomList) == 0:
                                thismom = rdata
                                ReadMom = True
                                thisoutputdict[thismom] = OrderedDict()
                            else:
                                ReadMom = False
                        elif ReadMom:
                            rdata = rdata.split()
                            if rdata[0] not in thisoutputdict[thismom].keys():
                                thisoutputdict[thismom][rdata[0]] = OrderedDict()
                            thisoutputdict[thismom][rdata[0]]['Avg'] = float(rdata[1])*renorm
                            thisoutputdict[thismom][rdata[0]]['Std'] = float(rdata[2])*renorm
                            thisoutputdict[thismom][rdata[0]]['Chi'] = float(rdata[3])
            return thisoutputdict
        outputdict = thisReadFile(outputdict)

        if os.path.isfile(bootfn):
            def thisReadBoot(thisoutputdict):
                ReadMom = False
                currMomList = []
                with open(bootfn,'r') as f:
                    for line in f:
                        rdata = line.strip()
                        if len(rdata) > 0:
                            if rdata.split()[0] == 'nboot':
                                if nboot != int(rdata.split()[1]): raise IOError("nboot missmatch")
                            elif rdata[0] == 'q':
                                if len(thisMomList) > 0 and all([imom in currMomList for imom in thisMomList]): 
                                    BootNdimDict(thisoutputdict)
                                    return thisoutputdict
                                if rdata in thisMomList or len(thisMomList) == 0:
                                    thismom = rdata
                                    currMomList.append(thismom)
                                    ReadMom = True
                                    if rdata not in thisoutputdict.keys():
                                        print 'WARNING: '+rdata+' found in:'
                                        print bootfn + ' but not in:'
                                        print filename
                                        thisoutputdict[thismom] = OrderedDict()
                                else:
                                    ReadMom = False
                            elif rdata[0] == 'c' and ReadMom:
                                rdata = rdata.split()
                                cutpar = rdata[0]
                                if cutpar not in thisoutputdict[thismom].keys():
                                    thisoutputdict[thismom][cutpar] = OrderedDict()
                                thisoutputdict[thismom][cutpar]['Boot'] = BootStrap1(nboot,0.9)
                            elif ReadMom:
                                try:
                                    thisoutputdict[thismom][cutpar]['Boot'].values[int(rdata.split()[0])] = float(rdata.split()[1])
                                except:
                                    thisoutputdict[thismom][cutpar]['Boot'].values[int(rdata[0][:3])] = float('NaN')
                BootNdimDict(thisoutputdict)
                return thisoutputdict
            outputdict = thisReadBoot(outputdict)
    return outputdict
    

##outputdict = { thismom , cutpar , tsinkrpar/tsinkval , Avg / Std / Chi / Boot (bs) }
def ReadSumFile(filename,bootfn='',thisMomList=[]):
    ReadMom = False
    renorm = 1.
    outputdict = OrderedDict()
    if os.path.isfile(filename):
        def thisReadFile(thisoutputdict):
            with open(filename,'r') as f:
                for line in f:
                    rdata = line.strip()
                    if len(rdata) > 0:
                        if rdata[0] == 'q':
                            if len(thisMomList) > 0 and all([imom in thisoutputdict.keys() for imom in thisMomList]): return thisoutputdict
                            if rdata in thisMomList or len(thisMomList) == 0:
                                thismom = rdata
                                ReadMom = True
                                try:
                                    thisoutputdict[thismom] = OrderedDict()
                                except:
                                    raise IOError("double momenta in file"+filename)
                            else:
                                ReadMom = False
                        elif rdata[0] == 'c' and ReadMom:
                            rdata = map(str,rdata.split())
                            cutpar = rdata[0]
                            thistsink = rdata[1][:-1]
                            if cutpar not in thisoutputdict[thismom].keys():
                                thisoutputdict[thismom][cutpar] = OrderedDict()
                            if thistsink not in thisoutputdict[thismom][cutpar].keys():
                                thisoutputdict[thismom][cutpar][thistsink] = OrderedDict()
                            try:
                                thisoutputdict[thismom][cutpar][thistsink]['Avg'] = float(rdata[2])*renorm
                                thisoutputdict[thismom][cutpar][thistsink]['Std'] = float(rdata[3])*renorm
                            except:
                                thisoutputdict[thismom][cutpar][thistsink]['Avg'] = float('NaN')
                                thisoutputdict[thismom][cutpar][thistsink]['Std'] = float('NaN')                            
                        elif rdata[0] == 'f' and ReadMom:
                            rdata = map(str,rdata.split())
                            cutpar = rdata[1]
                            tsinkrpar = 'fit '+rdata[2]+' '+rdata[3]+'-'+rdata[4][:-1]
                            if cutpar not in thisoutputdict[thismom].keys():
                                thisoutputdict[thismom][cutpar] = OrderedDict()
                            if tsinkrpar not in thisoutputdict[thismom][cutpar].keys():
                                thisoutputdict[thismom][cutpar][tsinkrpar] = OrderedDict()
                            try:
                                thisoutputdict[thismom][cutpar][tsinkrpar]['Avg'] = float(rdata[5])*renorm
                                thisoutputdict[thismom][cutpar][tsinkrpar]['Std'] = float(rdata[6])*renorm
                                thisoutputdict[thismom][cutpar][tsinkrpar]['Chi'] = float(rdata[7])
                            except:
                                thisoutputdict[thismom][cutpar][tsinkrpar]['Avg'] = float('NaN')
                                thisoutputdict[thismom][cutpar][tsinkrpar]['Std'] = float('NaN')
                                thisoutputdict[thismom][cutpar][tsinkrpar]['Chi'] = float('NaN')
            return thisoutputdict
        outputdict = thisReadFile(outputdict)

        if os.path.isfile(bootfn):
            def thisReadBoot(thisoutputdict):
                currMomList = []
                ReadMom = False
                with open(bootfn,'r') as f:
                    for line in f:
                        rdata = line.strip()
                        if len(rdata) > 0:
                            if rdata.split()[0] == 'nboot':
                                if nboot != int(rdata.split()[1]): raise IOError("nboot missmatch")
                            elif rdata[0] == 'q':
                                if len(thisMomList) > 0 and all([imom in currMomList for imom in thisMomList ]): 
                                    BootNdimDict(thisoutputdict)
                                    return thisoutputdict
                                if rdata in thisMomList or len(thisMomList) == 0:
                                    thismom = rdata
                                    currMomList.append(thismom)
                                    ReadMom = True
                                    if rdata not in thisoutputdict.keys():
                                        print 'WARNING: '+rdata+' found in:'
                                        print bootfn + ' but not in:'
                                        print filename
                                        thisoutputdict[thismom] = OrderedDict()
                                else:
                                    ReadMom = False
                            elif rdata[0] == 'c' and ReadMom:
                                rdata = map(str,rdata.split())
                                cutpar = rdata[0]
                                tsinkrpar = rdata[1][:-1]
                                if cutpar not in thisoutputdict[thismom].keys():
                                    thisoutputdict[thismom][cutpar] = OrderedDict()
                                if tsinkrpar not in thisoutputdict[thismom][cutpar].keys():
                                    thisoutputdict[thismom][cutpar][tsinkrpar] = OrderedDict()
                                thisoutputdict[thismom][cutpar][tsinkrpar]['Boot'] = BootStrap1(nboot,0.9)
                            elif rdata[0] == 'f' and ReadMom:
                                rdata = map(str,rdata.split())
                                cutpar = rdata[1]
                                tsinkrpar = 'fit '+rdata[2]+' '+rdata[3]+'-'+rdata[4][:-1]
                                if cutpar not in thisoutputdict[thismom].keys():
                                    thisoutputdict[thismom][cutpar] = OrderedDict()
                                if tsinkrpar not in thisoutputdict[thismom][cutpar].keys():
                                    thisoutputdict[thismom][cutpar][tsinkrpar] = OrderedDict()
                                thisoutputdict[thismom][cutpar][tsinkrpar]['Boot'] = BootStrap1(nboot,0.9)
                            elif ReadMom:
                                try:
                                    thisoutputdict[thismom][cutpar][tsinkrpar]['Boot'].values[int(rdata.split()[0])] = float(rdata.split()[1])*renorm
                                except:
                                    thisoutputdict[thismom][cutpar][tsinkrpar]['Boot'].values[int(rdata[0][:3])] = float('NaN')
                BootNdimDict(thisoutputdict)
                return thisoutputdict
            outputdict = thisReadBoot(outputdict)
    return outputdict

def ReadTSFFile(filename,bootfn='',thisMomList=[]):
    return ReadSFFile(filename,bootfn=bootfn,OneOrTwo='Two',thisMomList=thisMomList)

def ReadOSFFile(filename,bootfn='',thisMomList=[]):
    return ReadSFFile(filename,bootfn=bootfn,OneOrTwo='One',thisMomList=thisMomList)


##outputdict = { thismom , fitpar , 2corfitr , 3corcutr , Avg / Std / Chi / Boot (bs) }
## put ## as parameter
def ReadSFFile(filename,bootfn='',OneOrTwo='Two',thisMomList=[]):
    ReadMom = False
    renorm = 1.
    outputdict = OrderedDict()
    if 'twopt' in filename:
        twoptread = True
        thisCorr = 'C2'
    else:
        twoptread = False
        thisCorr = 'C3'
    for ipar in StateParList[OneOrTwo][thisCorr]:
        if ipar in ['m0','Dm'] and OneOrTwo == 'Two':
            filename = re.sub('sm.*twopt','twopt',filename)
            bootfn = re.sub('sm.*twopt','twopt',bootfn)
            filename = re.sub('state.*twopt','twopt',filename)
            bootfn = re.sub('state.*twopt','twopt',bootfn)
        # print filename.replace('##',ipar)
        if os.path.isfile(filename.replace('##',ipar)):
            def thisReadFile(thisoutputdict):
                currMomList = []
                with open(filename.replace('##',ipar),'r') as f:
                    for line in f:
                        rdata = line.strip()
                        if len(rdata) > 0:
                            if rdata[0] == 'q':
                                if len(thisMomList) > 0 and all([imom in currMomList for imom in thisMomList]): return thisoutputdict
                                if rdata in thisMomList or len(thisMomList) == 0:
                                    thismom = rdata
                                    currMomList.append(thismom)
                                    ReadMom = True
                                    if thismom not in thisoutputdict.keys():
                                        thisoutputdict[thismom] = OrderedDict()
                                    if ipar not in thisoutputdict[thismom].keys():
                                        thisoutputdict[thismom][ipar] = OrderedDict()
                                else:
                                    ReadMom = False
                            elif ReadMom and twoptread:
                                rdata = map(str,rdata.split())
                                dictfitr = rdata[0]+'-'+rdata[1]
                                if dictfitr not in thisoutputdict[thismom][ipar].keys():
                                    thisoutputdict[thismom][ipar][dictfitr] = OrderedDict()
                                try:
                                    thisoutputdict[thismom][ipar][dictfitr]['Avg'] = float(rdata[2])*renorm
                                    thisoutputdict[thismom][ipar][dictfitr]['Std'] = float(rdata[3])*renorm
                                    thisoutputdict[thismom][ipar][dictfitr]['Chi'] = float(rdata[4])
                                except:
                                    thisoutputdict[thismom][ipar][dictfitr]['Avg'] = float('NaN')
                                    thisoutputdict[thismom][ipar][dictfitr]['Std'] = float('NaN')
                                    thisoutputdict[thismom][ipar][dictfitr]['Chi'] = float('NaN')                                
                            elif ReadMom:
                                rdata = map(str,rdata.split())
                                dictfitr = rdata[0]+'-'+rdata[1]
                                cutr = 'cut'+rdata[2]
                                if dictfitr not in thisoutputdict[thismom][ipar].keys():
                                    thisoutputdict[thismom][ipar][dictfitr] = OrderedDict()
                                if cutr not in thisoutputdict[thismom][ipar][dictfitr].keys():
                                    thisoutputdict[thismom][ipar][dictfitr][cutr] = OrderedDict()
                                try:
                                    thisoutputdict[thismom][ipar][dictfitr][cutr]['Avg'] = float(rdata[3])*renorm
                                    thisoutputdict[thismom][ipar][dictfitr][cutr]['Std'] = float(rdata[4])*renorm
                                    thisoutputdict[thismom][ipar][dictfitr][cutr]['Chi'] = float(rdata[5])
                                except:
                                    thisoutputdict[thismom][ipar][dictfitr][cutr]['Avg'] = float('NaN')
                                    thisoutputdict[thismom][ipar][dictfitr][cutr]['Std'] = float('NaN')
                                    thisoutputdict[thismom][ipar][dictfitr][cutr]['Chi'] = float('NaN')
                return thisoutputdict
            outputdict = thisReadFile(outputdict)

            thisbootfn = bootfn.replace('##',ipar)
            # if twoptread: print thisbootfn
            if os.path.isfile(thisbootfn):
                def thisReadBoot(thisoutputdict):
                    currMomList = []
                    ReadMom = False
                    with open(thisbootfn,'r') as f:
                        for line in f:
                            rdata = line.strip()
                            if len(rdata) > 0:
                                if rdata.split()[0] == 'nboot':
                                    if nboot != int(rdata.split()[1]): raise IOError("nboot missmatch")
                                elif rdata[0] == 'q':
                                    if len(thisMomList) > 0 and all([imom in currMomList for imom in thisMomList]): 
                                        BootNdimDict(thisoutputdict)
                                        return thisoutputdict
                                    if rdata in thisMomList or len(thisMomList) == 0:
                                        thismom = rdata
                                        currMomList.append(thismom)
                                        ReadMom = True
                                        if thismom not in thisoutputdict.keys():
                                            if ipar not in thisoutputdict[thismom].keys():
                                                print 'WARNING: '+thismom + ' '+ ipar +' found in:'
                                                print thisbootfn + ' but not in:'
                                                print filename
                                                thisoutputdict[thismom] = OrderedDict()
                                                thisoutputdict[thismom][ipar] = OrderedDict()
                                    else:
                                        ReadMom = False
                                elif rdata[0] == 'f' and ReadMom and twoptread:
                                    rdata = map(str,rdata.split())
                                    dictfitr = rdata[1]+'-'+rdata[2]
                                    if dictfitr not in thisoutputdict[thismom][ipar].keys():
                                        thisoutputdict[thismom][ipar][dictfitr] = OrderedDict()
                                    thisoutputdict[thismom][ipar][dictfitr]['Boot'] = BootStrap1(nboot,0.9)
                                elif rdata[0] == 'f' and ReadMom:
                                    rdata = map(str,rdata.split())
                                    dictfitr = rdata[1]+'-'+rdata[2]
                                    cutr = 'cut'+rdata[3]
                                    if dictfitr not in thisoutputdict[thismom][ipar].keys():
                                        thisoutputdict[thismom][ipar][dictfitr] = OrderedDict()
                                    if cutr not in thisoutputdict[thismom][ipar][dictfitr].keys():
                                        thisoutputdict[thismom][ipar][dictfitr][cutr] = OrderedDict()
                                    thisoutputdict[thismom][ipar][dictfitr][cutr]['Boot'] = BootStrap1(nboot,0.9)
                                elif ReadMom:
                                    try:
                                        if twoptread:
                                            thisoutputdict[thismom][ipar][dictfitr]['Boot'].values[int(rdata.split()[0])] = float(rdata.split()[1])*renorm
                                        else:
                                            thisoutputdict[thismom][ipar][dictfitr][cutr]['Boot'].values[int(rdata.split()[0])] = float(rdata.split()[1])*renorm
                                    except:
                                        if twoptread:
                                            thisoutputdict[thismom][ipar][dictfitr]['Boot'].values[int(rdata[0][:3])] = float('NaN')
                                        else:
                                            thisoutputdict[thismom][ipar][dictfitr][cutr]['Boot'].values[int(rdata[0][:3])] = float('NaN')
                    BootNdimDict(thisoutputdict)
                    return thisoutputdict
                outputdict = thisReadBoot(outputdict)
    return outputdict


## dataout = { Mass:Set/Avg/Std/Chi/Boot , FF#:qsqrd:Avg/Std/Boot , Chi:qsqrd}
def ReadFFFile(filename,bootfn=''):
    outputdict = OrderedDict()
    if os.path.isfile(filename):
        with open(filename,'r') as f:
            for line in f:
                rdata = line.strip()
                if len(rdata) == 0: continue
                rdata = rdata.split()
                if 'Mass' in rdata[0]:
                    outputdict['Mass'] = {'Set':rdata[1],
                                          'Avg':float(rdata[2]),
                                          'Std':float(rdata[3]),
                                          'Chi':float(rdata[4])}
                elif 'q**2' in rdata[0]:
                    for idata in rdata:
                        if 'FF' in idata and 'Err' not in idata :
                            outputdict[idata] = OrderedDict()
                        elif 'Chi' in idata:
                            outputdict['Chi'] = OrderedDict()
                elif 'qsqrd' in rdata[0]:
                    thisikey = 0
                    for ikey,dictkey in enumerate(outputdict.keys()):
                        if 'Mass' in dictkey: continue
                        thisikey += 1
                        if 'FF' in dictkey:
                            outputdict[dictkey][rdata[0]] = {'Avg':float(rdata[thisikey]),
                                                             'Std':float(rdata[thisikey+1])}
                            thisikey +=1
                        else:
                            outputdict[dictkey][rdata[0]] = float(rdata[thisikey])                            
        if os.path.isfile(bootfn):
            RFF,qsqrdread = None,None
            with open(bootfn,'r') as f:
                for line in f:
                    rdata = line.strip()
                    if len(rdata) == 0: continue
                    rdata = rdata.split()
                    if 'nboot' in rdata[0]:
                        nboot = int(rdata[2])
                    elif 'Mass' in rdata[0] and 'Mass' in outputdict.keys():
                        outputdict['Mass']['Boot'] = BootStrap1(nboot,1)
                        RFF = 'Mass'
                    elif 'FF' in rdata[0]:
                        RFF = rdata[0]
                    elif 'qsqrd' in rdata[0] and RFF in outputdict.keys():
                        qsqrdread = rdata[0]
                        outputdict[RFF][qsqrdread]['Boot'] = BootStrap1(nboot,1)
                    else:
                        if 'Mass' in RFF:
                            outputdict['Mass']['Boot'].values[int(rdata[0])] = float(rdata[1])
                        else:
                            outputdict[RFF][qsqrdread]['Boot'].values[int(rdata[0])] = float(rdata[1])
            BootNdimDict(outputdict)
    return outputdict
    




# def SetSumDictToList(DictData):
#     dataoutSum = []
#     for igamma,gammadata in enumerate(DictData.itervalues()):
#         dataoutSum.append([])
#         for imom,momdata in enumerate(gammadata.itervalues()):
#             dataoutSum[igamma].append([])
#             for iset,setdata in enumerate(momdata['SumMeth'].itervalues()):
#                 dataoutSum[igamma][imom].append([])
#                 for idict1,thedict1 in enumerate(setdata.itervalues()):
#                     dataoutSum[igamma][imom][iset].append([])
#                     for idict2,thedict2 in enumerate(thedict1.itervalues()):
#                         if 'Boot' in thedict2.keys():
#                             dataoutSum[igamma][imom][iset][idict1].append(thedict2['Boot'])
#                             BorA = 'Boot'
#                         else:
#                             dataoutSum[iset][igamma][imom][idict1][iset].append([thedict2['Avg'],thedict2['Std']])
#                             BorA = 'Avg'
#     return dataoutRF,BorA

# def SetSFDictToList(DictData):
#     dataoutSF = []
#     for igamma,gammadata in enumerate(DictData.itervalues()):
#         dataoutSF.append([])
#         for imom,momdata in enumerate(gammadata.itervalues()):
#             dataoutSF[igamma].append([])
#             thisCol = 0
#             for iCol,(theCol,coldata) in enumerate(momdata.iteritems()):
#                 if 'SF' in theCol:
#                     thisCol += 1
#                     dataoutSF[igamma][imom].append([])
#                     for iset,setdata in enumerate(gammadata.itervalues()):
#                         dataoutSF[igamma][imom][thisCol].append([])
#                         for idict1,thedict1 in enumerate(setdata.itervalues()):
#                             dataoutSF[igamma][imom][thisCol][iset].append([])
#                             for idict2,thedict2 in enumerate(thedict1.itervalues()):
#                                 dataoutSF[igamma][imom][thisCol][iset][idict1].append([])
#                                 for idict3,thedict3 in enumerate(thedict2.itervalues()):
#                                     if 'Boot' in thedict3.keys():
#                                         dataoutSF[igamma][imom][thisCol][iset][idict1][idict2].append(thedict3['Boot'])
#                                         BorA = 'Boot'
#                                     else:
#                                         dataoutSF[iset][igamma][imom][iCol][idict1][idict2].append([thedict3['Avg'],thedict3['Std']])
#                                         BorA = 'Avg'
#     return dataoutRF,BorA










# def OppDicts(datad,datas,opp,RFList,NT,thisSetList):
#     datacomb = deepcopy(datad)
#     datacomb['Title'] = NT
#     if 'RFBoot' in RFList:
#         for thekey,setd in datad['Boot'].iteritems():
#             for it,td in enumerate(setd):
#                 datacomb['Boot'][thekey][it] = ops[opp](td,datas['Boot'][thekey][it])
#                 datacomb['Boot'][thekey][it].Stats()
#         if 'RF' in RFList:
#             for keyval,dataset in datacomb['Boot'].iteritems():
#                 datacomb['Vals'][keyval] = Pullflag(dataset,'Avg')
#                 datacomb['Valserr'][keyval] = Pullflag(dataset,'Std')
#                 # for val,err in zip(Pullflag(dataset,'Avg'),Pullflag(dataset,'Std')):
#                 #     print val,err
#     if 'FitBoot' in RFList:
#         for thekey,setd in datad['Fit'].iteritems():
#             datacomb['Fit'][thekey]['Boot'] = ops[opp](setd['Boot'],datas['Fit'][thekey]['Boot'])
#             datacomb['Fit'][thekey]['Boot'].Stats()
#             if 'Fit' in RFList:
#                 datacomb['Fit'][thekey]['Avg'] = datacomb['Fit'][thekey]['Boot'].Avg
#                 datacomb['Fit'][thekey]['Std'] = datacomb['Fit'][thekey]['Boot'].Std
#     if 'SumBoot' in RFList:
#         for key1,datak1 in datad['Sum'].iteritems():
#             for key2,datak2 in datak1.iteritems():
#                 for key3,datak3 in datak2.iteritems():
#                     datacomb['Sum'][key1][key2][key3]['Boot'] = ops[opp](datak3['Boot'],datas['Sum'][key1][key2][key3]['Boot'])
#                     datacomb['Sum'][key1][key2][key3]['Boot'].Stats()
#                     if 'Sum' in RFList:
#                         datacomb['Sum'][key1][key2][key3]['Avg'] = datacomb['Sum'][key1][key2][key3]['Boot'].Avg
#                         datacomb['Sum'][key1][key2][key3]['Std'] = datacomb['Sum'][key1][key2][key3]['Boot'].Std
#     if 'TwoStateFit' in RFList:
#         datacomb.update(ReadTSFCombs(NT,thisSetList))
#     return datacomb




# ## DataDict { { StateList } { gamma } { mom } { contents } }
# def ReadSetFF(thisGammaList,thisSet,ReadList,thisindir):
#     DataDict = OrderedDict()
#     for thisstate in thisSet:
#         DataDict[thisstate] = OrderedDict()
#         for thisgamma in thisGammaList:
#             gammadirin = thisindir+CreateOppDir(thisgamma)+'/fits/'
#             DataDict[thisstate][thisgamma] = ReadDirFitsFF(gammadirin,thisstate,thisgamma,ReadList)
#     return DataDict




# def ReadDirFitsFF(inputdir,State,Gamma,ReadList):
#     momdict = OrderedDict()
#     filename = State + Gamma+'.fit.txt'
#     if Gamma == 'Mass' : filename = State[7:]+Gamma+'.fit.txt'
#     PickFit = FitRFCutPicked[GetTSinkFromSet(State)]
#     PickMassFit = FitMassPicked[State[7:]]
#     if 'Fit' in ReadList:
#         with open(inputdir + filename, 'r') as f:
#             for line in f:
#                 rdata = line.strip()
#                 if rdata[0] == 'q':
#                     thismom = rdata
#                     if thismom not in momdict.keys():
#                         momdict[thismom] = {}
#                 else:
#                     rdata = rdata.split()
#                     if Gamma == 'Mass':
#                         if (int(rdata[1]),int(rdata[2])) == PickMassFit:
#                             momdict[thismom]['FitAvg'] = float(rdata[3])
#                             momdict[thismom]['FitStd'] = float(rdata[4])
#                             if rdata[4] == 'nan':
#                                 momdict[thismom]['FitChi'] = float(0)
#                             else:
#                                 momdict[thismom]['FitChi'] = float(rdata[5])
#                     else:
#                         if int(rdata[1]) == PickFit :
#                             momdict[thismom]['FitAvg'] = float(rdata[2])
#                             momdict[thismom]['FitStd'] = float(rdata[3])
#                             if rdata[4] == 'nan':
#                                 momdict[thismom]['FitChi'] = float(0)
#                             else:
#                                 momdict[thismom]['FitChi'] = float(rdata[4])

#     if 'FitBoot' in ReadList:
#         ReadThis,ReadMassThis = False,False
#         with open(inputdir +'boots/' +filename.replace('.fit.txt','.fit.boot.txt'), 'r') as f:
#             for line in f:
#                 rdata = line.strip()
#                 rdatas = rdata.split()
#                 if rdatas[0] == 'q':
#                     thismom = rdata
#                     if thismom not in momdict.keys():
#                         momdict[thismom] = {}
#                     momdict[thismom]['FitBoot'] = BootStrap1(nboot,0.9)
#                 elif rdatas[0] == 'cut':
#                     if Gamma == 'Mass':
#                         if (int(rdatas[1]),int(rdatas[2])) == PickMassFit:
#                             ReadThis = True
#                         else:
#                             ReadThis = False
#                     else:
#                         if rdatas[1] == str(PickFit):
#                             ReadThis = True
#                         else:
#                             ReadThis = False
#                 elif rdatas[0] == 'nboot':
#                     nboot = int(rdatas[1])
#                 else:
#                     if ReadThis:
#                         iboot = int(rdatas[0])
#                         momdict[thismom]['FitBoot'].values[iboot] = float(rdatas[1])
#         for idict in momdict.itervalues():
#             idict['FitBoot'].Stats()
#     return momdict


# def ReadTSFCombs(TitlePrefix,LegList):
#     dictout = OrderedDict()
#     inputdir = outputdir+'/'+TitlePrefix+'/'
#     filelist = [ifl +'.txt' for ifl in LegList]
#     TSFfilelist = ReduceTsink(filelist)
#     TSFPars = TwoStateParList['C3']
#     Cuts = True
#     for iFF in TSFFileFlags:
#         dictout['TSF'+iFF] = OrderedDict()
#         for ifile in TSFfilelist:
#             dictfile = ifile.replace('.txt','')
#             for ipar in TSFPars:
#                 fname = (inputdir.replace(TitlePrefix,'TwoStateFit/'+TitlePrefix)
#                          +iFF+ifile.replace('.txt',TitlePrefix+ipar+'.txt'))
#                 if os.path.isfile(fname):
#                     with open(fname, 'r') as f:
#                               ##Two State Fit only dealing with zero momenta##
#                         if dictfile not in dictout['TSF'+iFF].keys():
#                             dictout['TSF'+iFF][dictfile] = OrderedDict()
#                         dictout['TSF'+iFF][dictfile][ipar] = OrderedDict()
#                         for line in f:
#                             try:
#                                 rdata = line.strip()
#                                 if len(rdata) > 0:
#                                     rdata = map(str,rdata.split())
#                                     dictfitr = rdata[0]+'-'+rdata[1]
#                                     if Cuts:
#                                         cutr = rdata[2]
#                                         if dictfitr not in dictout['TSF'+iFF][dictfile][ipar].keys():
#                                             dictout['TSF'+iFF][dictfile][ipar][dictfitr] = OrderedDict()
#                                         if cutr not in dictout['TSF'+iFF][dictfile][ipar][dictfitr].keys():
#                                             dictout['TSF'+iFF][dictfile][ipar][dictfitr][cutr] = OrderedDict()
#                                             dictout['TSF'+iFF][dictfile][ipar][dictfitr][cutr]['Avg'] = float(rdata[3])
#                                             dictout['TSF'+iFF][dictfile][ipar][dictfitr][cutr]['Std'] = float(rdata[4])
#                                             dictout['TSF'+iFF][dictfile][ipar][dictfitr][cutr]['Chi'] = float(rdata[5])
#                                     else:
#                                         if dictfitr not in dictout['TSF'+iFF][dictfile][ipar].keys():
#                                             dictout['TSF'+iFF][dictfile][ipar][dictfitr] = OrderedDict()
#                                         dictout['TSF'+iFF][dictfile][ipar][dictfitr]['Avg'] = float(rdata[2])
#                                         dictout['TSF'+iFF][dictfile][ipar][dictfitr]['Std'] = float(rdata[3])
#                                         dictout['TSF'+iFF][dictfile][ipar][dictfitr]['Chi'] = float(rdata[4])
#                             except:continue
#     return dictout


# def ReadDirFilesDict(inputdir,LegList,TitlePrefix,ReadList):
#     momdict = OrderedDict()
#     filelist = [ifl + TitlePrefix+'.txt' for ifl in LegList]
#     # print 'reffile: ' + inputdir + filelist[0]
#     with open(inputdir + filelist[0],'r') as f:
#         for line in f:
#             thismom = line.strip()
#             if thismom[0] == 'q':
#                 momdict[thismom] = OrderedDict()
#                 if thismom == 'q = 0 0 0':  momdict[thismom]['Title'] = TitlePrefix
#                 else: momdict[thismom]['Title'] = TitlePrefix+thismom
#                 momdict[thismom]['xLable'] = r'$ \frac{\tau}{a} - \frac{t}{2a}$'
#                 momdict[thismom]['yLable'] = r'$ R \left(\tau,t\right) $'
#                 momdict[thismom]['tVals'] = OrderedDict()
#                 momdict[thismom]['Vals'] = OrderedDict()
#                 momdict[thismom]['Valserr'] = OrderedDict()
#                 momdict[thismom]['Boot'] = OrderedDict()
#                 if 'Fit' in ReadList:
#                     momdict[thismom]['Fit'] = OrderedDict()
#                 if 'Sum' in ReadList:
#                     momdict[thismom]['Sum'] = OrderedDict()
#                 if 'TwoStateFit' in ReadList and thismom == 'q = 0 0 0':
#                     for iFF in TSFFileFlags:
#                         momdict[thismom]['TSF'+iFF] = OrderedDict()

#     for ifile in filelist:
#         dictfile = ifile.replace(TitlePrefix+'.txt','')
#         if 'RF' in ReadList:
#             # print inputdir + ifile
#             with open(inputdir + ifile, 'r') as f:
#                 for line in f:
#                     rdata = line.strip()
#                     if rdata[0] == 'q':
#                         thismom = rdata
#                         momdict[thismom]['tVals'][dictfile] = []
#                         momdict[thismom]['Vals'][dictfile] = []
#                         momdict[thismom]['Valserr'][dictfile] = []
#                     else:
#                         rdata = rdata.split()
#                         momdict[thismom]['tVals'][dictfile].append(int(rdata[0]))
#                         momdict[thismom]['Vals'][dictfile].append(float(rdata[1]))
#                         momdict[thismom]['Valserr'][dictfile].append(float(rdata[2]))

#         if 'RFBoot' in ReadList:
#             # print inputdir +'boots/' +ifile.replace('.txt','.boot.txt')
#             with open(inputdir +'boots/' +ifile.replace('.txt','.boot.txt'), 'r') as f:
#                 for line in f:
#                     rdata = line.strip()
#                     if rdata[0] == 'q':
#                         it = -1
#                         thismom = rdata
#                         momdict[thismom]['tVals'][dictfile] = []
#                         momdict[thismom]['Boot'][dictfile] = []
#                     elif rdata[0] == 't':
#                         it += 1
#                         rdata = rdata.split()
#                         momdict[thismom]['tVals'][dictfile].append(int(rdata[1]))
#                         momdict[thismom]['Boot'][dictfile].append(BootStrap1(nboot,0.9))
#                     elif rdata[0:5] == 'nboot':
#                         nboot = int(rdata.split()[1])
#                     else:
#                         rdata = rdata.split()
#                         iboot = int(rdata[0])
#                         momdict[thismom]['Boot'][dictfile][it].values[iboot] = float(rdata[1])
#             for idict in momdict.itervalues():
#                 for tdict in idict['Boot'][dictfile]:
#                     tdict.Stats()

#         if 'Fit' in ReadList:
#             if 'state' in dictfile: PickFit = FitCutPickedCM
#             elif 'tsink26' in dictfile: PickFit = FitCutPicked-2
#             elif 'tsink29' in dictfile: PickFit = FitCutPicked-1
#             else: PickFit = FitCutPicked
#             if TitlePrefix == 'twopt':
#                 thisinputdir = inputdir.replace('cfuns/twopt','Mass')
#                 thisifile = ifile.replace('twopt','Mass')
#                 PickFit = FitMassPicked[dictfile]
#             # print inputdir +'boots/' +ifile.replace('.txt','.boot.txt')
#             with open(thisinputdir +'fits/' +thisifile.replace('.txt','.fit.txt'), 'r') as f:
#                 for line in f:
#                     rdata = line.strip()
#                     if rdata[0] == 'q':
#                         thismom = rdata
#                         momdict[thismom]['Fit'][dictfile] = {}
#                     else:
#                         rdata = rdata.split()
#                         if TitlePrefix == 'twopt':
#                             if (int(rdata[1]),int(rdata[2])) == PickFit:
#                                 momdict[thismom]['Fit'][dictfile]['Avg'] = float(rdata[3])
#                                 momdict[thismom]['Fit'][dictfile]['Std'] = float(rdata[4])
#                                 if rdata[4] == 'nan':
#                                     momdict[thismom]['Fit'][dictfile]['Chi'] = float(0)
#                                 else:
#                                     momdict[thismom]['Fit'][dictfile]['Chi'] = float(rdata[5])

#                         else:
#                             if int(rdata[1]) == PickFit:
#                                 momdict[thismom]['Fit'][dictfile]['Avg'] = float(rdata[2])
#                                 momdict[thismom]['Fit'][dictfile]['Std'] = float(rdata[3])
#                                 if rdata[4] == 'nan':
#                                     momdict[thismom]['Fit'][dictfile]['Chi'] = float(0)
#                                 else:
#                                     momdict[thismom]['Fit'][dictfile]['Chi'] = float(rdata[4])

#         if 'FitBoot' in ReadList:
#             # print inputdir +'boots/' +ifile.replace('.txt','.boot.txt')
#             if 'state' in dictfile: PickFit = FitCutPickedCM
#             elif 'tsink26' in dictfile: PickFit = FitCutPicked-2
#             elif 'tsink29' in dictfile: PickFit = FitCutPicked-1
#             else: PickFit = FitCutPicked
#             with open(inputdir +'fits/boots/' +ifile.replace('.txt','.fit.boot.txt'), 'r') as f:
#                 for line in f:
#                     rdata = line.strip()
#                     rdatas = rdata.split()
#                     if rdatas[0] == 'q':
#                         thismom = rdata
#                         momdict[thismom]['Fit'][dictfile]['Boot'] = BootStrap1(nboot,0.9)
#                     elif rdatas[0] == 'cut':
#                         if rdatas[1] == str(PickFit):
#                             ReadThis = True
#                         else:
#                             ReadThis = False
#                     elif rdatas[0] == 'nboot':
#                         nboot = int(rdatas[1])
#                     else:
#                         if ReadThis:
#                             iboot = int(rdatas[0])
#                             momdict[thismom]['Fit'][dictfile]['Boot'].values[iboot] = float(rdatas[1])
#             for idict in momdict.itervalues():
#                 idict['Fit'][dictfile]['Boot'].Stats()

#     if 'Sum' in ReadList:
#         Sumfilelist = ReduceTsink(filelist,NoCM=True)
#         SumDefed = True
#         for ifile in Sumfilelist:
#             if 'sm32' in ifile:
#                 dictfile = ifile.replace(TitlePrefix+'.txt','')
#                 with open(inputdir +'SumMeth/' +ifile, 'r') as f:
#                     for line in f:
#                         rdata = line.strip()
#                         if len(rdata) > 0:
#                             if rdata[0] == 'q':
#                                 thismom = rdata
#                                 momdict[thismom]['Sum'][dictfile] = {}
#                             elif rdata[0] == 'c':
#                                 rdata = map(str,rdata.split())
#                                 dicpar1 = rdata[0]
#                                 dicpar2 = rdata[1][:-1]
#                                 if dicpar1 not in momdict[thismom]['Sum'][dictfile].keys():
#                                     momdict[thismom]['Sum'][dictfile][dicpar1] = OrderedDict()
#                                 if dicpar2 not in momdict[thismom]['Sum'][dictfile][dicpar1].keys():
#                                     momdict[thismom]['Sum'][dictfile][dicpar1][dicpar2] = OrderedDict()
#                                 momdict[thismom]['Sum'][dictfile][dicpar1][dicpar2]['Avg'] = float(rdata[2])
#                                 momdict[thismom]['Sum'][dictfile][dicpar1][dicpar2]['Std'] = float(rdata[3])
#                             elif rdata[0] == 'f':
#                                 rdata = map(str,rdata.split())
#                                 dicpar1 = rdata[1]
#                                 dicpar2 = 'fit '+rdata[2]+' '+rdata[3]+'-'+rdata[4][:-1]
#                                 if dicpar1 not in momdict[thismom]['Sum'][dictfile].keys():
#                                     momdict[thismom]['Sum'][dictfile][dicpar1] = OrderedDict()
#                                 if dicpar2 not in momdict[thismom]['Sum'][dictfile][dicpar1].keys():
#                                     momdict[thismom]['Sum'][dictfile][dicpar1][dicpar2] = OrderedDict()
#                                 momdict[thismom]['Sum'][dictfile][dicpar1][dicpar2]['Avg'] = float(rdata[5])
#                                 momdict[thismom]['Sum'][dictfile][dicpar1][dicpar2]['Std'] = float(rdata[6])
#                                 momdict[thismom]['Sum'][dictfile][dicpar1][dicpar2]['Chi'] = float(rdata[7])

#     if 'SumBoot' in ReadList:
#         Sumfilelist = ReduceTsink(filelist,NoCM=True)
#         for ifile in Sumfilelist:
#             if 'sm32' in ifile:
#                 dictfile = ifile.replace(TitlePrefix+'.txt','')
#                 with open(inputdir +'SumMeth/boots/' +ifile.replace('.txt','.boot.txt'), 'r') as f:
#                     for line in f:
#                         rdata = line.strip()
#                         if len(rdata) > 0:
#                             if rdata[0:5] == 'nboot':
#                                 nboot = int(rdata.split()[1])
#                             elif rdata[0] == 'q':
#                                 thismom = rdata
#                                 if not SumDefed: momdict[thismom]['Sum'][dictfile] = {}
#                             elif rdata[0] == 'c':
#                                 rdata = map(str,rdata.split())
#                                 dicpar1 = rdata[0]
#                                 dicpar2 = rdata[1][:-1]
#                                 if dicpar1 not in momdict[thismom]['Sum'][dictfile].keys():
#                                     momdict[thismom]['Sum'][dictfile][dicpar1] = {}
#                                 if dicpar2 not in momdict[thismom]['Sum'][dictfile][dicpar1].keys():
#                                     momdict[thismom]['Sum'][dictfile][dicpar1][dicpar2] = {}
#                                 momdict[thismom]['Sum'][dictfile][dicpar1][dicpar2]['Boot'] = BootStrap1(nboot,0.9)
#                             elif rdata[0] == 'f':
#                                 rdata = map(str,rdata.split())
#                                 dicpar1 = rdata[1]
#                                 dicpar2 = 'fit '+rdata[2]+' '+rdata[3]+'-'+rdata[4][:-1]
#                                 if dicpar1 not in momdict[thismom]['Sum'][dictfile].keys():
#                                     momdict[thismom]['Sum'][dictfile][dicpar1] = {}
#                                 if dicpar2 not in momdict[thismom]['Sum'][dictfile][dicpar1].keys():
#                                     momdict[thismom]['Sum'][dictfile][dicpar1][dicpar2] = {}
#                                 momdict[thismom]['Sum'][dictfile][dicpar1][dicpar2]['Boot'] = BootStrap1(nboot,0.9)
#                             else:
#                                 rdata = rdata.split()
#                                 momdict[thismom]['Sum'][dictfile][dicpar1][dicpar2]['Boot'].values[int(rdata[0])] = float(rdata[1])
#                     for idict in momdict.itervalues():
#                         for idicpar in idict['Sum'][dictfile].itervalues():
#                             for jdicpar in idicpar.itervalues():
#                                 jdicpar['Boot'].Stats()

#     if 'TwoStateFit' in ReadList:
#         thismom = 'q = 0 0 0'
#         if TitlePrefix == 'Mass' or TitlePrefix == 'twopt':
#             TSFfilelist = filelist
#             TSFPars = TwoStateParList['C2']
#             Cuts = False
#         else:
#             TSFfilelist = ReduceTsink(filelist)
#             TSFPars = TwoStateParList['C3']
#             Cuts = True
#         thisTSFFF = TSFFileFlags
#         for ifile in TSFfilelist:
#             dictfile = ifile.replace(TitlePrefix+'.txt','')
#             for ipar in TSFPars:
#                 for iFF in thisTSFFF:
#                     fname = (inputdir.replace(TitlePrefix,'TwoStateFit/'+TitlePrefix)
#                              +iFF+ifile.replace('.txt',ipar+'.txt'))
#                     if TitlePrefix == 'twopt':
#                         fname = fname.replace('twopt','Mass').replace('cfuns/','')
#                         if ipar == 'm0' or ipar == 'Dm':
#                             fname = fname.replace(dictfile,'')
#                     if os.path.isfile(fname):
#                         # print fname
#                         with open(fname, 'r') as f:
#                               ##Two State Fit only dealing with zero momenta##
#                             if dictfile not in momdict[thismom]['TSF'+iFF].keys():
#                                 momdict[thismom]['TSF'+iFF][dictfile] = OrderedDict()
#                             momdict[thismom]['TSF'+iFF][dictfile][ipar] = OrderedDict()
#                             for line in f:
#                                 try:
#                                     rdata = line.strip()
#                                     if len(rdata) > 0:
#                                         rdata = map(str,rdata.split())
#                                         dictfitr = rdata[0]+'-'+rdata[1]
#                                         if Cuts:
#                                             cutr = rdata[2]
#                                             if dictfitr not in momdict[thismom]['TSF'+iFF][dictfile][ipar].keys():
#                                                 momdict[thismom]['TSF'+iFF][dictfile][ipar][dictfitr] = OrderedDict()
#                                             if cutr not in momdict[thismom]['TSF'+iFF][dictfile][ipar][dictfitr].keys():
#                                                 momdict[thismom]['TSF'+iFF][dictfile][ipar][dictfitr][cutr] = OrderedDict()
#                                                 momdict[thismom]['TSF'+iFF][dictfile][ipar][dictfitr][cutr]['Avg'] = float(rdata[3])
#                                                 momdict[thismom]['TSF'+iFF][dictfile][ipar][dictfitr][cutr]['Std'] = float(rdata[4])
#                                                 momdict[thismom]['TSF'+iFF][dictfile][ipar][dictfitr][cutr]['Chi'] = float(rdata[5])
#                                         else:
#                                             if dictfitr not in momdict[thismom]['TSF'+iFF][dictfile][ipar].keys():
#                                                 momdict[thismom]['TSF'+iFF][dictfile][ipar][dictfitr] = OrderedDict()
#                                             momdict[thismom]['TSF'+iFF][dictfile][ipar][dictfitr]['Avg'] = float(rdata[2])
#                                             momdict[thismom]['TSF'+iFF][dictfile][ipar][dictfitr]['Std'] = float(rdata[3])
#                                             momdict[thismom]['TSF'+iFF][dictfile][ipar][dictfitr]['Chi'] = float(rdata[4])
#                                 except: continue
#     return momdict


# ## DataDict { gamma } { mom } { `above` } { (thisCombStateList) }
# def ReadSetDict(thisGammaList,thisCombStateList,ReadList,thisindir):
#     DataDict = OrderedDict()
#     for thisgamma in thisGammaList:
#         gammadirin = thisindir+CreateOppDir(thisgamma)+'/'
#         RLout = ReadList
#         CLSout = thisCombStateList
#         if thisgamma in ['twopt'] :
#             DataDict[thisgamma] = ReadDirFilesDict(gammadirin.replace(thisgamma,'cfuns/'+thisgamma)
#                                                    ,CLSout,thisgamma,RLout)
#         else:
#             DataDict[thisgamma] = ReadDirFilesDict(gammadirin,CLSout,thisgamma,RLout)
#     return DataDict

# def ReadSetnp(thisGammaList,thisCombStateList):
#     dictdata = ReadSetDict(thisGammaList,thisCombStateList,['RF','RFBoot'],outputdir)
#     return ConvertDictTOnp(dictdata)

# def ConvertDictTOnp(datadict):
#     dataout = []
#     for igamma,(thisgamma,datagamma) in enumerate(datadict.iteritems()):
#         dataout.append([])
#         for iq,(thisq,dataq) in enumerate(datagamma.iteritems()):
#             dataout[igamma].append([])
#             for thisstate,datastate in dataq['Boot'].iteritems():
#                 dataout[igamma][iq].append(datadict[thisgamma][thisq]['Boot'][thisstate])
#     return [dataout,SmallestMomList(datadict)]

# def ConvertCfunTOnp(datadict):
#     out3pt = []
#     for igamma,(thisgamma,datagamma) in enumerate(datadict.iteritems()):
#         if thisgamma == 'twopt':
#             out2pt = []
#             for iq,(thisq,dataq) in enumerate(datagamma.iteritems()):
#                 out2pt.append([])
#                 for thisstate,datastate in dataq['Boot'].iteritems():
#                     out2pt[iq].append(datadict[thisgamma][thisq]['Boot'][thisstate])
#         else:
#             out3pt.append([])
#             for iq,(thisq,dataq) in enumerate(datagamma.iteritems()):
#                 out3pt[igamma].append([])
#                 for thisstate,datastate in dataq['Boot'].iteritems():
#                     out3pt[igamma][iq].append(datadict[thisgamma][thisq]['Boot'][thisstate])
#     return [out2pt,out3pt,SmallestMomList(datadict)]

# ## DataDict { gamma } { mom } { `above` } { (thisCombStateList) }
# def ReadSetcfunDict(thisGammaList,thisCombStateList,ReadList,thisindir):
#     DataDict = OrderedDict()
#     for thisgamma in thisGammaList:
#         gammadirin = thisindir+CreateOppDir(thisgamma)+'/'
#         RLout = ReadList
#         CLSout = thisCombStateList
#         if thisgamma in ['twopt'] :
#             CLSout = CutDupSet(thisCombStateList)
#             if 'Sum' in RLout: RLout.remove('Sum')
#             if 'SumBoot' in RLout: RLout.remove('SumBoot')
#             if 'TwoStateFit' in RLout: RLout.remove('TwoStateFit')
#         DataDict[thisgamma] = ReadDirFilesDict(gammadirin,CLSout,thisgamma,RLout)
#     return DataDict


# def ReadCfunsnp(thisGammaList,thisCombStateList):
#     dictdata = ReadSetcfunDict(thisGammaList + ['twopt'], thisCombStateList,['RF','RFBoot'],outputdir+'cfuns/')
#     return ConvertCfunTOnp(dictdata)

#  LocalWords:  iSF

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

def StripSrc(filelist):
    output = []
    for ifile in filelist:
        thisout = re.sub('_xsrc.*','',ifile)[-4:]
        if thisout not in output:
            output.append(thisout)
    return output
        

## tflow [ itflow ]
## topcharge [ itflow ]

def ReadTopCharge(thisfile):
    if 'weinopp' in thisfile:
        thisfile = thisfile.replace('q_flow','W_flow')
    if os.path.isfile(thisfile):
        if Debug: print 'Reading: ',thisfile
        data = np.rollaxis(np.loadtxt(thisfile),1)
        if data.shape[1] < tflowlist[-1]:
            print 'file: ',thisfile,' does not have all flows'
            return [],[]
        else:
            tflow,topcharge = data[:,tflowlist]
            return tflow,topcharge
    else:
        return [],[]

## cfglistout [ icfg ]
## tflow [ icfg , itflow ]
## topcharge [ icfg , itflow ]
def ReadTopList(thisdir,thiscfglist,OnlyCheck=False):
    cfglistout,tflow,topcharge = [],[],[]
    for root, thedir, thesefiles in os.walk(thisdir):
        for icfg in thiscfglist:
            for ifile in thesefiles:
                if Debug: 'checking ',icfg,ifile
                if icfg in ifile:
                    if OnlyCheck:
                        if Debug: ifile, ' found for config ', icfg
                        cfglistout.append(icfg)
                    else:
                        thistflow,thistcharge = ReadTopCharge(root+ifile)
                        if icfg in cfglistout:
                            print 'warning, duplicate file in directory for config', icfg
                        elif len(thistflow) != 0:
                            cfglistout.append(icfg)
                            tflow.append(thistflow)
                            topcharge.append(thistcharge)
    return cfglistout,topcharge,tflow


## cfglistout [ icfg ]
## tflow [ icfg , itflow ]
## topcharge [ icfg , itflow ]
def ReadTopAll(thisdir):
    cfglistout,tflow,topcharge = [],[],[]
    filelist = []
    for root, thedir, thesefiles in os.walk(thisdir):
        filelist = thesefiles
        thisroot = root
    if len(filelist) == 0: raise IOError('No files in directory ' + thisdir)
    for ifile in filesort(filelist,[filelist[0][:-8],'.out']):
        thistflow,thistcharge = ReadTopCharge(thisroot+ifile)
        if len(thistflow) != 0:
            cfglistout.append(ifile)
            tflow.append(thistflow)
            topcharge.append(thistcharge)
    return cfglistout,topcharge,tflow

## dictout = { iset , imom , Info/ Boots: itflow , tsink }
def ReadAlphaSet(thisSetList,thisMomList,Wein=False):
    dictout = OrderedDict()
    for iset in thisSetList:
        dictout[iset] = ReadAlphaFile(outputdir[0],iset,thisMomList=thisMomList,Wein=Wein)['RF']
    return dictout

## output = [iAvg,boot1,boot2,...,bootn]
def ReadAlphaList(thisset,Wein=False):
    twoptset = ReduceTooMassSet([thisset])[0]
    DictRead,thisfile = ReadAlphaFitFile(outputdir[0],twoptset,thisMomList=['q = 0 0 0'],givefile=True,Wein=Wein)
    outdict = [1.0]
    infodict = OrderedDict()
    WeinOrTop = 'Top'
    if Wein: WeinOrTop = 'Wein'
    if 'q = 0 0 0' not in DictRead.keys():
        raise IOError(' Alpha Fit File not read in properly for' + twoptset )
    for thisflow,flowdict in DictRead['q = 0 0 0']['Boots'].iteritems():
        if thisflow in thisset:
            if AlphaFitRPick[WeinOrTop] in flowdict.keys():
                outdict = [flowdict[AlphaFitRPick[WeinOrTop]].Avg] + np.array(flowdict[AlphaFitRPick[WeinOrTop]].values).tolist()
                thischi = DictRead['q = 0 0 0']['Values'][thisflow][AlphaFitRPick[WeinOrTop]]['Chi']
            else:
                raise IOError(AlphaFitRPick[WeinOrTop]+ ' not in alpha file ' + twoptset)
    if outdict == [1.0]:
        raise IOError(thisset + ' does not contain the flow time in file ' )
    infodict['Avg'] = outdict[0]
    infodict['Std'] = np.std(outdict[1:])
    infodict['fit_range'] = AlphaFitRPick[WeinOrTop]
    infodict['Chi'] = thischi
    infodict['File'] = thisfile
    return outdict,infodict
    
# R/L Evecs [ ip , istate , ival ]
# Emass [ ip , istate ]
def ReadLREM(todtval,thisMomList,filepref,NoWar=False):
    ##OLD##
    # filename = REvecDir+filepref+'to'+str(todtval[0])+'dt'+str(todtval[1])+'LREM.txt'
    LEvec,REvec,Emass = [],[],[]
    for imom in thisMomList:
        LEvec.append([])
        REvec.append([])
        Emass.append([])
        imomCond = ipTOqcond(imom,Avg=True)
        filename = REvecDir+filepref+'to'+str(todtval[0])+'dt'+str(todtval[1])+'LREM'+imomCond+'.xml'
        # filename = REvecDir+VarPref+'TestLREM.txt'
        if Debug: print filename
        if os.path.isfile(filename):
            data,dump = ReadXmlDict(filename,Boot=False)
            for istate in GetStateSet('PoF'):
                # for istate,idata in data[imomCond]['Values'].iteritems():
                thisstate = 'State'+str(istate)
                Emass[-1].append(float(data[imomCond]['Values'][thisstate]['Emass']))
                LEvec[-1].append([])
                REvec[-1].append([])
                # for leflag,levec in data['Left_Evec'].iteritems():

                if 'Source_Evec' in data[imomCond]['Values'][thisstate]:
                    for thissm in DefjSmList:
                        if thissm not in data[imomCond]['Values'][thisstate]['Source_Evec'].keys():
                            LEvec[-1][-1].append(float(data[imomCond]['Values'][thisstate]['Source_Evec'][thissm.replace('jsm','sm')]))
                        else:
                            LEvec[-1][-1].append(float(data[imomCond]['Values'][thisstate]['Source_Evec'][thissm]))
                    for thissm in DefiSmList:
                        if thissm not in data[imomCond]['Values'][thisstate]['Sink_Evec'].keys():
                            REvec[-1][-1].append(float(data[imomCond]['Values'][thisstate]['Sink_Evec'][thissm.replace('ism','sm')]))
                        else:
                            REvec[-1][-1].append(float(data[imomCond]['Values'][thisstate]['Sink_Evec'][thissm]))
                    
                
                ## Legacy, Left ==  jsm == Source
                ##         Right == ism == Sink
                else:
                    for thissm in DefjSmList:
                        if thissm not in data[imomCond]['Values'][thisstate]['Left_Evec'].keys():
                            LEvec[-1][-1].append(float(data[imomCond]['Values'][thisstate]['Left_Evec'][thissm.replace('jsm','sm')]))
                        else:
                            LEvec[-1][-1].append(float(data[imomCond]['Values'][thisstate]['Left_Evec'][thissm]))
                    for thissm in DefiSmList:
                        if thissm not in data[imomCond]['Values'][thisstate]['Right_Evec'].keys():
                            REvec[-1][-1].append(float(data[imomCond]['Values'][thisstate]['Right_Evec'][thissm.replace('ism','sm')]))
                        else:
                            REvec[-1][-1].append(float(data[imomCond]['Values'][thisstate]['Right_Evec'][thissm]))
        else:
            if not NoWar:
                print 'warning, weight file not found', filename
            return None,None,None        
    return np.array(LEvec),np.array(REvec),np.array(Emass)

## readdata { gamma } { mom } { method } { set }
## datadictout { collection } { gamma } { mom }
## datamassout { collection }
def ExtractValues(thisindir,thisGammaList,thisSetList,thisMethodList,thisMomList=qvecSet,thisPrintRead=PrintRead,TopRead=False):
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
                    for itvar in DefTvarPicked:
                        if itvar in ism: ism = PickedStateStr+itvar
                    if 'TSF' in iMeth:
                        if CheckDict(thismassdict,'m0',TSFfitr):
                            datamassout[ism+iMeth] = thismassdict['m0'][TSFfitr]
                            datamassout[ism+iMeth]['Info'] = thismassdict['m0']['Info']
                    elif 'OSF' in iMeth:
                        if CheckDict(thismassdict,'m0',OSFfitr[CreateOSFfitKey(ism)[0]]):
                            datamassout[ism+iMeth] = thismassdict['m0'][OSFfitr[CreateOSFfitKey(ism)[0]]]
                            datamassout[ism+iMeth]['Info'] = thismassdict['m0']['Info']
        else:
            for imom,momdata in readdata[igamma].iteritems():
                for iMeth,Methdata in momdata.iteritems():
                    for iSet,thisdict in Methdata.iteritems():
                        itsink,ism = SplitTSinkString(iSet)
                        fitsm = CreateOSFfitKey(ism)[0]
                        for itvar in DefTvarPicked:
                            if itvar in ism: ism = PickedStateStr+itvar
                        if 'TSF' in iMeth:
                            for icut in TSFCutList:
                                if CheckDict(thisdict,'B00',TSFfitr,icut):
                                    datadictout = SetupDict(datadictout,igamma,iSet+iMeth+icut)
                                    datadictout[iSet+iMeth+icut][igamma][imom] = thisdict['B00'][TSFfitr][icut]
                                    datadictout[iSet+iMeth+icut][igamma][imom]['Info'] = thisdict['B00']['Info']
                        elif 'OSF' in iMeth:
                            for icut in OSFCutList:
                                if CheckDict(thisdict,'B00',OSFfitr[fitsm],icut):
                                    datadictout = SetupDict(datadictout,igamma,iSet+iMeth+icut)
                                    datadictout[iSet+iMeth+icut][igamma][imom] = thisdict['B00'][OSFfitr[fitsm]][icut]
                                    datadictout[iSet+iMeth+icut][igamma][imom]['Info'] = thisdict['B00']['Info']
                        elif 'SumMeth' in iMeth:
                            for ifit in SumFitRList:
                                for icut in SumCutList:
                                    if CheckDict(thisdict,icut,ifit):
                                        fitdict = ifit.replace('fit sl ','fitr')
                                        datadictout = SetupDict(datadictout,igamma,iSet+iMeth+icut+fitdict)
                                        datadictout[iSet+iMeth+icut+fitdict][igamma][imom] = thisdict[icut][ifit]
                                        datadictout[iSet+iMeth+icut+fitdict][igamma][imom]['Info'] = thisdict['Info']
                        elif 'Fits' in iMeth:
                            if 'Top' in igamma or 'Wein' in igamma:
                                if 'Wein' in igamma: thisfa = WeinFlowArgs
                                else: thisfa = FlowArgs
                                for iflow in thisfa:
                                    if iflow not in thisdict.keys(): continue
                                    for icut in FitCutArgs:
                                        if icut not in thisdict[iflow].keys(): continue
                                        datadictout = SetupDict(datadictout,igamma,iSet+iMeth+icut+iflow)
                                        datadictout[iSet+iMeth+icut+iflow][igamma][imom] = thisdict[iflow][icut]
                                        datadictout[iSet+iMeth+icut+iflow][igamma][imom]['Info'] = thisdict['Info']
                            elif TopRead:
                                if 'Wein' == TopRead: thisfa = WeinFlowArgs
                                else: thisfa = FlowArgs
                                for iflow in thisfa:
                                    for icut in FitCutArgs:
                                        if icut not in thisdict.keys(): continue
                                        datadictout = SetupDict(datadictout,igamma,iSet+iMeth+icut+iflow)
                                        datadictout[iSet+iMeth+icut+iflow][igamma][imom] = thisdict[icut]
                                        datadictout[iSet+iMeth+icut+iflow][igamma][imom]['Info'] = thisdict['Info']
                            else:
                                for icut in FitCutArgs:
                                    if icut in thisdict.keys():
                                        datadictout = SetupDict(datadictout,igamma,iSet+iMeth+icut)
                                        datadictout[iSet+iMeth+icut][igamma][imom] = thisdict[icut]
                                        datadictout[iSet+iMeth+icut][igamma][imom]['Info'] = thisdict['Info']
    if thisPrintRead: print 'Extracting data took: ', str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s                  '
    return datadictout,datamassout




def Get2ptSetMoms(thisoutputdir,MomListIn,tvarlist=[],ismlist=[],jsmlist=[],tsrclist=[]):
    momlist = set([])
    xmlMomList = map(qstrTOqcond,MomListIn)
    for iflag in ['cfun/twopt','Mass']:
        for ip in xmlMomList:
            thisdir = thisoutputdir+iflag+MakeMomDir(ip)
            for itvar in tvarlist:
                if iflag == 'Mass':
                    ifile = thisdir+itvar+'LREM'+ip+'.xml'
                    if not CheckMomFile(ifile): momlist.add(qcondTOqstr(ip))
                for istate in GetStateSet(itvar):
                    ifile = thisdir+'state'+istate+itvar+iflag.replace('cfun/','')+ip+'.xml'
                    if not CheckMomFile(ifile): momlist.add(qcondTOqstr(ip))
            for its in tsrclist:
                for ism in ismlist:
                    for jsm in jsmlist:
                        ts,ismstr,jsmstr = str(its),str(ism),str(jsm)
                        if 'tsrc' not in ts: ts = 'tsrc'+ts
                        if 'ism' not in ismstr: ismstr = 'ism'+ismstr
                        if 'jsm' not in jsmstr: jsmstr = 'jsm'+jsmstr
                        ifile = thisdir+ts+ismstr+jsmstr+iflag.replace('cfun/','')+ip+'.xml'
                        if not CheckMomFile(ifile): momlist.add(qcondTOqstr(ip))
    return OrderMomList(momlist)

def Get2ptSetMoms(thisoutputdir,MomListIn,tvarlist=[],ismlist=[],jsmlist=[],tsrclist=[]):
    momlist = set([])
    xmlMomList = map(qstrTOqcond,MomListIn)
    for iflag in ['cfun/twopt','Mass']:
        for ip in xmlMomList:
            thisdir = thisoutputdir+iflag+MakeMomDir(ip)
            for itvar in tvarlist:
                if iflag == 'Mass':
                    ifile = thisdir+itvar+'LREM'+ip+'.xml'
                    if not CheckMomFile(ifile): momlist.add(qcondTOqstr(ip))
                for istate in GetStateSet(itvar):
                    ifile = thisdir+'state'+istate+itvar+iflag.replace('cfun/','')+ip+'.xml'
                    if not CheckMomFile(ifile): momlist.add(qcondTOqstr(ip))
            for its in tsrclist:
                for ism in ismlist:
                    for jsm in jsmlist:
                        ts,ismstr,jsmstr = str(its),str(ism),str(jsm)
                        if 'tsrc' not in ts: ts = 'tsrc'+ts
                        if 'ism' not in ismstr: ismstr = 'ism'+ismstr
                        if 'jsm' not in jsmstr: jsmstr = 'jsm'+jsmstr
                        ifile = thisdir+ts+ismstr+jsmstr+iflag.replace('cfun/','')+ip+'.xml'
                        if not CheckMomFile(ifile): momlist.add(qcondTOqstr(ip))
    return OrderMomList(momlist)


def GetTopSetMoms(thisoutputdir,MomListIn,tvarlist=[],ismlist=[],jsmlist=[],tsrclist=[]):
    momlist = set([])
    xmlMomList = map(qstrTOqcond,MomListIn)
    for iflag in ['cfun/twopt','cfun/NNQ','Alpha']:
        flagrep = iflag.replace('twopt','')
        for ip in xmlMomList:
            thisdir = thisoutputdir+iflag+MakeMomDir(ip)
            for itvar in tvarlist:
                for istate in GetStateSet(itvar):
                    ifile = thisdir+'state'+istate+itvar+iflag.replace(flagrep,'')+ip+'.xml'
                    if not CheckMomFile(ifile): momlist.add(qcondTOqstr(ip))
            for its in tsrclist:
                for ism in ismlist:
                    for jsm in jsmlist:
                        ts,ismstr,jsmstr = str(its),str(ism),str(jsm)
                        if 'tsrc' not in ts: ts = 'tsrc'+ts
                        if 'ism' not in ismstr: ismstr = 'ism'+ismstr
                        if 'jsm' not in jsmstr: jsmstr = 'jsm'+jsmstr
                        ifile = thisdir+ts+ismstr+jsmstr+iflag.replace(flagrep,'')+ip+'.xml'
                        if not CheckMomFile(ifile): momlist.add(qcondTOqstr(ip))
    return OrderMomList(momlist)



def Get3SM(thisoutputdir,thisGammaList,MomListIn,setlist):
    momlist = set([])
    xmlMomList = map(qstrTOqcond,MomListIn)
    for igamma in thisGammaList:
        thisdir = thisoutputdir+CreateOppDir(igamma)
        if 'cfun' in thisdir:
            thisdir = thisdir.replace('cfun/Top','Top/cfun')
        for iset in setlist:
            for ip in xmlMomList:
                ifile = thisdir+MakeMomDir(ip)+iset+igamma+ip+'.xml'
                if not CheckMomFile(ifile): momlist.add(qcondTOqstr(ip))
    return OrderMomList(momlist)
            

def Get3ptSetMoms(thisoutputdir,thisGammaList,MomListIn,setlist):
    return OrderMomList(set(Get3SM(thisoutputdir,thisGammaList,MomListIn,setlist)) |
                        set(Get3SM(thisoutputdir+'cfun/',thisGammaList,MomListIn,setlist)))


##NEW FUNCTIONS##

#data2pt = [ mom , set , it ] bs
#data3pt = [ gamma , mom , set , it ] bs
def ReadCfunsnp(thisGammaList,thisSetList,thisMomList=RunMomList,thisPrintRead=PrintRead):
    datadict = ReadCfunsDict(outputdir[0],thisSetList,thisGammaList,thisMomList=thisMomList,thisPrintRead=thisPrintRead)
    return SetRFDictToList(datadict,thisPrintRead=thisPrintRead)

def ReadRFnp(thisGammaList,thisSetList,thisMomList=RunMomList,thisPrintRead=PrintRead):
    datadict = ReadSetDict(outputdir[0],thisSetList,thisGammaList,['RF'],thisMomList=thisMomList,thisPrintRead=thisPrintRead)
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
    sflist = ['OSF'+iOSF for iOSF in OSFFileFlags]+['TSF'+iTSF for iTSF in TSFFileFlags]
    if any([imeth in sflist for imeth in thisMethodList]):
        data3ptdict = ReadCfunsDict(thisindir,thisSetList,thisGammaList,thisMomList=thisMomList)
    datadict = ReadSetDict(thisindir,thisSetList,thisGammaList,thisMethodList,thisMomList=thisMomList,thisPrintRead=thisPrintRead)
    zmomstr = 'q = 0 0 0'
    start = time.time()
    if ScaleByP4g4: datadict = RewriteRF(datadict,data3ptdict,'P4g4',zmomstr)
    for igamma in datadict.keys():
        if igamma == 'twopt': continue
        if thisPrintRead: print 'Constructing Fitted RF Values: ' , igamma , '     \r',
        for imom in datadict[igamma].iterkeys():
            twoptmom = GetAvgMom(imom)
            # if zmomstr not in datadict[igamma].keys(): continue
            if not CheckDict(datadict,igamma,imom,'RF'): continue
            for iset in datadict[igamma][imom]['RF'].keys():
                twoptiset = 'tsrc'+str(tsource)+RemoveTSink(iset)
                if not CheckDict(datadict,'twopt',twoptmom,'RF',twoptiset): continue
                # if thisPrintRead: print RemoveTSinkTsrc(iset)+' not in two point set list, not constructing RF'
                data3pt = data3ptdict[igamma][imom]['RF'][iset]['Boot']
                for iSF in sflist:
                    # print 'ONE ',datadict['twopt'][imom].keys(), iSF
                    # if CheckDict(datadict['twopt'][imom],iSF): print 'Two ',datadict['twopt'][imom][iSF].keys(), twoptiset, iset
                    if CheckDict(datadict['twopt'][twoptmom],iSF,twoptiset):
                        pars2pt = []
                        if 'OSF' in iSF:
                            if not all([iState in datadict['twopt'][twoptmom][iSF][twoptiset].keys() for iState in StateParList['One']['C2']]): continue
                            for ipar in StateParList['One']['C2']:
                                fitrkey = RemoveTSinkTsrc(iset)
                                for itvar in DefTvarPicked:
                                    if itvar in fitrkey: fitrkey = PickedStateStr+itvar
                                pars2pt.append(datadict['twopt'][twoptmom][iSF][twoptiset][ipar][OSFfitr[CreateOSFfitKey(fitrkey)[0]]]['Boot'])
                            data2ptZ =  ff.C2OneStateFitFunNoExp(np.array([GetintTSink(iset)-tsource]),pars2pt)
                        elif 'TSF' in iSF:
                            if not all([iState in datadict['twopt'][twoptmom][iSF][twoptiset].keys() for iState in StateParList['Two']['C2']]): continue
                            for ipar in StateParList['Two']['C2']:
                                pars2pt.append(datadict['twopt'][twoptmom][iSF][twoptiset][ipar][TSFfitr]['Boot'])
                            data2ptZ =  ff.C2TSFLineFun(GetintTSink(iset)-tsource,pars2pt)
                        if 'RF'+iSF not in datadict[igamma][imom].keys(): datadict[igamma][imom]['RF'+iSF] = OrderedDict()
                        if iset not in datadict[igamma][imom]['RF'+iSF].keys(): datadict[igamma][imom]['RF'+iSF][iset] = OrderedDict()
                        datadict[igamma][imom]['RF'+iSF][iset]['Boot'] = [idata3pt/data2ptZ for idata3pt in data3pt[tsource-1:GetintTSink(iset)]]
                        datadict[igamma][imom]['RF'+iSF][iset]['tVals'] = datadict[igamma][imom]['RF'][iset]['tVals']
                        GetBootStats(datadict[igamma][imom]['RF'+iSF][iset]['Boot'])
                        datadict[igamma][imom]['RF'+iSF][iset]['Vals'] = Pullflag(datadict[igamma][imom]['RF'+iSF][iset]['Boot'],'Avg')
                        datadict[igamma][imom]['RF'+iSF][iset]['Valserr'] = Pullflag(datadict[igamma][imom]['RF'+iSF][iset]['Boot'],'Std')
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
            thisfile = thisdir+CreateCurrCombFn(thisFF)+thisset+'.xml'
            mprint('FFread: ' ,thisfile)
            if '/' in thisFF:
                DataDict[thisFF][thisset] = ReadFFCombFile(thisfile)            
            else:
                DataDict[thisFF][thisset] = ReadFFFile(thisfile)                                
        if thisPrintRead: print 'Reading ',thisFF,'took : ' , str(datetime.timedelta(seconds=time.time()-startff)) , ' h:m:s                     '
    if thisPrintRead: print 'Reading all FFs took : ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s                     '
    return DataDict


## DataDict {kappa} { FormFactor } { Set } { Mass:Set/Avg/Std/Chi/Boot , FF#:qsqrd:Avg/Std/Boot , Chi:qsqrd}
## N.B. Set under Mass is different to Set, it is the set sellected for the mass extraction.
## thisFFDict = { keys=(scalar,vector etc) } {kappa } {values=(setlist for current) / Mass:Set/Avg/Std/Chi/Boot / FF#:qsqrd:Avg/Std/Boot / Chi:qsqrd}

def ReadMKFFDict(thisindir,thisFFDict,thisPrintRead=PrintRead):
    DataDict = OrderedDict()
    if thisPrintRead: print ''
    start = time.time()
    for iFF,(thisFF,FFKappaList) in enumerate(thisFFDict.iteritems()):
        DataDict[thisFF] = OrderedDict()
        for ikappa,FFSetList in FFKappaList.iteritems():
            startff = time.time()
            DataDict[thisFF][ikappa] = OrderedDict()
            for iset,thisset in enumerate(FFSetList):
                if thisPrintRead: print 'Reading ',thisFF,'at : ' ,int((iset*100)/float(len(FFSetList))),'%     \r',
                thisdir = thisindir.replace(str(kappa),ikappa)+'FormFactors/'+thisFF+'/'
                thisfile = thisdir+CreateCurrCombFn(thisFF)+thisset+'.xml'
                mprint('FFread: ' ,thisfile)
                if '/' in thisFF:
                    DataDict[thisFF][ikappa][thisset] = ReadFFCombFile(thisfile)            
                else:
                    DataDict[thisFF][ikappa][thisset] = ReadFFFile(thisfile)                                
            if thisPrintRead: print 'Reading ',thisFF,'took : ' , str(datetime.timedelta(seconds=time.time()-startff)) , ' h:m:s                     '
    if thisPrintRead: print 'Reading all MKFFs took : ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s                     '
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
            gammadirin = thisindir+'cfun/twopt/'
            massSetList = []
            for iSet in thisSetList:
                if 'sm' in iSet and 'state' not in iSet:
                    massSetList.append('tsrc'+str(tsource)+iSet)
                else:
                    massSetList.append(iSet)
                
            DataDict[thisgamma] = MakeMethodsDict(gammadirin,'twopt.xml',
                                                  massMethodList,thisSetList,thisMomList=thisMomList,
                                                  thisPrintRead=thisPrintRead)
        else:
            gammadirin = thisindir+CreateOppDir(thisgamma)+'/'
            DataDict[thisgamma] = MakeMethodsDict(gammadirin,thisgamma+'.xml',
                                                  thisMethodList,thisSetList,thisMomList=thisMomList,
                                                  thisPrintRead=thisPrintRead)
    # if 'RF' in thisMethodList: DataDict = CombSetBoot(DataDict,'-',thisPrintRead=thisPrintRead)
    if thisPrintRead: print 'Reading took: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s                     '
    return DataDict

def ReadCfunsDict(thisindir,thisSetList,thisGammaList,thisMomList=RunMomList,thisPrintRead=PrintRead):
    DataDict = OrderedDict()
    if thisPrintRead: print ''
    start = time.time()
    for ig,thisgamma in enumerate(thisGammaList):
        if thisPrintRead: print 'Reading Cfuns at: ' ,int((ig*100)/float(len(thisGammaList))),'%     \r',
        if thisgamma == 'twopt':
            gammadirin = thisindir+'cfun/twopt/'
            DataDict[thisgamma] = MakeMethodsDict(gammadirin,'twopt.xml',
                                                  ['RF'],thisSetList,thisMomList=thisMomList,
                                                  thisPrintRead=thisPrintRead)
        else:
            gammadirin = thisindir+'cfun/'+CreateOppDir(thisgamma)+'/'
            DataDict[thisgamma] = MakeMethodsDict(gammadirin,thisgamma+'.xml',
                                                  ['RF'],thisSetList,thisMomList=thisMomList,
                                                  thisPrintRead=thisPrintRead)
    if thisPrintRead: print 'Reading Cfuns took: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s                     '
    # DataDict = CombSetBoot(DataDict,'-',thisPrintRead=thisPrintRead)
    return DataDict


def SetRFDictToList(DictData,thisPrintRead=PrintRead):
    def BAswitch(BorA):
        if BorA == 'Avg': return 'Mixed'
        else: return 'Boot'

    dataoutRF,dataout2pt = [],[]
    infolistRF,infolist2pt = [],[]
    gammalistout = OrderedDict()
    flowlist = []
    BorA = ''
    start = time.time()
    for igamma,(thisgamma,gammadata) in enumerate(DictData.iteritems()):
        if thisPrintRead: print 'Converting DictToList ' + thisgamma+ '          \r',
        gammalistout[thisgamma] = []
        if 'twopt' in thisgamma:
            for imom,(thismom,momdata) in enumerate(gammadata.iteritems()):
                gammalistout[thisgamma].append(thismom)
                dataout2pt.append([])
                infolist2pt.append([])
                for iset,setdata in enumerate(momdata['RF'].itervalues()):
                    infolist2pt[imom].append(setdata['Info'])
                    if 'Boot' in setdata.keys():
                        dataout2pt[imom].append(setdata['Boot'])
                    elif 'Avg' in setdata.keys():
                        dataout2pt[imom].append([setdata['Avg'],setdata['Std']])
                    BorA = BAswitch(BorA)
                            
        else:
            dataoutRF.append([])
            infolistRF.append([])
            for imom,(thismom,momdata) in enumerate(gammadata.iteritems()):
                gammalistout[thisgamma].append(thismom)
                dataoutRF[-1].append([])
                infolistRF[-1].append([])
                for iset,(thisset,setdata) in enumerate(momdata['RF'].iteritems()):
                    infolistRF[-1][imom].append(setdata['Info'])
                    if 'Top' in thisgamma or 'Wein' in thisgamma:
                        flowlist = []
                        dataoutRF[-1][imom].append([])
                        for icf,(iflow,flowdata) in enumerate(setdata.iteritems()):
                            if 't_flow' not in iflow: continue
                            flowlist.append(iflow)
                            if 'Boot' in flowdata.keys():
                                dataoutRF[-1][imom][iset].append(flowdata['Boot'])
                            elif 'Avg' in flowdata.keys():
                                dataoutRF[-1][imom][iset].append([flowdata['Avg'],flowdata['Std']])
                    elif 'Boot' in setdata.keys():
                        dataoutRF[-1][imom].append(setdata['Boot'])
                    else:
                        dataoutRF[-1][imom].append([setdata['Avg'],setdata['Std']])
                    BorA = BAswitch(BorA)
    if 'twopt' in gammalistout.keys():
        if 'q = 0 0 0' in gammalistout['twopt']:
            dataout2pt.insert(0, dataout2pt.pop(gammalistout['twopt'].index('q = 0 0 0')))
            infolist2pt.insert(0, infolist2pt.pop(gammalistout['twopt'].index('q = 0 0 0')))
            gammalistout['twopt'].insert(0, gammalistout['twopt'].pop(gammalistout['twopt'].index('q = 0 0 0')))
    if thisPrintRead: print 'Converting DictToList took: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s                     '
    return dataoutRF,dataout2pt,gammalistout,BorA,infolistRF,infolist2pt,flowlist
##datadict.keys() is gammalist



# def CombBoot(udata,ddata,opp):
#     dataout = {'Boot':ops[opp](np.array(udata['Boot']),np.array(ddata['Boot']))}
#     GetBootStats(dataout['Boot'])
#     dataout['tVals'] = udata['tVals']
#     dataout['Info'] = min(udata['Info'],ddata['Info'])
#     dataout['Vals'] = Pullflag(dataout['Boot'],'Avg')
#     dataout['Valserr'] = Pullflag(dataout['Boot'],'Std')
#     return dataout

# def CombMethBoot(udata,ddata,opp,imom,igamma):
#     dataout = OrderedDict()
#     for iSet in udata.keys():
#         if 'Boot' in udata[iSet].keys() and 'Boot' in ddata[iSet].keys():
#             dataout[iSet] = CombBoot(udata[iSet],ddata[iSet],opp)
#         else:
#             raise IOError(iSet +' does not contain Boot for ' + igamma + ' ' + imom)
#     return dataout


# def CombSetBoot(data,opp,thisPrintRead=PrintRead):
#     dataout = deepcopy(data)
#     start = time.time()
#     for dgamma,dgammadata in data.iteritems():
#         if 'doub' in dgamma:
#             combgamma = dgamma.replace('doub','')
#             for sgamma,sgammadata in data.iteritems():
#                 if 'sing'+combgamma == sgamma:
#                     if thisPrintRead: print 'Combining ' , dgamma , ' ' , sgamma ,'             \r',
#                     if combgamma not in dataout.keys():
#                         dataout[combgamma] = OrderedDict()
#                     for imom in dgammadata.iterkeys():
#                         if imom in sgammadata.keys():
#                             if 'RF' in sgammadata[imom].keys():
#                                 if imom not in dataout[combgamma].keys():
#                                     dataout[combgamma][imom] = OrderedDict()
#                                 dataout[combgamma][imom]['RF'] = CombMethBoot(data[sgamma][imom]['RF'],
#                                                                               data[dgamma][imom]['RF'],opp,imom,combgamma)
#                             else:
#                                 print ''
#                                 print 'Warning: ', imom , ' in ' , dgamma , ' but not in ' , sgamma 
#                         else:
#                             print ''
#                             print 'Warning: ', imom , ' in ' , dgamma , ' but not in ' , sgamma 
#     if thisPrintRead: print 'Combining took: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s                     '
#     return dataout


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
            for iSet in ReduceTsink(loopSetList):
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



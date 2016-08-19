#!/usr/bin/env python

from Params import *
from FitParams import *
from BootTest import BootStrap1
from MiscFuns import *
from SetLists import GetTsinkSmLists
from OppFuns import *
from FormFactors import NoFFPars
import xmltodict
import cPickle as pickle
from collections import OrderedDict
from XmlFuns import *
from XmlFormatting import *
from ReadXml import *
import cPickle as pickle
import os


# def WriteXmlOutput(thisfile,outputdict):
#     with open(thisfile+'.xml','w') as f:
#         f.write( xmltodict.unparse(outputdict,pretty=True))

## listin = [ism] {'Info': {'nconfig':###}}
## dictout = {'Info': {'nconfig':###}}
def CombineSetInfo(listin):
    nconf = listin[0]['nconfig']
    if len(listin) > 1:
        for ismlist in listin[1:]:
            nconf = min(ismlist['nconfig'],nconf)
    return {'nconfig':nconf}

def WriteXml(thisfile,outputdict):
    if Debug: print 'Writing to: ' , thisfile
    with open(thisfile+'.xml','w') as f:
        f.write( xmltodict.unparse(outputdict,pretty=True))
    
def WriteXmlOutput(thisfile,outputdict):
    firstkey = outputdict.keys()[0]
    Vals = {firstkey:{'Values':outputdict[firstkey]['Values']}}
    if 'Info' in outputdict[firstkey].keys(): Vals[firstkey]['Info'] = outputdict[firstkey]['Info']
    Boots = outputdict[firstkey]['Boots']
    outdirlist = thisfile.split('/')
    bootdir = '/'.join(outdirlist[:-1]+['boots'])
    bootout = '/'.join(outdirlist[:-1]+['boots']+[outdirlist[-1]])
    mkdir_p(bootdir)
    Vals[firstkey]['Boots'] = bootout+'.boot.p'
    WriteXml(thisfile,Vals)
    with open( bootout+'.boot.p', "wb" ) as pfile:
        pickle.dump( Boots, pfile )
    
        
def MergeXmlOutput(thisfile,outputdict,CheckMom=True):
    if CheckMom:
        if CheckMomFile(thisfile+'.xml'):
            thisdict = ReadXmlAndPickle(thisfile+'.xml')[0]
            outputdict = merge_dicts(outputdict,thisdict)
    else:
        thisdict = ReadXmlAndPickle(thisfile+'.xml')[0]
        outputdict = merge_dicts(outputdict,thisdict)
    if len(outputdict.keys()) > 1:
        raise IOError('Xml main key not single: ' + ', '.join(outputdict.keys()))
    WriteXmlOutput(thisfile,outputdict)

def SetUpPDict(ip,filedir,filename):
    datadict = {ip:{'Values':OrderedDict(),'Boots':OrderedDict()}}
    outputfile = filedir+MakeMomDir(ip)
    mkdir_p(outputfile)
    outputfile = outputfile+filename+ip
    return datadict,outputfile
    

def PrintToFile(thisdata,filedir,filename,thisTList,thisMomList,AddDict={},frmtflag='f'):
    xmlMomList = map(ipTOqcond,thisMomList)
    tkeyList = map(tstr,thisTList)
    outputfilelist = []
    for ip,pdata in zip(xmlMomList,thisdata):
        datadict,outputfile = SetUpPDict(ip,filedir,filename)
        datadict[ip]['Info'] = AddDict
        datadict[ip]['Values'] = OrderedDict(zip(tkeyList,map(BootAvgStdToFormat,pdata,[frmtflag]*len(pdata))))
        datadict[ip]['Boots'] = OrderedDict()
        for itstr,tdata in zip(tkeyList,pdata):
            datadict[ip]['Boots'][itstr] = tdata.values
        WriteXmlOutput(outputfile,datadict)

# data = [ ip , icut ]
def PrintFitMassToFile(data,dataChi,iset,filedir,filename,thisMomList,FitRanges,mominfo2pt):
    xmlMomList = map(qstrTOqcond,thisMomList)
    xmlCutList = map(xmlcut,thisCutList)
    xmlFitRanges = map(xmlfitr,FitRanges)
    for ipc,(ip,pdata,pdataChi) in enumerate(zip(xmlMomList,data,dataChi)):        
        datadict,outputfile = SetUpPDict(ip,filedir,filename)
        datadict[ip]['Info'] = mominfo2pt[ipc]
        datadict[ip]['Values'] = OrderedDict()
        datadict[ip]['Boots'] = OrderedDict()
        for ifit,fitdata,fitdataChi in zip(xmlFitRanges,qdata,qdataChi):
            datadict[ip]['Values'][icutstr] = BootAvgStdChiToFormat(fitdata[iset],fitdataChi[iset])
            datadict[ip]['Boots'][icutstr] = cutdata[iset].values
        # MergeXmlOutput(outputfile,datadict)
        WriteXmlOutput(outputfile,datadict)
        
# data = [ ip , icut , iset ]
def PrintFitToFile(data,dataChi,iset,filedir,filename,thisMomList,thisCutList,mominfoRF):
    xmlMomList = map(qstrTOqcond,thisMomList)
    xmlCutList = map(xmlcut,thisCutList)
    for ipc,(ip,pdata,pdataChi) in enumerate(zip(xmlMomList,data,dataChi)):        
        datadict,outputfile = SetUpPDict(ip,filedir,filename)
        datadict[ip]['Info'] = mominfoRF[ipc]
        datadict[ip]['Values'] = OrderedDict()
        datadict[ip]['Boots'] = OrderedDict()
        for icutstr,cutdata,cutdataChi in zip(xmlCutList,pdata,pdataChi):
            datadict[ip]['Values'][icutstr] = BootAvgStdChiToFormat(cutdata[iset],cutdataChi[iset])
            datadict[ip]['Boots'][icutstr] = cutdata[iset].values
        # MergeXmlOutput(outputfile,datadict)
        WriteXmlOutput(outputfile,datadict)


def PrintLREvecMassToFile(thisLE,thisRE,thisEMass,thisMomList,thisTvar,AddDict={},DoPoF=True):
    xmlMomList = map(ipTOqcond,thisMomList)
    for ip,pLE,pRE,pEMass in zip(xmlMomList,thisLE,thisRE,thisEMass):
        mkdir_p(outputdir+'/Mass/')
        datadict,outputfile = SetUpPDict(ip,outputdir+'/Mass/',thisTvar+'LREM')
        datadict[ip]['Info'] = AddDict
        datadict[ip]['Values'] = OrderedDict()
        for istate,iLE,iRE,iEM in zip(GetStateSet(thisTvar),pLE,pRE,pEMass):
            datadict[ip]['Values']['State'+str(istate)] = LREVecToFormat(iLE,iRE,iEM,DoPoF)
        # MergeXmlOutput(outputfile,datadict)
        WriteXmlOutput(outputfile,datadict)



        
##data [ ip , icut , itsink ]
##datafit [ ip , icut , ifit , par ] bs1
def PrintSumToFile(data,datafit,datafitchi,filedir,filename,thisFitList,thisMomList,thisTSinkList,thisCutList,mominfoRF,frmtflag='f'):
    xmlMomList = map(qstrTOqcond,thisMomList)
    xmlCutList = map(xmlcut,thisCutList)
    xmlTSinkList = map(xmlTSink,thisTSinkList)
    for ipc,(ip,pdata,pdatafit,pdatafitchi,pfitlist) in enumerate(zip(xmlMomList,data,datafit,datafitchi,thisFitList)):
        xmlFitList = ParamsToFitFlag(pfitlist)
        datadict,outputfile = SetUpPDict(ip,filedir,filename)
        datadict[ip]['Info'] = CombineSetInfo(mominfoRF[ipc])
        datadict[ip]['Values'] = OrderedDict()
        datadict[ip]['Boots'] = OrderedDict()
        for icut,cutdata,cutdatafit,cutdatafitchi,cutfitlist in zip(xmlCutList,pdata,pdatafit,pdatafitchi,xmlFitList):
            datadict[ip]['Values'][icut] = OrderedDict()
            datadict[ip]['Values'][icut]['slope'] = OrderedDict()
            datadict[ip]['Values'][icut]['constant'] = OrderedDict()
            datadict[ip]['Boots'][icut] = OrderedDict()
            datadict[ip]['Boots'][icut]['slope'] = OrderedDict()
            datadict[ip]['Boots'][icut]['constant'] = OrderedDict()
            for itsink,tsinkdata in zip(xmlTSinkList,cutdata):
                datadict[ip]['Values'][icut][itsink] = BootAvgStdToFormat(tsinkdata)
                datadict[ip]['Boots'][icut][itsink] = tsinkdata.values
            for ifit,fitdata,fitdatachi in zip(cutfitlist,cutdatafit,cutdatafitchi):
                datadict[ip]['Values'][icut]['slope'][ifit] = BootAvgStdChiToFormat(fitdata[0],fitdatachi)
                datadict[ip]['Boots'][icut]['slope'][ifit] = fitdata[0].values
            for ifit,fitdata,fitdatachi in zip(cutfitlist,cutdatafit,cutdatafitchi):
                datadict[ip]['Values'][icut]['constant'][ifit] = BootAvgStdChiToFormat(fitdata[1],fitdatachi)
                datadict[ip]['Boots'][icut]['constant'][ifit] = fitdata[1].values
        WriteXmlOutput(outputfile,datadict)
        # MergeXmlOutput(outputfile,datadict)


def PrintFFSet(FFin,Set,Mass,SetMass,theCurr,infoFF):
    FFdir = outputdir +'/FormFactors/'+theCurr+'/'
    FFbootdir = FFdir + 'boots/'
    mkdir_p(FFbootdir)
    thisfile = FFdir +theCurr+Set
    datadict = {'Form_Factors':{'Values':OrderedDict(),'Boots':OrderedDict()}}
    if 'Chi' not in Mass.keys(): Mass['Chi'] = float('NaN')
    datadict['Form_Factors']['Values']['Mass'] = OrderedDict()
    datadict['Form_Factors']['Info'] = OrderedDict()
    datadict['Form_Factors']['Values']['Mass']['Set'] = SetMass
    datadict['Form_Factors']['Values']['Mass']['Avg'] = Mass['Avg']
    datadict['Form_Factors']['Values']['Mass']['Std'] = Mass['Std']
    datadict['Form_Factors']['Values']['Mass']['Chi'] = Mass['Chi']
    for iqsqrd,qdata in FFin.iteritems():
        if len(qdata.keys()) > 0:
            datadict['Form_Factors']['Info'][iqsqrd] = infoFF[iqsqrd]
            datadict['Form_Factors']['Values'][iqsqrd] = OrderedDict()
            datadict['Form_Factors']['Values'][iqsqrd]['Chi'] = qdata['Chi']
            for ic,iFF in enumerate(qdata['Boot']):
                datadict['Form_Factors']['Values'][iqsqrd]['FF'+str(ic+1)] = BootAvgStdToFormat(iFF)
    for iqsqrd,qdata in FFin.iteritems():
        if len(qdata.keys()) > 0:
            datadict['Form_Factors']['Boots'][iqsqrd] = OrderedDict()
            for ic,iFF in enumerate(qdata['Boot']):
                datadict['Form_Factors']['Boots'][iqsqrd]['FF'+str(ic+1)] = iFF.values
    # MergeXmlOutput(thisfile,datadict,CheckMom=False)
    WriteXmlOutput(thisfile,datadict)





def PickTF(iset,iA,setsize):
    return setsize*iA + iset


# #data2pt       = [ ifit2pt , ip , params ]
# #data2ptChi    = [ ifit2pt , ip ]

def PrintTSFMassToFile(data2pt,data2ptChi,thisSetList,thisFit2ptList,fileprefix,thisMomList,info2pt):
    thisTSinkList,thisSmList = GetTsinkSmLists(thisSetList)
    masspardir = outputdir + 'cfuns/twopt/TSF'+fileprefix+'/'
    xmlMomList = map(qstrTOqcond,thisMomList)
    xml2ptFitList = map(xmlfitr,thisFit2ptList)
    for im in [-2,-1]: #TwoStateParList m0 and dm
        filename = 'twopt'+TwoStateParList['C2'][im]
        for ipc,ip in enumerate(xmlMomList):        
            datadict,outputfile = SetUpPDict(ip,masspardir,filename)
            datadict[ip]['Info'] = CombineSetInfo(info2pt[ipc])
            datadict[ip]['Values'] = OrderedDict()
            datadict[ip]['Boots'] = OrderedDict()
            for icutstr,cutdata,cutdataChi in zip(xml2ptFitList,data2pt,data2ptChi):
                mcutdata = cutdata[ipc][im].exp(1)
                mcutdata.Stats()
                datadict[ip]['Values'][icutstr] = BootAvgStdChiToFormat(mcutdata,cutdataChi[ipc])
                datadict[ip]['Boots'][icutstr] = mcutdata.values
            WriteXmlOutput(outputfile,datadict)
            # MergeXmlOutput(outputfile,datadict)

    for iA,theA in enumerate(TwoStateParList['C2'][:-2]):
        for ism,thesm in enumerate(thisSmList):
            filename = thesm+'twopt'+ theA
            for ipc,ip in enumerate(xmlMomList):        
                datadict,outputfile = SetUpPDict(ip,masspardir,filename)
                datadict[ip]['Info'] = info2pt[ipc][ism]
                datadict[ip]['Values'] = OrderedDict()
                datadict[ip]['Boots'] = OrderedDict()
                for icutstr,cutdata,cutdataChi in zip(xml2ptFitList,data2pt,data2ptChi):
                    mcutdata = cutdata[ipc][PickTF(ism,iA,len(thisSmList))]
                    datadict[ip]['Values'][icutstr] = BootAvgStdChiToFormat(mcutdata,cutdataChi[ipc],frmtflag='e')
                    datadict[ip]['Boots'][icutstr] = mcutdata.values
                # MergeXmlOutput(outputfile,datadict)
                WriteXmlOutput(outputfile,datadict)
            

#data3pt       = [ ifit2pt , ip , igamma , istate , ifit3pt , params ] bs1
#data3ptChi    = [ ifit2pt , ip , igamma , istate , ifit3pt ]

def PrintTSFToFile(filedir,filename,thisMomList,xml2ptFitList,xmlTSFList,data3pt,data3ptChi,ipar,igamma,ism,mominfoRF):
    xmlMomList = map(qstrTOqcond,thisMomList)
    for ipc,ip in enumerate(xmlMomList):        
        datadict,outputfile = SetUpPDict(ip,filedir,filename)
        datadict[ip]['Info'] = CombineSetInfo(mominfoRF[igamma][ipc])
        datadict[ip]['Values'] = OrderedDict()
        datadict[ip]['Boots'] = OrderedDict()
        for ic2pt,icut2ptstr in enumerate(xml2ptFitList):
            datadict[ip]['Values'][icut2ptstr] = OrderedDict()
            datadict[ip]['Boots'][icut2ptstr] = OrderedDict()
            for icutstr,cutdata,cutdataChi in zip(xmlTSFList,data3pt[ic2pt][igamma][ipc][ism],data3ptChi[ic2pt][igamma][ipc][ism]):
                datadict[ip]['Values'][icut2ptstr][icutstr] = BootAvgStdChiToFormat(cutdata[ipar],cutdataChi)
                datadict[ip]['Boots'][icut2ptstr][icutstr] = cutdata[ipar].values
        # MergeXmlOutput(outputfile,datadict)
        WriteXmlOutput(outputfile,datadict)
    
                
                
#data3pt       = [ ifit2pt , ip , igamma , istate , ifit3pt , params ] bs1
#data3ptChi    = [ ifit2pt , ip , igamma , istate , ifit3pt ]
def PrintTSFSetToFile(data3pt,data3ptChi,thisGammaMomList,thisSetList,thisFit2ptList,fileprefix,infoRF):
    thisTSinkList,thisSmList = GetTsinkSmLists(thisSetList)
    xml2ptFitList = map(xmlfitr,thisFit2ptList)
    xmlTSFList = map(xmlcut,TSF3ptCutList)
    del thisGammaMomList['twopt']
    for ipar,thispar in enumerate(TwoStateParList['C3']):
        for igamma,(thisgamma,thismomlist) in enumerate(thisGammaMomList.iteritems()):
            mprint('Printing ' , thispar , ' ' , thisgamma , ' to file      \r',)
            gammapardir = outputdir+CreateOppDir(thisgamma)+'/TSF'+fileprefix+'/'
            mkdir_p(gammapardir)
            for ism,thesm in enumerate(thisSmList):                
                PrintTSFToFile(gammapardir,thesm+thisgamma+thispar,thismomlist,xml2ptFitList,xmlTSFList,data3pt,data3ptChi,ipar,igamma,ism,infoRF)


                
#OneFit2pt    = [ ifit2pt , ip , ism  , params ] bs1
#OneFit2ptChi    = [ ifit2pt , ip , ism ]
def PrintOSFMassToFile(data2pt,data2ptChi,thisSetList,thisFit2ptList,fileprefix,thisMomList,info2pt):
    thisTSinkList,thisSmList = GetTsinkSmLists(thisSetList)
    masspardir = outputdir + 'cfuns/twopt/OSF'+fileprefix+'/'
    for im in [1,0]:
        for ism,thesm in enumerate(thisSmList):
            filename = thesm+'twopt'+OneStateParList['C2'][im]
            xmlMomList = map(qstrTOqcond,thisMomList)
            xml2ptFitList = map(xmlfitr,thisFit2ptList)
            for ipc,ip in enumerate(xmlMomList):        
                datadict,outputfile = SetUpPDict(ip,masspardir,filename)
                datadict[ip]['Info'] = info2pt[ipc][ism]
                datadict[ip]['Values'] = OrderedDict()
                datadict[ip]['Boots'] = OrderedDict()
                for icutstr,cutdata,cutdataChi in zip(xml2ptFitList,data2pt,data2ptChi):
                    if im == 1:
                        mcutdata = cutdata[ipc][ism][im].exp(1)
                        mcutdata.Stats()
                        thisformat = 'f'
                    else:
                        mcutdata = cutdata[ipc][ism][im]
                        thisformat = 'e'
                    datadict[ip]['Values'][icutstr] = BootAvgStdChiToFormat(mcutdata,cutdataChi[ipc][ism],frmtflag=thisformat)
                    datadict[ip]['Boots'][icutstr] = mcutdata.values
                WriteXmlOutput(outputfile,datadict)
                # MergeXmlOutput(outputfile,datadict)


#OneFit3pt    = [ ifit2pt , igamma , ip , iset , ifit3pt , params ] bs1

def PrintOSFToFile(filedir,filename,thisMomList,xml2ptFitList,xmlOSFList,data3pt,data3ptChi,ipar,igamma,ism,infoRF):
    xmlMomList = map(qstrTOqcond,thisMomList)
    for ipc,ip in enumerate(xmlMomList):        
        datadict,outputfile = SetUpPDict(ip,filedir,filename)
        datadict[ip]['Info'] = infoRF[igamma][ipc][ism]
        datadict[ip]['Values'] = OrderedDict()
        datadict[ip]['Boots'] = OrderedDict()
        for ic2pt,icut2ptstr in enumerate(xml2ptFitList):
            datadict[ip]['Values'][icut2ptstr] = OrderedDict()
            datadict[ip]['Boots'][icut2ptstr] = OrderedDict()
            for icutstr,cutdata,cutdataChi in zip(xmlOSFList,data3pt[ic2pt][igamma][ipc][ism],data3ptChi[ic2pt][igamma][ipc][ism]):
                datadict[ip]['Values'][icut2ptstr][icutstr] = BootAvgStdChiToFormat(cutdata[ipar],cutdataChi)
                datadict[ip]['Boots'][icut2ptstr][icutstr] = cutdata[ipar].values
        # MergeXmlOutput(outputfile,datadict)
        WriteXmlOutput(outputfile,datadict)
        

def PrintOSFSetToFile(data3pt,data3ptChi,thisGammaMomList,thisSetList,thisFit2ptList,fileprefix,infoRF):
    xml2ptFitList = map(xmlfitr,thisFit2ptList)
    xmlOSFList = map(xmlcut,OSF3ptCutList)
    del thisGammaMomList['twopt']
    for ipar,thispar in enumerate(OneStateParList['C3']):
        for igamma,(thisgamma,thismomlist) in enumerate(thisGammaMomList.iteritems()):
            mprint('Printing ' , thispar , ' ' , thisgamma , ' to file      \r',)
            gammapardir = outputdir+CreateOppDir(thisgamma)+'/OSF'+fileprefix+'/'
            mkdir_p(gammapardir)
            for ism,thesm in enumerate(thisSetList):
                PrintOSFToFile(gammapardir,thesm+thisgamma+thispar,thismomlist,xml2ptFitList,xmlOSFList,data3pt,data3ptChi,ipar,igamma,ism,infoRF)









# #  LocalWords:  thisSetList

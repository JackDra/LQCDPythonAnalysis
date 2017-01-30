#!/usr/bin/env python

from Params import *
from FitParams import *
from FitFunctions import DPfitfunDer,DPfitfun2Der
from BootTest import BootStrap1
from MiscFuns import *
from SetLists import GetTsinkSmLists
from OppFuns import CreateOppDir
import xmltodict
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
        f.write( xmltodict.unparse(outputdict,pretty=True).replace('\t','    '))

def WriteChromaXml(thisfile,outputdict):
    if Debug: print 'Writing to: ' , thisfile
    with open(thisfile+'.xml','w') as f:
        output = xmltodict.unparse(outputdict,pretty=True).replace('\t','    ')
        output = re.sub('elem.>','elem>',output)
        output = re.sub('elem..>','elem>',output)
        output = re.sub('elem...>','elem>',output)
        output = re.sub('elem....>','elem>',output)
        output = re.sub('elem.....>','elem>',output)
        output = re.sub('elem......>','elem>',output)
        f.write(output )

        
def WriteXmlOutput(thisfile,outputdict):
    firstkey = outputdict.keys()[0]
    Vals = {firstkey:{'Values':outputdict[firstkey]['Values']}}
    if 'Info' in outputdict[firstkey].keys(): Vals[firstkey]['Info'] = outputdict[firstkey]['Info']
    Boots = outputdict[firstkey]['Boots']
    outdirlist = thisfile.split('/')
    thisdir = '/'.join(outdirlist[:-1])
    mkdir_p(thisdir)
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
    Is2pt = False
    if 'Mass' in filedir or 'twopt' in filedir: Is2pt = True
    xmlMomList = [ipTOqcond(imom,Avg=Is2pt) for imom in thisMomList]
    tkeyList = map(tstr,thisTList)
    for ip,pdata in zip(xmlMomList,thisdata):
        datadict,outputfile = SetUpPDict(ip,filedir,filename)
        datadict[ip]['Info'] = AddDict
        datadict[ip]['Values'] = OrderedDict(zip(tkeyList,map(BootAvgStdToFormat,pdata,[frmtflag]*len(pdata))))
        datadict[ip]['Boots'] = OrderedDict()
        for itstr,tdata in zip(tkeyList,pdata):
            datadict[ip]['Boots'][itstr] = tdata.values
        WriteXmlOutput(outputfile,datadict)







##pdata [  ip , iflow , it ] bs1
def Print3ptTopToFile(pdata,filedir,filename,thisTopList,thisTList,thisMomList,AddDict={},frmtflag='f'):
    xmlMomList = map(ipTOqcond,thisMomList)
    tkeyList = map(tstr,thisTList)
    for icp,(ip,topdata) in enumerate(zip(xmlMomList,pdata)):
        datadictTop,outputfileTop = SetUpPDict(ip,filedir,filename)
        datadictTop[ip]['Info'] = AddDict
        datadictTop[ip]['Values'] = OrderedDict()
        datadictTop[ip]['Boots'] = OrderedDict()
        for itflow,flowdata in zip(thisTopList,topdata):
            Top = []
            datadictTop[ip]['Boots'][tflowstr(itflow)] = OrderedDict()
            for itstr,tflowdata in zip(tkeyList,flowdata):
                Top.append(tflowdata)
                datadictTop[ip]['Boots'][tflowstr(itflow)][itstr] = tflowdata.values
            datadictTop[ip]['Values'][tflowstr(itflow)] = OrderedDict(zip(tkeyList,map(BootAvgStdToFormat,Top,[frmtflag]*len(Top))))
        WriteXmlOutput(outputfileTop,datadictTop)


        













        
def PrintTopToFile(topdata,thisdata,filedir,filename,thisTopList,thisTList,thisMomList,AddDict={},frmtflag='f'):
    xmlMomList = [ipTOqcond(imom,Avg=True) for imom in thisMomList]
    tkeyList = map(tstr,thisTList)
    for icp,(ip,pdata) in enumerate(zip(xmlMomList,thisdata)):
        datadictAlpha,outputfileAlpha = SetUpPDict(ip,filedir,filename)
        datadict,outputfile = SetUpPDict(ip,filedir.replace('Alpha/','cfun/NNQ/'),filename)
        datadict[ip]['Info'] = AddDict
        datadict[ip]['Values'] = OrderedDict()
        datadict[ip]['Boots'] = OrderedDict()
        datadictAlpha[ip]['Info'] = AddDict
        datadictAlpha[ip]['Values'] = OrderedDict()
        datadictAlpha[ip]['Boots'] = OrderedDict()
        for itflow,flowdata in zip(thisTopList,topdata):
            Alpha,NNQ = [],[]
            datadict[ip]['Boots'][tflowstr(itflow)] = OrderedDict()
            datadictAlpha[ip]['Boots'][tflowstr(itflow)] = OrderedDict()
            for itstr,tflowdata,tdata in zip(tkeyList,flowdata[icp],pdata):
                Alpha.append(tflowdata/tdata)
                NNQ.append(tflowdata)
                Alpha[-1].Stats()
                datadictAlpha[ip]['Boots'][tflowstr(itflow)][itstr] = Alpha[-1].values
                datadict[ip]['Boots'][tflowstr(itflow)][itstr] = NNQ[-1].values
            datadictAlpha[ip]['Values'][tflowstr(itflow)] = OrderedDict(zip(tkeyList,map(BootAvgStdToFormat,Alpha,[frmtflag]*len(pdata))))
            datadict[ip]['Values'][tflowstr(itflow)] = OrderedDict(zip(tkeyList,map(BootAvgStdToFormat,NNQ,['e']*len(pdata))))
        WriteXmlOutput(outputfileAlpha,datadictAlpha)
        WriteXmlOutput(outputfile,datadict)

## data = {ip , iflow , ifit , Boot/Avg/Chi }
def PrintAlphaFitFile(data,iset,filedir):
    for ip,pdata in data.iteritems():
        pcond = qstrTOqcond(ip)
        datadict,outputfile = SetUpPDict(pcond,filedir,iset)
        datadict[pcond]['Info'] = pdata['Info']
        datadict[pcond]['Boots'] = pdata['Boots']
        for (itflow,flowdata),chiflowdata in zip(pdata['Boots'].iteritems(),pdata['Chi'].itervalues()):
            datadict[pcond]['Values'][itflow] = OrderedDict()
            for (ifitr,fitdata),chifitdata in zip(flowdata.iteritems(),chiflowdata.itervalues()):
                datadict[pcond]['Values'][itflow][ifitr] = BootAvgStdChiToFormat(fitdata,chifitdata)
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
# OR
# data = [ ip , iflow , icut , iset ]
def PrintFitToFile(data,dataChi,iset,filedir,filename,thisMomList,thisCutList,mominfoRF,flowlist):
    xmlMomList = map(qstrTOqcond,thisMomList)
    xmlCutList = map(xmlcut,thisCutList)
    for ipc,(ip,pdata,pdataChi) in enumerate(zip(xmlMomList,data,dataChi)):        
        datadict,outputfile = SetUpPDict(ip,filedir,filename)
        datadict[ip]['Info'] = mominfoRF[ipc]
        if 'Top' in outputfile.replace('Top'+kappalist[0],''):
            for iflow,flowdata,flowdataChi in zip(flowlist,pdata,pdataChi):
                datadict[ip]['Boots'][iflow] = OrderedDict()
                datadict[ip]['Values'][iflow] = OrderedDict()
                for icutstr,cutdata,cutdataChi in zip(xmlCutList,flowdata,flowdataChi):
                    datadict[ip]['Values'][iflow][icutstr] = BootAvgStdChiToFormat(cutdata[iset],cutdataChi[iset])
                    datadict[ip]['Boots'][iflow][icutstr] = cutdata[iset].values
        else:   
            datadict[ip]['Values'] = OrderedDict()
            datadict[ip]['Boots'] = OrderedDict()
            for icutstr,cutdata,cutdataChi in zip(xmlCutList,pdata,pdataChi):
                datadict[ip]['Values'][icutstr] = BootAvgStdChiToFormat(cutdata[iset],cutdataChi[iset])
                datadict[ip]['Boots'][icutstr] = cutdata[iset].values
        # MergeXmlOutput(outputfile,datadict)
        WriteXmlOutput(outputfile,datadict)


def PrintLREvecMassToFile(thisLE,thisRE,thisEMass,thisMomList,thisTvar,AddDict={},DoPoF=True):
    xmlMomList = [ipTOqcond(imom,Avg=True) for imom in thisMomList]
    for ip,pLE,pRE,pEMass in zip(xmlMomList,thisLE,thisRE,thisEMass):
        mkdir_p(outputdir[0]+'/Mass/')
        datadict,outputfile = SetUpPDict(ip,outputdir[0]+'/Mass/',thisTvar+'LREM')
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
    FFdir = outputdir[0] +'/FormFactors/'+theCurr+'/'
    FFbootdir = FFdir + 'boots/'
    mkdir_p(FFbootdir)
    thisfile = FFdir +theCurr+Set
    datadict = {'Form_Factors':{'Values':OrderedDict(),'Boots':OrderedDict()}}
    if 'Chi' not in Mass.keys(): Mass['Chi'] = float('NaN')
    datadict['Form_Factors']['Values']['Mass'] = OrderedDict()
    datadict['Form_Factors']['Info'] = infoFF
    datadict['Form_Factors']['Values']['Mass']['Set'] = SetMass
    datadict['Form_Factors']['Values']['Mass']['Avg'] = Mass['Avg']
    datadict['Form_Factors']['Values']['Mass']['Std'] = Mass['Std']
    datadict['Form_Factors']['Values']['Mass']['Chi'] = Mass['Chi']
    for iqsqrd,qdata in FFin.iteritems():
        if len(qdata.keys()) > 0:
            # datadict['Form_Factors']['Info'][iqsqrd] = infoFF[iqsqrd]
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
    


## outputdict : {iSet : nFF : Boot/Avg/Chi : F(0)/mEM }
def PrintDPfit(iFF,outputdict,InfoDict):
    FFdir = outputdir[0] +'/FormFactors/DPfits/'
    FFbootdir = FFdir + 'boots/'
    mkdir_p(FFbootdir)
    thisfile = FFdir +iFF
    # datadict = {'DP_Fits':{'Values':OrderedDict(),'Boots':OrderedDict()}}
    datadict = {'DP_Fits':{'Values':OrderedDict(),'Boots':OrderedDict()}}
    datadict['DP_Fits']['Info'] = InfoDict
    for iFF in outputdict[outputdict.keys()[0]].iterkeys():
        datadict['DP_Fits']['Values'][iFF] = OrderedDict()
        datadict['DP_Fits']['Boots'][iFF] = OrderedDict()
        datadict['DP_Fits']['Values'][iFF]['Fzero'] = OrderedDict()
        datadict['DP_Fits']['Values'][iFF]['mEM'] = OrderedDict()
        # datadict['DP_Fits']['Values'][iFF]['zero_slope'] = OrderedDict()
        datadict['DP_Fits']['Values'][iFF]['Radius'] = OrderedDict() 
        datadict['DP_Fits']['Values'][iFF]['Chi2_pdf'] = OrderedDict() 
        for iSet,setdict in outputdict.iteritems():
            # if not CheckDict(outputdict,iSet,iFF,'Avg'): continue
            if not CheckDict(outputdict,iSet,iFF,'Boot'): continue
            datadict['DP_Fits']['Values'][iFF]['Fzero'][iSet] = BootAvgStdToFormat(outputdict[iSet][iFF]['Boot'][0])
            datadict['DP_Fits']['Values'][iFF]['mEM'][iSet] = BootAvgStdToFormat(outputdict[iSet][iFF]['Boot'][1])
            # datadict['DP_Fits']['Values'][iFF]['zero_slope'][iSet] = DPfitfun2Der(np.array([0.]),outputdict[iSet][iFF]['Avg'])[1]
            # datadict['DP_Fits']['Values'][iFF]['zero_slope'][iSet] = DPfitfunDer(np.array([0.]),outputdict[iSet][iFF]['Avg'])[1]
            thisrad = 12*hbarc**2 *outputdict[iSet][iFF]['Boot'][0]/(outputdict[iSet][iFF]['Boot'][1])
            # thisrad = -12 *outputdict[iSet][iFF]['Boot'][0]/outputdict[iSet][iFF]['Boot'][1]
            thisrad.Stats()
            datadict['DP_Fits']['Values'][iFF]['Radius'][iSet] = BootAvgStdToFormat(thisrad)
            if not CheckDict(outputdict,iSet,iFF,'Chi'): continue
            datadict['DP_Fits']['Values'][iFF]['Chi2_pdf'][iSet] = outputdict[iSet][iFF]['Chi'][0]
    # MergeXmlOutput(thisfile,datadict,CheckMom=False)
    WriteXmlOutput(thisfile,datadict)



def PickTF(iset,iA,setsize):
    return setsize*iA + iset


# #data2pt       = [ ifit2pt , ip , params ]
# #data2ptChi    = [ ifit2pt , ip ]

def PrintTSFMassToFile(data2pt,data2ptChi,thisSetList,thisFit2ptList,fileprefix,thisMomList,info2pt):
    thisTSinkList,thisSmList = GetTsinkSmLists(thisSetList)
    masspardir = outputdir[0] + 'cfun/twopt/TSF'+fileprefix+'/'
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
            gammapardir = outputdir[0]+CreateOppDir(thisgamma)+'/TSF'+fileprefix+'/'
            mkdir_p(gammapardir)
            for ism,thesm in enumerate(thisSmList):                
                PrintTSFToFile(gammapardir,thesm+thisgamma+thispar,thismomlist,xml2ptFitList,xmlTSFList,data3pt,data3ptChi,ipar,igamma,ism,infoRF)


                
#OneFit2pt    = [ ifit2pt , ip , ism  , params ] bs1
#OneFit2ptChi    = [ ifit2pt , ip , ism ]
def PrintOSFMassToFile(data2pt,data2ptChi,thisSetList,thisFit2ptList,fileprefix,thisMomList,info2pt):
    thisTSinkList,thisSmList = GetTsinkSmLists(thisSetList)
    masspardir = outputdir[0] + 'cfun/twopt/OSF'+fileprefix+'/'
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
            gammapardir = outputdir[0]+CreateOppDir(thisgamma)+'/OSF'+fileprefix+'/'
            mkdir_p(gammapardir)
            for ism,thesm in enumerate(thisSetList):
                PrintOSFToFile(gammapardir,thesm+thisgamma+thispar,thismomlist,xml2ptFitList,xmlOSFList,data3pt,data3ptChi,ipar,igamma,ism,infoRF)









# #  LocalWords:  thisSetList

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

def PrintToFile(thisdata,filename,thisTList,thisMomList,CalcFlag,frmtflag='f'):
    datadict = {CalcFlag:{'Values':OrderedDict(),'Boots':OrderedDict()}}
    xmlMomList = map(ipTOqcond,thisMomList)
    tkeyList = map(tstr,thisTList)
    for ip,pdata in zip(xmlMomList,thisdata):
        datadict[CalcFlag]['Values'][ip] = OrderedDict(zip(tkeyList,map(BootAvgStdToFormat,pdata,[frmtflag]*len(pdata))))
    for ip,pdata in zip(xmlMomList,thisdata):
        datadict[CalcFlag]['Boots'][ip] = OrderedDict()
        for itstr,tdata in zip(tkeyList,pdata):
            datadict[CalcFlag]['Boots'][ip][itstr] = tdata.values
    with open(filename+'.xml','w') as f:
        f.write( xmltodict.unparse(datadict,pretty=True))



# data = [ ip , icut , iset ]
def PrintFitToFile(data,dataChi,iset,filename,thisMomList,thisCutList):
    datadict = {'Fits':{'Values':OrderedDict(),'Boots':OrderedDict()}}
    xmlMomList = map(qstrTOqcond,thisMomList)
    xmlCutList = map(xmlcut,thisCutList)
    for ip,pdata,pdataChi in zip(xmlMomList,data,dataChi):        
        datadict['Fits']['Values'][ip] = OrderedDict()
        for icutstr,cutdata,cutdataChi in zip(xmlCutList,pdata,pdataChi):
            datadict['Fits']['Values'][ip][icutstr] = BootAvgStdChiToFormat(cutdata[iset],cutdataChi[iset])
    for ip,pdata in zip(xmlMomList,data):        
        datadict['Fits']['Boots'][ip] = OrderedDict()
        for icutstr,cutdata in zip(xmlCutList,pdata):
            datadict['Fits']['Boots'][ip][icutstr] = cutdata[iset].values
    with open(filename+'.xml','w') as f:
        f.write( xmltodict.unparse(datadict,pretty=True))


def PrintLREvecMassToFile(thisLE,thisRE,thisEMass,thisMomList,thisTvar,DoPoF=True):
    datadict = {'Eigen_solutions':{'Values':OrderedDict()}}
    xmlMomList = map(ipTOqcond,thisMomList)
    mkdir_p(outputdir+'/Mass/')
    for ip,pLE,pRE,pEMass in zip(xmlMomList,thisLE,thisRE,thisEMass):
        datadict['Eigen_solutions']['Values'][ip] = OrderedDict()
        for istate,iLE,iRE,iEM in zip(StateSet,pLE,pRE,pEMass):
            datadict['Eigen_solutions']['Values'][ip]['State'+str(istate)] = LREVecToFormat(iLE,iRE,iEM,DoPoF)
    with open(outputdir+'/Mass/'+thisTvar+'LREM.xml','w') as f:
        f.write( xmltodict.unparse(datadict,pretty=True))



        
##data [ ip , icut , itsink ]
##datafit [ ip , icut , ifit , par ] bs1
def PrintSumToFile(data,datafit,datafitchi,filename,thisFitList,thisMomList,thisTSinkList,thisCutList,frmtflag='f'):
    datadict = {'Sum':{'Values':OrderedDict(),'Boots':OrderedDict()}}
    xmlMomList = map(qstrTOqcond,thisMomList)
    xmlCutList = map(xmlcut,thisCutList)
    xmlTSinkList = map(xmlTSink,thisTSinkList)
    for ip,pdata,pdatafit,pdatafitchi,pfitlist in zip(xmlMomList,data,datafit,datafitchi,thisFitList):
        xmlFitList = ParamsToFitFlag(pfitlist)
        datadict['Sum']['Values'][ip] = OrderedDict()
        for icut,cutdata,cutdatafit,cutdatafitchi,cutfitlist in zip(xmlCutList,pdata,pdatafit,pdatafitchi,xmlFitList):
            datadict['Sum']['Values'][ip][icut] = OrderedDict()
            datadict['Sum']['Values'][ip][icut]['slope'] = OrderedDict()
            datadict['Sum']['Values'][ip][icut]['constant'] = OrderedDict()
            for itsink,tsinkdata in zip(xmlTSinkList,cutdata):
                datadict['Sum']['Values'][ip][icut][itsink] = BootAvgStdToFormat(tsinkdata)
            for ifit,fitdata,fitdatachi in zip(cutfitlist,cutdatafit,cutdatafitchi):
                datadict['Sum']['Values'][ip][icut]['slope'][ifit] = BootAvgStdChiToFormat(fitdata[0],fitdatachi)
            for ifit,fitdata,fitdatachi in zip(cutfitlist,cutdatafit,cutdatafitchi):
                datadict['Sum']['Values'][ip][icut]['constant'][ifit] = BootAvgStdChiToFormat(fitdata[1],fitdatachi)

    for ip,pdata,pdatafit,pdatafitchi,pfitlist in zip(xmlMomList,data,datafit,datafitchi,thisFitList):
        xmlFitList = ParamsToFitFlag(pfitlist)
        datadict['Sum']['Boots'][ip] = OrderedDict()
        for icut,cutdata,cutdatafit,cutdatafitchi,cutfitlist in zip(xmlCutList,pdata,pdatafit,pdatafitchi,xmlFitList):
            datadict['Sum']['Boots'][ip][icut] = OrderedDict()
            datadict['Sum']['Boots'][ip][icut]['slope'] = OrderedDict()
            datadict['Sum']['Boots'][ip][icut]['constant'] = OrderedDict()
            for itsink,tsinkdata in zip(xmlTSinkList,cutdata):
                datadict['Sum']['Boots'][ip][icut][itsink] = tsinkdata.values
            for ifit,fitdata,fitdatachi in zip(cutfitlist,cutdatafit,cutdatafitchi):
                datadict['Sum']['Boots'][ip][icut]['slope'][ifit] = fitdata[0].values
            for ifit,fitdata,fitdatachi in zip(cutfitlist,cutdatafit,cutdatafitchi):
                datadict['Sum']['Boots'][ip][icut]['constant'][ifit] = fitdata[1].values

    with open(filename+'.xml','w') as f:
        f.write( xmltodict.unparse(datadict,pretty=True))


def PrintFFSet(FFin,Set,Mass,SetMass,theCurr):
    FFdir = outputdir +'/FormFactors/'+theCurr+'/'
    FFbootdir = FFdir + 'boots/'
    mkdir_p(FFbootdir)
    thisfile = FFdir +theCurr+Set
    datadict = {'Form_Factors':{'Values':OrderedDict(),'Boots':OrderedDict()}}
    if 'Chi' not in Mass.keys(): Mass['Chi'] = float('NaN')
    datadict['Form_Factors']['Values']['Mass'] = OrderedDict()
    datadict['Form_Factors']['Values']['Mass']['Set'] = SetMass
    datadict['Form_Factors']['Values']['Mass']['Avg'] = Mass['Avg']
    datadict['Form_Factors']['Values']['Mass']['Std'] = Mass['Std']
    datadict['Form_Factors']['Values']['Mass']['Chi'] = Mass['Chi']
    for iqsqrd,qdata in FFin.iteritems():
        if len(qdata.keys()) > 0:
            datadict['Form_Factors']['Values'][iqsqrd] = OrderedDict()
            datadict['Form_Factors']['Values'][iqsqrd]['Chi'] = qdata['Chi']
            for ic,iFF in enumerate(qdata['Boot']):
                datadict['Form_Factors']['Values'][iqsqrd]['FF'+str(ic)] = BootAvgStdToFormat(iFF)
    for iqsqrd,qdata in FFin.iteritems():
        if len(qdata.keys()) > 0:
            datadict['Form_Factors']['Boots'][iqsqrd] = OrderedDict()
            for ic,iFF in enumerate(qdata['Boot']):
                datadict['Form_Factors']['Boots'][iqsqrd]['FF'+str(ic)] = iFF.values
    with open(thisfile+'.xml','w') as f:
        f.write( xmltodict.unparse(datadict,pretty=True))





def PickTF(iset,iA,setsize):
    return setsize*iA + iset


# #data2pt       = [ ifit2pt , ip , params ]
# #data2ptChi    = [ ifit2pt , ip ]

def PrintTSFMassToFile(data2pt,data2ptChi,thisSetList,thisFit2ptList,fileprefix,thisMomList):
    thisTSinkList,thisSmList = GetTsinkSmLists(thisSetList)
    masspardir = outputdir + 'cfuns/twopt/TSF'+fileprefix+'/'
    for im in [-2,-1]: #TwoStateParList m0 and dm
        filename = masspardir+'twopt'+TwoStateParList['C2'][im]
        datadict = {'TSFMass':{'Values':OrderedDict(),'Boots':OrderedDict()}}
        xmlMomList = map(qstrTOqcond,thisMomList)
        xml2ptFitList = map(xmlfitr,thisFit2ptList)
        for ipc,ip in enumerate(xmlMomList):        
            datadict['TSFMass']['Values'][ip] = OrderedDict()
            for icutstr,cutdata,cutdataChi in zip(xml2ptFitList,data2pt,data2ptChi):
                mcutdata = cutdata[ipc][im].exp(1)
                mcutdata.Stats()
                datadict['TSFMass']['Values'][ip][icutstr] = BootAvgStdChiToFormat(mcutdata,cutdataChi[ipc])
        for ipc,ip in enumerate(xmlMomList):        
            datadict['TSFMass']['Boots'][ip] = OrderedDict()
            for icutstr,cutdata in zip(xml2ptFitList,data2pt):
                mcutdata = cutdata[ipc][im].exp(1)
                mcutdata.Stats()
                datadict['TSFMass']['Boots'][ip][icutstr] = mcutdata.values
        with open(filename+'.xml','w') as f:
            f.write( xmltodict.unparse(datadict,pretty=True))

    for iA,theA in enumerate(TwoStateParList['C2'][:-2]):
        for ism,thesm in enumerate(thisSmList):
            filename = masspardir +thesm+'twopt'+ theA
            datadict = {'TSFMass':{'Values':OrderedDict(),'Boots':OrderedDict()}}
            xmlMomList = map(qstrTOqcond,thisMomList)
            xml2ptFitList = map(xmlfitr,thisFit2ptList)
            for ipc,ip in enumerate(xmlMomList):        
                datadict['TSFMass']['Values'][ip] = OrderedDict()
                for icutstr,cutdata,cutdataChi in zip(xml2ptFitList,data2pt,data2ptChi):
                    mcutdata = cutdata[ipc][PickTF(ism,iA,len(thisSmList))]
                    datadict['TSFMass']['Values'][ip][icutstr] = BootAvgStdChiToFormat(mcutdata,cutdataChi[ipc],frmtflag='e')
            for ipc,ip in enumerate(xmlMomList):        
                datadict['TSFMass']['Boots'][ip] = OrderedDict()
                for icutstr,cutdata in zip(xml2ptFitList,data2pt):
                    mcutdata = cutdata[ipc][PickTF(ism,iA,len(thisSmList))]
                    datadict['TSFMass']['Boots'][ip][icutstr] = mcutdata.values
            with open(filename+'.xml','w') as f:
                f.write( xmltodict.unparse(datadict,pretty=True))


#data3pt       = [ ifit2pt , ip , igamma , istate , ifit3pt , params ] bs1
#data3ptChi    = [ ifit2pt , ip , igamma , istate , ifit3pt ]

def PrintTSFToFile(filename,thisMomList,xml2ptFitList,xmlTSFList,data3pt,data3ptChi,ipar,igamma,ism):
    datadict = {'TSF':{'Values':OrderedDict(),'Boots':OrderedDict()}}
    xmlMomList = map(qstrTOqcond,thisMomList)
    for ipc,ip in enumerate(xmlMomList):        
        datadict['TSF']['Values'][ip] = OrderedDict()
        for ic2pt,icut2ptstr in enumerate(xml2ptFitList):
            datadict['TSF']['Values'][ip][icut2ptstr] = OrderedDict()
            for icutstr,cutdata,cutdataChi in zip(xmlTSFList,data3pt[ic2pt][ipc][igamma][ism],data3ptChi[ic2pt][ipc][igamma][ism]):
                datadict['TSF']['Values'][ip][icut2ptstr][icutstr] = BootAvgStdChiToFormat(cutdata[ipar],cutdataChi)
    for ipc,ip in enumerate(xmlMomList):        
        datadict['TSF']['Boots'][ip] = OrderedDict()
        for ic2pt,icut2ptstr in enumerate(xml2ptFitList):
            datadict['TSF']['Boots'][ip][icut2ptstr] = OrderedDict()
            for icutstr,cutdata in zip(xmlTSFList,data3pt[ic2pt][ipc][igamma][ism]):
                datadict['TSF']['Boots'][ip][icut2ptstr][icutstr] = cutdata[ipar].values
    print filename+'.xml'
    with open(filename+'.xml','w') as f:
        f.write( xmltodict.unparse(datadict,pretty=True))
    
                
                
#data3pt       = [ ifit2pt , ip , igamma , istate , ifit3pt , params ] bs1
#data3ptChi    = [ ifit2pt , ip , igamma , istate , ifit3pt ]
def PrintTSFSetToFile(data3pt,data3ptChi,thisGammaMomList,thisSetList,thisFit2ptList,fileprefix):
    thisTSinkList,thisSmList = GetTsinkSmLists(thisSetList)
    xml2ptFitList = map(xmlfitr,thisFit2ptList)
    xmlTSFList = map(xmlcut,TSF3ptCutList)
    del thisGammaMomList['twopt']
    for ipar,thispar in enumerate(TwoStateParList['C3']):
        for igamma,(thisgamma,thismomlist) in enumerate(thisGammaMomList.iteritems()):
            print 'Printing ' , thispar , ' ' , thisgamma , ' to file      \r',
            gammapardir = outputdir+CreateOppDir(thisgamma)+'/TSF'+fileprefix+'/'
            for ism,thesm in enumerate(thisSmList):
                filename = gammapardir+thesm+thisgamma+ thispar
                PrintTSFToFile(filename,thismomlist,xml2ptFitList,xmlTSFList,data3pt,data3ptChi,ipar,igamma,ism)


                
#OneFit2pt    = [ ifit2pt , ip , ism  , params ] bs1
#OneFit2ptChi    = [ ifit2pt , ip , ism ]
def PrintOSFMassToFile(data2pt,data2ptChi,thisSetList,thisFit2ptList,fileprefix,thisMomList):
    thisTSinkList,thisSmList = GetTsinkSmLists(thisSetList)
    masspardir = outputdir + 'cfuns/twopt/OSF'+fileprefix+'/'
    for im in [1,0]:
        for ism,thesm in enumerate(thisSmList):
            filename = masspardir+thesm+'twopt'+OneStateParList['C2'][im]
            datadict = {'OSFMass':{'Values':OrderedDict(),'Boots':OrderedDict()}}
            xmlMomList = map(qstrTOqcond,thisMomList)
            xml2ptFitList = map(xmlfitr,thisFit2ptList)
            for ipc,ip in enumerate(xmlMomList):        
                datadict['OSFMass']['Values'][ip] = OrderedDict()
                for icutstr,cutdata,cutdataChi in zip(xml2ptFitList,data2pt,data2ptChi):
                    if im == 1:
                        mcutdata = cutdata[ipc][im].exp(1)
                        mcutdata.Stats()
                        thisformat = 'f'
                    else:
                        mcutdata = cutdata[ipc][im]
                        thisformat = 'e'
                    datadict['OSFMass']['Values'][ip][icutstr] = BootAvgStdChiToFormat(mcutdata,cutdataChi[ipc],frmtflag=thisformat)
            for ipc,ip in enumerate(xmlMomList):        
                datadict['OSFMass']['Boots'][ip] = OrderedDict()
                for icutstr,cutdata in zip(xml2ptFitList,data2pt):
                    if im == 1:
                        mcutdata = cutdata[ipc][im].exp(1)
                        mcutdata.Stats()
                    else:
                        mcutdata = cutdata[ipc][im]
                    datadict['OSFMass']['Boots'][ip][icutstr] = mcutdata.values
            with open(filename+'.xml','w') as f:
                f.write( xmltodict.unparse(datadict,pretty=True))


#OneFit3pt    = [ ifit2pt , igamma , ip , iset , ifit3pt , params ] bs1

def PrintOSFToFile(filename,thisMomList,xml2ptFitList,xmlOSFList,data3pt,data3ptChi,ipar,igamma,ism):
    datadict = {'OSF':{'Values':OrderedDict(),'Boots':OrderedDict()}}
    xmlMomList = map(qstrTOqcond,thisMomList)
    for ipc,ip in enumerate(xmlMomList):        
        datadict['OSF']['Values'][ip] = OrderedDict()
        for ic2pt,icut2ptstr in enumerate(xml2ptFitList):
            datadict['OSF']['Values'][ip][icut2ptstr] = OrderedDict()
            for icutstr,cutdata,cutdataChi in zip(xmlOSFList,data3pt[ic2pt][igamma][ipc][ism],data3ptChi[ic2pt][igamma][ipc][ism]):
                datadict['OSF']['Values'][ip][icut2ptstr][icutstr] = BootAvgStdChiToFormat(cutdata[ipar],cutdataChi)
    for ipc,ip in enumerate(xmlMomList):        
        datadict['OSF']['Boots'][ip] = OrderedDict()
        for ic2pt,icut2ptstr in enumerate(xml2ptFitList):
            datadict['OSF']['Boots'][ip][icut2ptstr] = OrderedDict()
            for icutstr,cutdata in zip(xmlOSFList,data3pt[ic2pt][igamma][ipc][ism]):
                datadict['OSF']['Boots'][ip][icut2ptstr][icutstr] = cutdata[ipar].values
    with open(filename+'.xml','w') as f:
        f.write( xmltodict.unparse(datadict,pretty=True))
        

def PrintOSFSetToFile(data3pt,data3ptChi,thisGammaMomList,thisSetList,thisFit2ptList,fileprefix):
    xml2ptFitList = map(xmlfitr,thisFit2ptList)
    xmlOSFList = map(xmlcut,OSF3ptCutList)
    del thisGammaMomList['twopt']
    for ipar,thispar in enumerate(OneStateParList['C3']):
        for igamma,(thisgamma,thismomlist) in enumerate(thisGammaMomList.iteritems()):
            print 'Printing ' , thispar , ' ' , thisgamma , ' to file      \r',
            gammapardir = outputdir+CreateOppDir(thisgamma)+'/OSF'+fileprefix+'/'
            for ism,thesm in enumerate(thisSmList):
                filename = gammapardir+thesm+thisgamma+ thispar
                PrintOSFToFile(filename,thismomlist,xml2ptFitList,xmlOSFList,data3pt,data3ptChi,ipar,igamma,ism)








# # # data = [ ip , icut ]
# # def PrintFitMassToFile(data,dataChi,filename,thisMomList,FitRanges):
# #     with open(filename+'.fit.txt','a+') as fb:
# #         for theq,qdata,qdataChi in zip(thisMomList,data,dataChi):
# #             fb.write('      '+theq+'\n')
# #             for ifit,fitdata,fitdataChi in zip(FitRanges,qdata,qdataChi):
# #                 # print ifitmin , ifitmax , fitdata.Avg , fitdata.Std , fitdataChi/float(fitlen)
# #                 fb.write( 'cut{0:4}{1:4}{2:20.10f}{3:20.10f}{4:20.10f}'.format(ifit[0],ifit[1],
# #                                                                         fitdata.Avg,fitdata.Std,
# #                                                                         fitdataChi)+ '\n' )

# # # data = [ ip , icut ]
# # def PrintFitMassBootToFile(data,filename,thisMomList,FitRanges):
# #     with open(filename+'.fit.boot.txt','a+') as fb:
# #         fb.write('         nboot '+str(nboot) + '\n')
# #         for theq,qdata in zip(thisMomList,data):
# #             fb.write('      '+theq+'\n')
# #             for ifit,fitdata in zip(FitRanges,qdata):
# #                 fb.write('   cut{0:4} {1:4}\n'.format(ifit[0],ifit[1]))
# #                 for iboot,bootdata in zip(range(fitdata.nboot),fitdata.values):
# #                     fb.write( '{0}{1:20.10f}'.format(repr(iboot).rjust(4), bootdata)+ '\n' )

# #  LocalWords:  thisSetList

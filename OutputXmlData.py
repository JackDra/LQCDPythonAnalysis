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
    with open(filename+'.xml','w') as f:
        f.write( xmltodict.unparse(datadict,pretty=True))





# def PickTF(iset,iA,setsize):
#     return setsize*iA + iset


# #data2pt       = [ ifit2pt , ip , params ]
# #data2ptChi    = [ ifit2pt , ip ]

# def PrintTSFMassToFile(data2pt,data2ptChi,thisSetList,thisFit2ptList,fileprefix,thisMomList):
#     thisTSinkList,thisSmList = GetTsinkSmLists(thisSetList)
#     masspardir = outputdir + 'cfuns/twopt/TSF'+fileprefix+'/'
#     bootdir = masspardir + 'boots/'
#     mkdir_p(bootdir)
#     for im in [-2,-1]: #TwoStateParList m0 and dm
#         filename = masspardir+'twopt'+TwoStateParList['C2'][im]
#         with open(filename+'.txt','a+') as f:
#             for imom,thismom in enumerate(thisMomList):
#                 f.write(thismom+'\n')
#                 for icfit2pt,ifit2pt in enumerate(thisFit2ptList):
#                     TSF2pt = data2pt[icfit2pt][imom][im]
#                     TSF2ptChi = data2ptChi[icfit2pt][imom]
#                     TSF2ptout = TSF2pt.exp(1)
#                     TSF2ptout.Stats()
#                     frmstr = '{3:3}{4:3} {0:20.10f} {1:20.10f} {2:20.10f}'
#                     f.write(frmstr.format(float(TSF2ptout.Avg),float(TSF2ptout.Std),float(TSF2ptChi),ifit2pt[0],ifit2pt[1])+'\n')
#         bootfn = bootdir+'twopt'+TwoStateParList['C2'][im]
#         with open(bootfn+'.boot.txt','a+') as f:
#             f.write(' nboot ' + str(nboot)+'\n')
#             for imom,thismom in enumerate(thisMomList):
#                 f.write(thismom+'\n')
#                 for icfit2pt,ifit2pt in enumerate(thisFit2ptList):
#                     f.write(' fitr {0:3}{1:3}'.format(ifit2pt[0],ifit2pt[1])+'\n')
#                     TSF2pt = data2pt[icfit2pt][imom][im]
#                     for iboot,bootdata in enumerate(TSF2pt.values):
#                         frmstr = '{0:3}{1:20.10f}'
#                         f.write(frmstr.format(iboot,float(np.exp(bootdata)))+'\n')
#     for iA,theA in enumerate(TwoStateParList['C2'][:-2]):
#         for ism,thesm in enumerate(thisSmList):
#             filename = masspardir +thesm+'twopt'+ theA
#             with open(filename+'.txt','a+') as f:
#                 for imom,thismom in enumerate(thisMomList):
#                     f.write(thismom+'\n')
#                     for icfit2pt,ifit2pt in enumerate(thisFit2ptList):
#                         TSF2pt = data2pt[icfit2pt][imom][PickTF(ism,iA,len(thisSmList))]
#                         TSF2ptChi = data2ptChi[icfit2pt][imom]
#                         frmstr = '{3:3}{4:3} {0:20.10e} {1:20.10e} {2:20.10e}'
#                         f.write(frmstr.format(float(TSF2pt.Avg),float(TSF2pt.Std),float(TSF2ptChi),ifit2pt[0],ifit2pt[1] )+'\n')
#             bootfn = bootdir+thesm+'twopt'+theA
#             with open(bootfn+'.boot.txt','a+') as f:
#                 f.write(' nboot ' + str(nboot)+'\n')
#                 for imom,thismom in enumerate(thisMomList):
#                     f.write(thismom+'\n')
#                     for icfit2pt,ifit2pt in enumerate(thisFit2ptList):
#                         f.write(' fitr {0:3}{1:3}'.format(ifit2pt[0],ifit2pt[1])+'\n')
#                         TSF2pt = data2pt[icfit2pt][imom][PickTF(ism,iA,len(thisSmList))]                        
#                         for iboot,bootdata in enumerate(TSF2pt.values):
#                             frmstr = '{0:3} {1:20.10e}'
#                             f.write(frmstr.format(iboot,float(bootdata))+'\n')

# #data3pt       = [ ifit2pt , ip , igamma , istate , ifit3pt , params ] bs1
# #data3ptChi    = [ ifit2pt , ip , igamma , istate , ifit3pt ]
# def PrintTSFToFile(data3pt,data3ptChi,thisGammaMomList,thisSetList,thisFit2ptList,fileprefix):
#     thisTSinkList,thisSmList = GetTsinkSmLists(thisSetList)
#     del thisGammaMomList['twopt']
#     for ipar,thispar in enumerate(TwoStateParList['C3']):
#         for igamma,(thisgamma,thismomlist) in enumerate(thisGammaMomList.iteritems()):
#             print 'Printing ' , thispar , ' ' , thisgamma , ' to file      \r',
#             gammapardir = outputdir+CreateOppDir(thisgamma)+'/TSF'+fileprefix+'/'
#             bootdir = gammapardir+'boots/'
#             mkdir_p(bootdir)
#             for ism,thesm in enumerate(thisSmList):
#                 filename = gammapardir+thesm+thisgamma+ thispar
#                 with open(filename+'.txt','a+') as f:
#                     for imom,thismom in enumerate(thismomlist):
#                         f.write(thismom+'\n')
#                         for icfit2pt,ifit2pt in enumerate(thisFit2ptList):
#                             f.write('\n')
#                             for icfit3pt,ifit3pt in enumerate(TSF3ptCutList):
#                                 fit3TF = data3pt[icfit2pt][igamma][imom][ism][icfit3pt][ipar]
#                                 fit3TFChi = data3ptChi[icfit2pt][igamma][imom][ism][icfit3pt]
#                                 f.write('{3:3}{4:3}{5:4} {0:20.10f} {1:20.10f} {2:20.10f}'
#                                         .format(fit3TF.Avg,fit3TF.Std,fit3TFChi,ifit2pt[0],ifit2pt[1],ifit3pt)+'\n')
#                 bootfn = bootdir+thesm+thisgamma+ thispar
#                 with open(bootfn+'.boot.txt','a+') as f:
#                     f.write(' nboot ' + str(nboot)+'\n')
#                     for imom,thismom in enumerate(thismomlist):
#                         f.write(thismom+'\n')
#                         for icfit2pt,ifit2pt in enumerate(thisFit2ptList):
#                             f.write('\n')
#                             for icfit3pt,ifit3pt in enumerate(TSF3ptCutList):
#                                 f.write('fitr {0:3}{1:3}{2:4}'.format(ifit2pt[0],ifit2pt[1],ifit3pt)+'\n')
#                                 fit3TF = data3pt[icfit2pt][igamma][imom][ism][icfit3pt][ipar]
#                                 for iboot,bootdata in enumerate(fit3TF.values):
#                                     frmstr = '{0:3} {1:20.10f}'
#                                     f.write(frmstr.format(iboot,float(bootdata))+'\n')
#     print '                                     '

# #OneFit2pt    = [ ifit2pt , ip , ism  , params ] bs1
# #OneFit2ptChi    = [ ifit2pt , ip , ism ]
# def PrintOSFMassToFile(data2pt,data2ptChi,thisSetList,thisFit2ptList,fileprefix,thisMomList):
#     thisTSinkList,thisSmList = GetTsinkSmLists(thisSetList)
#     masspardir = outputdir + 'cfuns/twopt/OSF'+fileprefix+'/'
#     bootdir = masspardir + 'boots/'
#     mkdir_p(bootdir)    
#     for im in [1,0]:
#         for ism,thesm in enumerate(thisSmList):
#             filename = masspardir+thesm+'twopt'+OneStateParList['C2'][im]
#             frmstr = '{3:3}{4:3} {0:20.10f} {1:20.10f} {2:20.10f}'
#             if im == 0: frmstr = '{3:3}{4:3} {0:20.10e} {1:20.10e} {2:20.10f}'
#             with open(filename+'.txt','a+') as f:
#                 for imom,thismom in enumerate(thisMomList):
#                     f.write(thismom+'\n')
#                     for icfit2pt,ifit2pt in enumerate(thisFit2ptList):
#                         OSF2pt = data2pt[icfit2pt][imom][ism][im]
#                         OSF2ptChi = data2ptChi[icfit2pt][imom][ism]
#                         if im == 1: OSF2pt = OSF2pt.exp(1)
#                         OSF2pt.Stats()
#                         f.write(frmstr.format(float(OSF2pt.Avg),float(OSF2pt.Std),float(OSF2ptChi),ifit2pt[0],ifit2pt[1])+'\n')
#             bootfn = bootdir+thesm+'twopt'+OneStateParList['C2'][im]
#             frmstr = '{0:3} {1:20.10f}'
#             if im == 0: frmstr = '{0:3} {1:20.10e}'
#             with open(bootfn+'.boot.txt','a+') as f:
#                 f.write(' nboot ' + str(nboot)+'\n')
#                 for imom,thismom in enumerate(thisMomList):
#                     f.write(thismom+'\n')
#                     for icfit2pt,ifit2pt in enumerate(thisFit2ptList):
#                         f.write(' fitr {0:3}{1:3}'.format(ifit2pt[0],ifit2pt[1])+'\n')
#                         OSF2pt = data2pt[icfit2pt][imom][ism][im]                    
#                         for iboot,bootdata in enumerate(OSF2pt.values):
#                             if im == 1: bootdata = np.exp(bootdata)
#                             f.write(frmstr.format(iboot,float(bootdata))+'\n')

#     # for theset,setTF,setTFChi in zip(thisSetList,OneFit2pt[iA],OneFit2ptChi):
#     #     filename = masspardir+theset+'twopt'+OneStateParList['C2'][iA]
#     #     with open(filename+'.txt','a+') as f:
#     #         for momTF,momTFChi,imom in zip(setTF,setTFChi,thisMomList):
#     #             f.write(imom+'\n')
#     #             for ifit2pt,OSF2pt,OSF2ptChi in zip(thisFit2ptList,momTF,momTFChi):
#     #                 f.write(frmstr.format(float(OSF2pt.Avg),float(OSF2pt.Std),float(OSF2ptChi),ifit2pt[0],ifit2pt[1])+'\n')
#     #     bootfn = bootdir+theset+'twopt'+OneStateParList['C2'][iA]
#     #     with open(bootfn+'.boot.txt','a+') as f:
#     #         f.write(' nboot ' + str(nboot)+'\n')
#     #         for momTF,imom in zip(setTF,thisMomList):
#     #             f.write(imom+'\n')
#     #             for ifit2pt,OSF2pt in zip(thisFit2ptList,momTF):
#     #                 f.write(' fitr {0:3}{1:3}'.format(ifit2pt[0],ifit2pt[1])+'\n')
#     #                 for iboot,bootdata in enumerate(OSF2pt.values):
#     #                     frmstr = '{0:3}{1:20.10e}'
#     #                     f.write(frmstr.format(iboot,float(bootdata))+'\n')


# #OneFit3pt    = [ ifit2pt , igamma , ip , iset , ifit3pt , params ] bs1
# def PrintOSFToFile(data3pt,data3ptChi,thisGammaMomList,thisSetList,thisFit2ptList,fileprefix):
#     del thisGammaMomList['twopt']
#     for ipar,thispar in enumerate(OneStateParList['C3']):
#         for igamma,(thisgamma,thismomlist) in enumerate(thisGammaMomList.iteritems()):
#             print 'Printing ' , thispar , ' ' , thisgamma , ' to file      \r',
#             gammapardir = outputdir+CreateOppDir(thisgamma)+'/OSF'+fileprefix+'/'
#             bootdir = gammapardir+'boots/'
#             mkdir_p(bootdir)
#             for iset,theset in enumerate(thisSetList):
#                 filename = gammapardir +theset+thisgamma+thispar
#                 with open(filename+'.txt','a+') as f:
#                     for imom,thismom in enumerate(thismomlist):
#                         f.write(thismom+'\n')
#                         for icfit2pt,ifit2pt in enumerate(thisFit2ptList):
#                             f.write('\n')
#                             for icfit3pt,ifit3pt in enumerate(OSF3ptCutList):
#                                 fit3TF = data3pt[icfit2pt][igamma][imom][iset][icfit3pt][ipar]
#                                 fit3TFChi = data3ptChi[icfit2pt][igamma][imom][iset][icfit3pt]
#                                 f.write('{3:3}{4:3}{5:4} {0:20.10f} {1:20.10f} {2:20.10f}'
#                                         .format(fit3TF.Avg,fit3TF.Std,fit3TFChi,ifit2pt[0],ifit2pt[1],ifit3pt)+'\n')
#                 bootfn = bootdir+theset+thisgamma+ thispar
#                 with open(bootfn+'.boot.txt','a+') as f:
#                     f.write(' nboot ' + str(nboot)+'\n')
#                     for imom,thismom in enumerate(thismomlist):
#                         f.write(thismom+'\n')
#                         for icfit2pt,ifit2pt in enumerate(thisFit2ptList):
#                             f.write('\n')
#                             for icfit3pt,ifit3pt in enumerate(OSF3ptCutList):
#                                 f.write('fitr {0:3}{1:3}{2:4}'.format(ifit2pt[0],ifit2pt[1],ifit3pt)+'\n')
#                                 fit3TF = data3pt[icfit2pt][igamma][imom][iset][icfit3pt][ipar]
#                                 for iboot,bootdata in enumerate(fit3TF.values):
#                                     frmstr = '{0:3} {1:20.10f}'
#                                     f.write(frmstr.format(iboot,float(bootdata))+'\n')
#     print '                                             '








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

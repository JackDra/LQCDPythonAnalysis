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
        datadict[CalcFlag]['Values'][ip] = OrderedDict(zip(tkeyList,map(BootAvgStdToFormat,pdata)))
    for ip,pdata in zip(xmlMomList,thisdata):
        datadict[CalcFlag]['Boots'][ip] = OrderedDict()
        for itstr,tdata in zip(tkeyList,pdata):
            datadict[CalcFlag]['Boots'][ip][itstr] = tdata.values
    with open(filename+'.xml','w') as f:
        f.write( xmltodict.unparse(datadict,pretty=True))



# data = [ ip , icut , iset ]
def PrintFitToFile(data,dataChi,iset,filename,thisMomList,thisCutList):
    datadict = {'Fits':{'Values':OrderedDict(),'Boots':OrderedDict()}}
    xmlMomList = map(ipTOqcond,thisMomList)
    xmlCutList = map(xmlcut,thisCutList)
    for ip,pdata,pdataChi in zip(xmlMomList,data,dataChi):        
        datadict['Fits']['Values'][ip] = OrderedDict()
        for icutstr,cutdata,cutdataChi in zip(xmlCutList,pdata,pdataChi):
            datadict['Fits']['Values'][ip][icutstr] = OrderedDict(zip(xmlCutList,map(BootAvgStdChiToFormat,cutdata,cutdataChi)))
    for ip,pdata in zip(xmlMomList,data):        
        datadict['Fits']['Boots'][ip] = OrderedDict()
        for icutstr,cutdata in zip(xmlCutList,pdata):
            datadict['Fits']['Boots'][ip][icutstr] = cutdata.values
    with open(filename+'.xml','w') as f:
        f.write( xmltodict.unparse(datadict,pretty=True))

    # with open(filename+'.xml','a+') as fb:
    #     # print thisMomList
    #     for iq,theq in enumerate(thisMomList):
    #         fb.write('      '+theq+'\n')
    #         for icut,fitdata,fitdataChi in zip(thisCutList,data[iq],dataChi[iq]):
    #             # print ifitmin , ifitmax , fitdata.Avg , fitdata.Std , fitdataChi/float(fitlen)
    #             fb.write( frmtstr.format(icut,fitdata[iset].Avg,fitdata[iset].Std,fitdataChi[iset])+ '\n' )
                
# # data = [ ip , icut , iset ]
# def PrintFitBootToFile(data,filename,iset,thisMomList,thisCutList):
#     with open(filename+'.boot.txt','a+') as fb:
#         fb.write('         nboot '+str(nboot) + '\n')
#         for iq,theq in enumerate(thisMomList):
#             fb.write('      '+theq+'\n')
#             for icut,fitdata in zip(thisCutList,data[iq]):
#                 fb.write('   cut{0}\n'.format(icut))
#                 for iboot,bootdata in zip(range(fitdata[iset].nboot),fitdata[iset].values):
#                     fb.write( '{0} {1:20.10f}'.format(repr(iboot).rjust(4), bootdata)+ '\n' )

# ##data [ ip , icut , itsink ]
# ##datafit [ ip , icut , ifit , par ] bs1
# def PrintSumToFile(data,datafit,datafitchi,filename,thisFitList,thisMomList,thisTSinkList,thisCutList,frmtflag='f'):
#     frmtstr = '    cut{0} tsink{3}: {1:20.10'+frmtflag+'} {2:20.10'+frmtflag+'}'
#     frmtfitstr = ' {0:20.10'+frmtflag+'} {1:20.10'+frmtflag+'} {2:20.10'+frmtflag+'}'
#     with open(filename+'.txt','a+') as f:
#         for theq,qdata,qdatafit,qdatafitchi,qfitlist in zip(thisMomList,data,datafit,datafitchi,thisFitList):
#             f.write(theq+'\n')
#             for icut,cutdata,cutdatafit,cutdatafitchi,cutfitlist in zip(thisCutList,qdata,qdatafit,qdatafitchi,qfitlist):
#                 for itsink,tsinkdata in zip(thisTSinkList,cutdata):
#                     f.write( frmtstr.format(icut, tsinkdata.Avg, tsinkdata.Std,itsink) + '\n' )
#                 f.write('\n')
#                 for ifit,fitdata,fitdatachi in zip(cutfitlist,cutdatafit,cutdatafitchi):
#                     f.write('fit cut{0} sl {1:2}{2:2}:'.format(icut,ifit[0],ifit[1]) +
#                             frmtfitstr.format(fitdata[0].Avg,fitdata[0].Std,fitdatachi)+'\n')
#                 for ifit,fitdata,fitdatachi in zip(cutfitlist,cutdatafit,cutdatafitchi):
#                     f.write('fit cut{0} con{1:2}{2:2}:'.format(icut,ifit[0],ifit[1]) +
#                             frmtfitstr.format(fitdata[1].Avg,fitdata[1].Std,fitdatachi)+'\n')
#                 f.write('\n')


# def PrintSumBootToFile(data,datafit,filename,thisFitList,thisMomList,thisTSinkList,thisCutList,frmtflag='f'):
#     frmtstr = '{0:3} {1:20.10'+frmtflag+'}'
#     with open(filename+'.boot.txt','a+') as fb:
#         fb.write('         nboot '+str(nboot) + '\n')
#         for theq,qdata,qdatafit,qfitlist in zip(thisMomList,data,datafit,thisFitList):
#             fb.write('      '+theq+'\n')
#             for icut,cutdata,cutdatafit,cutfitlist in zip(thisCutList,qdata,qdatafit,qfitlist):
#                 for itsink,tsinkdata in zip(thisTSinkList,cutdata):
#                     fb.write('    cut{0} tsink{1}:'.format(icut,itsink)+ '\n' )
#                     for iboot,bootdata in zip(range(nboot),tsinkdata.values):
#                         fb.write( frmtstr.format(iboot, bootdata)+ '\n' )
#                 fb.write('\n')
#                 for ifit,fitdata in zip(cutfitlist,cutdatafit):
#                     fb.write('fit cut{0} sl {1:2}{2:2}:'.format(icut,ifit[0],ifit[1])+'\n')
#                     for iboot,bootdata in enumerate(fitdata[0].values):
#                         fb.write(frmtstr.format(iboot,bootdata)+'\n')
#                 for ifit,fitdata in zip(cutfitlist,cutdatafit):
#                     fb.write('fit cut{0} con{1:2}{2:2}:'.format(icut,ifit[0],ifit[1])+'\n')
#                     for iboot,bootdata in enumerate(fitdata[1].values):
#                         fb.write(frmtstr.format(iboot,bootdata)+'\n')
#                 fb.write('\n')



# def PrintLREvecMassToFile(thisLE,thisRE,thisEMass,thisMomList,thisTvar,DoPoF=True):
#     mkdir_p(outputdir+'/Mass/')
#     with open(outputdir+'/Mass/'+thisTvar+'LREM.txt','a+') as f:
#         for ip,pLE,pRE,pEMass in zip(thisMomList,thisLE,thisRE,thisEMass):
#             f.write(ipTOqstr(ip)+'\n')
#             if DoPoF:
#                 f.write('         '+''.join([' {0:>20}'.format('sm'+DefSmearList[i]) for i in range(len(pLE[0])//(PoFShifts+1))])+' {0:>20}'.format('E-Mass')+'\n')
#                 for istate,iLE,iRE,iEM in zip(StateSet,pLE,pRE,pEMass):
#                     for iPoF in range(PoFShifts+1):
#                         thisnsmears = len(iLE)//(PoFShifts+1)
#                         f.write( 'L/R'+istate+' PoF'+str(iPoF)+''.join(' {0:20.10f}'.format(k.real) for k in iLE.tolist()[iPoF*thisnsmears:(iPoF+1)*thisnsmears])+
#                                  ' {0:20.10f}'.format(iEM) + '\n' )
#                     f.write('\n')
#             else:
#                 f.write('     '+''.join([' {0:>20}'.format('sm'+DefSmearList[i]) for i in range(len(pLE[0]))])+' {0:>20}'.format('E-Mass')+'\n')
#                 for istate,iLE,iRE,iEM in zip(StateSet,pLE,pRE,pEMass):
#                     f.write( 'L/R'+istate+' '+''.join(' {0:20.10f}'.format(k.real) for k in iLE.tolist())+' {0:20.10f}'.format(iEM) + '\n' )

#                 # for iPoF in range(PoFShifts+1):
#                 #     thisnsmears = len(iRE)//(PoFShifts+1)
#                 #     f.write( 'R'+istate + ' PoF'+str(iPoF)+' '+ ''.join(' {0:20.10f}'.format(k.real) for k in iRE.tolist()[iPoF:iPoF+thisnsmears])+
#                 #              ' {0:20.10f}'.format(iEM) + '\n' )

# ##C3set [ igamma , iset , it ] bs1

# def PrintCfunToFile(C3set,thisSetList,thisMomList, thisGammaList):
#     cfundir = outputdir + 'cfuns/'
#     for thegamma,gammadata in zip(thisGammaList,C3set):
#         gammadir = cfundir+CreateOppDir(thegamma)+'/'
#         bootgammadir = gammadir + 'boots/'
#         mkdir_p(bootgammadir)
#         for iset,setdata in zip(thisSetList,gammadata):
#             print 'Printing : ' , thegamma , iset , '                \r',
#             filename = (gammadir + iset+thegamma)
#             PrintToFile(np.array(setdata),filename,range(64),thisMomList,frmtflag='e')
#             filename = (bootgammadir + iset+thegamma)
#             PrintBootToFile(np.array(setdata),filename,range(64),thisMomList,frmtflag='e')

# ##dataset [ igamma , iset , ip , it ] bs1

# def PrintSetToFile(dataset,thisSetList,thisMomList, thisGammaList,tsink):
#     for thegamma,gammadata in zip(thisGammaList,dataset):
#         gammadir = outputdir+CreateOppDir(thegamma)+'/'
#         bootgammadir = gammadir + 'boots/'
#         mkdir_p(bootgammadir)
#         for iset,setdata in zip(thisSetList,gammadata):
#             print 'Printing : ' , thegamma , iset , '                \r',
#             if thegamma == 'Mass':
#                 setdata = cfunTOmass(setdata)
#                 tlist = range(64)
#             else:
#                 tlist = range(tsource,int(tsink)+1)
#             filename = (gammadir +iset+thegamma)
#             PrintToFile(setdata,filename,tlist,thisMomList)
#             filename = (bootgammadir +iset+thegamma)
#             PrintBootToFile(setdata,filename,tlist,thisMomList)


# ##sumdata [ igamma , ip , icut , itsink ] bs1
# ##sumfits [ igamma , ip , icut , fitvals , par ] bs1
# ##sumfitsChi [ igamma , ip , icut , fitvals ]

# def PrintSumSetToFile(sumdata,sumfits,sumfitschi,thisFitList,thissm, thisGammaMomList,thisTSinkList,thisCutList):
#     for (thegamma,thisMomList),gammadata,gammafitdata,gfdchi,gfitlist in zip(thisGammaMomList.iteritems(),sumdata,sumfits,sumfitschi,thisFitList):
#         print 'Printing : ' , thegamma , '                \r',
#         gammadir = outputdir+CreateOppDir(thegamma)+'/SumMeth/'
#         bootgammadir = gammadir + 'boots/'
#         mkdir_p(bootgammadir)
#         filename = gammadir +thissm+thegamma
#         PrintSumToFile(gammadata,gammafitdata,gfdchi,filename,gfitlist,thisMomList,thisTSinkList,thisCutList)
#         filename = bootgammadir +thissm+thegamma
#         PrintSumBootToFile(gammadata,gammafitdata,filename,gfitlist,thisMomList,thisTSinkList,thisCutList)



# #FitData = [ igamma , ip , icut , iset ]

# def PrintFitSetToFile(dataset,datasetChi,thisGammaMomList,thisSetList,thisCutList):
#     for igamma,(thisgamma,thismomlist) in enumerate(thisGammaMomList.iteritems()):
#         gammadir = outputdir+CreateOppDir(thisgamma)+'/Fits/'
#         bootgammadir = gammadir + 'boots/'
#         mkdir_p(bootgammadir)
#         for iset,thisset in enumerate(thisSetList):
#             print 'Printing : ' , thisgamma , thisset , '                \r',
#             filename = gammadir +thisset+thisgamma
#             bootfilename = bootgammadir +thisset+thisgamma
#             PrintFitToFile(dataset[igamma],datasetChi[igamma],iset,filename,thismomlist,thisCutList)
#             PrintFitBootToFile(dataset[igamma],bootfilename,iset,thismomlist,thisCutList)


# #dataset    = [ cuts , ip , istate ] bs1
# #datasetChi = [ cuts , ip , istate ]
# #dataset (after roll)   = [ istate , ip , cuts ] bs1
# #datasetChi (after roll)= [ istate , ip , cuts ]
# #statedata = [ ip , icut ]

# def PrintFitMassSetToFile(dataset,datasetChi,thisMomList,thisStateList,thisFitR):
#     dataset = np.rollaxis(np.rollaxis(dataset,1),2)
#     datasetChi = np.rollaxis(np.rollaxis(datasetChi,1),2)
#     gammadir = outputdir+'Mass/fits/'
#     bootgammadir = gammadir + 'boots/'
#     mkdir_p(gammadir)
#     mkdir_p(bootgammadir)
#     for thisstate,statedata,statedataChi in zip(thisStateList,dataset,datasetChi):
#         filename = (gammadir +thisstate+'Mass')
#         bootfilename = (bootgammadir +thisstate+'Mass')
#         PrintFitMassToFile(statedata,statedataChi,filename,thisMomList,thisFitR)
#         PrintFitMassBootToFile(statedata,bootfilename,thisMomList,thisFitR)

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

# def PrintFFSet(FFin,Set,Mass,SetMass,theCurr):
#     FFdir = outputdir +'/FormFactors/'+theCurr+'/'
#     FFbootdir = FFdir + 'boots/'
#     mkdir_p(FFbootdir)
#     thisfile = FFdir +theCurr+Set
#     thisfileboot = FFbootdir +theCurr+Set
#     with open(thisfile+'.txt','w') as f:
#         if 'Chi' not in Mass.keys(): Mass['Chi'] = float('NaN')
#         f.write('Mass {0} {1:20.10f} {2:20.10f} {3:20.10f}'
#                 .format(SetMass,Mass['Avg'],Mass['Std'],Mass['Chi'])+' \n')
#         f.write('\n')
#         f.write('{0:>6}'.format('q**2'))
#         for iff in range(NoFFPars[theCurr]):
#             f.write('{0:>20}{1:>20}'.format('FF'+str(iff+1),'FF'+str(iff+1)+'Err'))
#         f.write('{0:>20}'.format('Chi**2 PDF\n'))
#         for iqsqrd,qdata in FFin.iteritems():
#             if 'Boot' not in qdata.keys(): continue
#             f.write('{0:>6}'.format(iqsqrd))
#             for iFF in qdata['Boot']:
#                 f.write(' {0:20.10f} {1:20.10f}'.format(iFF.Avg,iFF.Std))
#             f.write(' {0:20.10f} \n'.format(qdata['Chi']))
#     with open(thisfileboot+'.boot.txt','w') as f:
#         if 'Chi' not in Mass.keys(): Mass['Chi'] = float('NaN')
#         f.write('nboot = ' + str(nboot)+ '\n')
#         f.write('Mass {0} {1:20.10f} {2:20.10f} {3:20.10f}'
#                 .format(SetMass,Mass['Avg'],Mass['Std'],Mass['Chi'])+' \n')
#         for iFF in range(NoFFPars[theCurr]):
#             f.write('{0:>12}'.format('FF'+str(iFF+1))+'\n')
#             for iqsqrd,qdata in FFin.iteritems():
#                 if 'Boot' not in qdata.keys(): continue
#                 f.write('{0:>6}'.format(iqsqrd)+'\n')
#                 for iboot,bootval in enumerate(qdata['Boot'][iFF].values):
#                     f.write('{0:3} {1:20.10f}'.format(iboot,bootval)+'\n')







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

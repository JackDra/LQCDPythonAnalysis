#!/usr/bin/env python

import matplotlib
matplotlib.use('PDF') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as pl
from Params import *
from FitFunctions import DPfitfun,DPfitfun2
from FitFunctions import DPfitfunOnePar
import numpy as np
from BootTest import BootStrap1
# import pylab as pl
from collections import OrderedDict
from ReadTxt import *
from MiscFuns import *
from FitParams import *
import FitFunctions as ff
from copy import deepcopy
from StringFix import *
import itertools
from FormFactors import NoFFPars
from SetLists import SortMySet
from CombParams import *
from FFParams import Qtcut,GetCharRad

##FORCE TITLE PARAMETER, SET TO FALSE TO USE NORMAL TITLES#

ForceTitle = False
# ForceTitle = '$g_{A}$ Variational Comparison'
# ForceTitle = '$g_{A}$ Variational Parameter $t_{0}$ Comparison'
# ForceTitle = '$g_{A}$ Summed Ratio Function Comparison'

colourset8 = [ '#000080','#B22222','#00B800','#8B008B', '#0000FF','#FF0000','#000000','#32CD32','#FF0066']
markerset = ['o','s','^','v','d']
shiftmax = 3
shiftper = 0.05 ##R function shift
shiftperff = 0.01 ##Form Factor shift
shiftset = [0]
for ish in np.arange(1,shiftmax+1): shiftset += [-ish*shiftper,ish*shiftper]
shiftsetff = [0]
for ish in np.arange(1,shiftmax+1): shiftsetff += [-ish*shiftperff,ish*shiftperff]
incr = 0.01
thisalpha = 0.3

# MassTVals = 16,33
# MassTVals = 3,25
MassTVals = 3,25
# Massyrange = 0.44,0.54
# Massyrange = 0.40,0.60
Massyrange = 0.40,0.61
Qsqrdxlim = -0.03,1

ylimDict = {'VectorP4giDi':[-0.4,-0.15],
            'VectorP4g4':[0,4.0],
            'ProtonP3g2cmplx':[0,0.5],
            # 'IsoVectorP3g3g5':[1.0,1.2], ## Reg
            # 'IsoVectorP3g3g5':[0.9,1.3], ## PoF
            'IsoVectorP3g3g5':[1.0,1.1], ## ExStExample
            # 'IsoVectorP3g3g5':[.95,1.25], ## Tsink
            # 'IsoVectorP4I':[0.4,1.2], ## Tsink?
            'IsoVectorP4giDi':[0.32,0.18], ## Reg
            # 'IsoVectorP4giDi':[0.35,0.15], ## Tsink
            # 'IsoVectorP3g1g2':[1.2,1.04], ## Reg
            # 'IsoVectorP3g1g2':[1.15,.98], ## Tsink
            # 'IsoVectorP3g1g2':[1.23,1.], ## OSFTsink and TSFTsink
            # 'IsoVectorP3g1g2':[1.18,1.02], ## TSFCM 
            'IsoVectorP3g1g2':[1.22,1.02], ## Tsink Var 
            'IsoVectorP4g4':[0.8,1.1]}

ylimFFDict = {'ProtonGeGmFF1/F1divF2':[0.3,0.5],
              # 'NeutronGeGmFF2':[-1.2,-0.4],
              'ProtonVectorFF2':[0.4,1.2],
              'NeutronVectorFF2':[-1.2,-0.4],
              'IsoVectorPsVectorFF2':[10,2],
              'ProtonTensorFF2':[-0.5,-3],
              'ProtonTensorFF3':[0.9,0.15],
              'NeutronTensorFF2':[2.0,0.45],
              'NeutronTensorFF3':[-0.2,-1],
              'NeutronGeGmFF1/F1divF2':[0.0,0.06]}

leglocFFDict = {'NeutronVectorFF2':'upper left',
                'NeutronGeGmFF1':'upper left',
                'NeutronGeGmFF2':'upper right',
                'IsoVectorPsScalarFF1':'lower right',
                'IsoVectorPsVectorFF2':'upper left',
                'NeutronGeGmFF1/F1divF2':'upper left',
                'ProtonTensorFF3':'upper left',
                'NeutronTensorFF1':'upper left',
                'NeutronTensorFF2':'upper left'}

params = {'legend.fontsize': 10,
          'legend.numpoints': 1,
          'axes.labelsize' : 20,
          'axes.titlesize' : 20,
          'figure.autolayout': True,
          'axes.grid': True,
          'axes.xmargin':0.01,
          'axes.ymargin':0.01}

symcyc = itertools.cycle(markerset)
colcyc = itertools.cycle(colourset8)
shiftcyc = itertools.cycle(shiftset)
shiftcycff = itertools.cycle(shiftsetff)

pl.rcParams.update(params)
RFylab = r'$ R \left(\tau,t\right) $'
RFxlab = r'$ \tau - t/2$'

SUMylab = r'$ S(t)$'
SUMxlab = r'$ t$'


FFylab = r'$ F(Q^{2}) $'
FFxlab = r'$ Q^{2} (GeV)^{2}$'

DatFile = False

def AppendFFDat(xdata,ydata,yerr):
    global DatFile
    if DatFile == False: raise IOError("DatFile not defined yet")
    with open(DatFile,'a') as f:
        for ix,iy,iyerr in zip(xdata,ydata,yerr):
            f.write('{0:>3}  {1:10} \n'.format(ix,MakeValAndErr(iy,iyerr)))
            # f.write(' {1:10} \n'.format(ix,MakeValAndErr(iy,iyerr)))
            
def GetPlotIters():
    return itertools.cycle(markerset),itertools.cycle(colourset8),itertools.cycle(shiftset)

def GetPlotItersff():
    return itertools.cycle(markerset),itertools.cycle(colourset8),itertools.cycle(shiftsetff)

def SetxTicks(thisfig=False):
    xmin,xmax = pl.xlim()
    thisinc = 1
    if xmax-xmin > 16:
        thisinc = 2
    if thisfig == False:
        pl.xticks(np.arange(int(xmin),int(xmax)+thisinc,thisinc))
    else:
        thisinc = 3
        thisfig.axes[-1].set_xticks(np.arange(int(xmin),int(xmax)+thisinc,thisinc))
        

def getsmindex(thissm):return DefSmearList.index(thissm)

def PullThisDicts(datadict,igamma,imom):
    thisdatadict = {}
    for iset in datadict.keys():
        if igamma in datadict[iset]:
            thisdatadict[iset] = datadict[iset][igamma][imom]
    return thisdatadict

def SiftAndSort(thisSetList,comp,nocm=True):
    SetListOut = []
    for iset in thisSetList:
        if any([icomp in iset for icomp in comp]) and ('state' not in iset or not nocm): 
            SetListOut.append(iset)
    # SetListOut.sort()
    return SetListOut


def CreateFFFile(thisCol,thisCurr,thisFF):
    thistitle = thisCol + TitleFixFF(thisCurr.replace('/',' '),thisFF)
    if ForceTitle == False:
        pl.title(thistitle)
    else:
        pl.title(ForceTitle)
    thisdir = outputdir + 'graphs/FormFactors/'+thisCurr + '/'
    mkdir_p(thisdir)
    thisfile = thisCol+thisCurr.replace('/','') + thisFF
    return thisdir + thisfile

def CreateFile(thisflag,thisGamma,thisMom,TitlePref,thisfig=False):
    if 'twopt' in thisGamma:
        if 'Dt' in thisflag:
            thistitle = thisGamma+' '+TitlePref+' $' + thisflag.replace('Dt','\Delta t=') + '$'
        else:
            thistitle = thisGamma+' '+TitlePref+' '+thisflag
        thistitle = thistitle.replace('twopt ','')
    else:
        thistitle = thisGamma+' '+TitlePref+' '+thisflag
    if 'q = 0 0 0' not in thisMom: thistitle += ' '+thisMom
    if ForceTitle == False:
        if thisfig == False:
            pl.title(thistitle)
        else:
            thisfig.suptitle(thistitle, fontsize=20)
    else:
        # pl.title(ForceTitle+'$' + thisflag.replace('Dt','\Delta t') + '$')
        if thisfig == False:
            pl.title(ForceTitle)
        else:
            thisfig.suptitle(ForceTitle, fontsize=20)
    thisdir = outputdir + 'graphs/'+CreateOppDir(thisGamma)
    thisfile = TitlePref.replace(' ','')+thisflag
    thisdir += MakeMomDir(thisMom)
    mkdir_p(thisdir)
    return thisdir+thisfile

def SetSumFunAxies(DoY):
    pl.ylabel(SUMylab)
    if DoY:pl.xlabel(SUMxlab)

def SetRFAxies(thisGamma):
    pl.xlabel(RFxlab)
    pl.ylabel(RFylab)
    if Debug: print 'Hardcoding yaxis limits',thisGamma, ylimDict.keys()
    if thisGamma not in ylimDict.keys():
        pl.ylim(max(pl.ylim()[0],-5),min(pl.ylim()[1],5))
    else:
        # pl.ylim(max(ylimDict[thisGamma][0],pl.ylim()[0]),min(ylimDict[thisGamma][1],pl.ylim()[1]))
        pl.ylim(*ylimDict[thisGamma])
    SetxTicks()
    pl.legend()
    pl.tight_layout()
    
def SetFFAxies(thisCurr):
    pl.xlabel(FFxlab)
    pl.ylabel(FFylab)
    if thisCurr in ylimFFDict.keys():
        pl.ylim(ylimFFDict[thisCurr])
    # else:
        # pl.ylim(max(pl.ylim()[0],-2),min(pl.ylim()[1],2))
    pl.xlim(*Qsqrdxlim)
    if thisCurr not in leglocFFDict.keys():
        pl.legend()
    else:
        pl.legend(loc=leglocFFDict[thisCurr])
    pl.tight_layout()


def SetMassAxies():
    pl.xlabel(r'$t$')
    pl.ylabel(r'$aM_{N}$')
    pl.xlim(MassTVals)
    pl.ylim(Massyrange)
    SetxTicks()
    pl.legend()
    pl.tight_layout()
    
def SetLogAxies():
    pl.xlabel(r'$t$')
    pl.ylabel(r'$log(G_{2})$')
    pl.xlim(6,16)
    pl.ylim(-4,-10)
    SetxTicks()
    pl.legend()
    pl.tight_layout()


def PlotTSinkData(data,thisSetList,thisGamma,thisMom,FT,thissm='sm32'):
    global ForceTitle
    ForceTitle = FT
    PlotCol(data,thisSetList,[thissm],thisGamma,thisMom,'TSink Comparison ')

def PlotTSinkSumData(data,thisSetList,thisGamma,thisMom,FT,thissm='sm32'):
    global ForceTitle
    ForceTitle = FT
    for ifitr in SumFitRList:    
        PlotColSum(data,thisSetList,[thissm],thisGamma,thisMom,'Sum TSink Comparison ',thisTsinkR=ifitr)    
    # pl.rcParams['figure.autolayout'] = False
    thisfig = pl.figure()
    for ic,ifitr in enumerate(SumFitRList):    
        # pl.subplot(1,len(SumFitRList),ic+1)
        pl.subplot(len(SumFitRList),1,ic+1)
        PlotColSumFun(data,thisSetList,[thissm],thisGamma,thisMom,'Sum TSink Comparison ',thisfig,thisTsinkR=ifitr)
        SetSumFunAxies(ic==len(SumFitRList)-1)
    # thisfig.subplots_adjust(top=0.95)
    thisfig.savefig(CreateFile(thissm,thisGamma,thisMom,'Sum TSink Comparison ',thisfig=thisfig)+'Sfun.pdf')
    thisfig.clf()
    # pl.rcParams['figure.autolayout'] = True
        
def PlotTSinkSFData(data,data2pt,thisSetList,thisGamma,thisMom,FT,thisSF='TSFTsink',thissm='sm32'):
    global ForceTitle
    ForceTitle = FT
    if 'TSF' in thisSF:
        for icut in TSFCutList:
            PlotColTSF(data,data2pt,thisSetList,[thissm],thisGamma,thisMom,thisSF+' Comparison ',icut,thisSF)
    if 'OSF' in thisSF:
        for icut in OSFCutList:
            PlotColOSF(data,data2pt,thisSetList,[thissm],thisGamma,thisMom,thisSF+' Comparison ',icut,thisSF)
        
def PlotCMData(data,thisSetList,thisGamma,thisMom,FT,thistsink='tsink29'):
    global ForceTitle
    ForceTitle = FT
    PlotCol(data,thisSetList,[thistsink,'PoF'],thisGamma,thisMom,'Variational Comparison ')


def PlotMassData(data,thisSetList,thisMom,FT,TitleFlag=''):
    global ForceTitle
    ForceTitle = FT
    for thisDt in MassDtList:
        PlotRFSet(data,thisSetList,MassDt=thisDt)
        SetMassAxies()        
        pl.savefig(CreateFile('Dt'+str(thisDt),'twopt',thisMom,TitleFlag+' Mass Comparison')+'.pdf')
        pl.clf()
    PlotLogSet(data,thisSetList)
    SetLogAxies()
    pl.savefig(CreateFile('','twopt',thisMom,TitleFlag+' Log Comparison')+'.pdf')
    pl.clf()

def PlotMassSFData(data,thisSetList,thisMom,FT,thisSF='SFCM'):
    global ForceTitle
    ForceTitle = FT
    for thisDt in MassDtList:
        PlotMassSetOSF(data,thisSetList,thisDt,thisSF)
        SetMassAxies()
        pl.savefig(CreateFile('Dt'+str(thisDt),'twopt',thisMom,'O'+thisSF+' Mass Comparison')+'.pdf')
        pl.clf()
        PlotMassSetTSF(data,thisSetList,thisDt,thisSF)
        SetMassAxies()
        pl.savefig(CreateFile('Dt'+str(thisDt),'twopt',thisMom,'T'+thisSF+' Mass Comparison')+'.pdf')
        pl.clf()
    PlotLogSetOSF(data,thisSetList,thisSF)
    SetLogAxies()
    pl.savefig(CreateFile('','twopt',thisMom,'O'+thisSF+' Log Comparison')+'.pdf')
    pl.clf()
    PlotLogSetTSF(data,thisSetList,thisSF)
    SetLogAxies()
    pl.savefig(CreateFile('','twopt',thisMom,'T'+thisSF+' Log Comparison')+'.pdf')
    pl.clf()

def PlotCMOSFData(data,data2pt,thisSetList,thisGamma,thisMom,FT,thistsink='tsink29',thisSF='OSFCM'):
    global ForceTitle
    ForceTitle = FT
    for icut in OSFCutList:
        PlotColOSF(data,data2pt,thisSetList,[thistsink,'PoF','CM'],thisGamma,thisMom,thisSF+' Comparison ',icut,thisSF)
        
def PlotCMTSFData(data,data2pt,thisSetList,thisGamma,thisMom,FT,thistsink='tsink29',thisSF='TSFCM'):
    global ForceTitle
    ForceTitle = FT
    for icut in TSFCutList:
        PlotColTSF(data,data2pt,thisSetList,[thistsink],thisGamma,thisMom,thisSF+' Comparison ',icut,thisSF)



def PlotCol(data,thisSetList,thisflag,thisGamma,thisMom,TitlePref):
    if thisflag[0] == 's':
        thislegrem = 'state1'
    else:
        thislegrem = thisflag[0]
    PlotRFSet(data,SiftAndSort(thisSetList,thisflag,nocm=False),legrem=thislegrem)
    SetRFAxies(thisGamma)
    pl.savefig(CreateFile(thisflag[0],thisGamma,thisMom,TitlePref)+'.pdf')
    pl.clf()

def PlotColTSF(data,data2pt,thisSetList,thisflag,thisGamma,thisMom,TitlePref,TSFcut,thisTSF):
    PlotRFSetTSF(data,data2pt,SiftAndSort(thisSetList,thisflag),TSFcut,thisTSF,legrem=thisflag[0])
    SetRFAxies(thisGamma)
    pl.savefig(CreateFile(thisflag[0],thisGamma,thisMom,TitlePref)+str(TSFcut)+'.pdf')
    pl.clf()

def PlotColOSF(data,data2pt,thisSetList,thisflag,thisGamma,thisMom,TitlePref,OSFcut,thisOSF):
    PlotRFSetOSF(data,data2pt,SiftAndSort(thisSetList,thisflag,nocm=False),OSFcut,thisOSF,legrem=thisflag[0])
    SetRFAxies(thisGamma)
    pl.savefig(CreateFile(thisflag[0],thisGamma,thisMom,TitlePref)+str(OSFcut)+'.pdf')
    pl.clf()

def PlotColSum(data,thisSetList,thissm,thisGamma,thisMom,TitlePref,thisTsinkR='fit sl 0-4'):
    PlotRFSetSum(data,SiftAndSort(thisSetList,thissm),thisTsinkR,legrem=thissm[0])
    outTR = thisTsinkR.replace('-','_')
    SetRFAxies(thisGamma)
    pl.savefig(CreateFile(thissm[0],thisGamma,thisMom,TitlePref+outTR)+'.pdf')
    pl.clf()

def PlotColSumFun(data,thisSetList,thissm,thisGamma,thisMom,TitlePref,thisfig,thisTsinkR='fit sl 0-4'):
    if CheckDict(data,'SumMeth',thissm[0]): PlotSummedRF(data['SumMeth'][thissm[0]],thisTsinkR,thisfig)

def PlotRFSetTSF(data,data2pt,thisSetList,TSFcut,thisTSF,legrem=''):
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    for iset in SortMySet(thisSetList)[0]:
        thistsink,thissm = SplitTSinkString(iset)
        thiscol = thiscolcyc.next()
        thissym = thissymcyc.next()
        thisshift = thisshiftcyc.next()
        if not CheckDict(data,'RF'+thisTSF,iset): continue
        PlotRF(data['RF'+thisTSF][iset],thiscol,thissym,thisshift,LegLab(iset.replace(legrem,'')))
        if not CheckDict(data,thisTSF,thissm) or not CheckDict(data2pt,thisTSF,thissm): continue
        PlotTSFLine(data[thisTSF][thissm],data2pt[thisTSF][thissm],thistsink.replace('tsink',''),thiscol,thisshift,TSFcut,thissm)
        if 'CM' in thisTSF: PlotTSFValue(data[thisTSF][thissm],thiscol,thisshift,TSFcut,thissm,thistsink.replace('tsink','') )
    if 'CM' not in thisTSF: PlotTSFValue(data[thisTSF][thissm],thiscol,thisshift,TSFcut,thissm,thistsink.replace('tsink','') )

def PlotRFSetOSF(data,data2pt,thisSetList,OSFcut,thisOSF,legrem=''):
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    for iset in SortMySet(thisSetList)[0]:
        thisOSFcut = OSFcut
        intcut = int(OSFcut.replace('cut',''))
        if 'PoF' in iset: thisOSFcut = 'cut'+str(intcut-1)
        if 'sm32' in iset: thisOSFcut = 'cut'+str(intcut+1)
        thistsink,thissm = SplitTSinkString(iset)
        thiscol = thiscolcyc.next()
        thissym = thissymcyc.next()
        thisshift = thisshiftcyc.next()
        if not CheckDict(data,'RF'+thisOSF,iset): continue
        PlotRF(data['RF'+thisOSF][iset],thiscol,thissym,thisshift,LegLab(iset.replace(legrem,'')))
        # PlotOSFLine(data[thisOSF][thissm],data2pt[thisOSF][thissm],thistsink.replace('tsink',''),thiscol,thisOSFcut,thissm)
        OSFset = iset
        if not CheckDict(data,thisOSF,OSFset):
            if not CheckDict(data,thisOSF,thissm):
                continue
            else:
                OSFset = thissm
        PlotOSFValue(data[thisOSF][OSFset],thiscol,thisshift,thisOSFcut,thissm,thistsink.replace('tsink',''))



def PlotRFSetSum(data,thisSetList,thisTsinkR,legrem=''):
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    for iset in SortMySet(thisSetList)[0]:
        thistsink,thissm = SplitTSinkString(iset)
        if not CheckDict(data,'RF',iset): continue
        PlotRF(data['RF'][iset],thiscolcyc.next(),thissymcyc.next(),thisshiftcyc.next(),LegLab(iset.replace(legrem,'')))
    if CheckDict(data,'SumMeth',thissm):
        PlotSumMeth(data['SumMeth'][thissm],thiscolcyc.next(),'Sum '+SumCutPar,thisTsinkR)



def PlotRFSet(data,thisSetList,legrem='',MassDt = False):
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    if MassDt == False:
        iterSetList = SortMySet(thisSetList)[0]
    else:
        iterSetList = SortMySet(ReduceTooMassSet(thisSetList))[0]
    for iset in iterSetList:
        if not CheckDict(data,'RF',iset): continue
        if MassDt == False:
            thiscol,thisshift = thiscolcyc.next(),thisshiftcyc.next()
            thistsink = data['RF'][iset]['tVals'][-1]
            PlotRF(data['RF'][iset],thiscol,thissymcyc.next(),thisshift,LegLab(iset.replace(legrem,'')))
            if CheckDict(data,'Fits',iset):
                PlotFit(data['Fits'][iset],thiscol,thisshift,iset,thistsink)
        else:
            dataplot = deepcopy(data['RF'][iset])
            dataplot['Boot'] = MassFun(dataplot['Boot'],MassDt)
            # dataplot['tVals'] = dataplot['tVals'][:-MassDt] 
            dataplot['tVals'] = dataplot['tVals'][MassDt:] 
            PlotRF(dataplot,thiscolcyc.next(),thissymcyc.next(),thisshiftcyc.next(),LegLab(iset),MP=True)

def PlotLogSet(data,thisSetList,legrem=''):
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    iterSetList = SortMySet(ReduceTooMassSet(thisSetList))[0]
    for iset in iterSetList:
        if not CheckDict(data,'RF',iset): continue
        dataplot = deepcopy(data['RF'][iset])
        dataplot['Boot'] = np.log([tboot/dataplot['Boot'][tsource-1] for tboot in dataplot['Boot']])
        dataplot['Boot'] = GetBootStats(dataplot['Boot'])
        PlotRF(dataplot,thiscolcyc.next(),thissymcyc.next(),thisshiftcyc.next(),LegLab(iset.replace(legrem,'')),MP=True,Log=True)

def PlotLogSetOSF(data,thisSetList,thisSF,legrem=''):
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    iterSetList = SortMySet(ReduceTooMassSet(thisSetList))[0]
    for iset in iterSetList:
        if not CheckDict(data,'RF',iset): continue
        thiscol = thiscolcyc.next()
        thissym = thissymcyc.next()
        thisshift = thisshiftcyc.next()
        dataplot = deepcopy(data['RF'][iset])
        norm = dataplot['Boot'][tsource-1]
        dataplot['Boot'] = np.log([tboot/norm for tboot in dataplot['Boot']])
        dataplot['Boot'] = GetBootStats(dataplot['Boot'])
        PlotRF(dataplot,thiscol,thissym,thisshift,LegLab(iset.replace(legrem,'')),MP=True,Log=True)
        if not CheckDict(data,'O'+thisSF,iset): continue
        PlotOSFLog(data['O'+thisSF][iset],thiscol,iset,norm)



def PlotLogSetTSF(data,thisSetList,thisSF,legrem=''):
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    iterSetList = SortMySet(ReduceTooMassSet(thisSetList))[0]
    for iset in iterSetList:
        if not CheckDict(data,'RF',iset): continue
        thiscol = thiscolcyc.next()
        thissym = thissymcyc.next()
        thisshift = thisshiftcyc.next()
        dataplot = deepcopy(data['RF'][iset])
        norm = dataplot['Boot'][tsource-1]
        dataplot['Boot'] = np.log([tboot/norm for tboot in dataplot['Boot']])
        dataplot['Boot'] = GetBootStats(dataplot['Boot'])
        PlotRF(dataplot,thiscol,thissym,thisshift,LegLab(iset.replace(legrem,'')),MP=True,Log=True)
        if not CheckDict(data,'T'+thisSF,iset): continue
        PlotTSFLog(data['T'+thisSF][iset],thiscol,iset,norm)

def PlotMassSetOSF(data2pt,thisSetList,MassDt,thisSF):
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    iterSetList = SortMySet(ReduceTsink(thisSetList))[0]
    
    for iset in iterSetList:
        if not CheckDict(data2pt,'RF',iset): continue
        thiscol = thiscolcyc.next()
        thissym = thissymcyc.next()
        thisshift = thisshiftcyc.next()
        dataplot = deepcopy(data2pt['RF'][iset])
        dataplot['Boot'] = MassFun(dataplot['Boot'],MassDt)
        dataplot['tVals'] = dataplot['tVals'][MassDt:]
        PlotRF(dataplot,thiscol,thissym,thisshift,LegLab(iset),MP=True)
        if not CheckDict(data2pt,'O'+thisSF,iset): continue
        PlotOSFMassValue(data2pt['O'+thisSF][iset],thiscol,iset,MassDt)


def PlotMassSetTSF(data2pt,thisSetList,MassDt,thisSF):
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    iterSetList = SortMySet(ReduceTsink(thisSetList))[0]
    for iset in iterSetList:
        if not CheckDict(data2pt,'RF',iset): continue
        thiscol = thiscolcyc.next()
        thissym = thissymcyc.next()
        thisshift = thisshiftcyc.next()
        dataplot = deepcopy(data2pt['RF'][iset])
        dataplot['Boot'] = MassFun(dataplot['Boot'],MassDt)
        dataplot['tVals'] = dataplot['tVals'][MassDt:]
        PlotRF(dataplot,thiscol,thissym,thisshift,LegLab(iset),MP=True)
        if 'state1' not in iset and CheckDict(data2pt,'T'+thisSF,iset):
            PlotTSFMassLine(data2pt['T'+thisSF][iset],thiscol,iset,MassDt)
    if CheckDict(data2pt,'T'+thisSF,iterSetList[0]):
        PlotTSFMassValue(data2pt['T'+thisSF][iterSetList[0]],MassDt)

def PlotFFs(data,DSCurr,thisSetList,CollName,FT):
    global ForceTitle
    global DatFile
    ForceTitle = FT
    if len(thisSetList) == 0: return
    thisDS,thisCurr,thisFFComb = SplitDSCurr(DSCurr)
    for iFF in range(1,NoFFPars[thisCurr]+1):
        thisFF = 'FF'+str(iFF)
        if len(thisFFComb) > 1: thisFFComb = '/'+thisFFComb
        DatFile = CreateFFFile(CollName,DSCurr,thisFF)+'.dat'
        WipeFile(DatFile)
        PlotFFSet(data,thisFF,thisSetList,thisCurr,DSCurr.replace('/',''))
        SetFFAxies(thisDS+thisCurr+thisFF+thisFFComb)
        pl.savefig(DatFile.replace('.dat','.pdf'))
        pl.clf()
        

def SkipZeroFF(thisFF,thisset,thisCurr):
    skipzero = False
    if Debug: print thisset+thisFF+thisCurr
    if '2' in thisFF or '3' in thisFF:
        # if any([icheck in thisset for icheck in CheckFFZeroList]):
        skipzero = True
    if 'GedivGm' in thisset+thisFF:
        skipzero = True
    if 'PsVector' in thisset+thisFF+thisCurr or 'Tensor' in thisset+thisFF+thisCurr:
        return skipzero,True
    else:
        return skipzero,False
        
def PlotFFSet(dataset,thisFF,thisSetFlag,thisCurr,thisDSCurr):
    thissymcyc,thiscolcyc,thisshiftcycff = GetPlotItersff()
    collist = []
    if 'Ge' in thisDSCurr or ('Vector' in thisDSCurr and ('PsVector' not in thisDSCurr.replace('IsoVector',''))):
        if 'IsoVector' in thisDSCurr or 'Proton' in thisDSCurr or 'sing' in thisDSCurr:
            FixZ= 1
        elif 'Neutron' in thisDSCurr:
            FixZ= 0
        elif 'doub' in thisDSCurr:
            FixZ= 2
        else:
            FixZ=False
    else:
        FixZ=False    
    for thisset in SortMySet(thisSetFlag)[0]:
        ##make legend formatting function
        if not CheckDict(dataset,thisset,thisFF): continue
        if dataset[thisset][thisFF] == False: continue        
        thiscol = thiscolcyc.next()
        # thisshift = thisshiftcycff.next()
        thisshift = 0.0
        collist.append(thiscol)
        skipzero,flipsign = SkipZeroFF(thisFF,thisset,thisCurr)
        qrange = PlotFF(dataset[thisset][thisFF],thiscol,thissymcyc.next(),thisshift,LegLabFF(thisset),skipzero,flipsign,FixZ=FixZ)
        PlotDPFit(thisset,thisFF,thisDSCurr,thiscol,qrange,thisshift,flipsign)
    return collist

def PlotDPFit(thisset,thisFF,thisCurr,thiscol,qrange,thisshift,flipsign):
    if Debug: print thisset,thisFF
    Avg,Err = GetDPFitValue(thisset,thisFF,thisCurr)
    if len(Avg) == 0: return
    Avg,Err = np.array(Avg),np.array(Err)
    fitqdata = np.arange(qrange[0]-thisshift,qrange[-1]+incr-thisshift,incr)
    fitydataAvg = DPfitfun([fitqdata],Avg)
    fitydataup = np.array([max(iy1,iy2) for iy1,iy2 in zip(DPfitfun([fitqdata],Avg+Err),DPfitfun([fitqdata],np.array([Avg[0]-Err[0],Avg[1]+Err[1]])))])
    Avg,Err = np.array(Avg),np.array(Err)
    fitydatadown = np.array([min(iy1,iy2) for iy1,iy2 in zip(DPfitfun([fitqdata],Avg-Err),DPfitfun([fitqdata],np.array([Avg[0]+Err[0],Avg[1]-Err[1]])))])
    # fitydatadown = np.array([np.min(iy1,iy2) for iy1,iy2 in zip(DPfitfun([fitqdata],Avg-Err),DPfitfun([fitqdata],np.array([Avg[0]+Err[0],Avg[1]-Err[1]])))])
    if Debug: print Avg[1], Err[1]
    # if GetCharRad(Avg[1]) > 10 or Err[1]> 10: return
    ## Displays charge radius for FF1, and magnetic moment for FF2
    if 'FF1' in thisFF:
        LegVal = '$\\langle r^2 \\rangle='+MakeValAndErr(GetCharRad(Avg[1]),Err[1])+'\ fm^{2}$'        
    elif 'FF2' in thisFF:
        LegVal = '$\\mu='+MakeValAndErr(Avg[0],Err[0])+'$'        
    else:
        LegVal = 'nothing'        
    if flipsign:
        pl.plot(fitqdata+thisshift,-np.array(fitydataAvg),label=LegVal,color=thiscol)
        pl.fill_between(fitqdata+thisshift,-np.array(fitydataup),-np.array(fitydatadown),color=thiscol,alpha=thisalpha,edgecolor='none')
    else:
        pl.plot(fitqdata+thisshift,fitydataAvg,label=LegVal,color=thiscol)
        pl.fill_between(fitqdata+thisshift,fitydataup,fitydatadown,color=thiscol,alpha=thisalpha,edgecolor='none')

    
def PlotFF(data,col,sym,shift,lab,SkipZero,FlipSign,FixZ=False):
    qsqrdvals,dataavg,dataerr = [],[],[]
    for iqsqrd,(qsqrd,values) in enumerate(data.iteritems()):
        Qsqrd = GetQsqrd(float(qsqrd.replace('qsqrd','')),Phys=PhysicalUnits)
        thisshift = shift
        if PhysicalUnits: thisshift *= hbarcdivlat**2
        qsqrdvals.append(Qsqrd+shift)
        if 'Boot' in values.keys():
            dataavg.append(values['Boot'].Avg)
            dataerr.append(values['Boot'].Std)
        else:
            dataavg.append(values['Avg'])
            dataerr.append(values['Std'])
    if len(qsqrdvals) > 0:
        if len(qsqrdvals) > Qtcut:
            qsqrdvals,dataavg,dataerr = qsqrdvals[:Qtcut],dataavg[:Qtcut],dataerr[:Qtcut]
        if FlipSign: dataavg = -1*np.array(dataavg)
        if ForcePos: dataavg = np.abs(dataavg)
        AppendFFDat(qsqrdvals,dataavg,dataerr)
        if SkipZero and len(qsqrdvals) > 1:
            pl.errorbar(qsqrdvals[1:],dataavg[1:],dataerr[1:],color=col,fmt=sym,label=lab)
        elif FixZ != False:
            pl.plot([0],[FixZ],sym,color=col)            
            pl.errorbar(qsqrdvals[1:],dataavg[1:],dataerr[1:],color=col,fmt=sym,label=lab)            
        else:
            pl.errorbar(qsqrdvals,dataavg,dataerr,color=col,fmt=sym,label=lab)
    return qsqrdvals

def PlotRF(data,col,sym,shift,lab,MP=False,Log=False):
    if MP:
        thistsource = tsource +1
        if 'PoF' in lab and not Log:
            tvals = np.array(data['tVals'])+1+(2*PoFShifts) + shift
        else:
            tvals = np.array(data['tVals'])+1 + shift
    else:
        tvals = np.array(data['tVals'])
        tvals = tvals-(tvals[-1]+tvals[0])/2.0 + shift
        thistsource = 0
    dataavg = Pullflag(data['Boot'],'Avg')
    dataerr = Pullflag(data['Boot'],'Std')
    # if Debug:
    #     print lab
    #     for it,val,valerr in zip(tvals,dataavg,dataerr):
    #         print tvals,dataavg,dataerr
    if ForcePos: dataavg = np.abs(dataavg)
    pl.errorbar(tvals[thistsource:]-thistsource,dataavg[thistsource:],dataerr[thistsource:],color=col,fmt=sym,label=lab)


def PlotFit(data,col,shift,iset,thistsink):
    thiscut = GetCut(iset,FitCutPicked)
    if thiscut == False:
        if Debug:
            print 'warning', iset, 'not in FitCutPicked'
            print FitCutPicked.keys()
    else:
        thiscutint = int(thiscut.replace('cut',''))
        LRM = max((int(thistsink)-tsource)/2.-thiscutint,0)
        tvals = np.array([-LRM+shift,LRM+shift])
        # if Debug:
        #     print data.keys(), thiscut
        #     print data[thiscut].keys(), 'Boot'
        dataavg = Pullflag(data[thiscut]['Boot'],'Avg')
        dataerr = Pullflag(data[thiscut]['Boot'],'Std')
        if ForcePos: dataavg = np.abs(dataavg)
        dataup,datadown = dataavg+dataerr,dataavg-dataerr
        if Debug:
            print iset, thiscut
            print tvals[0],tvals[1],dataavg, dataerr
        pl.plot(tvals,[dataavg,dataavg],color=col)
        pl.fill_between(tvals,dataup,datadown,color=col,alpha=thisalpha,edgecolor='none')

def PlotSumMeth(data,col,lab,thisTsinkR):
    if not CheckDict(data,SumCutPar,thisTsinkR,'Avg'): return
    if not CheckDict(data,SumCutPar,thisTsinkR,'Std'): return
    dataavg = data[SumCutPar][thisTsinkR]['Avg']
    dataerr = data[SumCutPar][thisTsinkR]['Std']
    if ForcePos: dataavg = np.abs(dataavg)
    dataup,datadown = dataavg+dataerr,dataavg-dataerr
    trange = AllTSinkListVar[-1]-tsource - 2*int(SumCutPar.replace('cut',''))
    tvals = [-trange/2.,trange/2.]
    # pl.axhline(dataavg,color=col,label=lab)
    # pl.axhspan(datadown,dataup,facecolor=col,edgecolor='none',alpha=thisalpha)
    pl.plot(tvals,[dataavg,dataavg],color=col,label=lab)
    pl.fill_between(tvals,dataup,datadown,color=col,alpha=thisalpha,edgecolor='none')
    
def PlotSummedRF(data,thisfitr,thisfig):
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    thisfitr = thisfitr.split()[-1]
    thisfitmin,thisfitmax = thisfitr.split('-')
    for icut,cutdata in data.iteritems():
        if 'nconf' in cutdata.keys()[0]: continue
        tdata,dataplot,dataploterr = [],[],[]
        thiscol,thissym = thiscolcyc.next(),thissymcyc.next()
        
        for itsink,tsinkdata in cutdata.iteritems():
            if 'fit' not in itsink and 'nconf' not in itsink:
                tdata.append(int(itsink.replace('tsink',''))-tsource)
                dataplot.append(tsinkdata['Avg'])
                dataploterr.append(tsinkdata['Std'])
        tdata,dataplot,dataploterr = zip(*sorted(zip(tdata,dataplot,dataploterr)))
        if ForcePos:
            thisfig.axes[-1].errorbar(tdata,np.abs(dataplot),dataploterr,color=thiscol,fmt=thissym,label=icut)
        else:
            thisfig.axes[-1].errorbar(tdata,dataplot,dataploterr,color=thiscol,fmt=thissym,label=icut)
        if not CheckDict(cutdata,'fit con '+thisfitr,'Boot'): continue
        if not CheckDict(cutdata,'fit sl '+thisfitr,'Boot'): continue        
        parcon,parsl = cutdata['fit con '+thisfitr]['Boot'],cutdata['fit sl '+thisfitr]['Boot']
        fittdata = np.arange(tdata[int(thisfitmin)],tdata[int(thisfitmax)]+incr,incr)
        fittdashed = np.arange(tdata[0],tdata[int(thisfitmin)]+incr,incr)
        fitbootdata = [(parsl* (it+tsource)) + parcon for it in fittdata]
        fitbootdashed = [(parsl* (it+tsource)) + parcon for it in fittdashed]
        GetBootStats(fitbootdata),GetBootStats(fitbootdashed)
        plotup = Pullflag(fitbootdata,'Avg')+Pullflag(fitbootdata,'Std')
        plotdown = Pullflag(fitbootdata,'Avg')-Pullflag(fitbootdata,'Std')
        if len(fitbootdashed) == 0:
            plotdashedup,plotdasheddown = [],[]
        else:
            plotdashedup = Pullflag(fitbootdashed,'Avg')+Pullflag(fitbootdashed,'Std')
            plotdasheddown = Pullflag(fitbootdashed,'Avg')-Pullflag(fitbootdashed,'Std')
        if ForcePos:
            thisfig.axes[-1].plot(fittdata,np.abs(Pullflag(fitbootdata,'Avg')),
                                  label='slope='+MakeValAndErr(parsl.Avg,parsl.Std),color=thiscol)
            thisfig.axes[-1].fill_between(fittdata,np.abs(plotup),np.abs(plotdown),color=thiscol,alpha=thisalpha,edgecolor='none')
            thisfig.axes[-1].plot(fittdashed,np.abs(plotdashedup),color=thiscol,ls='--')
            thisfig.axes[-1].plot(fittdashed,np.abs(plotdasheddown),color=thiscol,ls='--')
        else:
            thisfig.axes[-1].plot(fittdata,Pullflag(fitbootdata,'Avg'),
                                  label='slope='+MakeValAndErr(parsl.Avg,parsl.Std),color=thiscol)
            thisfig.axes[-1].fill_between(fittdata,plotup,plotdown,color=thiscol,alpha=thisalpha,edgecolor='none')
            thisfig.axes[-1].plot(fittdashed,plotdashedup,color=thiscol,ls='--')
            thisfig.axes[-1].plot(fittdashed,plotdasheddown,color=thiscol,ls='--')
    # if Pullflag(fitbootdata,'Avg')[0] > 0:
    #     thisfig.legend(loc='upper left')
    # else:
    #     thisfig.legend(loc='upper right')
    SetxTicks(thisfig=thisfig)



# # PlotO/T SFLine only for zero momenta (no point otherwise)
# def PlotOSFLine(data,data2pt,thistsink,col,OSFcut,smear):
#     pars2pt,pars3pt = [],[]
#     thistsink = int(thistsink)-tsource + 1
#     cutint = int(OSFcut.replace('cut',''))
#     for ipar in StateParList['One']['C2']:
#         pars2pt.append(data2pt[ipar][OSFfitr[smear]]['Avg'])
#     for ipar in StateParList['One']['C3']:
#         pars3pt.append(data[ipar][OSFfitr[smear]][OSFcut]['Avg'])
#     fit3ptfun,f3fder = ff.CreateC3OSFitFun(pars2pt+pars2pt,False,NoExp=True)
#     tdata = np.arange(cutint,thistsink-cutint+incr,incr)
#     tplotdata = tdata - thistsink/2.0 
#     dataline = (fit3ptfun(np.array([tdata,[float(thistsink)]*len(tdata)]),pars3pt)/
#                 ff.C2OneStateFitFunNoExp(float(thistsink),pars2pt))
#     pl.plot(tplotdata,map(abs,dataline),color=col)


def PlotTSFLine(data,data2pt,thistsink,col,thisshift,TSFcut,smear):
    pars2pt,pars3pt = [],[]
    thistsink = int(thistsink)-tsource 
    cutint = int(TSFcut.replace('cut',''))
    LRM = max((thistsink/2.)-cutint,0)
    for ipar in StateParList['Two']['C2']:
        if not CheckDict(data2pt,ipar,TSFfitr,'Avg'): return
        pars2pt.append(data2pt[ipar][TSFfitr]['Avg'])
    for ipar in StateParList['Two']['C3']:
        if not CheckDict(data,ipar,TSFfitr,TSFcut,'Avg'): return
        pars3pt.append(data[ipar][TSFfitr][TSFcut]['Avg'])
    tdata = np.arange(cutint,max(thistsink-cutint+incr,cutint),incr)
    tplotdata = tdata - tdata[0] - LRM + thisshift
    # tplotdata = np.arange(-LRM,LRM+incr,incr) 
    dataline = (np.array(ff.C3TSFLineFun(pars2pt+pars3pt,tdata,thistsink))/
                ff.C2TSFLineFun(thistsink,pars2pt))
    if ForcePos: dataline = np.abs(dataline)
    pl.plot(tplotdata,dataline,color=col)


def PlotOSFValue(data,col,thisshift,OSFcut,smear,thistsink):
    LRM = max((int(thistsink)-tsource)/2.-int(OSFcut.replace('cut','')),0)
    if LRM == 0: return
    tvals = np.array([-LRM+thisshift,LRM+thisshift])
    osfsmear,dump = CreateOSFfitKey(smear)
    if not CheckDict(data,'B00',OSFfitr[osfsmear],OSFcut,'Avg'): return
    if not CheckDict(data,'B00',OSFfitr[osfsmear],OSFcut,'Std'): return
    dataval = data['B00'][OSFfitr[osfsmear]][OSFcut]['Avg']
    dataerr = data['B00'][OSFfitr[osfsmear]][OSFcut]['Std']
    if ForcePos: dataval = np.abs(dataval)
    dataup,datadown = dataval+dataerr,dataval-dataerr
    pl.fill_between(tvals,[datadown,datadown],[dataup,dataup],facecolor=col,edgecolor='none',alpha=thisalpha)
    pl.plot(tvals,[dataval,dataval],color = col)

def PlotTSFValue(data,col,thisshift,TSFcut,smear,thistsink):
    LRM = max((int(thistsink)-tsource)/2.-int(TSFcut.replace('cut','')),0)
    tvals = np.array([-LRM+thisshift,LRM+thisshift])
    if not CheckDict(data,'B00',TSFfitr,TSFcut,'Avg'): return
    dataAvg,dataErr = data['B00'][TSFfitr][TSFcut]['Avg'],data['B00'][TSFfitr][TSFcut]['Std']
    if ForcePos: dataAvg = np.abs(dataAvg)
    dataup,datadown = dataAvg-dataErr,dataAvg+dataErr
    pl.fill_between(tvals,[datadown,datadown],[dataup,dataup],facecolor=col,edgecolor='none',alpha=thisalpha)



# PlotO/T SFMassLine only for zero momenta (no point otherwise)

def PlotTSFMassLine(data2pt,col,smear,thisdt):
    pars2pt = []
    for ipar in StateParList['Two']['C2']:
        if not CheckDict(data2pt,ipar,TSFfitr,'Avg'): return
        pars2pt.append(data2pt[ipar][TSFfitr]['Avg'])
    tdata = np.arange(TSFfitvals[0]-tsource+1,TSFfitvals[1]-thisdt+incr-tsource+1,incr)
    fit2pt = ff.C2TSFLineFun(tdata,pars2pt)
    fit2ptdt = ff.C2TSFLineFun(tdata+thisdt,pars2pt)
    tplotdata = tdata + thisdt - 1
    pl.plot(tplotdata,map(abs,np.log(fit2ptdt/fit2pt)/thisdt),color=col)


def PlotOSFMassValue(data,col,smear,thisdt):
    smearindex,deltashift = CreateOSFfitKey(smear)
    if Debug:
        print data.keys(), 'm0'
        print data['m0'].keys(), OSFfitr[smearindex]
        print data['m0'][OSFfitr[smearindex]].keys(),'Boot','Avg','Std'
    if CheckDict(data,'m0',OSFfitr[smearindex],'Boot'): 
        databoot = data['m0'][OSFfitr[smearindex]]['Boot']
        dataval = abs(databoot.Avg)
        dataup,datadown = dataval+databoot.Std,dataval-databoot.Std
    elif CheckDict(data,'m0',OSFfitr[smearindex],'Avg') and CheckDict(data,'m0',OSFfitr[smearindex],'Std'):
        dataval = data['m0'][OSFfitr[smearindex]]['Avg']
        dataup,datadown = dataval+data['m0'][OSFfitr[smearindex]]['Std'],dataval-data['m0'][OSFfitr[smearindex]]['Std']
    else:
        if Debug: print 'OSF',smear, 'not found'
        return
    tdata = [OSFfitvals[smearindex][0]+thisdt+deltashift-tsource,OSFfitvals[smearindex][1]+deltashift-tsource]
    pl.fill_between(tdata,[datadown,datadown],[dataup,dataup],facecolor=col,edgecolor='none',alpha=thisalpha)
    pl.plot(tdata,[dataval,dataval],color=col)

def PlotTSFMassValue(data,thisdt):
    databoot = data['m0'][TSFfitr]['Boot']
    dataval = abs(databoot.Avg)
    dataup,datadown = dataval+databoot.Std,dataval-databoot.Std
    pl.fill_between([TSFfitvals[0]-tsource+thisdt ,TSFfitvals[1]-tsource],[datadown,datadown],[dataup,dataup],facecolor='k',edgecolor='none',alpha=thisalpha)


def PlotOSFLog(data,col,smear,norm):
    smearindex,deltashift = CreateOSFfitKey(smear)
    if 'sum' in smear: smearindex = PickedStateStr+'sum'    
    if not CheckDict(data,'m0',OSFfitr[smearindex],'Boot'): return
    if not CheckDict(data,'Am',OSFfitr[smearindex],'Boot'): return
    parm0 = data['m0'][OSFfitr[smearindex]]['Boot']
    parAm = data['Am'][OSFfitr[smearindex]]['Boot']
    tdata = np.arange(OSFfitvals[smearindex][0]+1,OSFfitvals[smearindex][1]+1+incr,incr)-tsource
    linedata = []
    for it in tdata:
        linedata.append(np.log((parAm*(parm0*(-it)).exp(1))/norm))
    GetBootStats(linedata)
    dataAvg,dataErr = np.array(Pullflag(linedata,'Avg')),np.array(Pullflag(linedata,'Std'))
    dataup = dataAvg+dataErr
    datadown = dataAvg-dataErr
    pl.fill_between(tdata-1,datadown,dataup,facecolor=col,edgecolor='none',alpha=thisalpha)
    pl.plot([tdata[0]-1,tdata[-1]-1],[dataAvg[0],dataAvg[-1]],color=col)


def PlotTSFLog(data,col,smear,norm):
    if not CheckDict(data,'m0',TSFfitr,'Boot'): return
    if not CheckDict(data,'Am',TSFfitr,'Boot'): return
    if not CheckDict(data,'Dm',TSFfitr,'Boot'): return
    if not CheckDict(data,'Amp',TSFfitr,'Boot'): return
    parm0 = data['m0'][TSFfitr]['Boot']
    parAm = data['Am'][TSFfitr]['Boot']
    parDm = data['Dm'][TSFfitr]['Boot']
    parAmp = data['Amp'][TSFfitr]['Boot']
    tdata = np.arange(TSFfitvals[0]+1,TSFfitvals[1]+1+incr,incr)-tsource
    linedata = []
    for it in tdata:
        thisit = it
        # thisit = it-tsource
        linedata.append(np.log((parAm*(((parm0*(-thisit)).exp(1)) + parAmp*((parm0+parDm)*(-thisit)).exp(1)))/norm))
    GetBootStats(linedata)
    dataAvg,dataErr = np.array(Pullflag(linedata,'Avg')),np.array(Pullflag(linedata,'Std'))
    dataup = dataAvg+dataErr
    datadown = dataAvg-dataErr
    pl.fill_between(tdata-1,datadown,dataup,facecolor=col,edgecolor='none',alpha=thisalpha)
    pl.plot(tdata-1,dataAvg,color=col)


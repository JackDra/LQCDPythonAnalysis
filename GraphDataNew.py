#!/usr/bin/env python

import matplotlib
matplotlib.use('PDF') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as pl
from Params import *
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

##FORCE TITLE PARAMETER, SET TO FALSE TO USE NORMAL TITLES#

# ForceTitle = False
ForceTitle = '$g_{A}$ Summation Comparison'

colourset8 = [ '#000080','#B22222','#00B800','#8B008B', '#0000FF','#FF0000','#000000','#32CD32','#FF0066']
markerset = ['o','s','^','v','d']
shiftmax = 3
shiftper = 0.05
shiftset = [0]
for ish in np.arange(1,shiftmax+1): shiftset += [-ish*shiftper,ish*shiftper]
incr = 0.01
thisalpha = 0.3

# MassTVals = 16,33
MassTVals = 5,19
Massyrange = 0.44,0.54
# Massyrange = 0.40,0.60

ylimDict = {'IsoVectorP4giDi':[-0.15,-0.05],
            'VectorP4giDi':[-0.4,-0.15],
            'VectorP4g4':[0,4.0],
            'ProtonP3g2cmplx':[0,0.5],
            'IsoVectorP3g3g5':[-1.0,-1.2],
            'IsoVectorP4g4':[0.8,1.1]}

ylimFFDict = {'ProtonGeGmFF1/F1divF2':[0.3,0.5],
              'NeutronGeGmFF1/F1divF2':[0.0,0.1]}

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

pl.rcParams.update(params)
RFylab = r'$ R \left(\tau,t\right) $'
RFxlab = r'$ \frac{\tau}{a} - \frac{t}{2a}$'

SUMylab = r'$ S(t)$'
SUMxlab = r'$ \frac{t}{a}$'


FFylab = r'$ FF $'
FFxlab = r'$ q^{2} $'

def GetPlotIters():
    return itertools.cycle(markerset),itertools.cycle(colourset8),itertools.cycle(shiftset)

def SetxTicks():
    xmin,xmax = pl.xlim()
    thisinc = 1
    if xmax-xmin > 15:
        thisinc = 2
    pl.xticks(np.arange(int(xmin),int(xmax)+thisinc,thisinc))

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
    thistitle = thisCol + ' ' + thisCurr.replace('/',' ') + ' ' + thisFF
    pl.title(thistitle)
    thisdir = outputdir + 'graphs/FormFactors/'+thisCurr + '/'
    mkdir_p(thisdir)
    thisfile = thisCol+thisCurr.replace('/','') + thisFF
    return thisdir + thisfile

def CreateFile(thisflag,thisGamma,thisMom,TitlePref,suptitle=False):
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
        if suptitle:
            pl.suptitle(thistitle)
        else:
            pl.title(thistitle)
    else:
        # pl.title(ForceTitle+'$' + thisflag.replace('Dt','\Delta t') + '$')
        if suptitle:
            pl.suptitle(ForceTitle)
        else:
            pl.title(ForceTitle)
    thisdir = outputdir + 'graphs/'+CreateOppDir(thisGamma)
    thisfile = TitlePref.replace(' ','')+thisflag
    thisdir += MakeMomDir(thisMom)
    mkdir_p(thisdir)
    return thisdir+thisfile

def SetSumFunAxies(DoY):
    if DoY:pl.ylabel(SUMylab)
    pl.xlabel(SUMxlab)

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
    if thisCurr not in ylimFFDict.keys():
        pl.ylim(max(pl.ylim()[0],-1),min(pl.ylim()[1],2))
    else:
        pl.ylim(ylimFFDict[thisCurr])
    pl.legend()
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


def PlotTSinkData(data,thisSetList,thisGamma,thisMom,thissm='sm32'):
    PlotCol(data,thisSetList,[thissm],thisGamma,thisMom,'TSink Comparison ')

def PlotTSinkSumData(data,thisSetList,thisGamma,thisMom,thissm='sm32'):
    for ifitr in SumFitRList:    
        PlotColSum(data,thisSetList,[thissm],thisGamma,thisMom,'Sum TSink Comparison ',thisTsinkR=ifitr)
    
    for ic,ifitr in enumerate(SumFitRList):    
        pl.subplot(1,len(SumFitRList),ic+1)
        PlotColSumFun(data,thisSetList,[thissm],thisGamma,thisMom,'Sum TSink Comparison ',thisTsinkR=ifitr)
        SetSumFunAxies(ic==0)
    pl.subplots_adjust(top=0.85)
    pl.savefig(CreateFile(thissm,thisGamma,thisMom,'Sum TSink Comparison ',suptitle=True)+'Sfun.pdf')
    pl.clf()
        
def PlotTSinkSFData(data,data2pt,thisSetList,thisGamma,thisMom,thisSF='TSFTsink',thissm='sm32'):
    if 'TSF' in thisSF:
        for icut in TSFCutList:
            PlotColTSF(data,data2pt,thisSetList,[thissm],thisGamma,thisMom,thisSF+' Comparison ',icut,thisSF)
    if 'OSF' in thisSF:
        for icut in OSFCutList:
            PlotColOSF(data,data2pt,thisSetList,[thissm],thisGamma,thisMom,thisSF+' Comparison ',icut,thisSF)
        
def PlotCMData(data,thisSetList,thisGamma,thisMom,thistsink='tsink29'):
    PlotCol(data,thisSetList,[thistsink,'PoF'],thisGamma,thisMom,'Variational Comparison ')


def PlotMassData(data,thisSetList,thisMom,TitleFlag=''):
    for thisDt in MassDtList:
        PlotRFSet(data,thisSetList,MassDt=thisDt)
        SetMassAxies()        
        pl.savefig(CreateFile('Dt'+str(thisDt),'twopt',thisMom,TitleFlag+' Mass Comparison')+'.pdf')
        pl.clf()
    PlotLogSet(data,thisSetList)
    SetLogAxies()
    pl.savefig(CreateFile('','twopt',thisMom,TitleFlag+' Log Comparison')+'.pdf')
    pl.clf()

def PlotMassSFData(data,thisSetList,thisMom,thisSF='SFCM'):
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

def PlotCMOSFData(data,data2pt,thisSetList,thisGamma,thisMom,thistsink='tsink29',thisSF='OSFCM'):
    for icut in OSFCutList:
        PlotColOSF(data,data2pt,thisSetList,[thistsink,'PoF','CM'],thisGamma,thisMom,thisSF+' Comparison ',icut,thisSF)
        
def PlotCMTSFData(data,data2pt,thisSetList,thisGamma,thisMom,thistsink='tsink29',thisSF='TSFCM'):
    for icut in TSFCutList:
        PlotColTSF(data,data2pt,thisSetList,[thistsink],thisGamma,thisMom,thisSF+' Comparison ',icut,thisSF)



def PlotCol(data,thisSetList,thisflag,thisGamma,thisMom,TitlePref):
    PlotRFSet(data,SiftAndSort(thisSetList,thisflag,nocm=False),legrem=thisflag[0])
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

def PlotColSumFun(data,thisSetList,thissm,thisGamma,thisMom,TitlePref,thisTsinkR='fit sl 0-4'):
    if CheckDict(data,'SumMeth',thissm[0]): PlotSummedRF(data['SumMeth'][thissm[0]],thisTsinkR)

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
        PlotTSFValue(data[thisTSF][thissm],thiscol,thisshift,TSFcut,thissm,thistsink.replace('tsink','') )

def PlotRFSetOSF(data,data2pt,thisSetList,OSFcut,thisOSF,legrem=''):
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    for iset in SortMySet(thisSetList)[0]:
        thistsink,thissm = SplitTSinkString(iset)
        thiscol = thiscolcyc.next()
        thissym = thissymcyc.next()
        thisshift = thisshiftcyc.next()
        if not CheckDict(data,'RF'+thisOSF,iset): continue
        PlotRF(data['RF'+thisOSF][iset],thiscol,thissym,thisshift,LegLab(iset.replace(legrem,'')))
        # PlotOSFLine(data[thisOSF][thissm],data2pt[thisOSF][thissm],thistsink.replace('tsink',''),thiscol,OSFcut,thissm)
        OSFset = iset
        if not CheckDict(data,thisOSF,OSFset):
            if not CheckDict(data,thisOSF,thissm):
                continue
            else:
                OSFset = thissm
        PlotOSFValue(data[thisOSF][OSFset],thiscol,thisshift,OSFcut,thissm,thistsink.replace('tsink',''))



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
            if 'Fits' in data.keys():
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

def PlotFFs(data,DSCurr,thisSetList,CollName):
    if len(thisSetList) == 0: return
    thisDS,thisCurr,thisFFComb = SplitDSCurr(DSCurr)
    for iFF in range(1,NoFFPars[thisCurr]+1):
        if len(thisFFComb) > 1: thisFFComb = '/'+thisFFComb
        thisFF = 'FF'+str(iFF)
        PlotFFSet(data,thisFF,thisSetList)
        SetFFAxies(thisDS+thisCurr+thisFF+thisFFComb)
        pl.savefig(CreateFFFile(CollName,DSCurr,thisFF)+'.pdf')
        pl.clf()
        

def PlotFFSet(dataset,thisFF,thisSetFlag):
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    collist = []
    for thisset in SortMySet(thisSetFlag)[0]:
        ##make legend formatting function
        if not CheckDict(dataset,thisset,thisFF): continue
        if dataset[thisset][thisFF] == False: continue        
        thiscol = thiscolcyc.next()
        collist.append(thiscol)
        PlotFF(dataset[thisset][thisFF],thiscol,thissymcyc.next(),thisshiftcyc.next(),LegLab(thisset))
    return collist

def PlotFF(data,col,sym,shift,lab):
    tvals,dataavg,dataerr = [],[],[]
    for iqsqrd,(qsqrd,values) in enumerate(data.iteritems()):
        ##here is where we can put qsqrd to physical or lattice units##
        tvals.append(int(qsqrd.replace('qsqrd',''))+shift)
        if 'Boot' in values.keys():
            dataavg.append(values['Boot'].Avg)
            dataerr.append(values['Boot'].Std)
        else:
            dataavg.append(values['Avg'])
            dataerr.append(values['Std'])
    if len(tvals) > 0:
        pl.errorbar(tvals,map(abs,dataavg),dataerr,color=col,fmt=sym,label=lab)

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
    dataup,datadown = dataavg+dataerr,dataavg-dataerr
    pl.axhline(dataavg,color=col,label=lab)
    pl.axhspan(datadown,dataup,facecolor=col,edgecolor='none',alpha=thisalpha)

def PlotSummedRF(data,thisfitr):
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
        pl.errorbar(tdata,dataplot,dataploterr,color=thiscol,fmt=thissym,label=icut)
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
        pl.plot(fittdata,Pullflag(fitbootdata,'Avg'),
                label='slope='+MakeValAndErr(parsl.Avg,parsl.Std),color=thiscol)
        pl.fill_between(fittdata,plotup,plotdown,color=thiscol,alpha=thisalpha,edgecolor='none')
        pl.plot(fittdashed,plotdashedup,color=thiscol,ls='--')
        pl.plot(fittdashed,plotdasheddown,color=thiscol,ls='--')
    # if Pullflag(fitbootdata,'Avg')[0] > 0:
    #     pl.legend(loc='upper left')
    # else:
    #     pl.legend(loc='upper right')
    SetxTicks()



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
    pl.plot(tplotdata,dataline,color=col)


def PlotOSFValue(data,col,thisshift,OSFcut,smear,thistsink):
    LRM = max((int(thistsink)-tsource)/2.-int(OSFcut.replace('cut','')),0)
    tvals = np.array([-LRM+thisshift,LRM+thisshift])
    osfsmear,dump = CreateOSFfitKey(smear)
    if not CheckDict(data,'B00',OSFfitr[osfsmear],OSFcut,'Avg'): return
    if not CheckDict(data,'B00',OSFfitr[osfsmear],OSFcut,'Std'): return
    dataval = data['B00'][OSFfitr[osfsmear]][OSFcut]['Avg']
    dataerr = data['B00'][OSFfitr[osfsmear]][OSFcut]['Std']
    dataup,datadown = dataval+dataerr,dataval-dataerr
    pl.fill_between(tvals,[datadown,datadown],[dataup,dataup],facecolor=col,edgecolor='none',alpha=thisalpha)
    pl.plot(tvals,[dataval,dataval],color = col)

def PlotTSFValue(data,col,thisshift,TSFcut,smear,thistsink):
    LRM = max((int(thistsink)-tsource)/2.-int(TSFcut.replace('cut','')),0)
    tvals = np.array([-LRM+thisshift,LRM+thisshift])
    if not CheckDict(data,'B00',TSFfitr,TSFcut,'Avg'): return
    dataAvg,dataErr = data['B00'][TSFfitr][TSFcut]['Avg'],data['B00'][TSFfitr][TSFcut]['Std']
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


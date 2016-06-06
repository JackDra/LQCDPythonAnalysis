#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
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

colourset8 = [ '#000080','#B22222','#00B800','#8B008B', '#0000FF','#FF0000','#000000','#32CD32','#FF0066']
markerset = ['o','s','^','v','d']
shiftmax = 3
shiftper = 0.05
shiftset = [0]
for ish in np.arange(1,shiftmax+1): shiftset += [-ish*shiftper,ish*shiftper]
incr = 0.01
thisalpha = 0.3

MassTVals = 16,34
Massyrange = 0.35,1.0

ylimDict = {'P4giDi':[0.05,0.15]}

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

FFylab = r'$ FF $'
FFxlab = r'$ q^{2} $'

def GetPlotIters():
    return itertools.cycle(markerset),itertools.cycle(colourset8),itertools.cycle(shiftset)

def SetxTicks():
    xmin,xmax = pl.xlim()
    thisinc = 1
    if xmax-xmin > 16:
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
    thistitle = thisCol + ' ' + thisCurr + ' ' + thisFF
    pl.title(thistitle)
    thisdir = outputdir + 'graphs/FormFactors/'+thisCurr + '/'
    mkdir_p(thisdir)
    thisfile = thisCol+thisCurr + thisFF
    return thisdir + thisfile

def CreateFile(thisflag,thisGamma,thisMom,TitlePref,subdir=''):
    thistitle = thisGamma+' '+TitlePref+' '+thisflag
    if 'q = 0 0 0' not in thisMom: thistitle += ' '+thisMom
    pl.title(thistitle)
    thisdir = outputdir + 'graphs/'+CreateOppDir(thisGamma)
    thisfile = TitlePref.replace(' ','')+thisflag
    if 'q = 0 0 0' not in thisMom: thisdir += '/'+thisMom.replace(' ','').replace('=','')+'/'
    mkdir_p(thisdir)
    return thisdir+thisfile

def SetRFAxies(thisGamma):
    pl.xlabel(RFxlab)
    pl.ylabel(RFylab)
    if thisGamma not in ylimDict.keys():
        pl.ylim(max(pl.ylim()[0],0),min(pl.ylim()[1],2))
    else:
        pl.ylim(ylimDict[thisGamma])
    SetxTicks()
    pl.legend()

def SetFFAxies():
    pl.xlabel(FFxlab)
    pl.ylabel(FFylab)
    pl.legend()


def SetMassAxies():
    pl.xlabel('time')
    pl.ylabel('Mass')
    pl.xlim(MassTVals)
    pl.ylim(Massyrange)
    SetxTicks()
    pl.legend()

def SetLogAxies():
    pl.xlabel('time')
    pl.ylabel('LogC2')
    pl.xlim(15,40)
    pl.ylim(0,20)
    SetxTicks()
    pl.legend()


def PlotTSinkData(data,thisSetList,thisGamma,thisMom,thissm='sm32'):
    PlotCol(data,thisSetList,[thissm],thisGamma,thisMom,'TSink Comparison ')

def PlotTSinkSumData(data,thisSetList,thisGamma,thisMom,thissm='sm32'):
    for ifitr in SumFitRList:
        PlotColSum(data,thisSetList,[thissm],thisGamma,thisMom,'Sum TSink Comparison ',thisTsinkR=ifitr)
        
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
    for thisDt in [2,4,6,8]:
        PlotRFSet(data,thisSetList,MassDt=thisDt)
        pl.savefig(CreateFile('Dt'+str(thisDt),'twopt',thisMom,TitleFlag+' Mass Comparison')+'.pdf')
        pl.clf()
    PlotLogSet(data,thisSetList)
    SetMassAxies()
    pl.savefig(CreateFile('','twopt',thisMom,TitleFlag+' Log Comparison')+'.pdf')
    pl.clf()

def PlotMassSFData(data,thisSetList,thisMom,thisSF='SFCM'):
    for thisDt in [2,4,6,8]:
        PlotMassSetOSF(data,thisSetList,thisDt,thisSF)
        pl.savefig(CreateFile('Dt'+str(thisDt),'twopt',thisMom,'O'+thisSF+' Mass Comparison')+'.pdf')
        pl.clf()
        PlotMassSetTSF(data,thisSetList,thisDt,thisSF)
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
        PlotColOSF(data,data2pt,thisSetList,[thistsink,'PoF'],thisGamma,thisMom,thisSF+' Comparison ',icut,thisSF)
        
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
    if CheckDict(data,'SumMeth',thissm[0]):
        PlotSummedRF(data['SumMeth'][thissm[0]],thisTsinkR)
        pl.savefig(CreateFile(thissm[0],thisGamma,thisMom,TitlePref+outTR)+'Sfun'+'.pdf')
        pl.clf()

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
            PlotRF(data['RF'][iset],thiscolcyc.next(),thissymcyc.next(),thisshiftcyc.next(),LegLab(iset.replace(legrem,'')))
        else:
            dataplot = deepcopy(data['RF'][iset])
            dataplot['Boot'] = MassFun(dataplot['Boot'],MassDt)
            dataplot['tVals'] = dataplot['tVals'][:-MassDt]
            PlotRF(dataplot,thiscolcyc.next(),thissymcyc.next(),thisshiftcyc.next(),LegLab(redset),MP=True)

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
        dataplot['tVals'] = dataplot['tVals'][:-MassDt]
        PlotRF(dataplot,thiscol,thissym,thisshift,LegLab(iset),MP=True,)
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
        dataplot['tVals'] = dataplot['tVals'][:-MassDt]
        PlotRF(dataplot,thiscol,thissym,thisshift,LegLab(iset),MP=True)
        if 'state1' not in iset and CheckDict(data2pt,'T'+thisSF,iset):
            PlotTSFMassLine(data2pt['T'+thisSF][iset],thiscol,iset,MassDt)
    if CheckDict(data2pt,'T'+thisSF,iterSetList[0]):
        PlotTSFMassValue(data2pt['T'+thisSF][iterSetList[0]],MassDt)

def PlotFFs(data,thisCurr,thisSetList,CollName):
    if len(thisSetList) == 0: return
    for iFF in range(1,NoFFPars[thisCurr]+1):
        thisFF = 'FF'+str(iFF)
        PlotFFSet(data,thisFF,thisSetList)
        SetFFAxies()
        pl.savefig(CreateFFFile(CollName,thisCurr,thisFF)+'.pdf')
        pl.clf()
        

def PlotFFSet(dataset,thisFF,thisSetFlag):
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    collist = []
    for thisset in SortMySet(thisSetFlag)[0]:
        ##make legend formatting function
        if not CheckDict(dataset,thisset,thisFF): continue
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
        if 'PoF' in lab and not Log:
            tvals = np.array(data['tVals'])+1+(2*PoFShifts) + shift
        else:
            tvals = np.array(data['tVals'])+1 + shift            
    else:
        tvals = np.array(data['tVals'])
        tvals = tvals-(tvals[-1]+tvals[0])/2.0 + shift
    dataavg = Pullflag(data['Boot'],'Avg')
    dataerr = Pullflag(data['Boot'],'Std')
    pl.errorbar(tvals,map(abs,dataavg),dataerr,color=col,fmt=sym,label=lab)

def PlotSumMeth(data,col,lab,thisTsinkR):
    if not CheckDict(data,SumCutPar,thisTsinkR,'Avg'): return
    if not CheckDict(data,SumCutPar,thisTsinkR,'Std'): return
    dataavg = abs(data[SumCutPar][thisTsinkR]['Avg'])
    dataerr = data[SumCutPar][thisTsinkR]['Std']
    dataup,datadown = dataavg+dataerr,dataavg-dataerr
    pl.axhline(dataavg,color=col,label=lab)
    pl.axhspan(datadown,dataup,facecolor=col,edgecolor='none',alpha=thisalpha)

def PlotSummedRF(data,thisfitr):
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    thisfitr = thisfitr.split()[-1]
    thisfitmin,thisfitmax = thisfitr.split('-')
    for icut,cutdata in data.iteritems():
        if cutdata.keys()[0] == 'nconf': continue
        tdata,dataplot,dataploterr = [],[],[]
        thiscol,thissym = thiscolcyc.next(),thissymcyc.next()
        
        for itsink,tsinkdata in cutdata.iteritems():
            if 'fit' not in itsink and 'nconf' not in itsink:
                tdata.append(int(itsink.replace('tsink',''))-tsource)
                dataplot.append(abs(tsinkdata['Avg']))
                dataploterr.append(tsinkdata['Std'])
        tdata,dataplot,dataploterr = zip(*sorted(zip(tdata,dataplot,dataploterr)))
        pl.errorbar(tdata,dataplot,dataploterr,color=thiscol,fmt=thissym,label=icut)
        if not CheckDict(cutdata,'fit con '+thisfitr,'Boot'): continue
        if not CheckDict(cutdata,'fit sl '+thisfitr,'Boot'): continue        
        parcon,parsl = cutdata['fit con '+thisfitr]['Boot'],cutdata['fit sl '+thisfitr]['Boot']
        fittdata = np.arange(tdata[int(thisfitmin)],tdata[int(thisfitmax)]+incr,incr)
        fittdashed = np.arange(tdata[0],tdata[int(thisfitmin)]+incr,incr)
        fitbootdata = [abs((parsl* (it+tsource)) + parcon) for it in fittdata]
        fitbootdashed = [abs((parsl* (it+tsource)) + parcon) for it in fittdashed]
        GetBootStats(fitbootdata),GetBootStats(fitbootdashed)
        plotup = Pullflag(fitbootdata,'Avg')+Pullflag(fitbootdata,'Std')
        plotdown = Pullflag(fitbootdata,'Avg')-Pullflag(fitbootdata,'Std')
        if len(fitbootdashed) == 0:
            plotdashedup,plotdasheddown = [],[]
        else:
            plotdashedup = Pullflag(fitbootdashed,'Avg')+Pullflag(fitbootdashed,'Std')
            plotdasheddown = Pullflag(fitbootdashed,'Avg')-Pullflag(fitbootdashed,'Std')
        pl.plot(fittdata,map(abs,Pullflag(fitbootdata,'Avg')),
                label='slope='+MakeValAndErr(abs(parsl.Avg),parsl.Std),color=thiscol)
        pl.fill_between(fittdata,plotup,plotdown,color=thiscol,alpha=thisalpha,edgecolor='none')
        pl.plot(fittdashed,plotdashedup,color=thiscol,ls='--')
        pl.plot(fittdashed,plotdasheddown,color=thiscol,ls='--')
    pl.legend(loc='upper left')
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
    pl.plot(tplotdata,map(abs,dataline),color=col)


def PlotOSFValue(data,col,thisshift,OSFcut,smear,thistsink):
    LRM = max((int(thistsink)-tsource)/2.-int(OSFcut.replace('cut','')),0)
    tvals = np.array([-LRM+thisshift,LRM+thisshift])
    osfsmear = RemoveToDt(smear)
    if not CheckDict(data,'B00',OSFfitr[osfsmear],OSFcut,'Avg'): return
    if not CheckDict(data,'B00',OSFfitr[osfsmear],OSFcut,'Std'): return
    dataval = abs(data['B00'][OSFfitr[osfsmear]][OSFcut]['Avg'])
    dataerr = data['B00'][OSFfitr[osfsmear]][OSFcut]['Std']
    dataup,datadown = dataval+dataerr,dataval-dataerr
    pl.fill_between(tvals,[datadown,datadown],[dataup,dataup],facecolor=col,edgecolor='none',alpha=thisalpha)
    pl.plot(tvals,[dataval,dataval],color = col)

def PlotTSFValue(data,col,thisshift,TSFcut,smear,thistsink):
    LRM = max((int(thistsink)-tsource)/2.-int(TSFcut.replace('cut','')),0)
    tvals = np.array([-LRM+thisshift,LRM+thisshift])
    if not CheckDict(data,'B00',TSFfitr,TSFcut,'Avg'): return
    dataAvg,dataErr = abs(data['B00'][TSFfitr][TSFcut]['Avg']),data['B00'][TSFfitr][TSFcut]['Std']
    dataup,datadown = dataAvg-dataErr,dataAvg+dataErr
    pl.fill_between(tvals,[datadown,datadown],[dataup,dataup],facecolor=col,edgecolor='none',alpha=thisalpha)



# PlotO/T SFMassLine only for zero momenta (no point otherwise)

def PlotTSFMassLine(data2pt,col,smear,thisdt):
    pars2pt = []
    for ipar in StateParList['Two']['C2']:
        if not CheckDict(data2pt,ipar,TSFfitr,'Avg'): return
        pars2pt.append(data2pt[ipar][TSFfitr]['Avg'])
    tdata = np.arange(TSFfitvals[0]-tsource,TSFfitvals[1]-thisdt+incr-tsource,incr)
    fit2pt = ff.C2TSFLineFun(tdata,pars2pt)
    fit2ptdt = ff.C2TSFLineFun(tdata+thisdt,pars2pt)
    pl.plot(np.array(tdata)+tsource,map(abs,np.log(fit2ptdt/fit2pt)/thisdt),color=col)


def PlotOSFMassValue(data,col,smear,thisdt):
    smearindex,deltashift = RemoveToDt(smear),0
    if 'sm' not in smear:
        if 'PoF' in smear:
            deltashift = PoFShifts*2
            smearindex = PickedStateStr+'PoF'
        elif 'REvec' in smear:
            deltashift = 0
            smearindex = PickedStateStr+'REvec'
        else:
            deltashift = 0
            smearindex = PickedStateStr+'CM'
    if CheckDict(data,'m0',OSFfitr[smearindex],'Boot'): 
        databoot = data['m0'][OSFfitr[smearindex]]['Boot']
        dataval = abs(databoot.Avg)
        dataup,datadown = dataval+databoot.Std,dataval-databoot.Std
    elif CheckDict(data,'m0',OSFfitr[smearindex],'Avg') and CheckDict(data,'m0',OSFfitr[smearindex],'Std'):
        dataval = data['m0'][OSFfitr[smearindex]]['Avg']
        dataup,datadown = dataval+data['m0'][OSFfitr[smearindex]]['Std'],dataval-data['m0'][OSFfitr[smearindex]]['Std']
    else:
        return
    pl.fill_between([OSFfitvals[smearindex][0]+deltashift,OSFfitvals[smearindex][1]+deltashift-thisdt],
                    [datadown,datadown],[dataup,dataup],facecolor=col,edgecolor='none',alpha=thisalpha)
    pl.plot([OSFfitvals[smearindex][0]+deltashift,OSFfitvals[smearindex][1]+deltashift-thisdt],[dataval,dataval],color=col)

def PlotTSFMassValue(data,thisdt):
    databoot = data['m0'][TSFfitr]['Boot']
    dataval = abs(databoot.Avg)
    dataup,datadown = dataval+databoot.Std,dataval-databoot.Std
    pl.fill_between([TSFfitvals[0],TSFfitvals[1]-thisdt],[datadown,datadown],[dataup,dataup],facecolor='k',edgecolor='none',alpha=thisalpha)


def PlotOSFLog(data,col,smear,norm):
    smearindex = RemoveToDt(smear)
    if 'sum' in smear: smearindex = PickedStateStr+'sum'    
    if not CheckDict(data,'m0',OSFfitr[smearindex],'Boot'): return
    if not CheckDict(data,'Am',OSFfitr[smearindex],'Boot'): return
    parm0 = data['m0'][OSFfitr[smearindex]]['Boot']
    parAm = data['Am'][OSFfitr[smearindex]]['Boot']
    tdata = np.arange(OSFfitvals[smearindex][0],OSFfitvals[smearindex][1]+incr,incr)
    linedata = []
    for it in tdata:
        linedata.append(np.log((parAm*(parm0*(-it+tsource)).exp(1))/norm))
    GetBootStats(linedata)
    dataAvg,dataErr = abs(np.array(Pullflag(linedata,'Avg'))),np.array(Pullflag(linedata,'Std'))
    dataup = dataAvg+dataErr
    datadown = dataAvg-dataErr
    pl.fill_between(tdata,datadown,dataup,facecolor=col,edgecolor='none',alpha=thisalpha)
    pl.plot(OSFfitvals[smearindex],[dataAvg[0],dataAvg[-1]],color=col)


def PlotTSFLog(data,col,smear,norm):
    if not CheckDict(data,'m0',TSFfitr,'Boot'): return
    if not CheckDict(data,'Am',TSFfitr,'Boot'): return
    if not CheckDict(data,'Dm',TSFfitr,'Boot'): return
    if not CheckDict(data,'Amp',TSFfitr,'Boot'): return
    parm0 = data['m0'][TSFfitr]['Boot']
    parAm = data['Am'][TSFfitr]['Boot']
    parDm = data['Dm'][TSFfitr]['Boot']
    parAmp = data['Amp'][TSFfitr]['Boot']
    tdata = np.arange(TSFfitvals[0],TSFfitvals[1]+incr,incr)
    linedata = []
    for it in tdata:
        thisit = it-tsource
        linedata.append(np.log((parAm*(((parm0*(-thisit)).exp(1)) + parAmp*((parm0+parDm)*(-thisit)).exp(1)))/norm))
    GetBootStats(linedata)
    dataAvg,dataErr = abs(np.array(Pullflag(linedata,'Avg'))),np.array(Pullflag(linedata,'Std'))
    dataup = dataAvg+dataErr
    datadown = dataAvg-dataErr
    pl.fill_between(tdata,datadown,dataup,facecolor=col,edgecolor='none',alpha=thisalpha)
    pl.plot(tdata,dataAvg,color=col)

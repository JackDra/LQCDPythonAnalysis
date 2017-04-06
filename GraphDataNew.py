#!/usr/bin/env python

import matplotlib
matplotlib.use('PDF') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as pl
from Params import *
from FitFunctions import *
from LLSBoot import FitBoots
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
import BootTest as bt
from AutoCorr import uWerrMine


##FORCE TITLE PARAMETER, SET TO FALSE TO USE NORMAL TITLES#

ForceTitle = False
TitleShift = 1.04
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
MassTVals = 1,21
Massyrange = .5,0.9

# if kappa == 1375400:
#     Massyrange = .5,0.7
# elif kappa == 1370000:
#     Massyrange = .6,0.9
# else:
#     Massyrange = 0,3
    
# Massyrange = 0.40,0.60
# Massyrange = 0.40,0.61
# Qsqrdxlim = -0.03,1
Qsqrdxlim = -0.03,0.8
    
ErrTVals = 1,21
Erryrange = 0,3

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
            'ProtonP3g2':[0.1,0.5], ## Tsink Var 
            'NeutronP3g2':[-0.4,0.0], ## Tsink Var 
            'IsoVectorP4g4':[0.8,1.1]}

leglocDict = {'ProtonP4g4':'upper left'}


# ylimFFDict = {'ProtonGeGmFF1/F1divF2':[0.3,0.5],
#               # 'NeutronGeGmFF2':[-1.2,-0.4],
#               'ProtonVectorFF2':[0.4,1.2],
#               'NeutronVectorFF2':[-1.2,-0.4],
#               'IsoVectorPsVectorFF2':[2,10],
#               'ProtonTensorFF2':[-0.5,-3],
#               'ProtonTensorFF3':[0.9,0.15],
#               'NeutronTensorFF2':[2.0,0.45],
#               'NeutronTensorFF3':[-0.2,-1],
#               'NeutronGeGmFF1/F1divF2':[0.0,0.06]}
# ylimFFDict = {'ProtonVectorTopFF3':[-0.15,0.25],'NeutronVectorTopFF3':[-0.15,0.25],'PandNVectorTopFF3':[-0.15,0.25],
#               'ProtonVectorWeinFF3':[-20,15]}
# ylimFFDict = {'ProtonVectorTopFF3':[-0.05,0.35],'NeutronVectorTopFF3':[-0.25,0],'PandNVectorTopFF3':[-.2,.3],
#               'ProtonVectorWeinFF3':[-20,15]}
ylimFFDict = {'ProtonGeGmFF1':[0.4,1.1],
              'ProtonGeGmFF2':[1.0,2.9],
              'NeutronGeGmFF1':[-0.002,0.018],
              'NeutronGeGmFF2':[-2.0,-0.7]}

leglengthFFDict = {'ProtonVectorTopFF3':2,'NeutronVectorTopFF3':2,'PandNVectorTopFF3':2,
                   'ProtonVectorWeinFF3':2,'NeutronVectorWeinFF3':2,'PandNVectorWeinFF3':1}

leglocFFDict = {'NeutronVectorFF2':'upper left',
                'NeutronGeGmFF1':'upper left',
                'NeutronGeGmFF2':'lower right',
                'IsoVectorPsScalarFF1':'lower right',
                'IsoVectorPsVectorFF2':'upper right',
                'NeutronGeGmFF1/F1divF2':'upper left',
                'ProtonTensorFF3':'upper left',
                'NeutronVectorTopFF3':'lower right',
                # 'NeutronVectorTopFF3':'upper right',
                'NeutronTensorFF1':'upper left',
                'NeutronTensorFF2':'upper left'}

params = {'legend.fontsize': 10,
          'legend.numpoints': 1,
          'axes.labelsize' : 20,
          'figure.autolayout': True,
          'axes.grid': True,
          'axes.xmargin':0.01,
          'axes.titlesize' : 20,
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


    # if 'FF1' in thisFF :
    #     if 'ProtonGeGm' in thisDSCurr:
    #         PlotExp('ChargeRadius',thiscolcyc)
    #     elif 'ProtonIsoVector' in thisDSCurr:
    #         PlotExp('AxialChargeRadius',thiscolcyc)
    # elif 'FF2' in thisFF and 'GeGm' in thisDSCurr:
    #     if 'Neutron' in thisDSCurr:
    #         PlotExp('NeutronMagMom',thiscolcyc)
    #     elif 'Proton' in thisDSCurr:
    #         PlotExp('ProtonMagMom',thiscolcyc)
    # if 'FF3' in thisFF and 'Top' in thisDSCurr:
    #     if 'Neutron' in thisDSCurr:
    #         PlotExp('ProtonEDM',thiscolcyc)
    #     elif 'Proton' in thisDSCurr:
    #         PlotExp('NeutronEDM',thiscolcyc)

CRadConv = 1/(6*hbarc**2)
ExpxvalsCrad = np.array([0,0.05])
    
def PlotExp(flag,thiscolcyc):
    if 'ChargeRadius' in flag:
        yvals = 1-CRadConv*ExpValues['CRad_muon'][0]*ExpxvalsCrad
        yup,ydown = 1-CRadConv*np.sum(ExpValues['CRad_muon'])*ExpxvalsCrad,1-CRadConv*(ExpValues['CRad_muon'][0]-ExpValues['CRad_muon'][1])*ExpxvalsCrad
        thiscol = thiscolcyc.next()
        pl.plot(ExpxvalsCrad,yvals,color=thiscol,label=r'ANTOGNINI 13, $\mu p-atom$ Lamb shift $\langle r^2 \rangle='+MakeValAndErr(*np.abs(ExpValues['CRad_muon']))+'\ fm^{2}$')
        pl.fill_between(ExpxvalsCrad,yup,ydown,color=thiscol,alpha=thisalpha,edgecolor='none')
        thiscol = thiscolcyc.next()
        yvals = 1-CRadConv*ExpValues['CRad_electron'][0]*ExpxvalsCrad
        yup,ydown = 1-CRadConv*np.sum(ExpValues['CRad_electron'])*ExpxvalsCrad,1-CRadConv*(ExpValues['CRad_electron'][0]-ExpValues['CRad_electron'][1])*ExpxvalsCrad
        pl.plot(ExpxvalsCrad,yvals,color=thiscol,label=r'MOHR 12, 2010 CODATA $e p$ data $\langle r^2 \rangle='+MakeValAndErr(*np.abs(ExpValues['CRad_electron']))+'\ fm^{2}$')
        pl.fill_between(ExpxvalsCrad,yup,ydown,color=thiscol,alpha=thisalpha,edgecolor='none')
    # if 'NeutronCRad' in flag:
    #     ExpxvalsCrad = np.array([0,0.1])
    #     yvals = CRadConv*ExpValues['NeutronCRad'][0]*ExpxvalsCrad
    #     yup,ydown = CRadConv*np.sum(ExpValues['NeutronCRad'])*ExpxvalsCrad,CRadConv*(ExpValues['NeutronCRad'][0]-ExpValues['NeutronCRad'][1])*ExpxvalsCrad
    #     thiscol = thiscolcyc.next()
    #     pl.plot(ExpxvalsCrad,yvals,color=thiscol,label=r'$ne$ scattering , pdg average $\langle r_{n}^2 \rangle='+MakeValAndErr(*np.abs(ExpValues['NeutronCRad']))+'\ fm^{2}$')
    #     pl.fill_between(ExpxvalsCrad,yup,ydown,color=thiscol,alpha=thisalpha,edgecolor='none')
    if 'NeutronMagMom' in flag:
        pl.errorbar([0.0],[ExpValues['MagMomNeutron'][0]],[ExpValues['MagMomNeutron'][1]],fmt='x',color=thiscolcyc.next(),label='2010 CODATA, $\mu_{n}='+MakeValAndErr(*ExpValues['MagMomNeutron'])+'$')
    if 'ProtonMagMom' in flag:
        pl.errorbar([0.0],[ExpValues['MagMomProton'][0]],[ExpValues['MagMomProton'][1]],fmt='x',color=thiscolcyc.next(),label='2010 CODATA, $\mu_{p}='+MakeValAndErr(*ExpValues['MagMomProton'])+'$')
    if 'MRadProton' in flag:
        yvals = ExpValues['MagMomProton'][0]-CRadConv*ExpValues['MRadProton'][0]*ExpxvalsCrad*ExpValues['MagMomProton'][0]
        yup,ydown = (ExpValues['MagMomProton'][0]-CRadConv*np.sum(ExpValues['MRadProton'])*ExpxvalsCrad*ExpValues['MagMomProton'][0],
                     ExpValues['MagMomProton'][0]-CRadConv*(ExpValues['MRadProton'][0]-ExpValues['MRadProton'][1])*ExpxvalsCrad*ExpValues['MagMomProton'][0])
        thiscol = thiscolcyc.next()
        pl.plot(ExpxvalsCrad,yvals,color=thiscol,label=r'BELUSHKIN 07 Dispersion Analysis $\langle r_{\mu}^2 \rangle='+MakeValAndErr(*np.abs(ExpValues['MRadProton']))+'\ fm^{2}$')
        pl.fill_between(ExpxvalsCrad,yup,ydown,color=thiscol,alpha=thisalpha,edgecolor='none')
    if 'MRadNeutron' in flag:
        yvals = ExpValues['MagMomNeutron'][0]-CRadConv*ExpValues['MRadNeutron'][0]*ExpxvalsCrad*np.abs(ExpValues['MagMomNeutron'][0])
        yup,ydown = (ExpValues['MagMomNeutron'][0]-CRadConv*np.sum(ExpValues['MRadNeutron'])*ExpxvalsCrad*np.abs(ExpValues['MagMomNeutron'][0]),
                     ExpValues['MagMomNeutron'][0]-CRadConv*(ExpValues['MRadNeutron'][0]-ExpValues['MRadNeutron'][1])*ExpxvalsCrad*np.abs(ExpValues['MagMomNeutron'][0]))
        thiscol = thiscolcyc.next()
        thisleg = r'BELUSHKIN 07 Dispersion Analysis $\langle r_{\mu}^2 \rangle='+MakeValAndErr(*np.abs(ExpValues['MRadNeutron']))+'\ fm^{2}$'
        pl.plot(ExpxvalsCrad,yvals,color=thiscol,label=thisleg)
        pl.fill_between(ExpxvalsCrad,yup,ydown,color=thiscol,alpha=thisalpha,edgecolor='none')
    # if 'ProtonEDM' in flag or 'PandNEDM' in flag:
    #     thiscol = thiscolcyc.next()
    #     thislab = r'$ d_{p} < '+r'{:0.2e} \theta fm$'.format(ExpValues['ProtonEDMtfm'])
    #     yval = np.array([ExpValues['ProtonEDMtfm'],ExpValues['ProtonEDMtfm']])
    #     pl.plot([0.0],[0.0],color=thiscol,label=thislab)
    #     pl.fill_between([-100,100],yval,-yval,color=thiscol,alpha=thisalpha,edgecolor='none')
    if 'NeutronEDM' in flag or 'PandNEDM' in flag:
        thiscol = thiscolcyc.next()
        thislab = r'$ d_{n} < '+r'{:0.2e}\ \theta fm$'.format(ExpValues['NeutronEDMtfm'])
        yval = np.array([ExpValues['NeutronEDMtfm'],ExpValues['NeutronEDMtfm']])
        pl.plot([0.0],[0.0],color=thiscol,label=thislab)
        pl.fill_between([-0.05,0.05],yval,-yval,color=thiscol,alpha=thisalpha,edgecolor='none')




        
def TflowToPhys(tflowlist):
    if isinstance(tflowlist, list):
        return [np.sqrt(8*iflow)*latspace for iflow in tflowlist] ## in fermi
    elif isinstance(tflowlist,str):
        return np.sqrt(8*float(tflowlist))*latspace
    else:
        return np.sqrt(8*tflowlist)*latspace
        
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
    return [itertools.cycle(markerset),itertools.cycle(colourset8),itertools.cycle(shiftsetff)]

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
    # thistitle = thisCol
    thistitle = TitleFixFF(thisCurr.replace('/',' '),thisFF)
    if ForceTitle == False:
        pl.title(thistitle,y=TitleShift)
    else:
        pl.title(ForceTitle,y=TitleShift)
    thisdir = outputdir[0] + 'graphs/FormFactors/'+thisCurr + '/'
    mkdir_p(thisdir)
    thisfile = thisCol+thisCurr.replace('/','') + thisFF
    return thisdir + thisfile


def CreateMKFFFile(thisCol,thisCurr,thisFF):
    thistitle = thisCol + TitleFixFF(thisCurr.replace('/',' '),thisFF)
    if ForceTitle == False:
        pl.title(thistitle,y=TitleShift)
    else:
        pl.title(ForceTitle,y=TitleShift)
    thisdir = outputdir[0].replace(str(kappa),'Comb') + 'graphs/FormFactors/'+thisCurr + '/'
    mkdir_p(thisdir)
    thisfile = thisCol+thisCurr.replace('/','') + thisFF
    return thisdir + thisfile



def CreateFile(thisflag,thisGamma,thisMom,TitlePref,thisfig=False,subdir=False,NoMpi=False):
    thisMpi = '$\ '+GetMpi(kappa,Phys=True)+'$'
    if NoMpi: thisMpi = ''
    if 'twopt' in thisGamma:
        if 'Dt' in thisflag:
            thistitle = thisGamma+' '+TitlePref+' $' + thisflag.replace('Dt','\Delta t=') + r'$'
        else:
            thistitle = thisGamma+' '+TitlePref+' '+thisflag
        thistitle = thistitle.replace('twopt ','')
    else:
        thistitle = thisGamma+' '+TitlePref+' '+thisflag
    if 'q = 0 0 0' not in thisMom: thistitle += ' '+thisMom
    if ForceTitle == False:
        if thisfig == False:
            pl.title(thistitle+thisMpi,y=TitleShift)
        else:
            thisfig.suptitle(thistitle, fontsize=20)
    else:
        # pl.title(ForceTitle+'$' + thisflag.replace('Dt','\Delta t') + '$')
        if thisfig == False:
            if Debug: print 'title: ' , ForceTitle
            pl.title(ForceTitle,y=TitleShift)
        else:
            thisfig.suptitle(ForceTitle, fontsize=20)
    thisdir = outputdir[0] + 'graphs/'+CreateOppDir(thisGamma)
    thisfile = TitlePref.replace(' ','')+thisflag
    thisdir += MakeMomDir(thisMom)
    if subdir !=False: thisdir += '/'+subdir+'/'
    mkdir_p(thisdir)
    return thisdir+thisfile

def SetSumFunAxies(DoY):
    pl.ylabel(SUMylab)
    if DoY:pl.xlabel(SUMxlab)

def SetRFAxies(thisGamma):
    pl.xlabel(RFxlab)
    pl.ylabel(RFylab)
    # if Debug: print 'Hardcoding yaxis limits',thisGamma, ylimDict.keys()
    if thisGamma not in ylimDict.keys():
        if 'Wein' not in thisGamma:
            pl.ylim(max(pl.ylim()[0],-5),min(pl.ylim()[1],5))
    else:
        # pl.ylim(max(ylimDict[thisGamma][0],pl.ylim()[0]),min(ylimDict[thisGamma][1],pl.ylim()[1]))
        pl.ylim(*ylimDict[thisGamma])
    SetxTicks()
    if thisGamma not in leglocDict.keys():
        pl.legend()
    else:
        pl.legend(loc=leglocDict[thisGamma])
    pl.tight_layout()
    
def SetFFAxies(thisCurr):
    currnumb = thisCurr[-1]
    pl.xlabel(FFxlab)
    if currnumb == '3':
        pl.ylabel(FFylab.replace('F','\\frac{F_{'+currnumb+'}').replace(') $',')}{2m_{N,phys}} (\theta fm) $'))
    else:
        if 'GeGm' in thisCurr:
            if currnumb == '1':
                currnumb = 'E'
            elif currnumb == '2':
                currnumb = 'M'
            pl.ylabel(FFylab.replace('F','G_{'+currnumb+'}'))
        else:
            pl.ylabel(FFylab.replace('F','F_{'+currnumb+'}'))
    if thisCurr in ylimFFDict.keys():
        pl.ylim(ylimFFDict[thisCurr])
    else:
        pl.ylim(pl.ylim()[0],pl.ylim()[1]*1.2)
    pl.xlim(*Qsqrdxlim)
    nlegline = 1
    if thisCurr in leglengthFFDict.keys(): nlegline = leglengthFFDict[thisCurr]
    if thisCurr not in leglocFFDict.keys():
        pl.legend(ncol=nlegline)
    else:
        pl.legend(loc=leglocFFDict[thisCurr],ncol=nlegline)
    # pl.tight_layout()


def SetMassAxies():
    pl.xlabel(r'$t$')
    if PhysMassPlot:
        pl.ylabel(r'$M_{N} GeV$')
        pl.ylim(Massyrange[0]*hbarcdivlat,Massyrange[1]*hbarcdivlat)
    else:
        pl.ylabel(r'$aM_{N}$')
        pl.ylim(Massyrange)
    pl.xlim(MassTVals)
    SetxTicks()
    pl.legend()
    pl.tight_layout()

def SetDispAxies():
    pl.xlabel(r'$\vec{p}^2$')
    pl.ylabel(r'$E_{\vec{p}}$')
    pl.legend(loc='upper left')
    pl.tight_layout()


AlphaTflowList = np.arange(0.01,10,1)
AlphaTlist = np.arange(3,15)
def SetTopAxies(torflow,NNQ=False,Dt=2,Wein=False,MpiTitle=True):
    if torflow == 't':
        pl.xlabel(r'$t$')
        pl.xlim(MassTVals[0],AlphaTlist[-1])
        SetxTicks()
    elif torflow == 'flow':
        pl.xlabel(r'$\sqrt{8t_{f}} fm$')        
        # pl.xlim(-0.1,10)
    elif torflow == 'Mpi':
        pl.xlabel(r'$m_{\pi} (GeV)$')
        pl.xlim(0,pl.xlim()[-1])
    pl.legend()
    pl.tight_layout()
    if MpiTitle: thisMpi = '\ '+GetMpi(kappa,Phys=True)
    else: thisMpi = ''
    if NNQ:
        pl.ylabel(r'$ Eff\ Mass $')
        if Wein:
            pl.title(r'$ Eff\ Mass\ of\ \gamma_{5}P_{+}NNW\ and\ NN\ Dt='+str(Dt)+thisMpi+r'$',y=TitleShift)
        else:
            pl.title(r'$ Eff\ Mass\ of\ \gamma_{5}P_{+}NNQ\ and\ NN\ Dt='+str(Dt)+thisMpi+r'$',y=TitleShift)
    else:
        if Wein:
            pl.ylabel(r'$ \alpha_{W} $')
            pl.title(r'$\frac{\langle NNW \rangle }{\langle NN\rangle} = \alpha_{W} '+thisMpi+r'$',y=TitleShift)
        else:
            pl.ylabel(r'$ \alpha_{Q} $')
            pl.title(r'$\frac{\langle NNQ \rangle }{\langle NN\rangle} = \alpha_{Q} '+thisMpi+r'$',y=TitleShift)
    
    

def SetErrAxies():
    pl.xlabel(r'$t$')
    pl.ylabel(r'$\Delta C/ \Delta C$')
    pl.xlim(ErrTVals)
    pl.ylim(Erryrange)
    SetxTicks()
    pl.legend()
    pl.tight_layout()

def SetLogAxies():
    pl.xlabel(r'$t$')
    pl.ylabel(r'$log(G_{2})$')
    pl.xlim(1.9,15)
    pl.ylim(-8,1.5)
    # pl.xlim(0,6)
    # pl.ylim(-6,0)
    # pl.xlim(11,20)
    # pl.ylim(-12,-8)
    SetxTicks()
    pl.legend()
    pl.tight_layout()


VsChi = False
def SetChiAxies():
    if VsChi:
        pl.xlabel(r'$\chi^{2}$')
        pl.xlim(0,10)
        pl.ylim(0.55,0.65)
    else:
        pl.xlabel(r'$t_{min}$')
        # pl.xlim(3,10)
    pl.ylabel(r'$aM_{N}$')
    # pl.xlim(0,6)
    # pl.ylim(-6,0)
    # pl.xlim(11,20)
    # pl.ylim(-12,-8)
    # SetxTicks()
    pl.legend(loc='lower right')
    pl.tight_layout()

Stacked = False
Normed = False
# HistType = 'bar'
# HistType = 'barstacked'
HistType = 'step'
# HistType = 'stepfilled'
FitMaxCutoff = 30
FitMinCutoff = 0
ChiThreshList = np.arange(0.1,2,0.1)
# ChiThreshList = [0.4]
binwidth = 0.0005
Histx = {'q = 0 0 0' : (0.54,0.62),
         'q = 1 0 0' : (0.56,0.66),
         'q = 1 1 0' : (0.6,0.7),
         'q = 1 1 1' : (0.64,0.74),
         'q = 2 0 0' : (0.64,0.76)}
         

BinList = np.arange(0.5,1.0+binwidth,binwidth)
def SetHistAxies(imom):
    pl.xlabel(r'$aM_{N}$')
    pl.ylabel(r'$Frequency$')
    pl.xlim(Histx[imom][0],Histx[imom][1])
    # SetxTicks()
    pl.ylim(0,pl.ylim()[1])
    pl.legend()
    pl.tight_layout()

def SetAlphaAxies():
    pl.xlabel(r'$\chi^{2}/DoF$')
    pl.ylabel(r'$\alpha$')
    pl.ylim(0,1)
    # pl.ylim(0,0.25)
    # pl.xlim(0,6)
    # pl.ylim(-6,0)
    # pl.xlim(11,20)
    # pl.ylim(-12,-8)
    # SetxTicks()
    pl.legend(loc='lower right')
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
        
def PlotCMData(data,thisSetList,thisGamma,thisMom,FT,thistsink=CMTSinkStrList[0]):
    global ForceTitle
    ForceTitle = FT
    # PlotCol(data,thisSetList,[thistsink,'PoF'],thisGamma,thisMom,'Variational Comparison ')
    PlotCol(data,thisSetList,[thistsink,'PoF'],thisGamma,thisMom,'Var Comp ')


def PlotMassData(data,massdata,thisSetList,thisMom,FT,TitleFlag=''):
    global ForceTitle
    ForceTitle = FT
    for thisDt in MassDtList:
        PlotMassSet(data,massdata,thisSetList,thisMom,MassDt=thisDt)
        SetMassAxies()        
        pl.savefig(CreateFile('Dt'+str(thisDt),'twopt',thisMom,TitleFlag+' Mass Comparison')+'.pdf')
        pl.clf()
    PlotLogSet(data,thisSetList)
    SetLogAxies()
    pl.savefig(CreateFile('','twopt',thisMom,TitleFlag+' Log Comparison')+'.pdf')
    pl.clf()

def PlotMassErrComp(dataList,kappaflags,thisSetList,thisMom,FT,TitleFlag=''):
    global ForceTitle
    ForceTitle = FT
    if 'xsrc1k'+str(kappa) in kappaflags:
        PlotErrComp(dataList[0],dataList[kappaflags.index('xsrc1k'+str(kappa))],thisSetList)
        SetErrAxies()
        outfile = CreateFile('','twopt',thisMom,TitleFlag+' C5b1 div C1')+'.pdf'
        if Debug: print 'creating: ' , outfile
        pl.savefig(outfile)
        pl.clf()
    if 'nboot1kk'+str(kappa) in kappaflags:
        PlotErrComp(dataList[kappaflags.index('nboot1kk'+str(kappa))],dataList[0],thisSetList)
        SetErrAxies()
        outfile = CreateFile('','twopt',thisMom,TitleFlag+' C5b2 div C5b1')+'.pdf'
        if Debug: print 'creating: ' , outfile
        pl.savefig(outfile)
        pl.clf()
    if 'XAvgk'+str(kappa) in kappaflags and 'xsrc1k'+str(kappa) in kappaflags:
        PlotErrComp(dataList[kappaflags.index('XAvgk'+str(kappa))],dataList[kappaflags.index('xsrc1k'+str(kappa))],thisSetList)
        SetErrAxies()
        outfile = CreateFile('','twopt',thisMom,TitleFlag+' C5av div C1')+'.pdf'
        if Debug: print 'creating: ' , outfile
        pl.savefig(outfile)
        pl.clf()
    if 'nboot1kk'+str(kappa) in kappaflags and 'XAvgk'+str(kappa) in kappaflags :
        PlotErrComp(dataList[kappaflags.index('nboot1kk'+str(kappa))],dataList[kappaflags.index('XAvgk'+str(kappa))],thisSetList)
        SetErrAxies()
        outfile = CreateFile('','twopt',thisMom,TitleFlag+' C5b2 div C5av')+'.pdf'
        if Debug: print 'creating: ' , outfile
        pl.savefig(outfile)
        pl.clf()

def PlotDispersion(data,thisset,FT):
    global ForceTitle
    ForceTitle = FT
    plotdata,dispdata,momlist = [],[],[]
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    for iMom,momdata in data.iteritems():
        dispdata.append([])
        OSFKey,dump = CreateOSFfitKey(thisset)
        # print momdata.keys(), 'OSFCM'
        # print momdata['OSFCM'].keys(),thisset
        # print momdata['OSFCM'][thisset].keys(),'m0'
        # print momdata['OSFCM'][thisset]['m0'].keys(),OSFfitrMom[iMom][OSFKey]
        # print momdata['OSFCM'][thisset]['m0'][OSFfitrMom[iMom][OSFKey]].keys(),'Boot'

        # print data['q = 0 0 0'].keys(), 'OSFCM'
        # print data['q = 0 0 0']['OSFCM'].keys(),thisset
        # print data['q = 0 0 0']['OSFCM'][thisset].keys(),'m0'
        # print data['q = 0 0 0']['OSFCM'][thisset]['m0'].keys(),OSFfitrMom[iMom][OSFKey]
        # print data['q = 0 0 0']['OSFCM'][thisset]['m0'][OSFfitrMom[iMom][OSFKey]].keys(),'Boot'
        # print 
        if (CheckDict(momdata,'OSFCM',thisset,'m0',OSFfitrMom[iMom][OSFKey],'Boot') and
            CheckDict(data,'q = 0 0 0','OSFCM',thisset,'m0',OSFfitrMom['q = 0 0 0'][OSFKey],'Boot')):
            momlist.append(qsqrdstr(iMom)*(qunit**2))
            plotdata.append(momdata['OSFCM'][thisset]['m0'][OSFfitrMom[iMom][OSFKey]]['Boot'])
            for iDisp in  LatDispList:                
                dispdata[-1].append(ScaledEffMass(iMom,[data['q = 0 0 0']['OSFCM'][thisset]['m0'][OSFfitrMom['q = 0 0 0'][OSFKey]]['Boot']],DispIn=iDisp)[0])
    if len(momlist) == 0: return
    momlist,dispdata,plotdata  = zip(*sorted(zip(np.array(momlist),dispdata,plotdata)))
    dispdata = np.swapaxes(np.array(dispdata),0,1)
    pl.errorbar(np.array(momlist)+0.02*thisshiftcyc.next(),Pullflag(plotdata,'Avg'),Pullflag(plotdata,'Std'),
                color=thiscolcyc.next(),fmt=thissymcyc.next(),label='Lat Results')

    for keydisp,idispdata in zip(DispKeyList,dispdata):
        pl.errorbar(np.array(momlist)+0.02*thisshiftcyc.next(),Pullflag(idispdata,'Avg'),Pullflag(idispdata,'Std'),
                    color=thiscolcyc.next(),fmt=thissymcyc.next(),label=LegLab(keydisp))
    SetDispAxies()
    pl.savefig(CreateFile('','twopt','q = 0 0 0','OSF Dispersion '+thisset,subdir='Disp')+'.pdf')
    pl.clf()
    
        
def PlotMassSFData(data,massdata,thisSetList,thisMom,FT,thisSF='SFCM'):
    global ForceTitle
    ForceTitle = FT
    for thisDt in MassDtList:
        PlotMassSetOSF(data,massdata,thisSetList,thisDt,thisSF,thisMom)
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
    for iFitMax in FitMaxList:
        PlotSetChiDist(data,thisSetList,'O'+thisSF,iFitMax)
        SetChiAxies()
        pl.savefig(CreateFile('','twopt',thisMom,'O'+thisSF+' ChiDist tmax'+str(iFitMax),subdir='ChiDists')+'.pdf')
        pl.clf()
        PlotSetChiDist(data,thisSetList,'T'+thisSF,iFitMax)
        SetChiAxies()
        pl.savefig(CreateFile('','twopt',thisMom,'T'+thisSF+' ChiDist tmax'+str(iFitMax),subdir='ChiDists')+'.pdf')
        pl.clf()

    for iThresh in ChiThreshList:
        if thisMom != 'q = 0 0 0':
            PlotSetHistDist(data,thisSetList,'O'+thisSF,iThresh,thisMom,Zmom=massdata)
        else:
            PlotSetHistDist(data,thisSetList,'O'+thisSF,iThresh,thisMom)
        SetHistAxies(thisMom)
        pl.savefig(CreateFile('','twopt',thisMom,'O'+thisSF+' HistDist Thresh'+str(iThresh),subdir='Hists')+'.pdf')
        pl.clf()
        PlotSetHistDist(data,thisSetList,'T'+thisSF,iThresh,thisMom)
        SetHistAxies(thisMom)
        pl.savefig(CreateFile('','twopt',thisMom,'T'+thisSF+' HistDist Thresh'+str(iThresh),subdir='Hists')+'.pdf')
        pl.clf()

    
    # iterSetList = SortMySet(ReduceTooMassSet(thisSetList))[0]
    # for iset in iterSetList:
    #     if not CheckDict(data,'O'+thisSF,iset):
    #         if Debug:
    #             print data.keys(), thisSF
    #             print data[thisSF].keys(), iset
    #         continue
    #     PlotAlphaDist(data['O'+thisSF][iset])
    #     SetAlphaAxies()
    #     pl.savefig(CreateFile('','twopt',thisMom,'O'+thisSF+' AlphaDist Set '+iset,subdir='AlphaDists')+'.pdf')
    #     pl.clf()
    #     return 
        

        
def PlotCMOSFData(data,data2pt,thisSetList,thisGamma,thisMom,FT,thistsink=CMTSinkStrList[0],thisSF='OSFCM'):
    global ForceTitle
    ForceTitle = FT
    for icut in OSFCutList:
        PlotColOSF(data,data2pt,thisSetList,[thistsink,'PoF','CM'],thisGamma,thisMom,thisSF+' Comparison ',icut,thisSF)
        
def PlotCMTSFData(data,data2pt,thisSetList,thisGamma,thisMom,FT,thistsink=CMTSinkStrList[0],thisSF='TSFCM'):
    global ForceTitle
    ForceTitle = FT
    for icut in TSFCutList:
        PlotColTSF(data,data2pt,thisSetList,[thistsink],thisGamma,thisMom,thisSF+' Comparison ',icut,thisSF)



def PlotCol(data,thisSetList,thisflag,thisGamma,thisMom,TitlePref):
    if thisflag[0] == 's':
        thislegrem = 'state1'
    else:
        thislegrem = thisflag[0]
    Top = 'Top' in thisGamma
    if 'Wein' in thisGamma: Top = 'Wein'
    PlotRFSet(data,SiftAndSort(thisSetList,thisflag,nocm=False),legrem=thislegrem,Top=Top )
    SetRFAxies(thisGamma)
    thisfile = CreateFile(thisflag[0],thisGamma,thisMom,TitlePref)
    if Debug: print 'printing: ' +thisfile+'.pdf'
    pl.savefig(thisfile+'.pdf')
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
        twoptsm = 'tsrc'+str(tsource)+thissm
        thiscol = thiscolcyc.next()
        thissym = thissymcyc.next()
        thisshift = thisshiftcyc.next()
        if not CheckDict(data,'RF'+thisTSF,iset): continue
        PlotRF(data['RF'+thisTSF][iset],thiscol,thissym,thisshift,LegLab(iset.replace(legrem,'')))
        # if Debug: print data2pt.keys(), thisTSF
        # if  CheckDict(data2pt,thisTSF) and Debug: print data2pt[thisTSF].keys(), thissm
        if not CheckDict(data2pt,thisTSF,twoptsm): continue
        PlotTSFLine(data[thisTSF][thissm],data2pt[thisTSF][twoptsm],thistsink.replace('tsink',''),thiscol,thisshift,TSFcut,thissm)
        if 'CM' in thisTSF: PlotTSFValue(data[thisTSF][thissm],thiscol,thisshift,TSFcut,twoptsm,thistsink.replace('tsink','') )
    if 'CM' not in thisTSF:
        if CheckDict(data,thisTSF,thissm) and CheckDict(data2pt,thisTSF,twoptsm):
            PlotTSFValue(data[thisTSF][thissm],thiscol,thisshift,TSFcut,thissm,thistsink.replace('tsink','') )

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
        if not CheckDict(data,'RF'+thisOSF,iset):
            # if Debug: print data.keys(), 'RF'+thisOSF
            # if Debug: print iset
            continue
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



def PlotRFSet(data,thisSetList,legrem='',MassDt = False,Top=False):
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    if MassDt == False:
        iterSetList = SortMySet(thisSetList)[0]
    else:
        iterSetList = SortMySet(ReduceTooMassSet(thisSetList))[0]
    for iset in iterSetList:
        if not CheckDict(data,'RF',iset): continue
        if MassDt == False:
            thiscol,thisshift = thiscolcyc.next(),thisshiftcyc.next()
            if Top:                
                if 'Wein' == Top:
                    thistflow = 't_flow8.0'
                else:
                    thistflow = 't_flow4.01'
                thistsink = data['RF'][iset][thistflow]['tVals'][-1]
                PlotRF(data['RF'][iset][thistflow],thiscol,thissymcyc.next(),thisshift,LegLab(iset+thistflow.replace(legrem,'')))
                if CheckDict(data,'Fits',iset,thistflow):
                    PlotFit(data['Fits'][iset][thistflow],thiscol,thisshift,iset+thistflow,thistsink)
            else:
                thistsink = data['RF'][iset]['tVals'][-1]
                PlotRF(data['RF'][iset],thiscol,thissymcyc.next(),thisshift,LegLab(iset.replace(legrem,'')))
                if CheckDict(data,'Fits',iset):
                    PlotFit(data['Fits'][iset],thiscol,thisshift,iset,thistsink)
        else:
            dataplot = deepcopy(data['RF'][iset])
            dataplot['Boot'] = MassFun(dataplot['Boot'],MassDt)
            # dataplot['tVals'] = dataplot['tVals'][:-MassDt] 
            # dataplot['tVals'] = dataplot['tVals'][MassDt:] 
            PlotRF(dataplot,thiscolcyc.next(),thissymcyc.next(),thisshiftcyc.next(),LegLab(iset),MP=True)

def PlotMassSet(data,massdata,thisSetList,thisMom,legrem='',MassDt = 1):
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    iterSetList = SortMySet(ReduceTooMassSet(thisSetList))[0]
    for iset in iterSetList:
        if not CheckDict(data,'RF',iset): continue
        dataplot = deepcopy(data['RF'][iset])
        dataplot['Boot'] = MassFun(dataplot['Boot'],MassDt)
        thiscol,thissym,thisshift = thiscolcyc.next(),thissymcyc.next(),thisshiftcyc.next()
        # dataplot['tVals'] = dataplot['tVals'][:-MassDt] 
        # dataplot['tVals'] = dataplot['tVals'][MassDt:] 
        PlotRF(dataplot,thiscol,thissym,thisshift,LegLab(iset),MP=True)
        if thisMom != 'q = 0 0 0':
            massdataplot = deepcopy(massdata['RF'][iset])
            if not CheckDict(massdata,'RF',iset): continue
            massdataplot['Boot'] = MassFun(massdataplot['Boot'],MassDt)
            massdataplot['Boot'] = ScaledEffMass(thisMom,massdataplot['Boot'])
            PlotRF(massdataplot,thiscol,thissym,thisshift,LegLab(iset+'\ Dis'),MP=True,alpha=0.5)
        

    
def PlotLogSet(data,thisSetList,legrem=''):
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    iterSetList = SortMySet(ReduceTooMassSet(thisSetList))[0]
    for iset in iterSetList:
        if not CheckDict(data,'RF',iset): continue
        dataplot = deepcopy(data['RF'][iset])
        norm,DoPoFS = GetNorm(dataplot['Boot'],tsource,iset)
        dataplot['Boot'] = np.log([tboot/norm for tboot in np.roll(dataplot['Boot'],-DoPoFS,0)])
        dataplot['Boot'] = GetBootStats(dataplot['Boot'])
        PlotRF(dataplot,thiscolcyc.next(),thissymcyc.next(),thisshiftcyc.next(),LegLab(iset.replace(legrem,'')),MP=True,Log=True)


            
def PlotErrComp(dataTop,dataBot,thisSetList,legrem=''):
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    iterSetList = SortMySet(ReduceTooMassSet(thisSetList))[0]
    for iset in iterSetList:
        if not CheckDict(dataBot,'RF',iset): continue
        if not CheckDict(dataTop,'RF',iset): continue        
        dataplot = deepcopy(dataTop['RF'][iset])
        dataplot['Values'] = Pullflag(dataTop['RF'][iset]['Boot'],'Std')/ Pullflag(dataBot['RF'][iset]['Boot'],'Std')
        PlotRFNoBoot(dataplot,thiscolcyc.next(),thissymcyc.next(),thisshiftcyc.next(),LegLab(iset.replace(legrem,'')),MP=True,Log=True)

def PlotLogSetOSF(data,thisSetList,thisSF,legrem=''):
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    iterSetList = SortMySet(ReduceTooMassSet(thisSetList))[0]
    for iset in iterSetList:
        if not CheckDict(data,'RF',iset): continue
        thiscol = thiscolcyc.next()
        thissym = thissymcyc.next()
        thisshift = thisshiftcyc.next()
        dataplot = deepcopy(data['RF'][iset])
        norm,DoPoFS = GetNorm(dataplot['Boot'],tsource,iset)
        dataplot['Boot'] = np.log([tboot/norm for tboot in np.roll(dataplot['Boot'],-DoPoFS,0)])
        # norm = dataplot['Boot'][tsource-1]
        # dataplot['Boot'] = np.log([tboot/norm for tboot in dataplot['Boot']])
        dataplot['Boot'] = GetBootStats(dataplot['Boot'])
        PlotRF(dataplot,thiscol,thissym,thisshift,LegLab(iset.replace(legrem,'')),MP=True,Log=True)
        if not CheckDict(data,'O'+thisSF,iset): continue
        PlotOSFLog(data['O'+thisSF][iset],thiscol,iset,norm,thisshift,DoPoFS)



def PlotLogSetTSF(data,thisSetList,thisSF,legrem=''):
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    iterSetList = SortMySet(ReduceTooMassSet(thisSetList))[0]
    for iset in iterSetList:
        if not CheckDict(data,'RF',iset): continue
        thiscol = thiscolcyc.next()
        thissym = thissymcyc.next()
        thisshift = thisshiftcyc.next()
        dataplot = deepcopy(data['RF'][iset])
        norm,DoPoFS = GetNorm(dataplot['Boot'],tsource,iset)
        dataplot['Boot'] = np.log([tboot/norm for tboot in np.roll(dataplot['Boot'],-DoPoFS,0)])
        # norm = dataplot['Boot'][tsource-1]
        # dataplot['Boot'] = np.log([tboot/norm for tboot in dataplot['Boot']])
        dataplot['Boot'] = GetBootStats(dataplot['Boot'])
        PlotRF(dataplot,thiscol,thissym,thisshift,LegLab(iset.replace(legrem,'')),MP=True,Log=True)
        if not CheckDict(data,'T'+thisSF,iset): continue
        PlotTSFLog(data['T'+thisSF][iset],thiscol,iset,norm,thisshift,DoPoFS)

def PlotMassSetOSF(data2pt,massdata,thisSetList,MassDt,thisSF,thisMom):
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    iterSetList = SortMySet(ReduceTooMassSet(thisSetList))[0]    
    for iset in iterSetList:
        # if Debug: print iset
        if not CheckDict(data2pt,'RF',iset): continue
        # if Debug: print 'present ' ,iset
        thiscol = thiscolcyc.next()
        thissym = thissymcyc.next()
        thisshift = thisshiftcyc.next()
        dataplot = deepcopy(data2pt['RF'][iset])
        dataplot['Boot'] = MassFun(dataplot['Boot'],MassDt)
        # dataplot['tVals'] = dataplot['tVals'][MassDt:]
        PlotRF(dataplot,thiscol,thissym,thisshift,LegLab(iset),MP=True)
        if not CheckDict(data2pt,'O'+thisSF,iset): continue
        PlotOSFMassValue(data2pt['O'+thisSF][iset],thiscol,iset,MassDt)
        if thisMom != 'q = 0 0 0':
            massdataplot = deepcopy(massdata['RF'][iset])
            if not CheckDict(massdata,'RF',iset): continue
            massdataplot['Boot'] = MassFun(massdataplot['Boot'],MassDt)
            massdataplot['Boot'] = ScaledEffMass(thisMom,massdataplot['Boot'])
            PlotRF(massdataplot,thiscol,thissym,thisshift,LegLab(iset+'\ Dis'),MP=True,alpha=0.5)
            if not CheckDict(massdata,'O'+thisSF,iset): continue            
            PlotOSFMassValue(massdata['O'+thisSF][iset],thiscol,iset,MassDt,MomScale=thisMom)


def PlotMassSetTSF(data2pt,thisSetList,MassDt,thisSF):
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    iterSetList = SortMySet(ReduceTooMassSet(thisSetList))[0]
    for iset in iterSetList:
        if not CheckDict(data2pt,'RF',iset): continue
        thiscol = thiscolcyc.next()
        thissym = thissymcyc.next()
        thisshift = thisshiftcyc.next()
        dataplot = deepcopy(data2pt['RF'][iset])
        dataplot['Boot'] = MassFun(dataplot['Boot'],MassDt)
        # dataplot['tVals'] = dataplot['tVals'][MassDt:]
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
    if 'Top' in DSCurr or 'Wein' in DSCurr: 
        thisSetList = FlowOrderSetList(thisSetList)
    for iFF in xrange(1,NoFFPars[thisCurr]+1):
        thisFF = 'FF'+str(iFF)
        if len(thisFFComb) > 1: thisFFComb = '/'+thisFFComb
        DatFile = CreateFFFile(CollName,DSCurr,thisFF)+'.dat'
        WipeFile(DatFile)
        graphparams = GetPlotItersff()
        if Debug: print 'Running PlotFFSet for ', iFF
        PlotFFSet(data,thisFF,thisSetList,thisCurr,DSCurr.replace('/',''),graphparams)
        SetFFAxies(thisDS+thisCurr+thisFF+thisFFComb)
        pl.savefig(DatFile.replace('.dat','.pdf'))
        pl.clf()

def PlotFFsPN(data,dataPN,DSCurr,thisSetList,CollName,FT):
    global ForceTitle
    global DatFile
    ForceTitle = FT
    if len(thisSetList) == 0: return
    thisDS,thisCurr,thisFFComb = SplitDSCurr(DSCurr)
    thisDS = thisDS.replace('Proton','PandN')
    DSCurr = DSCurr.replace('Proton','PandN')
    for iFF in xrange(1,NoFFPars[thisCurr]+1):
        thisFF = 'FF'+str(iFF)
        if len(thisFFComb) > 1: thisFFComb = '/'+thisFFComb
        DatFile = CreateFFFile(CollName,DSCurr,thisFF)+'.dat'
        WipeFile(DatFile)
        graphparams = GetPlotItersff()
        PlotFFSet(data,thisFF,['Proton'+iset for iset in thisSetList],thisCurr,DSCurr.replace('/',''),graphparams,PandNshift=0.0)
        graphparams[2] = GetPlotItersff()[2]
        PlotFFSet(dataPN,thisFF,['Neutron'+iset for iset in thisSetList],thisCurr.replace('Proton','Neutron'),DSCurr.replace('/',''),graphparams,PandNshift=0.005)
        SetFFAxies(thisDS+thisCurr+thisFF+thisFFComb)
        pl.savefig(DatFile.replace('.dat','.pdf'))
        pl.clf()

        
def PlotMKFFs(kdata,DSCurr,thisSetList,CollName,FT):
    global ForceTitle
    global DatFile
    ForceTitle = FT
    if len(thisSetList) == 0: return
    thisDS,thisCurr,thisFFComb = SplitDSCurr(DSCurr)
    if len(thisFFComb) > 1: thisFFComb = '/'+thisFFComb
    for iFF in xrange(1,NoFFPars[thisCurr]+1):
        thisFF = 'FF'+str(iFF)
        DatFile = CreateMKFFFile(CollName,DSCurr,thisFF)+'.dat'
        WipeFile(DatFile)
        graphparams = GetPlotItersff()
        for ikappa,data in kdata.iteritems():
            thiskSetList = []
            for iset in thisSetList:
                if ikappa in iset:
                    thiskSetList.append(iset)
            PlotFFSet(data,thisFF,thiskSetList,thisCurr,DSCurr.replace('/',''),graphparams)
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
        
def PlotFFSet(dataset,thisFF,thisSetFlag,thisCurr,thisDSCurr,graphparams,PandNshift=0):
    thissymcyc,thiscolcyc,thisshiftcycff = graphparams
    collist = []
    FixZ=False    
    # if 'Ge' in thisDSCurr or ('Vector' in thisDSCurr and ('PsVector' not in thisDSCurr.replace('IsoVector',''))):
    #     if 'IsoVector' in thisDSCurr or 'Proton' in thisDSCurr or 'sing' in thisDSCurr:
    #         FixZ= 1
    #     elif 'Neutron' in thisDSCurr:
    #         FixZ= 0
    #     elif 'doub' in thisDSCurr:
    #         FixZ= 2
    # DatDPFile = './'+thisFF+thisDSCurr+'.dat'
    datf = open(DatFile.replace('.dat','Params.dat'),'w')
    if 'FF1' in thisFF:
        if 'PsVector' in thisCurr:
            datf.write('\(\langle r_{A}^2 \\rangle\) \(fm^{2}\) \n')
        else:
            datf.write('\(\langle r^2 \\rangle\) \(fm^{2}\) \n')
    elif 'FF2' in thisFF:
        datf.write('\(\mu\)  \n')
    else:
        datf.write('F_{3}(0) \n')
    # for thisset in SortMySet(thisSetFlag)[0]:
    if Debug: print thisSetFlag
    for thisset in thisSetFlag:
        ##make legend formatting function
        thiskappa,keyset = SplitKappa(thisset)
        # if Debug:
        #     print dataset.keys(), keyset
        #     if keyset in dataset.keys(): dataset[keyset],thisFF
        #     print

        if ('Proton' in keyset or 'Neutron' in keyset):
            if FlowPlot[0] not in keyset and 'Top' in thisDSCurr: continue
            if WeinFlowPlot[0] not in keyset and 'Wein' in thisDSCurr: continue
        PorN = ''
        if 'Proton' in keyset: PorN = 'Proton'
        if 'Neutron' in keyset: PorN = 'Neutron'
        keyset = keyset.replace('Proton','')
        keyset = keyset.replace('Neutron','')
        if not CheckDict(dataset,keyset,thisFF): continue
        if dataset[keyset][thisFF] == False: continue        
        thiscol = thiscolcyc.next()
        collist.append(thiscol)
        skipzero,flipsign = SkipZeroFF(thisFF,keyset,thisCurr)
        if ('IsoVectorPsVector' in thisDSCurr and 'FF2' in thisFF) or ('NeutronGeGm' in thisDSCurr and 'FF1' in thisFF) or 'F1divF2' in thisDSCurr or 'FF3' in thisFF:
        # if ('IsoVectorPsVector' in thisDSCurr) or ('NeutronGeGm' in thisDSCurr and 'FF1' in thisFF):
            thisshift = thisshiftcycff.next()
            qrange = PlotFF(dataset[keyset][thisFF],thiscol,thissymcyc.next(),thisshift+PandNshift,LegLabFF(thisset),skipzero,flipsign,FixZ=FixZ)
            if 'FF3' in thisFF:
                PlotDPFit(keyset,thisFF,thisDSCurr,thiscol,qrange,thisshift+PandNshift,flipsign,datf,thiskappa,thisPN=PorN)
        else:
            if 'FF1' in thisFF and 'ProtonGeGm' in thisDSCurr:
                thisshift = 0.0
            else:
                thisshift = thisshiftcycff.next()*0.5
            qrange = PlotFF(dataset[keyset][thisFF],thiscol,thissymcyc.next(),thisshift,LegLabFF(thisset),skipzero,flipsign,FixZ=FixZ)
            # if 'sm32' in thisset or 'CM' in thisset or 'TSF' in thisset or '12104' in str(thiskappa):
            PlotDPFit(keyset,thisFF,thisDSCurr,thiscol,qrange,thisshift,flipsign,datf,thiskappa,thisPN=PorN)
    if 'FF1' in thisFF :
        if 'ProtonGeGm' in thisDSCurr:
            PlotExp('ChargeRadius',thiscolcyc)
        elif 'ProtonIsoVector' in thisDSCurr:
            PlotExp('AxialChargeRadius',thiscolcyc)
        elif 'NeutronGeGm' in thisDSCurr:
            PlotExp('NeutronCRad',thiscolcyc)            
    elif 'FF2' in thisFF and 'GeGm' in thisDSCurr:
        if 'Neutron' in thisDSCurr:
            PlotExp('NeutronMagMom',thiscolcyc)
            PlotExp('MRadNeutron',thiscolcyc)
        elif 'Proton' in thisDSCurr:
            PlotExp('ProtonMagMom',thiscolcyc)
            PlotExp('MRadProton',thiscolcyc)
    if 'FF3' in thisFF and 'Top' in thisDSCurr:
        if 'Neutron' in thisDSCurr:
            PlotExp('ProtonEDM',thiscolcyc)
        elif 'Proton' in thisDSCurr:
            PlotExp('NeutronEDM',thiscolcyc)
        elif 'PandN' in thisDSCurr and PandNshift == 0.0:
            PlotExp('PandNEDM',thiscolcyc)
    datf.close()
    return collist

def PlotDPFit(thisset,thisFF,thisCurr,thiscol,qrange,thisshift,flipsign,datf,thiskappa=kappa,thisPN=''):
    if Debug: print thisset,thisFF
    Avg,Err = GetDPFitValue(thisset,thisFF,thisCurr.replace('PandN',thisPN),thiskappa=thiskappa)
    if len(Avg) == 0: return

    Avg,Err = np.array(Avg),np.array(Err)
    fitqdata = np.arange(qrange[0],qrange[-1]+incr,incr)
    thisFitFun = DPfitfun
    if '3' in thisFF: thisFitFun = LinearFitFun
    fitydataAvg = thisFitFun([fitqdata],Avg)
    fitydataup = np.array([max(iy1,iy2) for iy1,iy2 in zip(thisFitFun([fitqdata],Avg+Err),thisFitFun([fitqdata],np.array([Avg[0]-Err[0],Avg[1]+Err[1]])))])
    fitydatadown = np.array([min(iy1,iy2) for iy1,iy2 in zip(thisFitFun([fitqdata],Avg-Err),thisFitFun([fitqdata],np.array([Avg[0]+Err[0],Avg[1]-Err[1]])))])
    # fitydatadown = np.array([np.min(iy1,iy2) for iy1,iy2 in zip(thisFitFun([fitqdata],Avg-Err),thisFitFun([fitqdata],np.array([Avg[0]+Err[0],Avg[1]-Err[1]])))])
    # if Debug: print Avg[1], Err[1]
    # if GetCharRad(Avg[1]) > 10 or Err[1]> 10: return
    ## Displays charge radius for FF1, and magnetic moment for FF2
            
    if 'FF1' in thisFF:
        datf.write(MakeValAndErr(Avg[2],Err[1])+' \n')
        if 'PsVector' in thisCurr:
            LegVal = '$\\langle r_{A}^2 \\rangle='+MakeValAndErr(Avg[2],Err[1])+'\ fm^{2}$'
        else:
            LegVal = '$\\langle r^2 \\rangle='+MakeValAndErr(Avg[2],Err[1])+'\ fm^{2}$'        
    elif 'FF2' in thisFF:
        LegVal = '$\\mu='+MakeValAndErr(Avg[0],Err[0])+r'$'        
        datf.write(MakeValAndErr(Avg[0],Err[0])+' \n')
    else:
        LegVal = '$\\frac{F_{3}(0)}{2m_{N}}='+MakeValAndErr(Avg[1],Err[1])+r'\ \theta fm$'        
        datf.write(MakeValAndErr(Avg[1],Err[1])+' \n')

    # print 'DPFit flip sign', flipsign
    # if flipsign:
    #     pl.plot(fitqdata+thisshift,-np.array(fitydataAvg),label=LegVal,color=thiscol)
    #     pl.fill_between(fitqdata+thisshift,-np.array(fitydataup),-np.array(fitydatadown),color=thiscol,alpha=thisalpha,edgecolor='none')
    # else:
    if flipsign:
        fitydataAvg = -np.array(fitydataAvg)
        Avg = -Avg
    pl.plot(fitqdata+thisshift,fitydataAvg,color=thiscol)
    if Err[0] < 1.0 :
        if 'FF1' in thisFF and 'Neutron' in thisCurr: return
        if 'FF3' in thisFF:
            pl.errorbar([0.0+thisshift],[Avg[1]],[Err[1]],fmt='-',color=thiscol,label=LegVal)
        else:
            pl.errorbar([0.0+thisshift],[Avg[0]],[Err[0]],fmt='-',color=thiscol,label=LegVal)
            # pl.fill_between(fitqdata+thisshift,fitydataup,fitydatadown,color=thiscol,alpha=thisalpha,edgecolor='none')
    
def PlotFF(data,col,sym,shift,lab,SkipZero,FlipSign,FixZ=False):
    qsqrdvals,dataavg,dataerr = [],[],[]
    for iqsqrd,(qsqrd,values) in enumerate(data.iteritems()):
        Qsqrd = GetQsqrd(float(qsqrd.replace('qsqrd','')),DefMass[str(kappa)],Phys=PhysicalUnits)
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
        # print 'FFvals flip sign', FlipSign
        if FlipSign: dataavg = -1*np.array(dataavg)
        if ForcePos: dataavg = np.abs(dataavg)
        AppendFFDat(qsqrdvals,dataavg,dataerr)
        if SkipZero and len(qsqrdvals) > 1:
            pl.errorbar(qsqrdvals[1:],dataavg[1:],dataerr[1:],color=col,fmt=sym,label=lab)
            if Debug: print 'Plotting ',lab
        elif FixZ != False:
            pl.plot([0],[FixZ],sym,color=col)            
            pl.errorbar(qsqrdvals[1:],dataavg[1:],dataerr[1:],color=col,fmt=sym,label=lab)            
        else:
            pl.errorbar(qsqrdvals,dataavg,dataerr,color=col,fmt=sym,label=lab)
    return np.array(qsqrdvals)-shift

def PlotRF(data,col,sym,shift,lab,MP=False,Log=False,alpha=1):
    if MP:
        if PhysMassPlot and not Log:
            data['Boot'] = [idata*hbarcdivlat for idata in data['Boot']]
            GetBootStats(data['Boot'])
        if 'PoF' in lab and not Log:
            tvals = np.array(data['tVals'])+1+(PoFShifts*PoFDelta) + shift
        else:
            tvals = np.array(data['tVals'])+1 + shift
    else:
        tvals = np.array(data['tVals'])
        tvals = tvals-(tvals[-1]+tvals[0])/2.0 + shift
    dataavg = Pullflag(data['Boot'],'Avg')
    dataerr = Pullflag(data['Boot'],'Std')
    if Debug and not Log:
        print lab
        # for it,val,valerr in zip(tvals,dataavg,dataerr):
        #     print it,val,valerr
    if ForcePos: dataavg = np.abs(dataavg)
    pl.errorbar(tvals[tsource:]-tsource,dataavg[tsource:],dataerr[tsource:],color=col,fmt=sym,label=lab,alpha=alpha)

def PlotRFNoBoot(data,col,sym,shift,lab,MP=False,Log=False):
    if MP:
        if 'PoF' in lab and not Log:
            tvals = np.array(data['tVals'])+1+(PoFShifts*PoFDelta) + shift
        else:
            tvals = np.array(data['tVals'])+1 + shift
    else:
        tvals = np.array(data['tVals'])
        tvals = tvals-(tvals[-1]+tvals[0])/2.0 + shift
    # dataavg = Pullflag(data['Boot'],'Avg')
    # dataerr = Pullflag(data['Boot'],'Std')
    if ForcePos: data['Values'] = np.abs(data['Values'])
    if Debug:
        print 
        print 'Error Ratios: ',lab
        # for it,idata in zip(tvals[tsource:]-tsource,data['Values'][tsource:]):
        #     print it, idata
    pl.plot(tvals[tsource:]-tsource,data['Values'][tsource:],sym,color=col,label=lab)


def PlotFit(data,col,shift,iset,thistsink):
    thiscut = GetCut(iset,FitCutPicked)
    if thiscut == False:
        if Debug:
            print 'warning', iset, 'not in FitCutPicked'
            print FitCutPicked.keys()
    else:
        cutlow,cuthigh = map(int,thiscut.replace('cut','').split('-'))
        LM = max((int(thistsink)-tsource)/2.-cutlow,0)
        RM = max((int(thistsink)-tsource)/2.-cuthigh,0)
        tvals = np.array([-LM+shift,RM+shift])
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
    tdata = np.arange(TSFfitvals[0]-tsource,TSFfitvals[1]+incr-tsource,incr)
    fit2pt = ff.C2TSFLineFun(tdata,pars2pt)
    fit2ptdt = ff.C2TSFLineFun(tdata+thisdt,pars2pt)
    effline = map(abs,np.log(fit2ptdt/fit2pt)/thisdt)
    if PhysMassPlot: effline = np.array(effline)*hbarcdivlat
    pl.plot(tdata,effline,color=col)


def PlotOSFMassValue(data,col,smear,thisdt,MomScale='q = 0 0 0'):
    smearindex,deltashift = CreateOSFfitKey(smear)
    units = '$'
    if PhysMassPlot: units = r'\ GeV$'
    if Debug:
        print data.keys(), 'm0'
        print OSFfitr.keys(), smearindex, smear
        print data['m0'].keys(), OSFfitr[smearindex]
        print data['m0'][OSFfitr[smearindex]].keys(),'Boot','Avg','Std'
    if CheckDict(data,'m0',OSFfitr[smearindex],'Boot'): 
        databoot = data['m0'][OSFfitr[smearindex]]['Boot']
        databoot = ScaledEffMass(MomScale,[databoot])[0]
        if PhysMassPlot:
            databoot = databoot*hbarcdivlat
            databoot.Stats()
        dataval = abs(databoot.Avg)
        datastd = databoot.Std
        dataup,datadown = dataval+databoot.Std,dataval-databoot.Std
    elif CheckDict(data,'m0',OSFfitr[smearindex],'Avg') and CheckDict(data,'m0',OSFfitr[smearindex],'Std'):
        dataval = data['m0'][OSFfitr[smearindex]]['Avg']
        datastd = data['m0'][OSFfitr[smearindex]]['Std']
        if PhysMassPlot:
            dataval,datastd = dataval*hbarcdivlat,datastd*hbarcdivlat
        dataup,datadown = dataval+datastd,dataval-datastd
    else:
        if Debug: print 'OSF',smear, 'not found'
        return
    tdata = [OSFfitvals[smearindex][0]-tsource,OSFfitvals[smearindex][1]-tsource]
    thislab = '$M_{N}='+MakeValAndErr(dataval,datastd)+units
    pl.fill_between(tdata,[datadown,datadown],[dataup,dataup],facecolor=col,edgecolor='none',alpha=thisalpha)
    pl.plot(tdata,[dataval,dataval],color=col,label=thislab)

def PlotTSFMassValue(data,thisdt):
    databoot = data['m0'][TSFfitr]['Boot']
    units = '$'
    if PhysMassPlot:
        units = r'\ GeV$'
        databoot = databoot*hbarcdivlat
        databoot.Stats()
    thislab = '$M_{N}='+MakeValAndErr(dataval,databoot.Std)+units
    dataval = abs(databoot.Avg)
    dataup,datadown = dataval+databoot.Std,dataval-databoot.Std
    pl.fill_between([TSFfitvals[0]-tsource ,TSFfitvals[1]-tsource],[datadown,datadown],[dataup,dataup],facecolor='k',edgecolor='none',alpha=thisalpha,label=thislab)


def PlotOSFLog(data,col,smear,norm,thisshift,DoPoFS):
    smearindex,deltashift = CreateOSFfitKey(smear)
    if 'sum' in smear: smearindex = PickedStateStr+'sum'    
    if not CheckDict(data,'m0',OSFfitr[smearindex],'Boot'): return
    if not CheckDict(data,'Am',OSFfitr[smearindex],'Boot'): return
    parm0 = data['m0'][OSFfitr[smearindex]]['Boot']
    parAm = data['Am'][OSFfitr[smearindex]]['Boot']
    tdata = np.arange(OSFfitvals[smearindex][0],OSFfitvals[smearindex][1]+incr,incr)-tsource
    linedata = []
    for it in tdata:
        linedata.append(np.log((parAm*(parm0*(-it)).exp(1))/norm))
    GetBootStats(linedata)
    dataAvg,dataErr = np.array(Pullflag(linedata,'Avg')),np.array(Pullflag(linedata,'Std'))
    units = '$'
    if PhysMassPlot:
        units = r'\ GeV$'
        dataAvg,dataErr = dataAvg*hbarcdivlat,dataErr*hbarcdivlat
    thislab = '$M_{N}='+MakeValAndErr(dataAvg,dataErr)+units
    dataup = dataAvg+dataErr
    datadown = dataAvg-dataErr
    pl.fill_between(tdata-DoPoFS+thisshift,datadown,dataup,facecolor=col,edgecolor='none',alpha=thisalpha,label=thislab)
    pl.plot([tdata[0]-DoPoFS+thisshift,tdata[-1]-DoPoFS+thisshift],[dataAvg[0],dataAvg[-1]],color=col)


def PlotTSFLog(data,col,smear,norm,thisshift,DoPoFS):
    if not CheckDict(data,'m0',TSFfitr,'Boot'): return
    if not CheckDict(data,'Am',TSFfitr,'Boot'): return
    if not CheckDict(data,'Dm',TSFfitr,'Boot'): return
    if not CheckDict(data,'Amp',TSFfitr,'Boot'): return
    parm0 = data['m0'][TSFfitr]['Boot']
    parAm = data['Am'][TSFfitr]['Boot']
    parDm = data['Dm'][TSFfitr]['Boot']
    parAmp = data['Amp'][TSFfitr]['Boot']
    tdata = np.arange(TSFfitvals[0],TSFfitvals[1]+incr,incr)-tsource
    linedata = []
    for it in tdata:
        linedata.append(np.log((parAm*(((parm0*(-it)).exp(1)) + parAmp*((parm0+parDm)*(-it)).exp(1)))/norm))
    GetBootStats(linedata)
    dataAvg,dataErr = np.array(Pullflag(linedata,'Avg')),np.array(Pullflag(linedata,'Std'))
    units = '$'
    if PhysMassPlot:
        units = r'\ GeV$'
        dataAvg,dataErr = dataAvg*hbarcdivlat,dataErr*hbarcdivlat
    thislab = '$M_{N}='+MakeValAndErr(dataAvg,dataErr)+units
    dataup = dataAvg+dataErr
    datadown = dataAvg-dataErr
    pl.fill_between(tdata-DoPoFS+thisshift,datadown,dataup,facecolor=col,edgecolor='none',alpha=thisalpha,label=thislab)
    pl.plot(tdata-DoPoFS+thisshift,dataAvg,color=col)



def PlotChiDist(data,col,sym,lab,thisFitMax):
    fitpar = len(data.keys())
    if fitpar == 0: return
    xdata,yavg,yerr,yup,ydown = [],[],[],[],[]
    for thisfitr,chidata in data['m0'].iteritems():
        if CheckDict(chidata,'Chi') and  CheckDict(chidata,'Boot'):
            thisfitmin,thisfitmax = map(int,unxmlfitr(thisfitr))            
            fitlen = thisfitmax-thisfitmin
            dof = fitlen-fitpar
            if thisFitMax == thisfitmax:
                if VsChi:
                    xdata.append(chidata['Chi'])
                else:
                    xdata.append(thisfitmin)
                yavg.append(chidata['Boot'].Avg)
                yerr.append(chidata['Boot'].Std)
                yup.append(chidata['Boot'].Avg + chidata['Boot'].Std)
                ydown.append(chidata['Boot'].Avg - chidata['Boot'].Std)
    if len(xdata) == 0: return
    xdata,yavg,yerr,yup,ydown = zip(*sorted(zip(xdata,yavg,yerr,yup,ydown)))
    pl.errorbar(xdata,yavg,yerr,color=col,label=lab,fmt=sym)
    pl.plot(xdata,yavg,color=col)
    pl.fill_between(xdata,yup,ydown,color=col,alpha=thisalpha,edgecolor='none')


def PlotSetChiDist(data,thisSetList,thisSF,thisFitMax):
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    iterSetList = SortMySet(ReduceTooMassSet(thisSetList))[0]
    for iset in iterSetList:
        if not CheckDict(data,thisSF,iset):
            if Debug:
                print data.keys(), thisSF
                print data[thisSF].keys(), iset
            continue
        thiscol = thiscolcyc.next()
        thissym = thissymcyc.next()
        thisshift = thisshiftcyc.next()
        PlotChiDist(data[thisSF][iset],thiscol,thissym,LegLab(iset),thisFitMax)




def PlotAlphaDist(data):
    fitpar = len(data.keys())
    if fitpar == 0: return
    doflist = np.arange(3,17,4)
    # for thisfitr,chidata in data[data.keys()[0]].iteritems():
    #     if CheckDict(chidata,'Chi'):
    #         fitlen = map(int,unxmlfitr(thisfitr))
    #         fitlen = fitlen[1] - fitlen[0]
    #         if fitlen-fitpar not in doflist:
    #             doflist.append(fitlen-fitpar)
    # doflist = sorted(doflist)
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    for idof in doflist:
        thiscol = thiscolcyc.next()
        xdata,ydata = AlphaVsChiDOF(idof)
        pl.plot(xdata,ydata,color=thiscol,label='DoF'+str(idof))



def GetSetHist(data,iThresh,thisMom=False):
    yvals = []
    if CheckDict(data,'m0'):
        for thisfitr,chidata in data['m0'].iteritems():
            if CheckDict(chidata,'Chi') and  CheckDict(chidata,'Boot'):
                thisfitmin,thisfitmax = map(int,unxmlfitr(thisfitr))            
                if float(chidata['Chi']) < iThresh and thisfitmax < FitMaxCutoff and thisfitmin > FitMinCutoff:
                    if thisMom != False:
                        yvals += ScaledEffMassList(thisMom,RemoveNAN(np.array(chidata['Boot'].values))).tolist()
                    else:
                        yvals += RemoveNAN(np.array(chidata['Boot'].values).tolist())  
    yavg,ystd = np.mean(yvals),np.std(yvals)
    return yvals,yavg,ystd


def PlotSetHistDist(data,thisSetList,thisSF,iThresh,thisMom,Zmom=False):
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    iterSetList = SortMySet(ReduceTooMassSet(thisSetList))[0]
    collist,setlist,leglist = [],[],[]
    for iset in iterSetList:
        if not CheckDict(data,thisSF,iset):
            if Debug:
                print data.keys(), thisSF
                print data[thisSF].keys(), iset
            continue
        thishist,thisavg,thisstd = GetSetHist(data[thisSF][iset],iThresh)
        if len(thishist) > 0:
            setlist.append(thishist)
            collist.append(thiscolcyc.next())
            leglist.append(LegLab(iset+'\ M='+ MakeValAndErr(thisavg,thisstd)))
            # if MakeValAndErr(thisavg,thisstd) == 'Err':
            #     print iset
            #     for ihist in thishist:
            #         print ihist
            
    if len(setlist) > 0:
        pl.hist(setlist,bins=BinList,color=collist,label=leglist,stacked=Stacked,histtype=HistType,normed=Normed)
    if Zmom != False:
        thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
        collist,setlist,leglist = [],[],[]
        for iset in iterSetList:
            if not CheckDict(Zmom,thisSF,iset): continue
            thishist,thisavg,thisstd = GetSetHist(Zmom[thisSF][iset],iThresh,thisMom=thisMom)
            if len(thishist) > 0:
                setlist.append(thishist)
                collist.append(thiscolcyc.next())
                leglist.append(LegLab(iset+'\ M='+ MakeValAndErr(thisavg,thisstd)+'\ Disp'))
        if len(setlist) > 0:
            pl.hist(setlist,bins=BinList,color=collist,label=leglist,stacked=Stacked,histtype=HistType,normed=Normed,alpha=0.5)
                


def PlotTopChargeOverFlow(data,iSet,iMom,tsink,thiscol,thissym,thisshift,NNQ=False,Dt=1):
    tflowlist,tlist,plotAvg,plotStd = [],[],[],[]
    plotAvgNN,plotStdNN = [],[]
    DictFlag,thisValue = 'RF','Alapha'
    if NNQ:  DictFlag,thisValue = 'NNQ','NNQ'
    momdata = data[DictFlag][iMom]['Boots']
    momdataNN = data['cfun'][iMom]['Boots']
    for itflow,flowdata in momdata.iteritems():
        if 't_flow0.0' == itflow: continue
        for (it,tdata),(itdump,tdataNN) in zip(flowdata.iteritems(),momdataNN.iteritems()):
            if tsink == int(untstr(it)):
                tflowlist.append(untflowstr(itflow)-thisshift)
                # tlist.append(untstr(it))
                if NNQ:
                    if untstr(it) +Dt <= len(flowdata.keys()):
                        EffMass = np.log(np.abs(np.array(tdata)/np.array(flowdata['t'+str(untstr(it)+Dt)])))/float(Dt)
                        EffMassNN = np.log(np.abs(tdataNN/momdataNN['t'+str(untstr(it)+Dt)]))/float(Dt)
                    else:
                        EffMass = 1.0
                        EffMassNN = 1.0
                    plotAvg.append(np.mean(EffMass))
                    plotStd.append(np.std(EffMass))
                    plotAvgNN.append(np.mean(EffMassNN))
                    plotStdNN.append(np.std(EffMassNN))
                else:
                    plotAvg.append(np.mean(tdata))
                    plotStd.append(np.std(tdata))
    if len(plotAvg) == 0: return
    pl.errorbar(TflowToPhys(tflowlist),plotAvg,plotStd,color=thiscol,fmt=thissym,label=LegLab(iSet+'\ t='+str(tsink)))
    if NNQ:
        pl.errorbar(np.array(TflowToPhys(tflowlist))-0.25,plotAvgNN,plotStdNN,color=thiscol,fmt=thissym,label=LegLab(iSet+' NN'),alpha=0.6)
        
def PlotTopChargeOvert(data,fitdata,iSet,iMom,tflow,thiscol,thissym,thisshift,NNQ=False,Dt=1,Wein=False):
    tflowlist,tlist,plotAvg,plotStd = [],[],[],[]
    plotAvgNN,plotStdNN = [],[]
    DictFlag,thisValue = 'RF','Alapha'
    if NNQ:  DictFlag,thisValue = 'NNQ','NNQ'
    momdata = data[DictFlag][iMom]['Boots']
    momdataNN = data['cfun'][iMom]['Boots']
    for itflow,flowdata in momdata.iteritems():        
        if 't_flow0.0' == itflow: continue
        if untflowstr(itflow) == tflow:
            for (it,tdata),(itdump,tdataNN) in zip(flowdata.iteritems(),momdataNN.iteritems()):
                if untstr(it) < AlphaTlist[-1] and untstr(it) > MassTVals[0] :
                    # tflowlist.append(untflowstr(itflow))
                    tlist.append(untstr(it)-thisshift)
                    if NNQ:
                        if untstr(it) + Dt <= len(flowdata.keys()):
                            EffMass = np.abs(np.log(np.abs(np.array(tdata)/np.array(flowdata['t'+str(untstr(it)+Dt)])))/float(Dt))
                            EffMassNN = np.abs(np.log(np.abs(tdataNN/momdataNN['t'+str(untstr(it)+Dt)]))/float(Dt))
                        else:
                            EffMass = 1.0
                            EffMassNN = 1.0
                        plotAvg.append(np.mean(EffMass))
                        plotStd.append(np.std(EffMass))
                        plotAvgNN.append(np.mean(EffMassNN))
                        plotStdNN.append(np.std(EffMassNN))
                    else:
                        plotAvg.append(np.mean(tdata))
                        plotStd.append(np.std(tdata))
    if len(plotAvg) == 0: return
    pl.errorbar(tlist,plotAvg,plotStd,color=thiscol,fmt=thissym,label=LegLab(iSet+'\ \sqrt{8t_{f}}='+'{:.2f}'.format(TflowToPhys(tflow))))
    if NNQ:
        pl.errorbar(np.array(tlist)-0.5,plotAvgNN,plotStdNN,color=thiscol,fmt=thissym,label=LegLab(iSet+' NN'),alpha=0.6)
    WeinOrTop = 'Top'
    if Wein: WeinOrTop = 'Wein'
    if not NNQ and CheckDict(fitdata,iMom,'Boots'):
        momdataFit = fitdata[iMom]['Boots']
        tflowlist = [] 
        for itflow,flowfitdata in momdataFit.iteritems():
            if untflowstr(itflow) == tflow:
                if AlphaFitRPick[WeinOrTop] in flowfitdata.keys():
                    tvals = map(int,unxmlfitr(AlphaFitRPick[WeinOrTop]))
                    dataavg = flowfitdata[AlphaFitRPick[WeinOrTop]].Avg
                    datastd = flowfitdata[AlphaFitRPick[WeinOrTop]].Std
                    dataup = dataavg+datastd
                    datadown = dataavg-datastd
                    pl.plot(tvals,[dataavg,dataavg],color=thiscol)
                    pl.fill_between(tvals,[dataup,dataup],[datadown,datadown],color=thiscol,alpha=thisalpha,edgecolor='none')
                    
def PlotTopSetCharge(data,thisSetList,imom,FT,NNQ=False,Wein=False):
    global ForceTitle    
    DictFlag,thisValue = 'RF','Alpha'
    if NNQ:  DictFlag,thisValue = 'NNQ','NNQ'
    Dt=1
    ForceTitle = FT
    if Wein: thissubdir = 'Wein'
    else: thissubdir = 'Top'
    for itsink in AlphaTlist:
        thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
        for iset,setdata in data.iteritems():
            MultKappa = len(setdata.keys()) > 1
            for ikappa,kappadata in setdata.iteritems():
                if CheckDict(kappadata,DictFlag,imom,'Boots'):
                    thisset = iset
                    if MultKappa: thisset += '\ '+GetMpi(ikappa,Phys=True)
                    # print 'plotting ', iset, imom
                    PlotTopChargeOverFlow(kappadata,thisset,imom,itsink,thiscolcyc.next(),thissymcyc.next(),thisshiftcyc.next(),NNQ=NNQ,Dt=Dt)
                
        filename = CreateFile('','twopt',imom,thisValue+'Overt'+str(itsink),subdir=thissubdir,NoMpi=MultKappa)        
        SetTopAxies('flow',NNQ=NNQ,Dt=Dt,Wein=Wein,MpiTitle=not MultKappa)
        pl.savefig(filename+'.pdf')
        pl.clf()
    for itflow in AlphaTflowList:
        thisitflow = itflow
        if Wein: thisitflow = itflow - 0.01
        thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
        for iset,setdata in data.iteritems():
            MultKappa = len(setdata.keys()) > 1
            for ikappa,kappadata in setdata.iteritems():
                if CheckDict(kappadata,DictFlag,imom,'Boots'):
                    # print 'plotting ', iset, imom
                    thisset = iset
                    if MultKappa: thisset += '\ '+GetMpi(ikappa,Phys=True)
                    if CheckDict(setdata,'Fits',ikappa):
                        fitdata = setdata['Fits'][ikappa]
                    else:
                        fitdata = {}
                    PlotTopChargeOvert(kappadata,fitdata,thisset,imom,thisitflow,thiscolcyc.next(),thissymcyc.next(),thisshiftcyc.next(),NNQ=NNQ,Dt=Dt,Wein=Wein)
        filename = CreateFile('','twopt',imom,thisValue+'OverFlow'+str(thisitflow),subdir=thissubdir,NoMpi= MultKappa)
    
        SetTopAxies('t',NNQ=NNQ,Dt=Dt,Wein=Wein,MpiTitle=not MultKappa)
        pl.savefig(filename+'.pdf')
        pl.clf()


        
def PlotTopChargeOverMpi(data,iSet,iMom,tflow,thiscol,thissym,thisshift,Wein=False):
    MpiList,plotAvg,plotStd = [],[],[]
    DictFlag,thisValue = 'RF','Alapha'
    tflow = untflowstr(tflow)
    WeinOrTop = 'Top'
    if Wein: WeinOrTop = 'Wein'
    for ikappa,kappadata in data.iteritems():
        if not CheckDict(kappadata,iMom,'Boots'): continue
        momdata = kappadata[iMom]['Boots']
        for itflow,flowdata in momdata.iteritems():        
            if 't_flow0.0' == itflow: continue
            if untflowstr(itflow) == tflow:
                if AlphaFitRPick[WeinOrTop] in flowdata.keys():                
                    MpiList.append(GetMpiNoForm(ikappa))
                    plotAvg.append(flowdata[AlphaFitRPick[WeinOrTop]].Avg)
                    plotStd.append(flowdata[AlphaFitRPick[WeinOrTop]].Std)
    pl.errorbar(np.array(MpiList)+thisshift*0.1,plotAvg,plotStd,color=thiscol,fmt=thissym,label=LegLab(iSet+'\ '+tflowstr(tflow)))

                  
def PlotTopChargeOverKappa(data,thisSetList,imom,FT,Wein=False):
    global ForceTitle    
    DictFlag,thisValue = 'RF','Alpha'
    ForceTitle = FT
    if Wein:
        thissubdir = 'Wein'
        thisFlowPlot = WeinFlowPlot
    else:
        thissubdir = 'Top'
        thisFlowPlot = FlowPlot
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    for iflow in thisFlowPlot:
        for iset,setdata in data.iteritems():
            thisset = iset
            # print 'plotting ', iset, imom, iflow
            PlotTopChargeOverMpi(setdata['Fits'],thisset,imom,iflow,thiscolcyc.next(),thissymcyc.next(),thisshiftcyc.next())

    filename = CreateFile('','twopt',imom,thisValue+'OverKappa',subdir=thissubdir) 
    SetTopAxies('Mpi',Wein=Wein,MpiTitle=False)
    pl.savefig(filename+'.pdf')
    pl.clf()

        
def GraphQExp(Qlist,flowlist):
    pl.errorbar(TflowToPhys(flowlist),np.mean(Qlist,axis=0),np.std(Qlist,axis=0),fmt='o')
    pl.xlim(TflowToPhys([flowlist[0]-0.1])[0],TflowToPhys([flowlist[-1]+0.1])[0])
    pl.xlabel(r'$ \sqrt{8t_{f}} fm$')
    pl.ylabel(r'$ \langle Q \rangle $')
    thisdir = outputdir[0] + 'graphs/Qdata/'
    pl.title(r'$ \langle Q \rangle \ '+GetMpi(kappa,Phys=True)+r'$',y=TitleShift)
    mkdir_p(thisdir)
    pl.savefig(thisdir+'QExp.pdf')
    pl.clf()

def GraphWExp(Wlist,flowlist):
    if Debug:
        print 'Printing WExpt'
        print len(TflowToPhys(flowlist))
        print  len(np.mean(Wlist,axis=0))
        for itflow,iW,iWerr in zip(TflowToPhys(flowlist),np.mean(Wlist,axis=0),np.std(Wlist,axis=0)):
            print itflow,iW,iWerr
            
    pl.errorbar(TflowToPhys(flowlist),np.mean(Wlist,axis=0),np.std(Wlist,axis=0),fmt='o')
    pl.xlim(TflowToPhys([flowlist[0]-0.1])[0],TflowToPhys([flowlist[-1]+0.1])[0])
    pl.xlabel(r'$ \sqrt{8t_{f}} fm$')
    pl.ylabel(r'$ \langle W \rangle $')
    thisdir = outputdir[0] + 'graphs/Wdata/'
    pl.title(r'$ \langle W \rangle \ '+GetMpi(kappa,Phys=True)+r'$',y=TitleShift)
    mkdir_p(thisdir)
    pl.savefig(thisdir+'WExp.pdf')
    pl.clf()

    

def Graphchit(Qlist,flowlist):
    ## Hard coded here....
    FormatChit = False
    flowlist = np.array(flowlist)
    coeff = (hbarc/(latspace*nx**(0.75)*nt**(0.25)))
    thisshift = 0.05
    flowpick = 6.0

    # Qboot,dump = bt.CreateBoot(Qlist,nboot,0)
    # Q2boot = np.array(Qboot)**2
    # chit = coeff*np.array(Q2boot)**(0.25)
    # chit = GetBootStats(chit)
    # pl.errorbar(flowlist,Pullflag(chit,'Avg'),Pullflag(chit,'Std'),fmt='o',label=r'$Q Boot$')

    Q2boot,dump = bt.CreateBoot(np.array(Qlist)**2,nboot,0)
    chit = coeff*np.array(Q2boot)**(0.25)
    chit = GetBootStats(chit)
    flowpick = 4.01
    # pl.errorbar(flowlist-0.02,Pullflag(chit,'Avg'),Pullflag(chit,'Std'),fmt='o',label=r'$Q^{2} Boot$')
    # pl.errorbar(flowlist[1:] ,Pullflag(chit,'Avg')[1:],Pullflag(chit,'Std')[1:],fmt='o')
    thisdir = outputdir[0] + 'graphs/Qdata/'
    taulist,tauerrlist,alphaerr,meanlist = [],[],[],[]
    for icf,(iflow,idata) in enumerate(zip(flowlist,np.rollaxis(np.array(Qlist),1))):
        # auto_gamma,Cw,Gfun,Wpick,auto_error = Gamma1D_est(idata)
        # taulist.append( auto_gamma[Wpick])
        # tauerrlist.append(auto_error[Wpick])
        # alphaerr.append(Cw)
        # meanlist.append(np.mean(idata))
        if iflow == flowpick:
            # mean, err, tint, dtint, G, W = tauint([[idata**2]], 0, True,thisdir+'AutoCorrQ2Flow'+str(flowpick))
            mean, err, tint, dtint = uWerrMine([idata**2], ff.eyeFun,ff.eyeFunDer,plot=thisdir+'AutoCorrQ2Flow'+str(flowpick))
        else:
            mean, err, tint, dtint = uWerrMine([idata**2], ff.eyeFun,ff.eyeFunDer)
            
            # auto_gamma,Cw,Gfun,Wpick,auto_error = GammaAlpha_estimate(iNNQ,iNN)
        taulist.append( tint)
        tauerrlist.append(dtint)
        alphaerr.append(err)
        meanlist.append(mean)
        # auto_gamma,Cw,Gfun,Wpick,auto_error = Gamma1D_est(idata)
        # taulist.append( auto_gamma[Wpick])
        # tauerrlist.append(auto_error[Wpick])
        # alphaerr.append(Cw)
        # meanlist.append(np.mean(idata))


    

        
    pl.errorbar(range(len(taulist)),taulist,tauerrlist,fmt='.')
    pl.axhline(0.5, color='k', linestyle='--')
    pl.ylabel(r'$ \tau_{int}$')
    
    pl.title(r'$\tau (Q^{2})$')
    pl.xlabel(r'$t_{flow}$')
    mkdir_p(thisdir)
    pl.savefig(thisdir+'IntAutoCorrQ2.pdf')
    pl.clf()
    pl.errorbar(TflowToPhys(flowlist[1:]),meanlist[1:],alphaerr[1:],fmt='.',label='Autocorr')
    pl.errorbar(TflowToPhys(flowlist[1:]+thisshift),Pullflag(Q2boot,'Avg')[1:],Pullflag(Q2boot,'Std')[1:],fmt='.',label='Bootstrap')
    pl.ylabel(r'$<Q^{2}>$')
    pl.legend()
    pl.title(r'$<Q^{2}>$')
    pl.xlabel(r'$\sqrt{8t_{f}} fm$')
    pl.savefig(thisdir+'AutoAlphaQ2.pdf')
    pl.clf()


    if FormatChit:
        pl.plot(flowlist,Pullflag(chit,'Avg'),color='b')
        pl.errorbar(flowlist,Pullflag(chit,'Avg'),Pullflag(chit,'Std'),fmt='k.',ecolor='r')
    else:
        pl.errorbar(flowlist[1:],Pullflag(chit,'Avg')[1:],Pullflag(chit,'Std')[1:],fmt='o')
    # Qavg = np.mean(np.array(Qlist)**2,axis=0)
    # Qstd = np.std(np.array(Qlist)**2,axis=0,ddof=1)
    # chitAvg = coeff*Qavg**(0.25)
    # chitStd = coeff*0.25*Qstd*Qavg**(0.25-1)
    # pl.errorbar(flowlist+0.1,chitAvg,chitStd,fmt='o',label=r'$No Boot$')
    pl.xlim(flowlist[0]-0.02,flowlist[-1]+0.1)
    pl.xlabel(r'$ \sqrt{8t_{f}} fm$')
    pl.ylabel(r'$\chi_{t}^{1/4} GeV$')
    # pl.ylim(0.14,0.3)
    pl.legend()
    pl.title(r'$ \chi_{t}^{1/4} = \frac{1}{V^{1/4}} \langle Q^2 \rangle^{1/4} \ '+GetMpi(kappa,Phys=True)+r'$',y=TitleShift)
    mkdir_p(thisdir)
    pl.savefig(thisdir+'chit.pdf')
    pl.clf()

def GraphWchit(Wlist,flowlist):
    ## Hard coded here....
    flowlist = np.array(flowlist)
    thisdir = outputdir[0] + 'graphs/Wdata/'
    coeff = (hbarc/(latspace*nx**(0.75)*nt**(0.25))**(0.5))
    thisshift = 0.05
    flowpick = 6.0
    # Wboot,dump = bt.CreateBoot(Wlist,nboot,0)
    # W2boot = np.array(Wboot)**2
    # chit = coeff*np.array(W2boot)**(0.25)
    # chit = GetBootStats(chit)
    # pl.errorbar(flowlist,Pullflag(chit,'Avg'),Pullflag(chit,'Std'),fmt='o',label=r'$W Boot$')
    
    taulist,tauerrlist,alphaerr,meanlist = [],[],[],[]
    for iflow,idata in zip(flowlist,np.rollaxis(np.array(Wlist)**2,1)):
        if iflow == flowpick:
            mean, err, tint, dtint = uWerrMine([idata**2], ff.eyeFun,ff.eyeFunDer)
        else:
            mean, err, tint, dtint = uWerrMine([idata**2], ff.eyeFun,ff.eyeFunDer,plot=thisdir+'AutoCorrW2Flow'+str(flowpick))
            # auto_gamma,Cw,Gfun,Wpick,auto_error = GammaAlpha_estimate(iNNQ,iNN)
        taulist.append( tint)
        tauerrlist.append(dtint)
        alphaerr.append(err)
        meanlist.append(mean)
        # auto_gamma,Cw,Gfun,Wpick,auto_error = Gamma1D_est(idata)
        # taulist.append( auto_gamma[Wpick])
        # tauerrlist.append(auto_error[Wpick])
        # alphaerr.append(Cw)
        # meanlist.append(np.mean(idata))
    #         PGfun = 
    #         PWpick = Wpick
    #         Pag = auto_gamma
    #         Pagerr = auto_error
    #         PCw = Cw
    # # pl.plot(range(len(PGfun[:3*PWpick+1])),Gfun[:3*PWpick+1],'.-')
    # pl.axvline(PWpick, color='k', linestyle='-')
    # pl.axhline(0.0, color='k', linestyle='--')
    # pl.ylabel(r'$ \Gamma$')
    # pl.xlabel('W')
    # pl.title('Autocorrelation of $ \\alpha\\ t_{flow}='+str(flowpick)+'$' )
    # pl.savefig(+'.pdf')
    # pl.clf()
    # pl.errorbar(range(len(Pag[:3*Wpick+1])),Pag[:3*Wpick+1],Pagerr[:3*Wpick+1],label='Error={:.2g}'.format(PCw))
    # pl.axvline(PWpick, color='k', linestyle='-')
    # pl.axhline(0.5, color='k', linestyle='--')
    # pl.ylabel(r'$ \tau_{int}$')
    # pl.xlabel('W')
    # pl.legend()
    # pl.title('Integrated Autocorrelation function of $ \\alpha\\ t_{flow}='+str(flowpick)+'$' )
    # pl.savefig(thisdir+'IntAutoCorrW2Flow'+str(flowpick)+'.pdf')
    # pl.clf()

            
    W2boot,dump = bt.CreateBoot(np.array(Wlist)**2,nboot,0)
    pl.errorbar(range(len(taulist)),taulist,tauerrlist,fmt='.')
    pl.axhline(0.5, color='k', linestyle='--')
    pl.ylabel(r'$ \tau_{int}$')
    
    pl.title(r'$\tau (W^{2})$')
    pl.xlabel(r'$t_{flow}$')
    mkdir_p(thisdir)
    pl.savefig(thisdir+'IntAutoCorrW2.pdf')
    pl.clf()
    pl.errorbar(TflowToPhys(flowlist[1:]),meanlist[1:],alphaerr[1:],fmt='.',label='Autocorr')
    pl.errorbar(TflowToPhys(flowlist[1:]+thisshift),Pullflag(W2boot,'Avg')[1:],Pullflag(W2boot,'Std')[1:],fmt='.',label='Bootstrap')
    pl.ylabel(r'$<W^{2}>$')
    pl.legend()
    pl.title(r'$<W^{2}>$')
    pl.xlabel(r'$\sqrt{8t_{f}} fm$')
    pl.savefig(thisdir+'AutoAlphaW2.pdf')
    pl.clf()

    chit = coeff*np.array(W2boot)**(0.125)
    chit = GetBootStats(chit)
    # pl.errorbar(flowlist-0.02,Pullflag(chit,'Avg'),Pullflag(chit,'Std'),fmt='o',label=r'$W^{2} Boot$')
    pl.errorbar(TflowToPhys(flowlist[1:]),Pullflag(chit,'Avg')[1:],Pullflag(chit,'Std')[1:],fmt='o')

    # Wavg = np.mean(np.array(Wlist)**2,axis=0)
    # Wstd = np.std(np.array(Wlist)**2,axis=0,ddof=1)
    # chitAvg = coeff*Wavg**(0.25)
    # chitStd = coeff*0.25*Wstd*Wavg**(0.25-1)
    # pl.errorbar(flowlist+0.1,chitAvg,chitStd,fmt='o',label=r'$No Boot$')
    pl.xlim(TflowToPhys([flowlist[1]-0.02])[0],TflowToPhys([flowlist[-1]+0.1])[0])
    pl.xlabel(r'$ \sqrt{8t_{f}} fm$')
    pl.ylabel(r'$ \frac{1}{V^{1/8}} \langle W^2 \rangle^{1/8} GeV$')
    # pl.ylim(0,0.4)
    pl.legend()
    pl.title(r'$\frac{1}{V^{1/8}} \langle W^2 \rangle^{1/8} \ '+GetMpi(kappa,Phys=True)+r'$',y=TitleShift)
    mkdir_p(thisdir)
    pl.savefig(thisdir+'chit.pdf')
    pl.clf()

    

def GraphchitKappas(Qlist,flowlist):
    ## Hard coded here....
    flowlist = np.array(flowlist)
    # coeff = (hbarc/(latspace*nx**(0.75)*nt**(0.25)))
    coeff = (1/(nx**(0.75)*nt**(0.25)))
    # Qboot,dump = bt.CreateBoot(Qlist,nboot,0)
    # Q2boot = np.array(Qboot)**2
    # chit = coeff*np.array(Q2boot)**(0.25)
    # chit = GetBootStats(chit)
    # pl.errorbar(flowlist,Pullflag(chit,'Avg'),Pullflag(chit,'Std'),fmt='o',label=r'$Q Boot$')
    tflowindex = flowlist[0].tolist().index(4.01)
    chitKappa = []
    # MpiList = [DefPionMass['1370000']*hbarc/latspace,DefPionMass['1375400']*hbarc/latspace]
    MpiList = []
    for icQ,iQ in enumerate(Qlist):
        Q2boot,dump = bt.CreateBoot(np.array(iQ)**2,nboot,0)
        chit = coeff*np.array(Q2boot)*(r0Som[0]/latspace)**4
        # MpiList.append(GetMpiNoForm(kappalist[icQ]))
        MpiList.append(GetMpiSom(kappalist[icQ],r0Som)**2)
        chit = GetBootStats(chit)
        chitKappa.append(chit[tflowindex])
    # pl.errorbar(flowlist-0.02,Pullflag(chit,'Avg'),Pullflag(chit,'Std'),fmt='o',label=r'$Q^{2} Boot$')
    pl.errorbar(MpiList,Pullflag(chitKappa,'Avg'),Pullflag(chitKappa,'Std'),fmt='o')

    # Qavg = np.mean(np.array(Qlist)**2,axis=0)
    # Qstd = np.std(np.array(Qlist)**2,axis=0,ddof=1)
    # chitAvg = coeff*Qavg**(0.25)
    # chitStd = coeff*0.25*Qstd*Qavg**(0.25-1)
    # pl.errorbar(flowlist+0.1,chitAvg,chitStd,fmt='o',label=r'$No Boot$')
    pl.xlim(0,pl.xlim()[1])
    pl.ylim(0,0.22)
    # pl.xlabel(r'$ m_{\pi} GeV $')
    # pl.ylabel(r'$\chi_{t}^{1/4} GeV$')
    pl.xlabel(r'$ (m_{\pi} r_{0})^{2} $')
    pl.ylabel(r'$\chi_{t} r_{0}^{4} $')
    # pl.ylim(0,0.4)
    pl.legend()
    thisdir = outputdir[0] + 'graphs/Qdata/'
    # pl.title(r'$ \chi_{t}^{1/4} = \frac{1}{V^{1/4}} \langle Q^2 \rangle^{1/4} $',y=TitleShift)
    pl.title(r'$ \chi_{t} = \frac{1}{V} \langle Q^2 \rangle $',y=TitleShift)
    mkdir_p(thisdir)
    pl.savefig(thisdir+'chitKappa.pdf')
    pl.clf()

def GraphchitKappasOverFlow(Qlist,flowlist,thiskappalist):
    ## Hard coded here....
    flowlist = np.array(flowlist)
    coeff = (hbarc/(latspace*nx**(0.75)*nt**(0.25)))
    
    # Qboot,dump = bt.CreateBoot(Qlist,nboot,0)
    # Q2boot = np.array(Qboot)**2
    # chit = coeff*np.array(Q2boot)**(0.25)
    # chit = GetBootStats(chit)
    # pl.errorbar(flowlist,Pullflag(chit,'Avg'),Pullflag(chit,'Std'),fmt='o',label=r'$Q Boot$')
    chitKappa = []
    thisshiftlist = [0,0.05]
    for iflowlist,iQ,ishift,ikappa in zip(flowlist,Qlist,thisshiftlist,thiskappalist):
        thisncfg= np.array(iQ).shape[0]
        Q2boot,dump = bt.CreateBoot(np.array(iQ)**2,nboot,0)
        chit = coeff*np.array(Q2boot)**(0.25)
        chit = GetBootStats(chit)
        
        pl.errorbar(TflowToPhys(iflowlist+ishift),Pullflag(chit,'Avg'),Pullflag(chit,'Std'),fmt='o',label=r'$'+GetMpi(ikappa)+r'\quad ncfg='+str(thisncfg)+'$')

    # Qavg = np.mean(np.array(Qlist)**2,axis=0)
    # Qstd = np.std(np.array(Qlist)**2,axis=0,ddof=1)
    # chitAvg = coeff*Qavg**(0.25)
    # chitStd = coeff*0.25*Qstd*Qavg**(0.25-1)
    # pl.errorbar(flowlist+0.1,chitAvg,chitStd,fmt='o',label=r'$No Boot$')
    pl.xlim(0,pl.xlim()[1])
    pl.xlabel(r'$ \sqrt{8t_{f}} fm$')
    pl.ylabel(r'$\chi_{t}^{1/4} GeV$')
    # pl.ylim(0,0.4)
    pl.legend()
    thisdir = outputdir[0] + 'graphs/Qdata/'
    pl.title(r'$ \chi_{t}^{1/4} = \frac{1}{V^{1/4}} \langle Q^2 \rangle^{1/4} $',y=TitleShift)
    mkdir_p(thisdir)
    pl.savefig(thisdir+'chitKappaOverFlow.pdf')
    pl.clf()

    
def GraphWchitKappas(Wlist,flowlist):
    ## Hard coded here....
    flowlist = np.array(flowlist)
    coeff = (hbarc/(latspace*nx**(0.75)*nt**(0.25))**(0.5))

    # Wboot,dump = bt.CreateBoot(Wlist,nboot,0)
    # W2boot = np.array(Wboot)**2
    # chit = coeff*np.array(W2boot)**(0.25)
    # chit = GetBootStats(chit)
    # pl.errorbar(flowlist,Pullflag(chit,'Avg'),Pullflag(chit,'Std'),fmt='o',label=r'$W Boot$')
    tflowindex = flowlist[0].tolist().index(4.0)
    chitKappa = []
    # MpiList = [DefPionMass['1370000']*hbarc/latspace,DefPionMass['1375400']*hbarc/latspace]
    MpiList = []
    for icW,iW in enumerate(Wlist):
        W2boot,dump = bt.CreateBoot(np.array(iW)**2,nboot,0)
        chit = coeff*np.array(W2boot)**(0.125)
        MpiList.append(GetMpiNoForm(kappalist[icW]))
        chit = GetBootStats(chit)
        chitKappa.append(chit[tflowindex])
    # pl.errorbar(flowlist-0.02,Pullflag(chit,'Avg'),Pullflag(chit,'Std'),fmt='o',label=r'$W^{2} Boot$')
    pl.errorbar(MpiList,Pullflag(chitKappa,'Avg'),Pullflag(chitKappa,'Std'),fmt='o')

    # Wavg = np.mean(np.array(Wlist)**2,axis=0)
    # Wstd = np.std(np.array(Wlist)**2,axis=0,ddof=1)
    # chitAvg = coeff*Wavg**(0.25)
    # chitStd = coeff*0.25*Wstd*Wavg**(0.25-1)
    # pl.errorbar(flowlist+0.1,chitAvg,chitStd,fmt='o',label=r'$No Boot$')
    pl.xlim(0,pl.xlim()[1])
    pl.ylim(0,pl.ylim()[1]*1.2)
    pl.xlabel(r'$ m_{\pi} GeV $')
    pl.ylabel(r'$\chi_{t}^{1/8} GeV$')
    # pl.ylim(0,0.4)
    pl.legend()
    thisdir = outputdir[0] + 'graphs/Wdata/'
    pl.title(r'$ \chi_{t}^{1/8} = \frac{1}{V^{1/8}} \langle W^2 \rangle^{1/8} $',y=TitleShift)
    mkdir_p(thisdir)
    pl.savefig(thisdir+'chitKappa.pdf')
    pl.clf()

    
def GraphWchitKappasOverFlow(Wlist,flowlist,thiskappalist):
    ## Hard coded here....
    flowlist = np.array(flowlist)
    coeff = (hbarc/(latspace*nx**(0.75)*nt**(0.25))**(0.5))
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    # Qboot,dump = bt.CreateBoot(Wlist,nboot,0)
    # Q2boot = np.array(Qboot)**2
    # chit = coeff*np.array(Q2boot)**(0.25)
    # chit = GetBootStats(chit)
    # pl.errorbar(flowlist,Pullflag(chit,'Avg'),Pullflag(chit,'Std'),fmt='o',label=r'$Q Boot$')
    chitKappa = []
    thisshiftlist = [0,0.05]
    for iflowlist,iW,ishift,ikappa in zip(flowlist,Wlist,thisshiftlist,thiskappalist):
        thisshift = ishift*0.1
        thiscol = thiscolcyc.next()
        thisncfg = np.array(iW).shape[0]
        W2boot,dump = bt.CreateBoot(np.array(iW)**2,nboot,0)
        # SF = np.log(np.array(TflowToPhys(iflowlist)))
        # SF = np.log(np.array(TflowToPhys(iflowlist))+2)
        chit = coeff*np.array(W2boot)**(0.125)
        xdata = TflowToPhys(iflowlist)
        # ydataboot = chit*SF
        ydataboot = chit
        ydataboot = GetBootStats(ydataboot)
        Nearmin,Nearmax = find_nearest(xdata,FitWchiMinMax[0]),find_nearest(xdata,FitWchiMinMax[1])
        fitxdata,fitydata = np.array(xdata[Nearmin:Nearmax]),ydataboot[Nearmin:Nearmax]
        thisfitfun = ChitFitFun
        
        [fitBoot,fitAvg,fitChi] = FitBoots(np.array(fitydata),np.array([fitxdata]),thisfitfun)

        ydataAvg,ydataErr = Pullflag(ydataboot,'Avg'),Pullflag(ydataboot,'Std')        
        pl.errorbar(xdata+thisshift,ydataAvg,ydataErr,fmt='o',color=thiscol,label=r'$'+GetMpi(ikappa)+r'\quad ncfg='+str(thisncfg)+'$')

        fityplot = []
        for ibootlist in np.swapaxes(Pullflag(fitBoot,'values'),0,1):
            fityplot.append(thisfitfun(np.array([fitxdata]),ibootlist))
        fityAvg,fityStd = np.mean(fityplot,axis=0),np.std(fityplot,axis=0)
        fityup,fitydown = fityAvg+fityStd,fityAvg-fityStd
        # thislab = r'$'+MakeValAndErr(fitBoot[1].Avg,fitBoot[1].Std)+r' + \frac{'+MakeValAndErr(fitBoot[0].Avg,fitBoot[0].Std)+r'}{\sqrt{8t_{f}}}$'
        thislab = r'$'+'\ '.join([' A_{'+str(icb)+'} = '+MakeValAndErr(iboot.Avg,iboot.Std) for icb,iboot in enumerate(fitBoot)])+r'$'
        pl.plot(fitxdata+thisshift,fityAvg,color=thiscol,label=thislab)
        pl.fill_between(fitxdata+thisshift,fityup,fitydown,color=thiscol,alpha=thisalpha,edgecolor='none')

        
        # chitdivlog = Pullflag(chit,'Avg')/np.log(np.array(iflowlist))
        # chitdivlogErr = Pullflag(chit,'Std')/np.log(np.array(iflowlist))
        # chitdivlog = Pullflag(chit,'Avg')/np.log(TflowToPhys(iflowlist))
        # chitdivlogErr = Pullflag(chit,'Std')/np.log(TflowToPhys(iflowlist))
        # chitdivlog = Pullflag(chit,'Avg')/np.log(np.array(TflowToPhys(iflowlist))+2)
        # chitdivlogErr = Pullflag(chit,'Std')/np.log(np.array(TflowToPhys(iflowlist))+2)
        # chitdivlog = Pullflag(chit,'Avg')/np.log(np.array(TflowToPhys(iflowlist))+2)
        # chitdivlogErr = Pullflag(chit,'Std')/np.log(np.array(TflowToPhys(iflowlist))+2)
        # pl.errorbar(TflowToPhys(iflowlist+ishift),chitdivlog,chitdivlogErr,fmt='o',color=thiscol,alpha=0.5,label=r'$ \chi / log(\sqrt{8t_{f}})$')
        # chitdivlog = Pullflag(chit,'Avg')*np.log(np.array(TflowToPhys(iflowlist))+2)
        # chitdivlogErr = Pullflag(chit,'Std')*np.log(np.array(TflowToPhys(iflowlist))+2)
        # pl.errorbar(TflowToPhys(iflowlist+ishift),chitdivlog,chitdivlogErr,fmt='s',color=thiscol,alpha=0.5,label=r'$ \chi * log(\sqrt{8t_{f}})$')

    # Qavg = np.mean(np.array(Wlist)**2,axis=0)
    # Qstd = np.std(np.array(Wlist)**2,axis=0,ddof=1)
    # chitAvg = coeff*Qavg**(0.25)
    # chitStd = coeff*0.25*Qstd*Qavg**(0.25-1)
    # pl.errorbar(flowlist+0.1,chitAvg,chitStd,fmt='o',label=r'$No Boot$')
    pl.xlim(0,pl.xlim()[1])
    pl.xlabel(r'$ \sqrt{8t_{f}} fm$')
    # pl.ylabel(r'$\frac{\chi_{t}^{1/8}}{ log\left(100*\sqrt{8t_{f}}\right)} GeV$')
    # pl.ylabel(r'$\chi_{t}^{1/8} log\left(\sqrt{8t_{f}}\right) GeV$')
    pl.ylabel(r'$\chi_{t}^{1/8} GeV$')
    # pl.ylim(0,0.4)
    pl.legend()
    # pl.legend(loc='lower right')
    thisdir = outputdir[0] + 'graphs/Wdata/'
    FitFunLatex = r'\chi_{t,fit}^{1/8} = \frac{A_{0}\sqrt{8t_{f}} + A_{1}}{\sqrt{8t_{f}}\left(log(\sqrt{8t_{f}}) + A_{2}\right)}'
    pl.title(r'$ \chi_{t}^{1/8} = \frac{1}{V^{1/8}} \langle W^2 \rangle^{1/8} \quad '+FitFunLatex+'$',y=TitleShift)
    mkdir_p(thisdir)
    pl.savefig(thisdir+'chitKappaOverFlow.pdf')
    pl.clf()


    
    
def GraphQ2Hist(Qlist,thisflow):
    ## Hard coded here....
    Q2list = np.array(Qlist)[:,thisflow]**2

    Q2min,Q2max = 0,300
    pl.hist(Q2list,bins=np.arange(Q2min,Q2max,10),histtype=HistType,normed=False)
    pl.xlim(Q2min,Q2max)
    pl.xlabel(r'$\langle Q^2 \rangle $')
    pl.ylabel(r'$ N $')
    thisdir = outputdir[0] + 'graphs/Qdata/'
    pl.title(r'$ Histogram\ \langle Q^2 \rangle \ '+GetMpi(kappa,Phys=True)+r'$',y=TitleShift)
    mkdir_p(thisdir)
    pl.savefig(thisdir+'Q2Hist.pdf')
    pl.clf()

    
def GraphQLines(Qlist,flowlist,cfglist):
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    for icfg in cfglist:
        pl.plot(TflowToPhys(flowlist),Qlist[icfg],color=thiscolcyc.next(),label='Q cfg='+str(icfg))
    pl.xlim(TflowToPhys([flowlist[0]-0.1])[0],TflowToPhys([flowlist[-1]+0.1])[0])
    pl.xlabel(r'$ \sqrt{8t_{f}} fm$')
    pl.ylabel(r'$ Q $')
    pl.legend()
    thisdir = outputdir[0] + 'graphs/Qdata/'
    pl.title(r'$ Q \left( icfg \right) \ '+GetMpi(kappa,Phys=True)+r'$',y=TitleShift)
    pl.tight_layout()
    mkdir_p(thisdir)
    pl.savefig(thisdir+'QLines.pdf')
    pl.clf()


def GraphW2Hist(Wlist,thisflow):
    ## Hard coded here....
    W2list = np.array(Wlist)[:,thisflow]**2

    W2min,W2max = 0,300
    pl.hist(W2list,bins=np.arange(W2min,W2max,10),histtype=HistType,normed=False)
    pl.xlim(W2min,W2max)
    pl.xlabel(r'$\langle W^2 \rangle $')
    pl.ylabel(r'$ N $')
    thisdir = outputdir[0] + 'graphs/Wdata/'
    pl.title(r'$ Histogram\ \langle W^2 \rangle \ '+GetMpi(kappa,Phys=True)+r'$',y=TitleShift)
    mkdir_p(thisdir)
    pl.savefig(thisdir+'W2Hist.pdf')
    pl.clf()

    
def GraphWLines(Wlist,flowlist,cfglist):
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    for icfg in cfglist:
        pl.plot(TflowToPhys(flowlist),Wlist[icfg],color=thiscolcyc.next(),label='W cfg='+str(icfg))
    pl.xlim(TflowToPhys([flowlist[0]-0.1])[0],TflowToPhys([flowlist[-1]+0.1])[0])
    pl.xlabel(r'$ \sqrt{8t_{f}} fm$')
    pl.ylabel(r'$ W $')
    pl.legend()
    thisdir = outputdir[0] + 'graphs/Wdata/'
    pl.title(r'$ W \left( icfg \right) \ '+GetMpi(kappa,Phys=True)+r'$',y=TitleShift)
    pl.tight_layout()
    mkdir_p(thisdir)
    pl.savefig(thisdir+'WLines.pdf')
    pl.clf()

    
    

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
from SetLists import *
from FFParams import *
from CreateCombs import CreateDictOldCombs

colourset8 = [ '#000080','#B22222','#00B800','#8B008B', '#0000FF','#FF0000','#000000','#32CD32','#FF0066']
markerset = ['o','s','^','v','d']
shiftmax = 3
shiftper = 0.05
shiftset = [0]
for ish in np.arange(1,shiftmax+1): shiftset += [-ish*shiftper,ish*shiftper]
incr = 0.01
thisalpha = 0.3

xlabshift = 0.05
manylim = (0,2.5)
lowlim = manylim[0]
highlim = manylim[1]

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

def GetPlotIters():
    return itertools.cycle(markerset),itertools.cycle(colourset8),itertools.cycle(shiftset)


## Add special exclusions here: ##
def GraphCondit(iDS,igamma,iq,thisMeth,iSet):
    graphthis = True
    # if 'OSF' in thisMeth and 'tsink26sm32' in iSet: graphthis = False
    if 'Fits' in thisMeth and 'tsink32' in iSet and 'to18dt2' in iSet: graphthis = False
    if 'Fits' in thisMeth and 'tsink26' in iSet and 'to18dt2' in iSet: graphthis = False
    if 'SumMeth' in thisMeth and 'fitr2-4' in iSet: graphthis = False
    return graphthis

##FIX STRINGS tsink and sm stuff##
def PlotXlabs(method,thissetlist,xminmax):
    xmin,xmax = xminmax
    xmid = (xmax+xmin-2)/float(2)
    tsinklist,smlist = GetTsinkSmLists(thissetlist)
    if 'SumMeth' in method:
        line1 = 'Summation'
        line2 = 'Method ' + smlist[0]
        line3 = []
        for isub in thissetlist:
            if '0-4' in isub and 'All' not in line3:
                line3.append('All')
            elif '1-4' in isub and 'fr1' not in line3:
                line3.append('fr1')
            # elif '2-4' in isub and 'fr2' not in line3:
            #     line3.append('fr2')
    elif 'FitsTsink' in method:
        line1 = 'Fits'
        line2  = smlist[0]
        # line3 = tsinklist
        line3 = []
    elif 'FitsSm' in method:
        line1 = 'Fits'
        line2 = tsinklist[0]
        # line3 = smlist
        line3 = []
    elif 'FitsVar' in method:
        line1 = 'Var'
        line2 = SplitToDt(smlist[0])[0]
        line3 = [SplitToDt(smlist[0])[1].replace('dt','\Delta t')]
    elif 'TSF' in method:
        line1 = '2SF'
        if 'CM' in method:
            if tsinklist[0] == None:
                line2 = 'tsink29'
            else:
                line2 = tsinklist[0]
            line3 = smlist
        elif 'Tsink' in method:
            line2 = 'All'
            line3 = smlist
        elif 'test32' in method:
            line2 = 'tst'
            line3 = smlist
        elif 'Small' in method:
            line2 = 'Sml'
            line3 = smlist
        else:
            raise IOError('method not implemented' + method)
    elif 'OSF' in method:
        line1 = '1SF'
        if 'CM' in method:
            line2 = tsinklist[0]
            line3 = map(ReducedVar,smlist)
            if any(['PoF' in il for il in line3]):
                line3 = ['PoF']*(len(PoFTSinkList)-1) + line3
        elif 'Tsink' in method:
            line2 = 'All ' + smlist[0]
            line3 = tsinklist
        else:
            raise IOError('method not implemented' + method)
    # ysize = pl.ylim()[1]-pl.ylim()[0]
    # ylow = pl.ylim()[0]
    ysize = manylim[1]-manylim[0]
    ylow = manylim[0]
    line1box = (xmid,LegLab(line1))
    line2box = (xmid,LegLab(line2))
    line3box = []
    for il,iline in enumerate(line3):
        thisx = ((xmax-xmin+2)*(il+1))/(float(len(line3)+1)) + xmin -2
        line3box.append((thisx,LegLab(iline)))
    return [line1box],[line2box],line3box

    

def CreateMethodSetList(thisMethodList,setlist):
    mslout = OrderedDict()
    for iMeth in thisMethodList:
        if 'Fits' in iMeth :
            mslout[iMeth+'Tsink'] = []
            mslout[iMeth+'Sm'] = []
            mslout[iMeth+'Var'] = []
            for iset in setlist:
                if CheckCut(iset,FitCutPicked):
                    if iMeth in iset and 'sm32' in iset and 'cut2':
                        mslout[iMeth+'Tsink'].append(iset.replace(iMeth,''))
                    if iMeth in iset and 'tsink29' in iset and 'sm' in iset:
                        mslout[iMeth+'Sm'].append(iset.replace(iMeth,''))
                    if iMeth in iset and 'to' in iset:
                        mslout[iMeth+'Var'].append(iset.replace(iMeth,''))
        # elif 'OSF' in iMeth:
        #     mslout[iMeth] = []
        #     for iset in setlist:
        #         if OSFCutPicked in iset and iMeth in iset:
        #             mslout[iMeth].append(iset.replace(iMeth,''))                    
        # elif 'TSF' in iMeth:
        #     mslout[iMeth] = []
        #     for iset in setlist:
        #         if TSFCutPicked in iset and iMeth in iset:
        #             mslout[iMeth].append(iset.replace(iMeth,''))                    
        else:
            mslout[iMeth] = []
            for iset in setlist:
                if iMeth in iset:
                    mslout[iMeth].append(iset.replace(iMeth,''))
    return mslout
    
# data { collection , gamma , mom }
def PlotSummaryMethods(data,thisMethodSetList,iDS,igamma,iq,outputdir,dirpref=''):
    xvalues,Xlabs,Xbox1,Xbox2,Xbox3 = 0,[],[],[],[]
    thissymcyc,thiscolcyc,thisshiftcyc = GetPlotIters()
    for iMeth,thisSetList in thisMethodSetList.iteritems():
        if len(thisSetList) == 0: continue
        prevxval = xvalues
        nextxvalue,thisxlab = PlotSummarySet(iDS,igamma,iq,data,iMeth,thisSetList,prevxval,thiscolcyc.next(),thissymcyc.next())
        if prevxval == nextxvalue: continue
        xvalues = 1+nextxvalue
        Xlabs += thisxlab
        thisXbox1,thisXbox2,thisXbox3 = PlotXlabs(iMeth,thisSetList,(prevxval,xvalues))
        Xbox1 += thisXbox1
        Xbox2 += thisXbox2
        Xbox3 += thisXbox3
    if xvalues == 0 : return
    pl.xticks(range(0,xvalues),Xlabs,rotation=90)
    ylow,yhigh = pl.ylim()
    ysize = yhigh-ylow
    pl.ylim(ylow-(3*ysize*xlabshift),yhigh)
    pl.ylim(max(lowlim,pl.ylim()[0]),min(highlim,pl.ylim()[1]))
    ylow,yhigh = pl.ylim()
    ysize = yhigh-ylow
    for xmid,line1 in Xbox1:
        pl.text(xmid,ylow+(3*ysize*xlabshift),line1,horizontalalignment='center')
    for xmid,line2 in Xbox2:
        pl.text(xmid,ylow+(2*ysize*xlabshift),line2,horizontalalignment='center')
    for xmid,line3 in Xbox3:
        pl.text(xmid,ylow+(ysize*xlabshift),line3,horizontalalignment='center')
    pl.xlabel('Methods')
    pl.ylabel('Value')
    pl.title(TitleFix('SummaryPlot ' +iDS + ' ' + igamma + ' ' + iq))
    pl.xlim(-1,xvalues-1)
    pl.grid(False,axis='x')
    if 'FF' in igamma:
        thisgammadir = dirpref + '/'
        thisdir = outputdir + 'graphs/Summarys/'+thisgammadir+'/'+iq+'/'
        mkdir_p(thisdir)
        pl.savefig(thisdir+'SummaryPlot'+dirpref+iDS+igamma+iq+'.pdf')
    else:
        thisgammadir = CreateOppDir(iDS+igamma)
        thisdir = outputdir + 'graphs/Summarys/'+thisgammadir+'/qsqrd'+str(qsqrdstr(iq))+'/'
        mkdir_p(thisdir)
        pl.savefig(thisdir+'SummaryPlot'+dirpref+iDS+igamma+qstrTOqcond(iq)+'.pdf')
    pl.clf()


        
def PlotSummarySet(iDS,igamma,iq,data,thisMeth,thisSetList,xstart,col,sym):
    dataval,dataerr,Xlables = [],[],[]
    keylist = SortMySet(data.keys())[0]
    for datakey in keylist:
        if iDS == False:
            idata = data
        elif iDS in DefDSList:
            idata = CreateDictOldCombs(data[datakey],[]) 
        else:
            idata = CreateDictOldCombs(data[datakey],[iDS]) 
        for iSet in thisSetList:
            methcomp = thisMeth
            if 'Fits' in thisMeth: methcomp = 'Fits'
            if iSet in datakey.replace(methcomp,'') and methcomp in datakey: 
                if iDS == False:
                    if not CheckDict(idata,iDS,igamma,iq): continue
                    if not GraphCondit(iDS,igamma,iq,methcomp,iSet): continue
                else:
                    print idata.keys(),igamma
                    print idata[igamma].keys(), iq
                    if not CheckDict(idata,igamma,iq): continue
                    if not GraphCondit('',igamma,iq,methcomp,iSet): continue
                Xlables.append(LabToXaxis(iSet,thisMeth))
                if 'Avg' in idata[igamma][iq].keys() and 'Std' in idata[igamma][iq].keys():
                    dataval.append(abs(idata[igamma][iq]['Avg']))
                    dataerr.append(idata[igamma][iq]['Std'])
                else:
                    dataval.append(abs(idata[igamma][iq]['Boot'].Avg))
                    dataerr.append(idata[igamma][iq]['Boot'].Std)
    if len(dataval) == 0: return xstart,[]
    xdata = range(xstart,xstart+len(dataval))
    pl.errorbar(xdata,dataval,dataerr,color=col,fmt=sym)
    pl.axvline(xdata[-1]+1, color='k')
    return xstart+len(xdata),Xlables+['']


def ReadAndPlotSummary(thisMethodList,thisGammaList,thisSetList,thisMomList,thisCombList):
    data,massdata = ExtractValues(outputdir,thisGammaList,thisSetList,thisMethodList,thisMomList=thisMomList)
    thisMethodSetList = CreateMethodSetList(thisMethodList,data.keys())
    for igamma in thisGammaList:
        if 'twopt' in igamma: continue
        if 'sing' in igamma:
            thisDSList = ['sing']
        elif 'doub' in igamma:
            thisDSList = thisCombList+['doub']
        else:
            thisDSList = []
        for iDS in thisDSList:
            for imom in thisMomList:
                gammastrip = igamma.replace('doub','').replace('sing','') 
                print 'Plotting: ', iDS , gammastrip , imom , ' Complete  '
                PlotSummaryMethods(data,thisMethodSetList,iDS,gammastrip,imom,outputdir)


def PlotFFSummary(thisSL,thiscurr,currdata):
    if 'RF' in MethodList: MethodList.remove('RF')
    thisMethodSetList = CreateMethodSetList(MethodList,thisSL)
    for iFF in NoFFList[thiscurr]:
        for iqsqrd in QsqrdSet:
            PlotSummaryMethods(currdata,thisMethodSetList,'',iFF,iqsqrd,outputdir,dirpref=thiscurr) 

# def PlotSummary(thisdata,thisopfile):
#     pl.rcParams.update({'axes.labelsize' : 15})
#     XSubAxis,XSubSubAxis,XAxis = [],[],[]
#     colorcyc = itertools.cycle(colourset8)
#     YaxisPoints,YaxisErr,YaxisBoot,Xlables = [],[],[],[]
#     if 'Fit' in thisdata.keys():
#         for fitkey,ifit in thisdata['Fit'].iteritems():
#             if 'tsink29sm' in fitkey:
#                 YaxisPoints.append(ifit['Avg'])
#                 YaxisErr.append(ifit['Std'])
#                 Xlables.append(LabToXaxis(str(fitkey),'SS'))
#                 if 'FitBoot' in ifit.keys():
#                     YaxisBoot.append(ifit['Boot'])
#     midpoint = len(Xlables)-(len(YaxisPoints)-1)/2.0
#     XAxis.append(('$Fits$',midpoint))
#     XSubAxis.append(('$Smears$',midpoint))
#     XSubSubAxis.append(('$t=13$',midpoint))
#     pl.errorbar(range(1+len(Xlables)-len(YaxisPoints),len(Xlables)+1),np.abs(YaxisPoints),YaxisErr,fmt='o',color=next(colorcyc))
#     Xlables.append('')
#     pl.axvline(len(Xlables), color='k')

#     YaxisPoints,YaxisErr,YaxisBoot = [],[],[]
#     if 'Fit' in thisdata.keys():
#         for fitkey,ifit in thisdata['Fit'].iteritems():
#             if 'sm32' in fitkey:
#                 YaxisPoints.append(ifit['Avg'])
#                 YaxisErr.append(ifit['Std'])
#                 Xlables.append(LabToXaxis(str(fitkey),'TSink'))
#                 if 'FitBoot' in ifit.keys():
#                     YaxisBoot.append(ifit['Boot'])
#     pl.errorbar(range(1+len(Xlables)-len(YaxisPoints),len(Xlables)+1),np.abs(YaxisPoints),YaxisErr,fmt='o',color=next(colorcyc))
#     midpoint = len(Xlables)-(len(YaxisPoints)-1)/2.0
#     XAxis.append(('$Fits$',midpoint))
#     XSubAxis.append(('$Sink\ Times$',midpoint))
#     XSubSubAxis.append(('$sm1$',midpoint))
#     Xlables.append('')
#     pl.axvline(len(Xlables), color='k')

#     if 'Sum' in thisdata.keys():
#         # YaxisPoints,YaxisErr,YaxisBoot = [],[],[]
#         thiscol = next(colorcyc)
#         prevx = len(Xlables)
#         for sumkey,isum in thisdata['Sum'].iteritems():
#             # for dictpar1,idata1 in isum.iteritems():
#             imrk = 0
#             for dictpar2,idata2 in isum[isum.keys()[0]].iteritems():
#                 if 'sl' in dictpar2 and '-4' in dictpar2:
#                     YaxisPoints,YaxisErr,YaxisBoot = [],[],[]
#                     for dictpar1,idata1 in isum.iteritems():
#                         YaxisPoints.append(isum[dictpar1][dictpar2]['Avg'])
#                         YaxisErr.append(isum[dictpar1][dictpar2]['Std'])
#                         Xlables.append(LabToXaxis(str(sumkey)+' '+dictpar1+' '+dictpar2,'sum'))
#                         if 'SumBoot' in isum[dictpar1][dictpar2].keys():
#                             YaxisBoot.append(isum[dictpar1][dictpar2]['Boot'])
#                     submids = len(Xlables)-(len(YaxisPoints)-1)/2.0
#                     lastlen = len(YaxisPoints)
#                     XSubSubAxis.append((LabToXaxis(dictpar2,'sub'),submids))
#                     pl.errorbar(range(1+len(Xlables)-len(YaxisPoints),len(Xlables)+1),np.abs(YaxisPoints),YaxisErr,fmt=markerset[imrk],color=thiscol)
#                     imrk += 1
#         midpoint = (len(Xlables)+prevx+1)/2.0
#         XAxis.append(('$Summation$',midpoint))
#         XSubAxis.append(('$Method\ sm1$',midpoint))
#         Xlables.append('')
#         pl.axvline(len(Xlables), color='k')

#     for iFF in TSFFileFlags:
#         if 'TSF'+iFF in thisdata.keys():
#             if iFF not in ['CM','Tsink','test32']: continue
#             prevx = len(Xlables)
#             thiscol = next(colorcyc)
#             imrk = 0
#             for TSFkey,iTSF in thisdata['TSF'+iFF].iteritems():
#                 if 'sm' not in TSFkey: continue
#                 YaxisPoints,YaxisErr,YaxisBoot = [],[],[]
#                 for dictpar1,idata1 in iTSF['B00'][TSFfitr[TSFkey]].iteritems():
#                     if int(dictpar1) > 1 and int(dictpar1) < 5:
#                         YaxisPoints.append(idata1['Avg'])
#                         YaxisErr.append(idata1['Std'])
#                         Xlables.append(LabToXaxis(iFF+ ' '+str(TSFkey)+' cut'+dictpar1,'TSF'+iFF))
#                 if len(YaxisPoints) > 0:
#                     submids = len(Xlables)-(len(YaxisPoints)-1)/2.0
#                     XSubSubAxis.append((LabToXaxis(TSFkey,'sub'),submids))
#                     xvals = range(1+len(Xlables)-len(YaxisPoints),len(Xlables)+1)
#                     pl.errorbar(xvals,np.abs(YaxisPoints),YaxisErr,fmt=markerset[imrk],color=thiscol)
#                     imrk += 1
#             midpoint = (len(Xlables)+prevx+1)/2.0
#             if iFF == 'CM':
#                 XAxis.append(('$2SF$',midpoint))
#                 XSubAxis.append(('$t=13$',midpoint))
#             elif iFF == 'Tsink':
#                 XAxis.append(('$2SF$',midpoint))
#                 XSubAxis.append(('$All$',midpoint))
#                 XSubSubAxis.append(('$sm1$',midpoint))
#             elif iFF == 'test32':
#                 XAxis.append(('$2SF$',midpoint))
#                 XSubAxis.append(('$t\\neq 10,13$',midpoint))
#                 XSubSubAxis.append(('$sm1$',midpoint))
#             Xlables.append('')
#             pl.axvline(len(Xlables), color='k')

#     YaxisPoints,YaxisErr,YaxisBoot = [],[],[]
#     if 'Fit' in thisdata.keys():
#         for fitkey,ifit in thisdata['Fit'].iteritems():
#             if 'state' in fitkey:
#                 YaxisPoints.append(ifit['Avg'])
#                 YaxisErr.append(ifit['Std'])
#                 Xlables.append(LabToXaxis(str(fitkey),'CM'))
#                 if 'FitBoot' in ifit.keys():
#                     YaxisBoot.append(ifit['Boot'])
#                 if 'tsink29state1to18dt2' in fitkey:
#                     bandthis = (ifit['Avg'],ifit['Std'])
#     midpoint = len(Xlables)-(len(YaxisPoints)-1)/2.0
#     XAxis.append(('$CM$',midpoint))
#     XSubAxis.append(('$t_{0}=2$',midpoint))
#     XSubSubAxis.append(('$\Delta t=2$',midpoint))
#     thiscol = next(colorcyc)
#     pl.errorbar(range(1+len(Xlables)-len(YaxisPoints),len(Xlables)+1),np.abs(YaxisPoints),YaxisErr,fmt='o',color=thiscol)

#     bandup,banddown = np.abs(bandthis[0])+bandthis[1],np.abs(bandthis[0])-bandthis[1]
#     pl.fill_between([0,len(Xlables)+1],[bandup,bandup],[banddown,banddown],color=thiscol,alpha=thisalpha,edgecolor='none')

#     pl.xlabel('Methods')
#     pl.ylabel('Value')
#     pl.xlim(pl.xlim()[0],pl.xlim()[1]-1)
#     pl.xticks(range(1,len(Xlables)+1),Xlables,rotation=90)
#     pl.tight_layout()
#     ysize = pl.ylim()[1]-pl.ylim()[0]
#     pl.ylim(pl.ylim()[0]-ysize*2*ypergraph,pl.ylim()[1])
#     ysize = pl.ylim()[1]-pl.ylim()[0]
#     for (xstr,xmid) in XAxis:
#         pl.text(xmid,pl.ylim()[0]+ysize*(shiftpersummary+2*ypertext),xstr,horizontalalignment='center')
#     for (xstr,xmid) in XSubAxis:
#         pl.text(xmid,pl.ylim()[0]+ysize*(shiftpersummary+ypertext),xstr,horizontalalignment='center')
#     for (xstr,xmid) in XSubSubAxis:
#         pl.text(xmid,pl.ylim()[0]+ysize*shiftpersummary,xstr,horizontalalignment='center')
#     pl.title(TitleFix(thisdata['Title'].replace('_','\ ').replace('q','\ q')
#                       + '\ Summary'),y=titleloc)

#     # pl.grid(True,axis='y')
#     pl.savefig(thisopfile+'.pdf')
#     pl.clf()
    


# def ReadAndPlotSummary(thisGammaList,thisCombStateList,FilePrefix,SetTitle=False):
#     RFList = ['Fit','FitBoot','Sum','SumBoot','TwoStateFit']
#     SetData = ReadSetDict(thisGammaList,thisCombStateList,RFList,outputdir)
#     for thisgamma,gammaData in SetData.iteritems():
#         gammadirout = outputdir+'graphs/'+thisgamma+'/'
#         mkdir_p(gammadirout)
#         newkey = thisgamma.replace('doub','')
#         singkey = thisgamma.replace('doub','sing')
#         for mom,data in gammaData.iteritems():
#             if SetTitle: data['Title'] = FilePrefix
#             PlotSummary(deepcopy(data),gammadirout+FilePrefix+mom.replace(' ',''))
#             if 'doub' in thisgamma and singkey in SetData.keys() and mom == 'q = 0 0 0':
#                 print 'found ' + newkey
#                 gammaDatasing = SetData[singkey]
#                 newdict = OppDicts(data,gammaDatasing[mom],'-',RFList,newkey,thisCombStateList)
#                 mkdir_p(outputdir+'graphs/'+newkey+'/')
#                 if SetTitle: newdict['Title'] = FilePrefix
#                 PlotSummary(deepcopy(newdict),outputdir+'graphs/'+newkey+'/'+FilePrefix)

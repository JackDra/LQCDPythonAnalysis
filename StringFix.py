#!/usr/bin/env python

from sets import Set
import numpy as np
import os
# from MiscFuns import *
from Params import *
from CombParams import CombList
from XmlFormatting import *
import re

def GetCfgNumb(name):
    cfgnum = re.search('-00....',name).group()
    return int(cfgnum.replace('-00',''))
               
def MakeValAndErr(Avg,Err,Dec=2):
    try:
        int(np.floor(np.log10(Err)))-(Dec-1)
    except:
        return 'Err'
    ErrMag = int(np.floor(np.log10(Err)))-(Dec-1)
    sigErr = round(Err, -ErrMag)
    if ErrMag > 0:
        sigAvg = round(Avg,-ErrMag)        
        return ('{0:.'+str(ErrMag)+'f}({1})').format(sigAvg,sigErr)
    else:
        sigErr = round(Err, -ErrMag)
        strErr = int(sigErr*10**(-ErrMag))
        sigAvg = round(Avg,-ErrMag)
        return ('{0:.'+str(-ErrMag)+'f}({1})').format(sigAvg,strErr)

        
# def NoSm(thestring): return thestring
# def NoTSink(thestring): return thestring
def NoSm(thestring):
    for i in DefSmList:
        thestring = thestring.replace(i,'')
    return thestring

def NoTSink(thestring):
    for i in AllTSinkStrList + CMTSinkStrList:
        thestring = thestring.replace(i,'')
    for i in AllTSinkStrListVar:
        thestring = thestring.replace(i,'')
    return thestring

def NoTSource(thestring):
    for i in PoFTsrcstrList:
        thestring = thestring.replace(i,'')
    return thestring

def NoCM(thestring):
    for itvar in PoFTvarList:
        thestring = thestring.replace(itvar,'')
    for itvar in REvecTvarList:
        thestring = thestring.replace(itvar,'')
    for itvar in DefTvarList:
        thestring = thestring.replace(itvar,'')
    return thestring

def ReducedVar(thestring):
    for itvar in PoFTvarList:
        thestring = thestring.replace(itvar,'PoF')
    for itvar in REvecTvarList:
        thestring = thestring.replace(itvar,'CM')
    for itvar in DefTvarList:
        thestring = thestring.replace(itvar,'CM')
    return thestring


def ProperTsink(thestring):
    for tsinkstr,itsink in zip(AllTSinkStrListVar,AllTSinkListVar):
        thestring = thestring.replace(tsinkstr,'SPACEt'+str(itsink-tsource))
    return thestring

def ProperSmear(thestring):
    # for i,ism in enumerate(DefSmList):
    #     thestring = thestring.replace(ism,'sm'+str(i+1))
    return thestring.replace('ism','SPACEism').replace('jsm','SPACEjsm')

def SplitToDt(tvar):
    try:
        toval = re.search('to.*dt',tvar)
        toval = toval.group().replace('dt','')
    except:
        return None,None
    try:
        dtval = re.search('dt.*',tvar).group()
    except:
        return None,None
    return toval,dtval

def ProperCM(thestring):
    thestring = thestring.replace('state1','')
    if kappa == 1375400:
        if 'PoF0to1dt3' in thestring :
            return thestring.replace('PoF0to1dt3','ism16')
        if 'PoF0to2dt3' in thestring :
            return thestring.replace('PoF0to2dt3','ism32')
        if 'PoF0to4dt3' in thestring :
            return thestring.replace('PoF0to4dt3','ism64')
    elif kappa == 1370000:
        if 'PoF0to1dt3' in thestring :
            return thestring.replace('PoF0to1dt3','ism16SPACEjsm64')
        if 'PoF0to2dt3' in thestring :
            return thestring.replace('PoF0to2dt3','ism32SPACEjsm64')
        if 'PoF0to4dt3' in thestring :
            return thestring.replace('PoF0to4dt3','ism64SPACEjsm64')        
    # thestring = thestring.replace('CM','CMSPACE')
    # thestring = thestring.replace('REvec','CMSPACE')
    # thestring = thestring.replace('PoF','PoFSPACE')
    thestring = thestring.replace('CM','VarSPACE')
    thestring = thestring.replace('REvec','VarSPACE')
    thestring = thestring.replace('PoF0','VarSPACE')
    if OverDetRun:
        toval = re.search('to.-.',thestring.replace('Proton','').replace('Neutron',''))
    else:
        toval = re.search('to.',thestring.replace('Proton','').replace('Neutron',''))
    try:
        toval = toval.group().replace('to','').replace('dt','')
    except:
        return thestring
    if 'to'+toval in thestring :
        if OverDetRun:
            if '-' in toval:
                itosh = '-'.join(map(str,np.array(map(int,toval.split('-')))-tsource))
            else:
                itosh = str(int(toval)-tsource)
        else:
            itosh = str(int(toval)-tsource)
        thestring = thestring.replace('to'+toval,'to'+itosh)
    thestring = thestring.replace('dt',r'\Delta t')        
    thestring = thestring.replace('to','SPACEt_{0}')    
    return thestring

def ProperAll(thestring):
    return ProperTsink(ProperSmear(ProperCM(thestring)))

def FixTflow(thisstr):
    return thisstr.replace('t_flow','SPACEt_{flow}')

def LegLab(string,thisNoSm=False,thisNoTSink=False):
    thisstr = string
    if len(string) == 0: return ''
    thisstr = FixTflow(thisstr)
    if thisNoSm:
        thisstr = NoSm(thisstr)
    # if thisNoTSink:
    ## only single sourcesink sep here
    thisstr = NoTSink(thisstr)    
    return r'$'+ProperAll(thisstr).replace('SPACE','\ ')+'$'


def cutTOfitr(thiscutstr,thistsink):
    intcut = unxmlcut(thiscutstr)
    intcut[1] = int(thistsink) - intcut[1]
    return xmlfitr(intcut)

def LegLabFF(string,thisNoSm=False,thisNoTSink=False):
    if len(string) == 0: return ''
    thisstr = string.replace('state1','')
    thisstr = FixTflow(thisstr)
    thisstr = thisstr.replace('OSFCM','SPACE1SF')
    thisstr = thisstr.replace('TSFTsink','SPACE2SF')
    thisstr = thisstr.replace('CM','SPACEVar')
    thisstr = thisstr.replace('Fits','SPACEFitSPACE')
    if thisNoSm:
        thisstr = NoSm(thisstr)
    ## if thisNoTSink:
    ## only single sourcesink sep here
    thisstr = NoTSink(thisstr)    
    cutstr = re.search('cut.-.',thisstr)
    try:
        cutstr = cutstr.group()
        ## WARNING, cut is forced to tsink 13, CHANGE LATER
        fitstr = cutTOfitr(cutstr,'13')
        thisstr = thisstr.replace(cutstr,fitstr)
    except:
        if 'Proton' in thisstr:
            return r'$Proton\ '+ProperAll(thisstr.replace('Proton','')).replace('SPACE','\ ')+'$'
        elif 'Neutron' in thisstr:
            return r'$Neutron\ '+ProperAll(thisstr.replace('Neutron','')).replace('SPACE','\ ')+'$'
        else:
            return r'$'+ProperAll(thisstr).replace('SPACE','\ ')+'$'
    if 'Proton' in thisstr:
        return r'$Proton\ '+ProperAll(thisstr.replace('Proton','')).replace('SPACE','\ ')+'$'
    elif 'Neutron' in thisstr:
        return r'$Neutron\ '+ProperAll(thisstr.replace('Neutron','')).replace('SPACE','\ ')+'$'
    else:
        return r'$'+ProperAll(thisstr).replace('SPACE','\ ')+'$'

def TitleFix(string):
    return (r'$'+string.replace('P4giDi','\\langle x \\rangle ')
            .replace('giDi','\\langle x \\rangle ')
            .replace('P4I','g_{S}')
            .replace('Ge','G_{E}')
            .replace('Gm','G_{M}')
            .replace('F1divF2','G_{E}/G_{M}')
            # .replace('I','g_{S}')
            .replace('P3g3g5','g_{A}')
            .replace('g3g5','g_{A}')
            .replace('P3',' \\Gamma_{3} ')
            .replace('P4',' \\Gamma_{4} ')
            .replace('g1','\\gamma_{1}')
            .replace('g2','\\gamma_{2}')
            .replace('g3','\\gamma_{3}')
            .replace('g4','\\gamma_{4}')
            .replace('g5','\\gamma_{5}')
            .replace('gA','g_{A}')
            .replace('  ',' ')
            .replace(' ','\\ ')+GetMpi(kappa)+r'$')

def TitleFixFF(string,FF):
    mpi =  GetMpi(kappa)
    if 'PandN' in string:
        string = string.replace('PandN', 'Proton and Neutron' )
    for iDS in DefDSList + CombList:
        string = string.replace(iDS,iDS+' ')
    string = string.replace('IsoVector','Iso-vector')
    if 'PsVector' in string:
        if 'FF1' in FF and 'PsVector' in string:
            string = string.replace('PsVector','Axial')
            return string+ r' $G_{A}\ '+GetMpi(kappa) + '$'
        elif 'FF2' in FF:    
            string = string.replace('PsVector','Induced Pseudoscalar')
            return string+ ' $G_{P}\ '+GetMpi(kappa) + '$'
    string = string.replace('PsScalar','Pseudo-scalar')
    string = string.replace('IsoScalar','Iso-scalar')
    if 'F1divF2' in string:
        return string.replace('F1divF2','').replace(' GeGm','')+ ' $G_{E}/G_{M}\ '+GetMpi(kappa) + '$'        
    elif 'GeGm' in string:
        string = string.replace(' GeGm','')
        FF = FF.replace('FF1','G_{E}')
        FF = FF.replace('FF2','G_{M}')
        return string+ ' $'+FF +'\ '+GetMpi(kappa) + '$'
    elif 'Pseudo-scalar' in string:
        FF = FF.replace('FF1','G_{P}')
        return string+ ' $'+FF +'\ '+GetMpi(kappa) + '$'
    elif 'Scalar' in string:
        FF = FF.replace('FF1','G_{S}')
        return string+ ' $'+FF +'\ '+GetMpi(kappa) + '$'
    elif 'Tensor' in string:
        FF = FF.replace('FF1','H_{T}')
        FF = FF.replace('FF2','E_{T}')
        FF = FF.replace('FF3','E_{1T}')
        return string+ ' $'+FF +'\ '+GetMpi(kappa) + '$'
    elif 'Vector' in string:
        string = string.replace(' VectorTop','')
        string = string.replace(' Vector Top','')
        string = string.replace(' VectorWein','')
        string = string.replace(' Vector Wein','')
        string = string.replace(' Vector','')
        return string+ ' $'+FF.replace('FF','F_{')+'} \ '+GetMpi(kappa) + '$'
    else:        
        return string+ ' $'+FF.replace('FF','F_{')+'} \ '+GetMpi(kappa) + '$'
        


def LabToXaxis(thestring,col):
    stringout = thestring
    if col in 'FitsSm':
        stringout = NoTSink(re.sub('cut.-.','',stringout))
    elif col in 'FitsTsink':
        stringout = NoCM(NoSm(re.sub('cut.-.','',stringout)))
    elif col in 'FitsVar':        
        stringout = re.sub('dt.','',re.sub('cut.-.','',stringout))
    elif col in 'SumMeth':
        stringout = NoCM(NoSm(stringout))
        stringout = stringout.replace('fitr0-4','')
        stringout = stringout.replace('fitr1-4','')
        stringout = stringout.replace('fitr2-4','')
    elif col in 'OSFCM':
        stringout = NoTSink(re.sub('cut.-.','',stringout))
        stringout = re.sub('dt.','',stringout)
    elif col in 'OSFTsink':
        stringout = NoCM(NoSm(re.sub('cut.-.','',stringout)))
    elif col in 'TSFCM':
        stringout = NoTSink(NoCM(NoSm(stringout)))
    elif col in 'TSFTsink':
        stringout = NoCM(NoSm(stringout))
    elif col in 'TSFtest32':
        stringout = NoCM(NoSm(stringout))
    elif col in 'TSFSmall':
        stringout = NoCM(NoSm(stringout))
    stringout = ProperAll(stringout)
    stringout = stringout.replace('cut',r'SPACE\delta t=')
    return r'$'+stringout.replace('SPACE','\ ')+'$'

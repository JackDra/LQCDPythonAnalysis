#!/usr/bin/env python

from sets import Set
import numpy as np
import os
# from MiscFuns import *
from Params import *
from CombParams import CombList
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
    for i in AllTSinkStrList:
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
    thestring = thestring.replace('CM','CMSPACE')
    thestring = thestring.replace('REvec','CMSPACE')
    thestring = thestring.replace('PoF','PoFSPACE')
    toval = re.search('to.',thestring)
    try:
        toval = toval.group().replace('to','').replace('dt','')
    except:
        return thestring
    if 'to'+toval in thestring:
        itosh = str(int(toval)-tsource)
        thestring = thestring.replace('to'+toval,'to'+itosh)
    thestring = thestring.replace('dt',r'\Delta t')        
    thestring = thestring.replace('to','SPACEt_{0}')    
    return thestring

def ProperAll(thestring):
    return ProperTsink(ProperSmear(ProperCM(thestring)))
        
def LegLab(string,NoSm=False,NoTSink=False):
    thisstr = string
    if len(string) == 0: return ''
    if NoSm:
        thisstr = NoSm(thisstr)
    if NoTSink:
        thisstr = NoTSink(thisstr)    
    return r'$'+ProperAll(thisstr).replace('SPACE','\ ')+'$'


def LegLabFF(string,NoSm=False,NoTSink=False):
    if len(string) == 0: return ''
    thisstr = re.sub('cut.','',string)
    thisstr = thisstr.replace('OSFCM','SPACE1SF')
    thisstr = thisstr.replace('TSFTsink','SPACE2SF')
    thisstr = thisstr.replace('CM','SPACEVar')
    thisstr = thisstr.replace('Fits','SPACEFits')
    if NoSm:
        thisstr = NoSm(thisstr)
    if NoTSink:
        thisstr = NoTSink(thisstr)    
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
            .replace(' ','\\ ')+'$')

def TitleFixFF(string,FF):
    for iDS in DefDSList + CombList:
        string = string.replace(iDS,iDS+' ')
    string = string.replace('IsoVector','Iso-vector')
    if 'PsVector' in string:
        if 'FF1' in FF and 'PsVector' in string:
            string = string.replace('PsVector','Axial')
            return string+ ' $G_{A}$'
        elif 'FF2' in FF:    
            string = string.replace('PsVector','Induced Pseudoscalar')
            return string+ ' $G_{P}$'
    string = string.replace('PsScalar','Pseudo-scalar')
    string = string.replace('IsoScalar','Iso-scalar')
    if 'F1divF2' in string:
        return string.replace('F1divF2','').replace(' GeGm','')+ ' $G_{E}/G_{M}$'        
    elif 'GeGm' in string:
        string = string.replace(' GeGm',' Vector')
        FF = FF.replace('FF1','G_{E}')
        FF = FF.replace('FF2','G_{M}')
        return string+ ' $'+FF +'$'
    elif 'Pseudo-scalar' in string:
        FF = FF.replace('FF1','G_{P}')
        return string+ ' $'+FF +'$'
    elif 'Scalar' in string:
        FF = FF.replace('FF1','G_{S}')
        return string+ ' $'+FF +'$'
    elif 'Tensor' in string:
        FF = FF.replace('FF1','H_{T}')
        FF = FF.replace('FF2','E_{T}')
        FF = FF.replace('FF3','E_{1T}')
        return string+ ' $'+FF +'$'
    else:        
        return string+ ' $'+FF.replace('FF','F_{')+'} $'
        


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

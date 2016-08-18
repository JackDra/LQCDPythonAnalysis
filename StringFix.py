#!/usr/bin/env python

from sets import Set
import numpy as np
import os
# from MiscFuns import *
from Params import *
import re

def MakeValAndErr(Avg,Err):
    try:
        int(np.floor(np.log10(Err)))
    except:
        return 'Err'
    ErrMag = int(np.floor(np.log10(Err)))
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
        thestring = thestring.replace(itvar,'REvec')
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
    return thestring.replace('sm','SPACEsm')

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
    thestring = thestring.replace('REvec','REvecSPACE')
    thestring = thestring.replace('PoF','PoFSPACE')
    toval = re.search('to.*dt',thestring)
    try:
        toval = toval.group().replace('to','').replace('dt','')
    except:
        return thestring
    if 'to'+toval in thestring:
        itosh = str(int(toval)-tsource)
        thestring = thestring.replace('to'+toval+'dt','to'+itosh+'dt')
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

def TitleFix(string):
    return (r'$'+string.replace('P4giDi','\\langle x \\rangle ')
            .replace('giDi','\\langle x \\rangle ')
            .replace('P4I','g_{S}')
            .replace('Ge','G_{e}')
            .replace('Gm','G_{m}')
            .replace('F1divF2','G_{e}/G_{m}')
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


def LabToXaxis(thestring,col):
    stringout = thestring
    if col in 'FitsSm':
        stringout = NoTSink(re.sub('cut.','',stringout))
    elif col in 'FitsTsink':
        stringout = NoCM(NoSm(re.sub('cut.','',stringout)))
    elif col in 'FitsVar':        
        stringout = re.sub('cut.','',stringout)
    elif col in 'SumMeth':
        stringout = NoCM(NoSm(stringout))
        stringout = stringout.replace('fitr0-4','')
        stringout = stringout.replace('fitr1-4','')
        stringout = stringout.replace('fitr2-4','')
    elif col in 'OSFCM':
        stringout = NoCM(NoSm(re.sub('cut.','',stringout)))
    elif col in 'OSFTsink':
        stringout = NoCM(NoSm(re.sub('cut.','',stringout)))
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

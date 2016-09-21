#!/usr/bin/env python

from FitFunctions import *
from multiprocessing import Pool
from MultiWrap import *
import numpy as np
from BootTest import BootStrap1
from scipy.optimize import leastsq
from Params import *
from FitParams import *
from MiscFuns import *

##Fitting Routines##


def LSCreate(Fun):
    def LSFun(par,val):
        xval = val[:-2]
        # if len(xval) == 1:
        #     xval = xval[0]
        yval = val[-2]
        errval = val[-1]
        # print xval
        # print par
        # print Fun(xval,par)
        # print yval
        # print 
        return (np.array(Fun(xval,par))-yval)/errval
    return LSFun

def LSDerCreate(FunDer):
    def LSDerFun(par,val):
        xval = val[:-2]
        # if len(xval) == 1:
        #     xval = xval[0]
        yval = val[-2]
        errval = val[-1]
        print xval, par
        print FunDer(xval,par)
        print errval
        print 
        return FunDer(xval,par)/errval
        # return np.transpose(FunDer(xval,par)/errval)
    return LSDerFun

def DerOfFun(Fun,Len=1):
    if Fun.__name__ == 'ConstantFitFun':
        return ConstFFDer
    elif Fun.__name__ == 'LinearFitFun':
        return LinFFDer
    elif Fun.__name__ == 'C2OneStateFitFun':
        return C2OSFFDer
    elif Fun.__name__ == 'C3OneStateFitFun':
        return C3OSFFDer
    elif Fun.__name__ == 'C2TwoStateFitFun':
        return C2TSFFDer
    elif Fun.__name__ == 'C2TwoStateFitFunCM':
        return C2TSFFCMDer
    elif Fun.__name__ == 'C3MomTwoStateFitFun':
        return C3MomTSFFDer
    elif Fun.__name__ == 'TestTwoVarFitFun':
        return TestTwoVarFFDer
    elif Fun.__name__ == 'FormFactorO1':
        return FormFactorO1Der
    elif Fun.__name__ == 'FormFactorO':
        return CreateFFO(Len)[1]
    elif Fun.__name__ == 'FormFactorO2':
        return FormFactorO2Der
    elif Fun.__name__ == 'FormFactorO3':
        return FormFactorO3Der
    elif Fun.__name__ == 'DPfitfun':
        return DPfitfunDer
    elif Fun.__name__ == 'DPfitfutn2':
        return DPfitfun2Der

def GetLSFuns(fitfun,derfun,iGuess,parlen):
   if iGuess == None:
       iGuess = FitDefGuess(fitfun,Len=parlen)
   LSfitfun = LSCreate(fitfun)
   if derfun == None:
       LSDerfitfun = LSDerCreate(DerOfFun(fitfun,parlen))
   else:
       LSDerfitfun = LSDerCreate(derfun)
   # if 'C3TwoStateFitFun' not in fitfun.__name__ :
   #     LSDerfitfun = LSDerCreate(DerOfFun(fitfun))
   # else:
   #     if derfun == None:
   #         raise AttributeError('Must Pass Derivative as derfun= for C3FitFun')
   #     else:
   #         LSDerfitfun = LSDerCreate(derfun)
   return LSfitfun,LSDerfitfun,iGuess

def CreateArgs(xdata,ydata,yerr):
    data = []
    if isinstance(xdata[0],list) or isinstance(xdata[0], np.ndarray):
        for ix in xdata:
            data.append(np.array(ix))
    else:
        data.append(xdata)
    data.append(np.array(ydata))
    data.append(np.array(yerr))
    return data


def LSFit(parlen,xdata,yerr,fitfun,ydata):
    iGuess = None
    MI = MaxIters
    derfun = None
    data = CreateArgs(xdata,ydata,yerr)
    LSfitfun,LSDerfitfun,iGuess = GetLSFuns(fitfun,derfun,iGuess,parlen)
    # if isinstance(ydata[0],complex):iGuess = map(complex,iGuess)
    # if Debug:
    #     print LSDerfitfun.__name__
    #     print LSfitfun.__name__
    #     print iGuess
    #     print data
    if ForceNoDer:
        x,covar, infodict, mesg, ier=leastsq(LSfitfun,iGuess,args=data, maxfev=MI, full_output=1)
    else:
        x,covar, infodict, mesg, ier=leastsq(LSfitfun,iGuess,args=data, Dfun=LSDerfitfun, maxfev=MI, full_output=1)
    if float(len(ydata)-len(x)) == 0:
        chisqpdf = float('NaN')
    else:
        chisqpdf=sum(infodict["fvec"]*infodict["fvec"])/float(len(ydata)-len(x))
    # if ier != 1:
    #    print x
    #    print "WARNING: Optimal parameters not found: " + mesg
    #      raise ValueError, "Optimal parameters not found: " + mesg
    #   print x,covar
    return x,covar,chisqpdf

def FitBoots(ydatain,xdatain,FitFun,DoW='T',MI=MaxIters,parlen=1,tBooted=False,thisnboot=nboot):
    GetBootStats(ydatain)
    # print ydatain
    # print Pullflag(ydatain,'Avg')
    # for iyd in ydatain:
    #     print len(iyd.values)
    # print Pullflag(ydatain,'values')
    ydataAvg = Pullflag(ydatain,'Avg')
    ydatavals = np.rollaxis(Pullflag(ydatain,'values'),1)
    fitdata = []
    if DoW:
        ydataStd = Pullflag(ydatain,'Std')
    else:
        ydataStd = [1]*len(ydataAvg)
    if tBooted:
        [fitdataAvg,fitdataAvgErr,fitdataChi] = LSFit(parlen,xdatain[0],ydataStd,FitFun,ydataAvg)
    else:
        [fitdataAvg,fitdataAvgErr,fitdataChi] = LSFit(parlen,xdatain,ydataStd,FitFun,ydataAvg)

    if MultiCoreFitting:
        makeContextFunctions(LSFit)
        FitPool = Pool(processes=AnaProc)
        if tBooted:
            inputdata = [(parlen,ix,ydataStd,FitFun,iy) for ix,iy in zip(xdatain[1:],ydatavals)]
        else:
            inputdata = [(parlen,xdatain,ydataStd,FitFun,iy) for iy in ydatavals]
        fitdatavals = FitPool.map(LSFit.mapper,inputdata)
        FitPool.close()
        for iy in range(len(fitdatavals[0][0])):
            fitdata.append(BootStrap1(thisnboot,0))
            for iboot in range(thisnboot):
                fitdata[iy].values[iboot] = fitdatavals[iboot][0][iy]
    else:
        fitdatavals = []
        for iboot,bootdata in enumerate(ydatavals):
            if tBooted:
                tempboot = LSFit(parlen,xdatain[iboot+1],ydataStd,FitFun,bootdata)
            else:
                tempboot = LSFit(parlen,xdatain,ydataStd,FitFun,bootdata)
            fitdatavals.append([])
            for iy,iyd in enumerate(tempboot[0]):
                fitdatavals[iboot].append(iyd)
        for iy in range(len(fitdatavals[0])):
            fitdata.append(BootStrap1(thisnboot,0))
            for iboot in range(thisnboot):
                fitdata[iy].values[iboot] = fitdatavals[iboot][iy]
    GetBootStats(fitdata)
    # fitdataChi = CalcChiSqrdPDF(FitFun,fitdataAvg,xdatain,ydataAvg,ydataStd)
    return fitdata,fitdataAvg,[fitdataChi]*len(fitdata)





# def FitVarFunBoots(ydatain,xdatain,FunAvg,FunBootList,DoW='T',MI=MaxIters,parlen=1):
#     # C3FitFun,C3FitFunDer = FunGen(Pullflag(BootPars,'Avg'),tsvar)
#     # BPlist = np.rollaxis(Pullflag(BootPars,'values'),1)
#     GetBootStats(ydatain)
#     ydataAvg = Pullflag(ydatain,'Avg')
#     ydatavals = np.rollaxis(Pullflag(ydatain,'values'),1)
#     fitdata = []
#     if DoW:
#         ydataStd = Pullflag(ydatain,'Std')
#     else:
#         ydataStd = [1]*len(ydataAvg)
#     [fitdataAvg,fitdataAvgErr,fitdataChi] = LSFit(parlen,xdatain,ydataStd,[FunAvg,ydataAvg])
#     if DoMultiCore:
#         # thisFunWrap = MakeWrap(LSFit,(parlen,xdatain,ydataStd))
#         # makeContextFunctionsNo2(thisFunWrap)
#         FitPool = Pool(processes=AnaProc)
#         # inputargs = [[[FunBootList[iboot][0].mapper,FunBootList[iboot][1].mapper],ydatavals[iboot]] for iboot in range(nboot)]
#         inputargs = [[parlen,xdatain,ydataStd,FunBootList[iboot],ydatavals[iboot]] for iboot in range(nboot)]
#         # inputargs = [[FunAvg,ydatavals[iboot]] for iboot in range(nboot)]
#         # fitdatavals = FitPool.map(thisFunWrap.mapper,inputargs)
#         makeContextFunctions(LSFit)
#         fitdatavals = FitPool.map(LSFit.mapper,inputargs)
#         FitPool.close()
#         for iy in range(len(fitdatavals[0][0])):
#             fitdata.append(BootStrap1(nboot,0))
#             for iboot in range(nboot):
#                 fitdata[iy].values[iboot] = fitdatavals[iboot][0][iy]
#     else:
#         fitdatavals = []
#         for iboot,bootdata in enumerate(ydatavals):
#             tempboot = LSFit(parlen,xdatain,ydataStd,[FunBootList[iboot],bootdata])
#             fitdatavals.append([])
#             for iy,iyd in enumerate(tempboot[0]):
#                 fitdatavals[iboot].append(iyd)
#         for iy in range(len(fitdatavals[0])):
#             fitdata.append(BootStrap1(nboot,0))
#             for iboot in range(nboot):
#                 fitdata[iy].values[iboot] = fitdatavals[iboot][iy]
#     GetBootStats(fitdata)
#     # fitdataChi = CalcChiSqrdPDF(C3FitFun,fitdataAvg,xdatain,ydataAvg,ydataStd)
#     return fitdata,fitdataAvg,[fitdataChi]*len(fitdata)





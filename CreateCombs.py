#!/usr/bin/env python
from array import array
import numpy as np
from Params import *
from MiscFuns import *
from collections import OrderedDict
from copy import deepcopy
from ReadXml import ReadXmlAndPickle
from OppFuns import CreateOppDir
from OutputXmlData import MergeXmlOutput
import operator

giDiVecSet = ['P4g1D1','P4g2D2','P4g3D3']
##Proton: doublet is up quark, singlet is down quark
##Neutron: doublet is down quark, singlet is up quark
upCharge = 2.0/3.0
downCharge = -1.0/3.0
DSCombs = [('P4I',upCharge,downCharge),('P4g4',upCharge,downCharge),('P4g3g5',-1,1)]
ops = { "+": operator.add, "-": operator.sub } 

def IsoVector(val1,val2):
    return val1 - val2

def FFProton(val1,val2):
    return upCharge*val1 - downCharge*val2

def FFNeutron(val1,val2):
    return downCharge*val1 - upCharge*val2

# ## data = [ itsink , ism/istate , igamma , ip , it ]
# def MakeUmD(data,gammalist):
#     dataout = deepcopy(data)
#     newgammalist,pairindex = UDIndex(gammalist)
#     for tsink,tsinkdata in enumerate(data):
#         for state,statedata in enumerate(tsinkdata):
#             for idu,index in enumerate(pairindex):
#                 dataout[tsink][state].append([])
#                 dgamma,ugamma = statedata[index[0]],statedata[index[1]]
#                 for ip,(pd,pu) in enumerate(zip(dgamma,ugamma)):
#                     dataout[tsink][state][len(data[0][0])+idu].append([])
#                     for it,(td,tu) in enumerate(zip(pd,pu)):
#                         dataout[tsink][state][len(data[0][0])+idu][ip].append(td-tu)
#                         dataout[tsink][state][len(data[0][0])+idu][ip][it].Stats()
#     return dataout,newgammalist

def CreateDS(datadoub,datasing,thisGammaList):
    dataout = Swap3ptSS(datadoub)
    datahold = Swap3ptSS(datasing)
    for igamma, Lcoef,Rcoef in DSCombs:
        if igamma in thisGammaList:
            igloc = thisGammaList.index(igamma)
            dataout[igloc] = Lcoef*dataout[igloc]+Rcoef*datahold[igloc]
    return SwapBack3ptSS(dataout)

def CreategiDi(data3pt,thisGammaList,thisDSList):
    data3ptout = Swap3ptSS(data3pt)
    if len(data3ptout) != 8:
        print 'dimensions of correlator (ism, jsm, igamma)'
        for ic,idata in enumerate(data3pt):
            for jc,jdata in enumerate(idata):
                print ic , jc , len(jdata)
        raise IOError("Need giDi with i =1,2,3,4 for giDi combination")
    thisGammaListOut = np.array(thisGammaList)
    for iDS in thisDSList:
        g4D4i = thisGammaListOut.tolist().index(iDS+'P4g4D4')
        thisGammaListOut = np.where(thisGammaListOut==iDS+'P4g4D4',iDS+'P4giDi',thisGammaListOut)
        giDii = []
        for ivec,igamma in enumerate(giDiVecSet):
            giDii.append(thisGammaListOut.tolist().index(iDS+igamma))
            data3ptout[g4D4i] -= data3ptout[giDii[ivec]]/3.0
            
            data3ptout = np.delete(data3ptout,giDii,axis=0)
            thisGammaListOut = np.delete(thisGammaListOut,giDii)
    return [SwapBack3ptSS(data3ptout),thisGammaListOut.tolist()]
    

#cfunin [  ism , jsm , igamma , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
#cfunout [ igamma  , ism , jsm , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)

def Swap3ptSS(cfunin):
    return np.rollaxis(np.array(cfunin),2)


#cfunin [ igamma  , ism , jsm , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
#cfunout [  ism , jsm , igamma , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)

def SwapBack3ptSS(cfunin):
    return np.rollaxis(np.array(cfunin),0,3)
        

def FunctOfDicts(a, b,Funct):
    for key in b:
        if key in a:
            if isinstance(a[key], dict) and isinstance(b[key], dict):
                FunctOfDicts(a[key], b[key],Funct)
            elif hasattr(a[key],"__len__") and hasattr(b[key],"__len__"):
                if len(a[key]) == nboot and len(b[key]) == nboot:
                    a[key] = [Funct(ia,ib) for ia,ib in zip(a[key],b[key])]
            elif key == 'Chi':
                a[key] = a[key] + b[key]
            else:
                pass
        else:
            raise IOError('Dictionaries not equal in keys')
    return a

def XmlBootToAvg(datadict,BootDict=None):
    if BootDict == None: BootDict = datadict
    for key in datadict:
        if key == 'Boot': continue
        if isinstance(datadict[key], dict):
            if key == 'Values':
                XmlBootToAvg(datadict[key],BootDict = BootDict['Boots'])
            else:
                XmlBootToAvg(datadict[key],BootDict = BootDict[key])          
        elif key == 'Avg':
            datadict[key] = np.mean(BootDict)
        elif key == 'Std':
            datadict[key] = np.std(BootDict)
        else:
            pass
    return datadict
    

def CombTwoFiles(file1,file2,funct):
    data1,dump = ReadXmlAndPickle(file1)
    data2,dump = ReadXmlAndPickle(file2)
    return XmlBootToAvg(FunctOfDicts(data1,data2,funct))

def ReadAndComb(inputargs,Funct,funname):
    for igamma in inputargs['gamma']:
        if 'doub' in igamma or 'sing' in igamma: continue
        doubgamma = 'doub'+igamma
        singgamma = 'sing'+igamma
        doubgammadir = CreateOppDir(doubgamma)
        singgammadir = CreateOppDir(singgamma)
        gammadir = CreateOppDir(igamma)
        for imethod in inputargs['method']:
            if 'TSF' in imethod:
                preflist = TwoStateParList['C3']
            if 'OSF' in imethod:
                preflist = OneStateParList['C3']
            else:
                preflist = ['']
            for ipref in preflist:
                for imom in inputargs['mom']:
                    momdir = MakeMomDir(imom)
                    for iset in inputargs['set']:
                        filedoub = outputdir +'/'+ doubgammadir + '/' + imethod + '/'+momdir + '/' + iset+doubgamma+ipref+imom+'.xml'
                        filesing = outputdir +'/'+ singgammadir + '/' + imethod + '/'+momdir + '/' + iset+singgamma+ipref+imom+'.xml'
                        outdata = CombTwoFiles(filedoub,filesing,Funct)
                        mkdir_p( outputdir +'/'+ gammadir + '/'+funname+'/' + imethod + '/'+momdir + '/')
                        outfile = outputdir +'/'+ gammadir + '/'+funname+'/' + imethod + '/'+momdir + '/' + iset+funname+igamma+ipref+imom+'.xml'
                        MergeXmlOutput(outfile,outdata)

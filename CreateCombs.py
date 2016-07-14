#!/usr/bin/env python
from array import array
import numpy as np
from Params import *
from MiscFuns import *
from collections import OrderedDict
from copy import deepcopy
from ReadXml import ReadXmlAndPickle
from OppFuns import CreateOppDir
from OutputXmlData import MergeXmlOutput,WriteXmlOutput
import operator
from MomParams import *
from XmlFormatting import *
from FitParams import *
from SetLists import *
from CombParams import *
from ReadDir import CheckCurrentSets



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
                if 'Avg' not in datadict[key].keys():
                    XmlBootToAvg(datadict[key],BootDict = BootDict[key])          
                elif key in BootDict.keys():
                    datadict[key]['Avg'] = np.mean(BootDict[key])
                    datadict[key]['Std'] = np.std(BootDict[key])
                    if 'Chi' in datadict[key].keys():
                        datadict[key] = DictAvgStdChiToFormat(datadict[key],datadict[key]['Chi'])
                    else:
                        datadict[key] = DictAvgStdToFormat(datadict[key])
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
        if 'doub' in igamma or 'sing' in igamma or 'twopt' in igamma: continue
        doubgamma = 'doub'+igamma
        singgamma = 'sing'+igamma
        doubgammadir = CreateOppDir(doubgamma)
        singgammadir = CreateOppDir(singgamma)
        gammadir = CreateOppDir(funname+igamma)
        for imethod in inputargs['method']:
            if imethod == 'RF': methoddir = ''
            else: methoddir = '/'+imethod
            thisSetList = inputargs['set']
            if 'TSF' in imethod or 'SumMeth' in imethod:
                thisSetList = ReduceTsink(inputargs['set'])
            if 'TSF' in imethod:
                preflist = TwoStateParList['C3']
            elif 'OSF' in imethod:
                preflist = OneStateParList['C3']
            else:
                preflist = ['']                
            for ipref in preflist:
                for imom in inputargs['mom']:
                    momstr = qstrTOqcond(imom)
                    momdir = MakeMomDir(imom)
                    for iset in thisSetList:
                        filedoub = outputdir +'/'+ doubgammadir + methoddir + '/'+momdir + '/' + iset+doubgamma+ipref+momstr+'.xml'
                        filesing = outputdir +'/'+ singgammadir + methoddir + '/'+momdir + '/' + iset+singgamma+ipref+momstr+'.xml'
                        if Debug: print filedoub
                        if Debug: print filesing                        
                        outdata = CombTwoFiles(filedoub,filesing,Funct)
                        mkdir_p( outputdir +'/'+ gammadir +methoddir + '/'+momdir + '/')
                        outfile = outputdir +'/'+ gammadir +methoddir + '/'+momdir + '/' + iset+funname+igamma+ipref+momstr
                        if Debug: print outfile
                        # WriteXmlOutput(outfile,outdata)
                        MergeXmlOutput(outfile,outdata)


                        
def ReadAndCombFF(thisCurrDict,Funct,funname):
    newCurrDict = CheckCurrentSets(thisCurrDict)
    for icurr,isetlist in newCurrDict.iteritems():
        doubcurr = 'doub'+icurr
        singcurr = 'sing'+icurr
        for iset in isetlist:
            filedoub = outputdir+'FormFactors/'+doubcurr+'/' +doubcurr+iset+'.xml'
            filesing = outputdir+'FormFactors/'+singcurr+'/' +singcurr+iset+'.xml'
            if Debug: print filedoub
            if Debug: print filesing                        
            outdata = CombTwoFiles(filedoub,filesing,Funct)
            mkdir_p( outputdir+'FormFactors/'+funname+icurr+'/')
            outfile = outputdir+'FormFactors/'+funname+icurr+'/'+ funname+icurr+iset
            MergeXmlOutput(outfile,outdata)

                        

def FunctOfDictsOld(a, b,Funct):
    for key in b:
        if key in a:
            if isinstance(a[key], dict) and isinstance(b[key], dict):
                FunctOfDictsOld(a[key], b[key],Funct)
            elif hasattr(a[key],"values") and hasattr(b[key],"values"):
                if len(a[key].values) == nboot and len(b[key].values) == nboot:
                    a[key].values = np.array([Funct(ia,ib) for ia,ib in zip(a[key].values,b[key].values)])
                    a[key].Stats()
                else:
                    raise IOError('nboot missmatch, file1: ', len(a[key].values), 'file2: ', len(b[key].values), ' params: ',nboot)
            elif hasattr(a[key],"__len__") and hasattr(b[key],"__len__"):
                for j,(ja,jb) in enumerate(zip(a[key],b[key])):
                    if hasattr(ja,"values") and hasattr(jb,"values"):
                        if len(ja.values) == nboot and len(jb.values) == nboot:
                            a[key][j].values = np.array([Funct(ia,ib) for ia,ib in zip(ja.values,ja.values)])
                            a[key][j].Stats()                    
                        else:
                            raise IOError('nboot missmatch, file1: ', len(ja.values), 'file2: ', len(jb.values),' params: ',nboot)
            elif key == 'Chi':
                a[key] = a[key] + b[key]
            else:
                pass
        else:
            raise IOError('Dictionaries not equal in keys')
    return a

def XmlBootToAvgOld(datadict):
    for key in datadict.keys():
        if isinstance(datadict[key], dict):
            XmlBootToAvgOld(datadict[key])          
        elif key == 'Avg' or key == 'Vals':
            datadict[key] = Pullflag(datadict['Boot'],'Avg')
        elif key == 'Std' or key == 'Valserr':
            datadict[key] = Pullflag(datadict['Boot'],'Std')
        else:
            pass
    return datadict


def CreateDictOldCombs(datadict,thisCombList):
    datadictout = {'doub':{} , 'sing':{}}
    for igamma,gammadict in datadict.iteritems():
        if 'doub' in igamma:
            doubgamma,singgamma,gamma = igamma,igamma.replace('doub','sing'),igamma.replace('doub','')
        else:
            continue
        if singgamma not in datadict.keys(): continue
        datadictout['doub'][gamma] = datadict[doubgamma]
        datadictout['sing'][gamma] = datadict[singgamma]
        for funtype in thisCombList:
            if funtype not in datadictout.keys(): datadictout[funtype] = {}
            datadictout[funtype][gamma] = deepcopy(XmlBootToAvgOld( FunctOfDictsOld(datadict[doubgamma],datadict[singgamma],CombFunsDict[funtype])))
            if Debug:
                print funtype, gamma, datadictout[funtype][gamma]
    return datadictout

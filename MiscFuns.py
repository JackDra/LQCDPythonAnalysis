#!/usr/bin/env python

import numpy as np
from BootTest import BootStrap1
import operator as op
import functools as ft
import os, errno
import sys
from copy import deepcopy
from collections import OrderedDict
# import pylab as pl
import time,datetime
from copy import deepcopy

def SplitCmplxReal(array):
   realout,cmplxout = [],[]
   for ia in array:
      if isinstance(ia, complex):
         cmplxout.append(ia.imag)
         realout.append(ia.real)
      else:
         realout.append(ia)
         cmplxout.append(0.)
   return realout,cmplxout
   
         
def CheckZip(*data):
   for itd,idata in enumerate(data[1:]):
      if len(idata) != len(data[0]):
         raise IOError(str(itd) +' input is not same length as 0')
      
      
def merge_dicts(a, b, path=None):
   if path is None: path = []
   for key in b:
      if key in a:
         if isinstance(a[key], dict) and isinstance(b[key], dict):
            merge_dicts(a[key], b[key], path + [str(key)])
         else:
            pass
         # elif a[key] == b[key]:
         #    pass # same leaf value
         # else:
         #    raise Exception('Conflict at %s' % '.'.join(path + [str(key)]))
      else:
         a[key] = b[key]
   return a
   


def DelDubs(listin):
   listout = []
   for il in listin:
      if il not in listout:
         listout.append(il)
   return listout

def chunks(l, n):
   n = max(1, n)
   return [l[i:i + n] for i in range(0, len(l), n)]

def BootNdimDict(dictin):
   if isinstance(dictin,dict):
      for dictkey,dictval in dictin.iteritems():
         if 'Boot' in dictkey:
            if isinstance(dictval,BootStrap1):
               dictval.Stats()
            else:
               for idict in dictval:
                  idict.Stats()
         else:
            BootNdimDict(dictval)
      
def CheckDict(thisdict,*dictkeys):
   dictcheck = deepcopy(thisdict)
   for idict in dictkeys:
      if idict in dictcheck.keys():
         dictcheck = dictcheck[idict]
      else:
         return False
   return True

def GetPercent(counter,totlen):
   return str(int((counter*100)/float(totlen))) + '% '

def GetTimeStr(thistime):
   return str(datetime.timedelta(seconds=thistime)) + ' h:m:s '

def GetTimeLeft(counter,totlen,timedone):
   perdone = counter/float(totlen)
   if perdone < .01:
      return float(0)
   else:
      return timedone*((1-perdone)/perdone)

def GetTimeLeftStr(counter,totlen,timedone):
   return GetTimeStr(GetTimeLeft(counter,totlen,timedone))

def GetTimeForm(counter,totlen,timedone):
   return 'Current Time: ' + GetTimeStr(timedone) + 'at ' + GetPercent(counter,totlen) + 'Time Left: ' + GetTimeLeftStr(counter,totlen,timedone)
   
class Unbuffered(object):
   def __init__(self, stream):
       self.stream = stream
   def write(self, data):
       self.stream.write(data)
       self.stream.flush()
   def __getattr__(self, attr):
       return getattr(self.stream, attr)
sys.stdout = Unbuffered(sys.stdout)

def ParInput(constargs,pararg):
   return [constargs+(ipar,) for ipar in pararg]

def ParInputTwo(constargs,pararg,pararg2):
   return [constargs+(ipar,ipar2) for ipar,ipar2 in zip(pararg,pararg2)]

def ParInputTwo(constargs,pararg,pararg2,pararg3):
   return [constargs+(ipar,ipar2,ipar3) for ipar,ipar2,ipar3 in zip(pararg,pararg2,pararg3)]

def WipeFile(fname):
   open(fname,'w').close()

def touch(fname):
    try:
        os.utime(fname, None)
    except:
        open(fname, 'a').close()

def removekey(d, *key):
    r = dict(d)
    for ikey in key:
       del r[ikey]
    return r

def VecDelta(tvec,val):
    return map(int,val==np.array(tvec))

#list a is outer loop
def Elongate(lista,listb):
    c = []
    for ia in lista:
        for ib in listb:
            c.append([ia,ib])
    return c

def ElongateName(lista,listb):
    c = []
    for ia in lista:
        for ib in listb:
            c.append(ia+ib)
    return c

def FindWhichTSF(Title):
    if 'TSinkRed' in Title or 'Reduced' in Title:
        return 'test32'
    elif 'TSink' in Title or 'Tsink' in Title:
        return 'Tsink'
    elif 'MassCM' in Title:
        return 'CM'
    elif 'Mass' in Title:
        if 'Comparison' in Title:
            return 'CM'
        else:
            return 'Tsink'
    elif 'CM' in Title or 'Smearing' in Title:
        return 'CM'
    elif 'Small' in Title:
        return 'Small'
    else:
        return 'CM'

def DelConf(data,iconf):
    for delism,p in enumerate(data):
        for deljsm,q in enumerate(p):
            del data[delism][deljsm][iconf]
    return data

# def get_delcmap(thiscmap,points,delcol):
#     colorrange = pl.linspace(0,1,points)
#     midcr = (delcol[0]+delcol[1])/2.0
#     widcr = (delcol[1]-delcol[0])/2.0
#     for i,icr in enumerate(colorrange):
#         if icr < midcr: colorrange[i] = icr - (icr*widcr)/midcr
#         elif icr > midcr: colorrange[i] = (1-(widcr/(1-midcr)))*(icr-midcr) + midcr + widcr
#     return pl.get_cmap(thiscmap)(colorrange)


def DiagSmear(data2pt):
    data2ptout = []
    for ism,ismdata in enumerate(data2pt):
        data2ptout.append(ismdata[ism])
    return np.array(data2ptout)


def DiagSmearWithTsrc(data2pt):
    data2ptout = []
    for its,tsdata in enumerate(data2pt):
       for ism,ismdata in enumerate(tsdata):
          data2ptout.append(ismdata[ism])
    return np.array(data2ptout)

def Diag3ptSmear(data3pt):
    data3ptout = []
    for itsink,tsinkdata in enumerate(data3pt):
        data3ptout.append([])
        for ism,ismdata in enumerate(tsinkdata):
            data3ptout[itsink].append(ismdata[ism])
    return np.array(data3ptout)

def GetBootStats(data):
   flatdata = np.array(data).flatten()
   map(lambda x : x.Stats() ,flatdata)
   return np.reshape(flatdata,np.array(data).shape)

def MassFun(cfun,Dt=1):
    mass = []
    for it,tcfun in enumerate(cfun):
       if it+Dt < len(cfun):
          mass.append(np.abs(np.log(np.abs(cfun[it+Dt]/tcfun)))/Dt)          
    return GetBootStats(mass)

def cfunTOmass(cfun):
   if len(np.array(cfun).shape) > 1:
      out = np.array(NDimOpp(cfun,1,MassFun))
      return np.rollaxis(out,0,len(out.shape))
   else:
      return MassFun(cfun)
def flattenAllBut(data,dimleft):
    data = np.array(data)
    return data.reshape(reduce(op.mul,data.shape[:-dimleft]), *data.shape[-dimleft:])

def MultArray(data):
    return ft.reduce(op.mul, list(data), 1)

##applies funin datavariables over all - dimleft variables with extra variables funvars using funin
def NDimOpp(data,dimleft,funin,*funvars):
    npdata = np.array(data)
    dataout = np.array([])
    if dimleft > 0:flatdata = flattenAllBut(npdata,dimleft)
    else: flatdata = npdata.flatten()
    for idata in flatdata:
        dataapp = funin(idata,*funvars)
        datainshape = np.array(dataapp).shape
        dataout = np.append(dataout,dataapp)
    if dimleft > 0: outdim = npdata.shape[:-dimleft] + datainshape
    else: outdim = npdata.shape
    dataout = np.reshape(dataout,outdim)
    for roll in datainshape:
        dataout = np.rollaxis(dataout,len(dataout.shape)-1)
    return dataout

def NDimOpp2(data,data2,dimleft,funin,*funvars):
    npdata,npdata2 = np.array(data),np.array(data2)
    dataout = np.array([])
    if dimleft > 0:flatdata = zip(flattenAllBut(npdata,dimleft),flattenAllBut(npdata2,dimleft))
    else: flatdata = zip(npdata.flatten(),npdata2.flatten())
    for idata,idata2 in flatdata:
        dataapp = funin(idata,idata2,*funvars)
        datainshape = np.array(dataapp).shape
        dataout = np.append(dataout,dataapp)
    if dimleft > 0: outdim = npdata.shape[:-dimleft] + datainshape
    else: outdim = npdata.shape
    if MultArray(dataout.shape) > MultArray(outdim): extra = (MultArray(dataout.shape)/MultArray(outdim),)
    else : extra = ()
    dataout = np.reshape(dataout,outdim)
    for roll in datainshape:
        dataout = np.rollaxis(dataout,len(dataout.shape)-1)
    return dataout

def NDimOpp3(data,data2,data3,dimleft,funin,*funvars):
    npdata,npdata2,npdata3 = np.array(data),np.array(data2),np.array(data3)
    dataout = np.array([])
    if dimleft > 0:flatdata = zip(flattenAllBut(npdata,dimleft),flattenAllBut(npdata2,dimleft),flattenAllBut(npdata3,dimleft))
    else: flatdata = zip(npdata.flatten(),npdata2.flatten(),npdata3.flatten())
    for idata,idata2,idata3 in flatdata:
        dataapp = funin(idata,idata2,idata3,*funvars)
        datainshape = np.array(dataapp).shape
        dataout = np.append(dataout,dataapp)
    if dimleft > 0: outdim = npdata.shape[:-dimleft] + datainshape
    else: outdim = npdata.shape
    if MultArray(dataout.shape) > MultArray(outdim): extra = (MultArray(dataout.shape)/MultArray(outdim),)
    else : extra = ()
    dataout = np.reshape(dataout,outdim)
    for roll in datainshape:
        dataout = np.rollaxis(dataout,len(dataout.shape)-1)
    return dataout

#Pulls out a class "atribute" out of a numpy tensor "data" of classes ##
#e.g. return data[:,:,:].atribute
def Pullflag(data,atribute):
    npdata = np.array(data)
    dataout = np.array([])
    if len(npdata.flatten()) == 0: return data
    for idata in npdata.flatten():
        flagdim = np.array(getattr(idata,atribute)).shape
        dataout = np.append(dataout,getattr(idata,atribute))
    return np.reshape(dataout,npdata.shape + flagdim)


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def SplitList(thelist):
    return [ilist[0] for ilist in thelist],[itlist[1] for ilist in thelist]
def UDIndex(gammalist):
   pairindex,newgammalist = [],deepcopy(gammalist)
   for dindex,dgamma in enumerate(gammalist):
      if 'doub' in dgamma:
         for uindex,ugamma in enumerate(gammalist):
            if dgamma.replace('doub','sing') == ugamma:
               pairindex.append([dindex,uindex])
               newgammalist.append(dgamma.replace('doub',''))
   return newgammalist,pairindex

def DeCorrBoot(thisBoot):
   npdata = np.array(thisBoot)
   dataout = np.array([])
   if len(npdata.flatten()) == 0: return thisBoot
   for ic,idata in enumerate(npdata.flatten()):
      vals = idata.values
      flagdim = vals.shape
      dataout = np.append(dataout,BootStrap1(len(vals),5))
      dataout[-1].values = np.array(list(reversed(vals)))
      dataout[-1].Stats()
   return np.reshape(dataout,npdata.shape)
   
 
## USED FOR OLD OPPERATOR DIRECTORIES##
# if kappa == 12090:
#     def CreateOppDir(Opp):
#         return Opp
# else:

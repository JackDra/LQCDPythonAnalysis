#!/usr/bin/env python

import numpy as np
import os, errno


MomSqrdSet = [ '0','1','2','3','4','5','6','7','8','9' ]
QsqrdSet = ['qsqrd'+iq for iq in MomSqrdSet]
Maxqsqrd = np.max(np.array(MomSqrdSet).astype(int))
nxyz = 32
qunit = (2.0*np.pi)/float(nxyz)



def makeqlist():
    qlist = np.array([])
    for iq1 in range(-Maxqsqrd,Maxqsqrd+1):
        for iq2 in range(-Maxqsqrd,Maxqsqrd+1):
            for iq3 in range(-Maxqsqrd,Maxqsqrd+1):
                if iq1**2 + iq2**2 + iq3**2 > Maxqsqrd: continue
                qlist = np.append(qlist,'q = ' + str(iq1) + ' ' + str(iq2) + ' ' + str(iq3))
    return qlist

qvecSet = makeqlist()
nmom = len(qvecSet)
qhigh = nmom/2

def OrderMomList(momset):
    qout = []
    for iq in qvecSet:
        if iq in momset:
            qout.append(iq)
    return qout

## momlist in p indexing
def DragpZ(momlist):
    if iqTOip(0) in momlist:
        momlist.insert(0, momlist.pop(momlist.index(iqTOip(0))))
    return momlist

## momlist in str indexing
def DragpZstr(momlist):
    if 'q = 0 0 0' in momlist:
        momlist.insert(0, momlist.pop(momlist.index('q = 0 0 0')))
    return momlist

qvecSetZfirst = DragpZstr(qvecSet.tolist())

# ZeroMomIndex = int(open('./parfiles/ZeroMom.par','r').read())

### NOTE: ALWAYS MAKE q = 0 0 0 the FIRST MOMENTA!!! ###
ZeroMomIndex = 0


def PrintZMomI(val):
    open('./parfiles/ZeroMom.par','w').write(str(val))

# def PrintAllMom():
#     f = open('./parfiles/AllMom.par','w')
#     for line in qvecSet:
#         f.write(line+'\n')

#ip from 0 -> nmom, iq from -qhigh -> qhigh
def ipTOiq(ip):
    return ip - qhigh
def iqTOip(iq):
    return iq + qhigh


def ipTOqstr(ip):
    return qvecSet[int(ip)]

def ipTOqcond(ip):
    return qvecSet[int(ip)].replace('=','').replace(' ','')

def qcondTOqstr(iq):
    return ' '.join(iq).replace('q ','q = ')

def iqTOqstr(iq):
    return qvecSet[iqTOip(iq)]

def ipZshiftTOqstr(ip):
    return qvecSetZfirst[ip]
    
def qvecTOqstr(qvec):
    return qvecSet[qvecTOip(qvec)]

def qstrTOqvec(qstr):
    sqvec = qstr.split()
    return [int(sqvec[2]),int(sqvec[3]),int(sqvec[4])]


def qstrTOip(qstr):
    return qvecTOip(qstrTOqvec(qstr))

def qstrTOiq(qstr):
    return qvecTOiq(qstrTOqvec(qstr))

def ipTOqvec(ip):
    sqvec = qvecSet[ip].split()
    return [int(sqvec[2]),int(sqvec[3]),int(sqvec[4])]

def qvecTOip(qvec):
    strqvec = 'q = ' + str(qvec[0]) + ' ' + str(qvec[1]) + ' ' + str(qvec[2])
    return np.where(qvecSet==strqvec)[0][0]

def iqTOqvec(iq):
    return ipTOqvec(iqTOip(iq))

def qvecTOiq(qvec):
    return ipTOiq(qvecTOip(qvec))

def qsqrdstr(thisqvec):
    sqvec = thisqvec.split()
    return int(sqvec[2])**2+int(sqvec[3])**2+int(sqvec[4])**2

def qvecINqsqrd(thisqsqrd):
    qvecout = []
    for iq in qvecSet:
        if qsqrdstr(iq) == thisqsqrd:
            qvecout.append(iq)
    return qvecout

DefMomList = [iqTOip(0),iqTOip(1),qvecTOip([-1,0,0])]

def SmallestMomList(data):
    dataout = []
    for idata in data[data.keys()[0]]:
        if all(idata in testdata for testdata in data.itervalues()):
            dataout.append(idata)
    return dataout


def MomOrderLists(MLread,MLsort,*SortThese):
    IndexShuff = [MLread.index(iMLs) for iMLs in MLsort]
    STout = []
    for iST in SortThese:
        STout.append(np.array(iST)[IndexShuff])
    return STout

def ipTOE(ip,mass):
    return np.sqrt((qsqrdstr(ipTOqstr(ip))*qunit)**2 + mass**2)


def qstrTOE(ip,mass):
    return np.sqrt((qsqrdstr(ip)*qunit)**2 + mass**2)

def iqTOE(ip,mass):
    return np.sqrt((qsqrdstr(iqTOqstr(ip))*qunit)**2 + mass**2)

def qvecTOE(ip,mass):
    return np.sqrt((qsqrdstr(qvecTOqstr(ip))*qunit)**2 + mass**2)

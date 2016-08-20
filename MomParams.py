#!/usr/bin/env python

import numpy as np
import os, errno


MomSqrdSet = [ '0','1','2','3','4','5','6','7','8','9' ]
QsqrdSet = ['qsqrd'+iq for iq in MomSqrdSet]
Maxqsqrd = np.max(np.array(MomSqrdSet).astype(int))
nxyz = 32
qunit = (2.0*np.pi)/float(nxyz)
hbarc = 0.1973269718 ## In GeV * fermi
latspace = 0.0074 ## In fermi

hbarcdivlat = hbarc/latspace
qunitPhys = qunit*hbarcdivlat

DefMass = 0.4662535526 #Lat Units
DefMassPhys = 0.4662535526*hbarcdivlat #Lat Units


def GetQsqrd(nqsqrd,Phys=True):
    if Phys:
        qsqrd = nqsqrd*(qunitPhys**2)
        Ep = np.sqrt(DefMassPhys**2 + qsqrd)
        return qsqrd - (Ep-DefMassPhys)**2
    else:
        qsqrd = nqsqrd*(qunit**2)
        Ep = np.sqrt(DefMass**2 + qsqrd)
        return qsqrd - (Ep-DefMass)**2

def makeqlist(thisMaxqsqrd):
    qlist = np.array([])
    for iq1 in range(-thisMaxqsqrd,thisMaxqsqrd+1):
        for iq2 in range(-thisMaxqsqrd,thisMaxqsqrd+1):
            for iq3 in range(-thisMaxqsqrd,thisMaxqsqrd+1):
                if iq1**2 + iq2**2 + iq3**2 > thisMaxqsqrd: continue
                qlist = np.append(qlist,'q = ' + str(iq1) + ' ' + str(iq2) + ' ' + str(iq3))
    return qlist

qvecSet = makeqlist(Maxqsqrd)
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


def qcondTOqstr(iq):
    return ' '.join(iq).replace('q ','q = ').replace('- ','-')

def qstrTOqcond(iq):
    return ' '.join(iq).replace('=','').replace(' ','')

def ipTOqcond(ip):
    return qstrTOqcond(ipTOqstr(ip))

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
                
def MakeMomDir(ip):
    thisip = ip
    if ' ' in ip: thisip = qstrTOqcond(ip)    
    iqsqrd = qsqrdstr(qcondTOqstr(thisip))
    return '/qsqrd'+str(iqsqrd)+'/'+thisip+'/'

def GetqcondFromFilename(filename):
    for ip in qvecSet:
        if qstrTOqcond(ip) in filename:
            return qstrTOqcond(ip)
    raise IOError('No Momenta in filename :' + filename)

## plistout = [ plistindex , pabc combinations , backward/forward , pxyzt]
def CreateSMOMPairs(plist):
    plistout = []
    for ip in plist:
        if len(ip) != 4: raise IOError('plist must contain 4vectors: ' + ' '.join(map(str,ip)))
        pa = ip
        pb = [-ip[0]] + ip[1:]
        pc = [ip[0],-ip[1],-ip[2],ip[3]]
        plistout.append([[pa,pb],[pa,pc],[pb,pc]])
    return plistout



def OutputSMOMPairs(filename):
    def strneg(intin):
        nspace = abs(intin)/10
        if intin < 0: nspace += 1
        return ' '*(2-nspace)+str(intin)
    thisplist = [[1, 1, 1, 2], [2, 2, 2, 4], [3, 3, 3, 6], [4, 4, 4, 8], [5, 5, 5, 10], [6, 6, 6, 12], [7, 7, 7, 14], [8, 8, 8, 16]]
    plistout = CreateSMOMPairs(thisplist)
    plistout = np.rollaxis(np.rollaxis(np.array(plistout),3),3)
    plistout = [[iip.flatten() for iip in ip] for ip in plistout]
    with open(filename,'w') as f:
        for pbf in plistout:
            f.write('\n')
            for pxyzt in pbf:
                f.write('( '+' '.join(map(strneg,pxyzt)) + ' ) \n')


def makeq4list(thisMinqsqrd,thisMaxqsqrd):
    qlist = []
    for iq1 in range(-thisMaxqsqrd,thisMaxqsqrd+1):
        for iq2 in range(-thisMaxqsqrd,thisMaxqsqrd+1):
            for iq3 in range(-thisMaxqsqrd,thisMaxqsqrd+1):
                for iq4 in range(-thisMaxqsqrd,thisMaxqsqrd+1):
                    if iq1**2 + iq2**2 + iq3**2 +  iq4**2 > thisMaxqsqrd: continue
                    if iq1**2 + iq2**2 + iq3**2 +  iq4**2 < thisMinqsqrd : continue
                    qlist.append([iq1,iq2,iq3,iq4])
    return np.array(qlist)




def CreateSMOMNewPairs(thisMinqsqrd,thisMaxqsqrd):
    def Myqsqrd(iq):
        return sum(iq**2)
    Plist = PPlist = makeq4list(thisMinqsqrd,thisMaxqsqrd)
    outplist,outpplist,outqlist,outomegalist = [],[],[],[]
    outp2list,outpp2list = [],[]
    for ip in Plist:
        for ipp in PPlist:
            ip2 = Myqsqrd(ip)
            ipp2 = Myqsqrd(ipp)
            if ip2 != ipp2: continue
            iq = ipp - ip
            iq2 = Myqsqrd(iq)
            w = iq2/ip2
            if not (0 < w < 4): continue
            if w in outomegalist and ip2 in outp2list: continue
            outplist.append(ip)
            outpplist.append(ipp)
            outp2list.append(ip2)
            outpp2list.append(ipp2)
            outqlist.append(iq)
            outomegalist.append(w)
    return outplist,outpplist,outqlist,outomegalist
            

            

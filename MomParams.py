#!/usr/bin/env python

import numpy as np
import os, errno


CHROMA = True ## For using chroma output

# MomSqrdSet = [ '0','1','2','3','4','5','6','7','8','9' ]
MomSqrdSet = [ '0','1','2','3','4']
QsqrdSet = ['qsqrd'+iq for iq in MomSqrdSet]
Maxqsqrd = np.max(np.array(MomSqrdSet).astype(int))
nxyz = 32
qunit = (2.0*np.pi)/float(nxyz)
hbarc = 0.1973269718 ## In GeV * fermi
# latspace = 0.074 ## In fermi
latspace = 0.0907 ## In fermi

hbarcdivlat = hbarc/latspace
qunitPhys = qunit*hbarcdivlat

DefMass = {}
DefMass['1370000'] = 0.7277 #Lat Units
DefMass['1375400'] = 0.5554004119 #Lat Units

DefPionMass = {}
DefPionMass['1370000'] = 0.32242 #Lat Units 
DefPionMass['1375400'] = 0.18903 #Lat Units 

def GetMpi(kappa,Phys=True):
    if Phys:
        return ' {:.3f} GeV '.format(DefPionMass[str(kappa)]*hbarcdivlat)
    else:
        return ' {:.3f} '.format(DefPionMass[str(kappa)])

DiagPList = [[1, 1, 1, 2], [2, 2, 2, 4], [3, 3, 3, 6], [4, 4, 4, 8], [5, 5, 5, 10], [6, 6, 6, 12], [7, 7, 7, 14], [8, 8, 8, 16]]
DiagPListtdiv = [[1, 1, 1, 1], [2, 2, 2, 2], [3, 3, 3, 3], [4, 4, 4, 4], [5, 5, 5, 5], [6, 6, 6, 6], [7, 7, 7, 7], [8, 8, 8, 8]]

def GetQsqrd(nqsqrd,Mass,Phys=True):
    if Phys:
        MassPhys = Mass*hbarcdivlat
        qsqrd = nqsqrd*(qunitPhys**2)
        Ep = np.sqrt(MassPhys**2 + qsqrd)
        return qsqrd - (Ep-MassPhys)**2
    else:
        qsqrd = nqsqrd*(qunit**2)
        Ep = np.sqrt(Mass**2 + qsqrd)
        return qsqrd - (Ep-Mass)**2

def makeqlist(thisMaxqsqrd):
    qlist = np.array([])
    for iq1 in range(-thisMaxqsqrd,thisMaxqsqrd+1):
        for iq2 in range(-thisMaxqsqrd,thisMaxqsqrd+1):
            for iq3 in range(-thisMaxqsqrd,thisMaxqsqrd+1):
                if iq1**2 + iq2**2 + iq3**2 > thisMaxqsqrd: continue
                qlist = np.append(qlist,'q = ' + str(iq1) + ' ' + str(iq2) + ' ' + str(iq3))
    return qlist

def Avgqlist(thisMaxqsqrd):
    qlist = np.array([])
    qsqrdlist = []
    for iq1 in range(thisMaxqsqrd+1):
        for iq2 in range(thisMaxqsqrd+1):
            for iq3 in range(thisMaxqsqrd+1):
                if iq1**2 + iq2**2 + iq3**2 > thisMaxqsqrd: continue
                if iq1**2 + iq2**2 + iq3**2 in qsqrdlist: continue
                qsqrdlist.append(iq1**2 + iq2**2 + iq3**2)
                # qlist = np.append(qlist,'q = ' + str(iq3) + ' ' + str(iq2) + ' ' + str(iq1))
                qlist = np.append(qlist,'q = ' + str(iq1) + ' ' + str(iq2) + ' ' + str(iq3))
    return qlist
    

def Chromaqlist(thisMaxqsqrd):
    qlist = np.array([])
    for iq1 in range(-thisMaxqsqrd,thisMaxqsqrd+1):
        for iq2 in range(-thisMaxqsqrd,thisMaxqsqrd+1):
            for iq3 in range(-thisMaxqsqrd,thisMaxqsqrd+1):
                if iq1**2 + iq2**2 + iq3**2 > thisMaxqsqrd: continue
                qlist = np.append(qlist,'q = ' + str(iq3) + ' ' + str(iq2) + ' ' + str(iq1))
    return qlist

if CHROMA:
    qvecSet = Chromaqlist(Maxqsqrd)
    qvecAvgSet = Avgqlist(Maxqsqrd)
else:
    qvecSet = makeqlist(Maxqsqrd)
    qvecAvgSet = makeqlist(Maxqsqrd)

    
nmom = len(qvecSet)
qhigh = nmom/2

def OrderMomList(momset):
    qout = []
    for iq in qvecSet:
        if iq in momset:
            qout.append(iq)
    return qout

## momlist in p indexing
def DragpZ(momlist,Avg=False):
    if iqTOip(0) in momlist:
        momlist.insert(0, momlist.pop(momlist.index(iqTOip(0))))
    return momlist

## momlist in str indexing
def DragpZstr(momlist,Avg=False):
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


def ipTOqstr(ip,Avg=False):
    if Avg:
        return qvecAvgSet[int(ip)]
    else:
        return qvecSet[int(ip)]


def qcondTOqstr(iq):
    return ' '.join(iq).replace('q ','q = ').replace('- ','-')

def qstrTOqcond(iq):
    return ' '.join(iq).replace('=','').replace(' ','')

def ipTOqcond(ip,Avg=False):
    return qstrTOqcond(ipTOqstr(ip,Avg=Avg))

def iqTOqstr(iq):
    return qvecSet[iqTOip(iq)]

def ipZshiftTOqstr(ip):
    return qvecSetZfirst[ip]
    
def qvecTOqstr(qvec,Avg=False):
    return qvecSet[qvecTOip(qvec,Avg=Avg)]

def qstrTOqvec(qstr):
    sqvec = qstr.split()
    return [int(sqvec[2]),int(sqvec[3]),int(sqvec[4])]


def qstrTOip(qstr,Avg=False):
    return qvecTOip(qstrTOqvec(qstr),Avg=Avg)

def qstrTOiq(qstr):
    return qvecTOiq(qstrTOqvec(qstr))

def ipTOqvec(ip):
    sqvec = qvecSet[ip].split()
    return [int(sqvec[2]),int(sqvec[3]),int(sqvec[4])]

def qvecTOip(qvec,Avg=False):
    strqvec = 'q = ' + str(qvec[0]) + ' ' + str(qvec[1]) + ' ' + str(qvec[2])
    if Avg:
        return np.where(qvecAvgSet==strqvec)[0][0]
    else:        
        return np.where(qvecSet==strqvec)[0][0]

def iqTOqvec(iq):
    return ipTOqvec(iqTOip(iq))

def qvecTOiq(qvec,Avg=False):
    return ipTOiq(qvecTOip(qvec,Avg=Avg))

def qsqrdstr(thisqvec):
    sqvec = thisqvec.split()
    return int(sqvec[2])**2+int(sqvec[3])**2+int(sqvec[4])**2

def qvecINqsqrd(thisqsqrd):
    qvecout = []
    for iq in qvecSet:
        if qsqrdstr(iq) == thisqsqrd:
            qvecout.append(iq)
    return qvecout

DefMomList = [iqTOip(0),iqTOip(1),qvecTOip([-1,0,0],Avg=False)]

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


def ipTOE(ip,mass,Avg=False):
    return np.sqrt((qsqrdstr(ipTOqstr(ip,Avg=Avg))*(qunit**2)) + mass**2)


def qstrTOE(ip,mass):
    return np.sqrt(qsqrdstr(ip)*(qunit**2) + mass**2)

def qstrTOEBoot(ip,mass):
    value = ((qsqrdstr(ip)*(qunit**2)) + mass**2)
    value.sqrt()
    return value


LatDisDenominator = 0.
# LatDispList = [0.,1.,2.,4.]
LatDispList = [0.]

DispKeyList = []
for iDisp in  LatDispList:
    if iDisp == 0.:
        DispKeyList.append('E^{2} = p^{2} + m^{2}')
    else:
        DispKeyList.append('Lat\ Disp\ Den='+str(iDisp))
                                                

def qstrTOLatEBoot(ip,mass,Disp=LatDisDenominator):
    value = (mass/Disp)
    value.sinh()
    ipvec = qstrTOqvec(ip)
    psqrdval = np.sum((np.sin(np.array(ipvec)*qunit)/Disp)**2)
    value = (psqrdval + value**2)
    value.sqrt()
    value.arcsinh()
    return value*Disp


def qstrTOLatEList(ip,mass,Disp=LatDisDenominator):
    value = np.sinh(np.array(mass)/Disp)
    ipvec = qstrTOqvec(ip)
    psqrdval = np.sum((np.sin(np.array(ipvec)*qunit)/Disp)**2)
    value = np.arcsinh(np.sqrt(psqrdval + value**2))
    return value*Disp

def iqTOE(ip,mass):
    return np.sqrt((qsqrdstr(iqTOqstr(ip))*(qunit**2)) + mass**2)

def qvecTOE(ip,mass):
    return np.sqrt((qsqrdstr(qvecTOqstr(ip))*(qunit**2)) + mass**2)



## expecting qstr
## MassBoot [tsink] BS 
def ScaledEffMass(ip,MassBoot,DispIn=LatDisDenominator):
    outdict = []
    for tboot in MassBoot:
        if DispIn == 0.:
            outdict.append(qstrTOEBoot(ip,tboot))
        else:
            outdict.append(qstrTOLatEBoot(ip,tboot,Disp=DispIn))
        outdict[-1].Stats()
    return np.array(outdict)


## expecting qstr
## MassBoot [tsink, bootlist ] 
def ScaledEffMassList(ip,MassList,DispIn=LatDisDenominator):
    outdict = []
    for tboot in MassList:
        if DispIn == 0.:
            outdict.append(qstrTOE(ip,np.array(tboot)))
        else:
            outdict.append(qstrTOLatEList(ip,tboot,Disp=DispIn))
    return np.array(outdict)

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
    plistout = CreateSMOMPairs(DiagPlist)
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



def GetAvgMom(qstr):
    if CHROMA:
        for ip in qvecAvgSet:
            if qsqrdstr(ip) == qsqrdstr(qstr):
                return ip
        raise IOError('Mom Not Found in AvgList')
    else:
        return qstr        

def GetAvgMomip(ip):
    if CHROMA:
        for ipc,thep in enumerate(qvecAvgSet):
            if qsqrdstr(thep) == qsqrdstr(ipTOqstr(ip)):
                return ipc
        raise IOError('Mom Not Found in AvgList')
    else:
        return qstr        

def GetAvgMomTotal(qlist):    
    outlist = []
    for iq in qlist:
        outlist.append(GetAvgMom(iq))
    return outlist

def GetAvgMomTotalip(iplist):
    outlist = []
    for iq in iplist:
        outlist.append(GetAvgMomip(iq))
    return outlist

def SortAvgMomList(qlist):
    qout = []
    for imom in qvecAvgSet:
        if imom in qlist:
            qout.append(imom)
    return qout

def SortAvgMomListip(qlistip):
    qout = []
    avgip = [qstrTOip(ip,True) for ip in qvecAvgSet]
    for imom in avgip:
        if imom in qlistip:
            qout.append(imom)
    return qout

def GetAvgMomList(qlist,sort=True):
    outlist = []
    for iq in qlist:
        iqavg = GetAvgMom(iq)
        if iqavg not in outlist:
            outlist.append(iqavg)
    if sort:
        return SortAvgMomList(outlist)
    else:
        return outlist

def GetAvgMomListip(iplist,sort=True):
    outlist = []
    for iq in iplist:
        iqavg = GetAvgMomip(iq)
        if iqavg not in outlist:
            outlist.append(iqavg)
    if sort:
        return SortAvgMomListip(outlist)
    else:
        return outlist


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
            


def GetCharRad(msqVal,Fval=1.):
    return 12*Fval* hbarc**2/msqVal 
# return msqVal


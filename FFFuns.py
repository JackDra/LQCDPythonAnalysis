#!/usr/bin/env python

import numpy as np
from Params import *
from MiscFuns import *
from FFSympy import *
import re
from collections import OrderedDict

CreateGamma('sakurai')

##REMEBER to remove 'cmplx' from opperators

def PicknFF(CT):
    if CT in ['Scalar','PsScalar']:
        return 1
    elif CT in ['Vector','PsVector']:
        return 2
    elif CT == 'Tensor':
        return 3
    else:
        return -1

def CreateEs(momvec,thismass,curr=False):
    momvecsqrd = sum([iq**2 for iq in momvec])
    return np.sqrt(thismass**2+momvecsqrd)

def Create4Mom(thisqvec,thisppvec,thismass):
    thispvec = [ipp-iq for ipp,iq in zip(thisppvec,thisqvec)]
    thisp = [1.0j*CreateEs(thispvec,thismass)] + thispvec + [1.0j*CreateEs(thispvec,thismass)] 
    thispp = [1.0j*CreateEs(thisppvec,thismass)] + thisppvec + [1.0j*CreateEs(thisppvec,thismass)]
    thisq = [ipp-ip for ipp,ip in zip(thispp,thisp)]
    return thisp,thisq,thispp


def SubsMom(thisFFun,thisp,thispp,thismass):
    thisFFun = subsvector(Epp,pxp,pyp,pzp,thispp[4]/(1j),thispp[1],thispp[2],thispp[3],thisFFun)
    thisFFun = subsvector(Ep,px,py,pz,thisp[4]/(1j),thisp[1],thisp[2],thisp[3],thisFFun)
    return thisFFun.subs(m,thismass)

def ScalarFF(opp,thisqvec,thisppvec,thismass):
    thisp,thisq,thispp = Create4Mom(thisqvec,thisppvec,thismass)
    term1 = FFunCheck(opp,thisp,thispp,thismass)
    rcheck,ccheck = abs(term1.real)<myeps,abs(term1.imag)<myeps    
    return [term1],not rcheck, not ccheck

def VectorFF(opp,thisqvec,thisppvec,thismass):
    thisp,thisq,thispp = Create4Mom(thisqvec,thisppvec,thismass)
    term1 = FFunCheck(opp,thisp,thispp,thismass)
    term2 = 0.0j
    for i in [1,2,3,4]:
        if i != int(opp[-1]) :
            term2 += FFunCheck(opp+'g'+str(i),thisp,thispp,thismass)*thisq[i]
    term2 = term2/(2.0*thismass)
    rcheck,ccheck = abs(term1.real)<myeps and abs(term2.real)<myeps,abs(term1.imag)<myeps and abs(term2.imag)<myeps
    return [term1,term2],not rcheck, not ccheck

def PsVectorFF(opp,thisqvec,thisppvec,thismass):
    thisp,thisq,thispp = Create4Mom(thisqvec,thisppvec,thismass)
    index1 = opp[-1]
    if index1 == '5': index1 = opp[-3]
    term1 = FFunCheck(opp,thisp,thispp,thismass)
    term2 = (1.0j*FFunCheck(opp.replace('g'+index1,''),thisp,thispp,thismass)*thisq[int(index1)])/(2.0*thismass)
    rcheck,ccheck = abs(term1.real)<myeps and abs(term2.real)<myeps,abs(term1.imag)<myeps and abs(term2.imag)<myeps
    return [term1,term2],not rcheck, not ccheck


def TensorFF(opp,thisqvec,thisppvec,thismass):
    thisp,thisq,thispp = Create4Mom(thisqvec,thisppvec,thismass)
    ## P = pp + p = 2*pp -q
    thisP = [ip+ipp for ip,ipp in zip(thisp,thispp)]
    gammalist,Proj = re.findall('g.',opp),re.search('P.',opp).group()
    index1,index2 = int(gammalist[0][-1]),int(gammalist[1][-1])
    term1 = 1.0j*FFunCheck(opp,thisp,thispp,thismass)
    term2 = 1.0j*(FFunCheck(opp.replace(gammalist[1],''),thisp,thispp,thismass)*thisq[index2] - 
             FFunCheck(opp.replace(gammalist[0],''),thisp,thispp,thismass)*thisq[index1] )/(2.0*thismass)
    # print ''
    # print 'term2'
    # print Proj,gammalist
    # print FFunCheck(opp.replace(gammalist[1],''),thisp,thispp,thismass) , thisq[index2]
    # print FFunCheck(opp.replace(gammalist[1],''),thisp,thispp,thismass)*thisq[index2]
    # print FFunCheck(opp.replace(gammalist[0],''),thisp,thispp,thismass) , thisq[index1]
    # print FFunCheck(opp.replace(gammalist[0],''),thisp,thispp,thismass)*thisq[index1]
    # print ''
    term3 = FFunCheck(Proj+'I',thisp,thispp,thismass)*(thisP[index1]*thisq[index2] -
                               thisP[index2]*thisq[index1] )/(2.0*thismass**2)
    # print 'term3'
    # print FFunCheck(Proj+'I',thisp,thispp,thismass)
    # print thisP[index1]*thisq[index2]
    # print thisP[index2]*thisq[index1] 
    # print ''
    ##Discrepancy? i think
    rcheck,ccheck = (abs(term1.real)<myeps and abs(term2.real)<myeps and abs(term3.real)<myeps,
                     abs(term1.imag)<myeps and abs(term2.imag)<myeps and abs(term3.imag)<myeps)
    return [term1,term2,term3], not rcheck,not ccheck

def CombineVector(thisFF,thisMass):
    ## FF { { momsqrd } { Boot/Avg/Chi } }
    FFout = OrderedDict()
    for iq,qFF in thisFF.iteritems():
        if len(qFF.keys()) > 0:
            FFout[iq] = {}
            qvecsqrd = int(iq.replace('qsqrd',''))*(qunit**2)
            Ep = np.sqrt(thisMass['Avg']**2 + qvecsqrd**2)
            Qsqrd = qvecsqrd - (Ep-thisMass['Avg'])**2
            FFout[iq]['Chi'] = qFF['Chi']
            if 'Boot' in qFF.keys():
                FFout[iq]['Boot'] = []
                FFout[iq]['Avg'] = []

                FFout[iq]['Boot'].append(qFF['Boot'][0] - (Qsqrd/(4*thisMass['Avg']**2))*qFF['Boot'][1])
                FFout[iq]['Boot'][-1].Stats()
                FFout[iq]['Avg'].append(FFout[iq]['Boot'][-1].Avg)

                FFout[iq]['Boot'].append(qFF['Boot'][0] + qFF['Boot'][1])
                FFout[iq]['Boot'][-1].Stats()
                FFout[iq]['Avg'].append(FFout[iq]['Boot'][-1].Avg)
            else:
                FFout[iq]['Avg'] = []
                FFout[iq]['Avg'].append(qFF['Avg'][0] - (Qsqrd/(4*thisMass['Avg']**2))*qFF['Avg'][1])
                FFout[iq]['Avg'].append(qFF['Avg'][0] + qFF['Avg'][1])
    return FFout

##Same as above, but only checks for 0:

#!/usr/bin/env python
from MiscFuns import *
from FFFuns import *
from FitFunctions import *
from Params import DefProjDerList
import cPickle as pickle

gi = ['g'+str(i) for i in [1,2,3,4]]
gig5 = ['g'+str(i)+'g5' for i in [1,2,3,4]]
Pi = ['P'+str(i) for i in [4,3]]
GMAi = ['GMA'+str(i) for i in [4,3]]
gigj = []
for i,gii in enumerate(gi):
    for j,gij in enumerate(gi):
         if j>i:
             gigj.append(gii+gij)

PiI    = ElongateName(Pi,['I'])
Pig5   = ElongateName(Pi,['g5'])
Pigi   = ElongateName(Pi,gi)
Pigig5 = ElongateName(Pi,gig5)
Pigigj = ElongateName(Pi,gigj)

DictGMAiI = {iPi:['I'] for iPi in GMAi}
DictGMAig5 = {iPi:['g5'] for iPi in GMAi}
DictGMAigi = {iPi:gi for iPi in GMAi}
DictGMAigig5 = {iPi:gig5 for iPi in GMAi}
DictGMAigigj = {iPi:gigj for iPi in GMAi}
DictGMA4giDi = ReadProjDerList

CurrTypes = ['Scalar','Vector','PsScalar','PsVector','Tensor']
DerCurrTypes = ['giDi']
AllCurrTypes = CurrTypes + DerCurrTypes

DictCurrOpps = {'Scalar'   : DictGMAiI,
                'Vector'   : DictGMAigi,
                'PsScalar' : DictGMAig5,
                'PsVector' : DictGMAigig5,
                'Tensor'   : DictGMAigigj,
                'giDi'   : DictGMA4giDi,
                'Test'    : {'GMA4':['g4']}}


CurrOpps = {'Scalar'   : PiI,
            'Vector'   : Pigi,
            'PsScalar' : Pig5,
            'PsVector' : Pigig5,
            'Tensor'   : Pigigj}

CurrFFs = {'Scalar'   : ScalarFF,
           'Vector'   : VectorFF,
           'PsScalar' : ScalarFF,
           'PsVector' : PsVectorFF,
           'Tensor'   : TensorFF}


FFFitFuns = {'Scalar'   : FormFactorO1,
             'Vector'   : FormFactorO2,
             'PsScalar' : FormFactorO1,
             'PsVector' : FormFactorO2,
             'Tensor'   : FormFactorO3}

NoFFPars = {'Scalar'   : 1,
            'Vector'   : 2,
            'PsScalar' : 1,
            'PsVector' : 2,
            'Tensor'   : 3}

NoFFList = {'Scalar'   : ['FF1'],
            'Vector'   : ['FF1','FF2'],
            'PsScalar' : ['FF1'],
            'PsVector' : ['FF1','FF2'],
            'Tensor'   : ['FF1','FF2','FF3']}



def FindMomFromGamma(igamma,thisMomList=qvecSet):
    momout = []
    testgamma = igamma.replace('cmplx','').replace('doub','').replace('sing','')
    for thisFF,iCurrOpps in CurrOpps.iteritems():
        if testgamma in iCurrOpps:
            thisFFtype = thisFF
    for iq in thisMomList:
        iqvec = np.array(qstrTOqvec(iq))*qunit
        dump,rcheck,ccheck =  CurrFFs[thisFFtype](testgamma,iqvec.tolist(),[0,0,0],1.0)
        if rcheck and 'cmplx' not in igamma:
            momout.append(iq)
        elif ccheck and 'cmplx' in igamma:
            momout.append(iq)
    return momout
                                    
            
def DumpAllMomLists():
    mygammalist = ['P4'+igamma for igamma in AllGammaSet] +['P3'+igamma for igamma in AllGammaSet] 
    for igamma in mygammalist:
        print ' Dumping MomList: ' + igamma
        outfile = momlistdir + igamma + '.p'
        with open(outfile,'wb') as f:
            pickle.dump(FindMomFromGamma(igamma),f)

def GetMomFromGamma(igamma,thisMomList=qvecSet):
    infile = momlistdir + igamma.replace('doub','').replace('sing','') + '.p'
    with open(infile,'rb') as f:
        momlistout = pickle.load(f)
    return momlistout

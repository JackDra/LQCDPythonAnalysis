#!/usr/bin/env python

import numpy as np
from MomParams import *
from MiscFuns import *
import operator as op
import socket
import os
import re
try:
    import xmltodict
except:
    raise EnvironmentError("please have xmltodict installed (git clone \"https://github.com/martinblech/xmltodict.git\") and install")

##CHANGE StringFix.py FUNCTION##

##DEAULT TAKEN FROM OSF from: /raid/jdragos/data/results/ReadSetk12090/cfuns/twopt/OSFCM/state1sum22to23dt4twoptm0.txt
##date/time calculated: 00:29 12 Oct 2015
# DefMassVal = {'22-31' : {'Avg':DefMass,'Std': 0.0028592256 ,'Chi': 0.4851841332}}

DefMassVal = {}
DefMassVal['fitr6-13'] = OrderedDict()
DefMassVal['fitr6-13']['Avg'] = DefMass
DefMassVal['fitr6-13']['Std'] = 0.0071123989
DefMassVal['fitr6-13']['Chi'] = 0.2782012297
DefMassVal['fitr6-13']['File'] = '/home/jackdra/PHD/CHROMA/TestVar/scratch/results/ReadSetBarN17k1375400/cfun/twopt/OSFCM/qsqrd0/q000/state1PoF0to3dt3twoptm0q000.xml'

# THISMACHINE = socket.gethostname()

# if 'erwin' in THISMACHINE:
#     scriptdir = "/home/accounts/jdragos/scripts/LQCDPythonAnalysis/"
#     datadir = '/raid/jdragos/data/'
#     AnaProc = 60 # number of max processors (computer specific)
# elif 'henley' in THISMACHINE:
#     datadir = '/home/henley/jdragos/PHD/ErwinBackup/'
#     scriptdir = datadir + "LQCDPythonAnalysis/"
#     AnaProc = 10 # number of max processors (computer specific)
# elif 'JackLappyUbu' in THISMACHINE:
    # datadir = '/home/jackdra/PHD/DataAnalysis/'
    # scriptdir = datadir + "LQCDPythonAnalysis/"
    # AnaProc = 10 # number of max processors (computer specific)
# else:
#     raise EnvironmentError('Machine not recongnised: \n'+THISMACHINE+'\nplease add to Params.py')

# ## kappa values, pick your kappa
# kappa = 12090
# # kappa = 12104
# # ListOrSet = 'ReadList'
# ListOrSet = 'ReadSet'

# ## PoF is pencil of function method, you can set VarPref to what you want, and set PoFShifts to 0 if you want regular variational method
# if kappa == 12104:
#     PoFShifts = 1
#     PoFOrSum = 'PoF'
#     VarPref = PoFOrSum+str(PoFShifts)
# elif kappa == 12090:
#     VarPref = 'sum'
#     PoFOrSum = VarPref
#     PoFShifts = 0
# PoFShifts = 1
# PoFOrSum = 'PoF'
# VarPref = PoFOrSum+str(PoFShifts)


if not os.path.isfile('./setup.cfg'): raise EnvironmentError('Please Run Setup.py to create setup.cfg')
RunNconfs = False
with open('./setup.cfg','r') as f:
    thisread = ''
    for line in f:
        thisline = line.strip()
        if len(thisline) > 0:
            if thisline[-1] == ':':
                thisread = thisline[:-1]
            else:
                if 'scriptdir' in thisread:
                    scriptdir = thisline
                elif 'datadir' in thisread:
                    datadir = thisline
                elif 'AnaProc' in thisread:
                    AnaProc = int(thisline)
                elif 'ListOrSet' in thisread:
                    ListOrSet = thisline
                elif 'PoFShifts' in thisread:
                    PoFShifts = int(thisline)
                elif 'PoFDelta' in thisread:
                    PoFDelta = int(thisline)
                elif 'kappa' in thisread:
                    kappa = int(thisline)
                elif 'Nconfs' in thisread:
                    if 'False' in thisline:
                        RunNconfs = False
                    else:
                        RunNconfs = int(thisline)
                elif 'Debug' in thisread:
                    Debug = 'True' in thisline or 'true' in thisline
                        
PoFC2C3Dis = ''
NewFileFlag = PoFC2C3Dis


kappalist = ['k'+str(kappa),'xsrc1k'+str(kappa),'nboot1kk'+str(kappa),'XAvgk'+str(kappa)]
# kappalist = ['k'+str(kappa),'xsrc1k'+str(kappa)]
kappaflags = [ik.replace('k'+str(kappa),'') for ik in kappalist]
ScalarNorm = 1 # normalisation for Scalar current
PsScalarNorm = 1 # normalisation for Pseudo Scalar current
VectorNorm = 0.736862 # normalisation for Vector Current from arXiv:1006.1164v2 beta = 1.9, Iwasaki gauge action, clover term Csw=1.7
PsVectorNorm = 1 # normalisation for Pseudo Vector current
TensorNorm = 1 # normalisation for Tensor Current
MomFracNorm = 1 # normalisation for Momentum Fraction

# ScalarNorm = 0.6822 # normalisation for Scalar current
# PsScalarNorm = 0.4948 # normalisation for Pseudo Scalar current
# VectorNorm = 0.8574 # normalisation for Vector Current
# PsVectorNorm = 0.8728 # normalisation for Pseudo Vector current
# TensorNorm = 0.9945 # normalisation for Tensor Current
# MomFracNorm = -1.067/DefMass # normalisation for Momentum Fraction


myeps = np.finfo(0.0).eps
DeleteNanCfgs = False ## Deletes configs that have nans
ForceVecNorm = False ## Forces the vector current to be normalised to 2 for doublet and 1 for singlet at zero momentum transfer
ForcePos = False ## Forces all non-form factor graphs to be positive
MultiCoreFitting = False # Multicore for Boot Fitting, not needed in current build
DoMulticore = True # Runs multicore wherever implemented
DoContentsCheck = False # True makes sure the xml file has the correct momenta first field, turn off for more performance
OnlySelVar = True # Selects "ThePickedSumVar" (see below) variable for all the method calculations instead of all
DoNorm = False # normalises the 2 point function (see CMSTech.py)
DoSym = False # symmetrises the 2 point function (see CMSTech.py)
# VarMethodMethod = 'Regular' # for solving the Variational method, different ways of doing it/
# VarMethodMethod = 'Symmetriceigh' ## Symmetic matrix construction WITH symmetric solver
# VarMethodMethod = 'Symmetric' ## Symmetic matrix construction WITHOUT symmetric solver
VarMethodMethod = 'AxBxlSolve' ## solve Ax = Bxc system directly (generalised eigenvalue problem).
NoSFRfacScale = False # Turn on to only scale the R function by sqrt((Epp+m)(Ep+m)/EppEp) for form factor creation
ReadPoF2pt = True # Create PoF using already calculated eigenvectors. This is used if the statistics or solver method has changed.
DeCorrPoF = False ## used for debugging the pencil of function method (decorrelation problem) !!!!!DEPRECIATED, LEAVE FALSE!!!!!
TimeInv = False ## uses time invariance to calculate the Pencil of Function method/ Oposed to calculating [tsource,tsource-1,...,tsource-PoFShifts]
DoCM = False ## does correlation matrix result ( no PoF) 
PlotMonte = False ## Plots montecarlo time history of NNQ at time slice MonteTime and flowtime MonteFlow
PlotXSrcDep = False ## Plots montecarlo time history of NNQ at time slice MonteTime and flowtime MonteFlow
PlotAutocorr = False
MonteTime = 7
MonteFlow = 40

if 'XAvg' in ListOrSet:
    XAvg = True ## averages over source position locatinos before bootstrapping
else:
    XAvg = False
ExactXSrcNumber = False ## makes it so there are the same number of sources for each configuration (hardwired to the first configuration found)
ForceXSrcLen = False ## forces so that any more sorce locations after XSrcLen are ignored
XSrcLen = 15
NoXAvg = False ## Does each source separatly for each 

ForceMinXSrcLen = True ## only calculates with a minumum of MinXSRCLen of sources per gauge field
MinXSrcLen = 10
##DEBUG toggles (True/False):
# Debug = True # for debugging, toggles alot of print statements on
# DEBUGPoF = True ## Temp debug parameter to test PoF (REMOVE ONCE DONE)
ScaleByP4g4 = False ## scales out all operators by P4g4 instead of 2 point correlator at tsink for Ratio value (RF)
ShowConfNum = Debug # debugging, show number of configs during read
PrintRead = not DoMulticore # Screws up output if on and doing mulitcore reading
DoCmplx = True # reads complex opperator values as well as real values, should be on
DoCons = False # reads Conserved vector current NOT IMPLEMENTED YET, USING BELOW UNTILL I GET AROUDN TO IT. ONLY WORKS WITH CHROMA
RepWithCons = False # TEMPORARY, overrides vector current with Conserved vector current ONLY WORKS WITH CHROMA

DefWipe = False # Wipes sets before running RunMcorr, only doing if debugging, if working, should be False
MakeGraphDat = True # Creates a .dat file of the values plotted for the corresponding graph (where implemented)
PhysicalUnits = True # uses lattice momenta or physical momenta
ForceNoDer = False # Forces the fitting LLSBoot to use manual derivaive calculation
# MesOrBar = 'Meson' ## selects meson or baryon to read
MesOrBar = 'Baryon'

if MesOrBar == 'Meson':
    InterpFlag = 'gamma_value'
    InterpNumb = '15'
    INg5 = '15'
elif MesOrBar == 'Baryon':
    InterpFlag = 'baryon_num'
    ## InterpNumb = '0' ## Old
    InterpNumb = '9' ## Proton
    if 'BarN' in ListOrSet:
        INg5 = ListOrSet.split('BarN')[-1]
        print 'INg5 set to ' , INg5
    else:
        INg5 = '17' ## 0 is Pp*g5 , 1 is g5 , 17 is g5*Pp 
    
if DoCons: RepWithCons = False

## Currenlty, not using Time Invariance in the Pencil of Function Analysis is only properly implemented and tested for 2-point correlation analysis
if TimeInv:
    CfunConfigCheck = True # Checks all two and three point correlators before running (turn to False if doing 2-pt corr analysis)
else:
    CfunConfigCheck = False # Checks all two and three point correlators before running (turn to False if doing 2-pt corr analysis)
CfunConfigCheck = True # Checks all two and three point correlators before running (turn to False if doing 2-pt corr analysis)


if Debug: print 'nconfs saved is: ' , RunNconfs

VarMassCutoff = 0.4 # used in correlation matrix for cutting artifacts out of eigenmass sorting.
if kappa == 12:
    kappas = kappa
    dirread = datadir+'/cfun/Kud0'+str(kappa)+'Ks0'+str(kappas)
    nt = 8
    nx = 4
else:
    kappas = 1364000
    if 'JackLappy' in socket.gethostname():
        dirread = datadir+'/cfunsg5/Kud0'+str(kappa)+'Ks0'+str(kappas)
    else:
        dirread = datadir+'/cfunPChroma/Kud0'+str(kappa)+'Ks0'+str(kappas)
    nt = 64
    nx = 32


outputdir = [datadir+'/results/'+ListOrSet+ikappa+'/' for ikappa in kappalist]
TCDir = datadir+'/topcharge/Kud0'+str(kappa)+'Ks0'+str(kappas)+'/PerGF/'
logdir = scriptdir+'../logdir/k'+str(kappa)+'/'
momlistdir = datadir+'momdir/'
pickledir = datadir+"pickledir/"
REvecDir = scriptdir+'REvecSave/k'+str(kappa)+'/'
RunMomList = qvecSet 
RunAvgMomList = qvecAvgSet
# For Debuggin, only use zero momenta
# RunMomList = [qvecSet[iqTOip(0)],qvecSet[iqTOip(3)]]
# RunMomList = [qvecSet[iqTOip(0)]]
# RunAvgMomList = [qvecAvgSet[0]]
# map(mkdir_p,outputdir)
mkdir_p(outputdir[0])
mkdir_p(pickledir)
mkdir_p(logdir)
mkdir_p(REvecDir)
mkdir_p(momlistdir)
ndim = [nx,nx,nx,nt]
latspace = 0.074 ## in fm
ns = 4


if 'ReadList' in ListOrSet:
    nboot = 2
elif 'ReadSet' in ListOrSet:
    if 'nboot1k' in ListOrSet:
        nboot = 1000
    else:
        if 'JackLappy' in socket.gethostname():
            nboot = 20
        else:
            nboot = 200
##DEBUGGING:
print nboot , ListOrSet, XAvg
# tsource = 17
tsource = 0
if TimeInv:
    # PoFtsourceList = map(str,[tsource]*(PoFShifts+1))
    PoFtsourceList = [str(tsource)]
else:
    # PoFtsourceList = map(str,range(tsource-(PoFShifts*PoFDelta),tsource+1,PoFDelta))
    PoFtsourceList = map(str,range(tsource,tsource+1+(PoFShifts*PoFDelta),PoFDelta))

# note: dim of StateSet < dim of SmearSet
GammaSet = ['I','g1','g2','g3','g4','g1g2','g1g3','g1g4','g2g3','g2g4','g3g4','g1g5','g2g5','g3g5','g4g5','g5']
GammaConsSet = ['Consg1','Consg2','Consg3','Consg4']
tflowlist = map(float,np.arange(0,1000,100)) ## indicies of flows to read
# tflowlist = map(float,np.arange(0,1000))

if DoCmplx:
    AllGammaSet = GammaSet + [igamma+'cmplx' for igamma in GammaSet]
else:
    AllGammaSet = GammaSet
DerSet = ['D1','D2','D3','D4']
GBSize = len(GammaSet)
DGBSize = GBSize * len(DerSet)

DGSet = []
for ider in DerSet:
    for igamma in GammaSet:
        DGSet.append(igamma+ider)

#this part is for ReadList (used for analying only specific configurations#
if kappa == 12:
    FileStruct = "Testing.lime"
    conflist = ['21_xsrc1_']
else:
    FileStruct = "RC32x64_B1900Kud01375400Ks01364000C1715"
    conflist = ['-a-004310_xsrc117',
                '-a-004010_xsrc105',
                '-a-003310_xsrc102']
                
                
filelist = dirread+"/@/"+FileStruct+"*"
            
nsrc = 200
xsrcList = ['xsrc'+str(ix)+'_' for ix in range(1,nsrc+1)]

SourceList = ['']



DefProjGammaList = { 'GMA4' : AllGammaSet ,'GMA3' : AllGammaSet }
ReadProjDerList = {'GMA4' :['g4D4','g1D1','g2D2','g3D3']}
# ReadProjDerList = {'GMA4' :[]}
# DefProjDerList = {'GMA4' :['giDi']}
DefProjDerList = {}
# DefDSList = ['doub']
DefDSList = ['doub','sing']
DefGammaList = []
for iDS in DefDSList:
    for Proj,GL in DefProjGammaList.iteritems():
        DefGammaList += [iDS+'P'+Proj[3]+iGL for iGL in GL]
    for Proj,GL in DefProjDerList.iteritems():
        DefGammaList += [iDS+'P'+Proj[3]+iGL for iGL in GL]
        
ReadGammaList = []
for iDS in DefDSList:
    for Proj,GL in DefProjGammaList.iteritems():
        ReadGammaList += [iDS+'P'+Proj[3]+iGL for iGL in GL]
    for Proj,GL in ReadProjDerList.iteritems():
        ReadGammaList += [iDS+'P'+Proj[3]+iGL for iGL in GL]

DefNoDSGammaList = []
for Proj,GL in DefProjGammaList.iteritems():
    DefNoDSGammaList += ['P'+Proj[3]+iGL for iGL in GL]
for Proj,GL in DefProjDerList.iteritems():
    DefNoDSGammaList += ['P'+Proj[3]+iGL for iGL in GL]
DefCombGammaList = DefGammaList+DefNoDSGammaList

# DeftoList = [18]
# DeftoList = range(1,10)
# DeftoList = range(1,10)
DeftoList = range(1,5)
# DeftoList = [20,21,22,23]
# DeftoList = [17,18,19,20,21,22,23]
# DeftoList = [17,18,19,20,21,22,23,24,25,26,27]
# DeftoList = [18,19,20,21,22]
# DeftoList = [16,17,18,19,20]
# DefdtList = [1,2,3,4,5,6]
# DefdtList = range(1,7)
# DefdtList = range(1,10)
DefdtList = range(1,8)
# DeftodtPicked = (18,2)
##MUST BE IN SORTING ORDER##
# DeftodtPicked = [(18,2),(20,2)]
# DefTvarPicked = ['CMto'+str(iDeftodtPicked[0])+'dt'+str(iDeftodtPicked[1]) for iDeftodtPicked in DeftodtPicked]
DeftodtPicked = [(3,3)]
# DeftodtPicked = [(1,1)]
# DeftodtPicked = [(18,2)]
DefTvarPicked = ['CMto'+str(iDeftodtPicked[0])+'dt'+str(iDeftodtPicked[1]) for iDeftodtPicked in DeftodtPicked]

# DeftodtPicked = Elongate(DeftoList,DefdtList)
# DefTvarPicked = ['CMto'+str(iDeftodtPicked[0])+'dt'+str(iDeftodtPicked[1]) for iDeftodtPicked in DeftodtPicked]

DeftodtList = []
DefTvarList = []
TwoPtDefTvarList = []
TwoTotDefTvarList = []
DefTvarDt1 = []
DefTvarDt2 = []
DefTvarDt3 = []
DefTvarDt4 = []
for ito in DeftoList:
    DefTvarDt1.append('CMto'+str(ito)+'dt1')
    DefTvarDt2.append('CMto'+str(ito)+'dt2')
    DefTvarDt3.append('CMto'+str(ito)+'dt3')
    DefTvarDt4.append('CMto'+str(ito)+'dt4')
    for idt in DefdtList:
        DeftodtList.append((ito,idt))
        DefTvarList.append('CMto'+str(ito)+'dt'+str(idt))
        TwoPtDefTvarList.append('to'+str(ito)+'dt'+str(idt))
        TwoTotDefTvarList.append('PoF'+str(PoFShifts)+'to'+str(ito)+'dt'+str(idt))
        TwoTotDefTvarList.append('REvecto'+str(ito)+'dt'+str(idt))

TwoTotDefTvarList += DefTvarList
DefTvarto16 = []
DefTvarto17 = []
DefTvarto18 = []
DefTvarto19 = []
DefTvarto20 = []
for idt in DefdtList:
    DefTvarto16.append('CMto16dt'+str(idt))
    DefTvarto17.append('CMto17dt'+str(idt))
    DefTvarto18.append('CMto18dt'+str(idt))
    DefTvarto19.append('CMto19dt'+str(idt))
    DefTvarto20.append('CMto20dt'+str(idt))
            
if OnlySelVar:
    AnatodtList = DeftodtPicked
    AnaTvarList = DefTvarPicked
else:
    AnatodtList = DeftodtList
    AnaTvarList = DefTvarList
    
# DefSmearList = ['16','32','64']
# DefSmearList = ['16']
DefSmearList = ['64']
# DefSmearList = ['8','16','32','64','128','256']
# DefSmearList = ['64','128']
# DefSmearList = ['32','64']
# DefSmearList = ['32','128']
# DefSmearList = ['32']
# DefSmearList = ['64']
# DefSmearList = ['128']
# if kappa == 12:
SingSmearList = [DefSmearList[0]] ## must be first element of DefSmearList for now FIX!!
# DefSmearList = ['32']
# DefSmearList = ['8','16','32']
StateSet = map(str,range(1,(PoFShifts+1)*len(DefSmearList)+1))
CMStateSet = map(str,range(1,len(DefSmearList)+1))
def GetStateSet(keystring):
    if 'PoF' in keystring: return StateSet
    else: return CMStateSet
PickedState = 1
PickedStateStr = 'state'+str(PickedState)
DefInterpList = ['nucleon']
# DefInterpList = ['nucleon','nucleon2']
DefSmList = ['sm'+ism for ism in DefSmearList]
SingSmList = ['sm'+ism for ism in SingSmearList]
DefInterpSmearList = ElongateName(DefInterpList,DefSmList)
##THIS NEEDS TO BE SET TO GET SIGN RIGHT ON CORRELATOR
CMTSinkList = [13]
AllTSinkList = [11,14,17,20,23]
AllTSinkListNoCM = []
for its in AllTSinkList:
    if its not in CMTSinkList:
        AllTSinkListNoCM.append(its)
AllTSinkListVar = [11,12,14,17,20,23]
AllTSinkShift = [it-tsource for it in AllTSinkList]
AllTSinkStrList = ['tsink'+str(its) for its in AllTSinkList]
AllTSinkStrListVar = ['tsink'+str(its) for its in AllTSinkListVar]
CMTSinkStrList = ['tsink'+str(its) for its in CMTSinkList]
AllTSinkStrListNoCM = ['tsink'+str(its) for its in AllTSinkListNoCM]

# AllREvecTSinkList = {'12104':[29],'12090':[32]}
# AllREvecTSinkList = {'12104':[29],'12090':[26,32]}
AllREvecTSinkList = {str(kappa):[]}
REvecTSinkList = AllREvecTSinkList[str(kappa)]
REvecTSinkStrList = ['tsink'+str(its) for its in REvecTSinkList]
DefREvecVarList = [18,2]
REvecTvarList = ['REvecto'+str(DefREvecVarList[0])+'dt'+str(DefREvecVarList[1])]
REvecFlagList = [PickedStateStr+iREvec for iREvec in REvecTvarList]

# REvecTvarList = []

if TimeInv:
    # DefPoFVarList = [[1,1]]
    DefPoFVarPicked = [[6,1]]
else:
    if kappa == 12:
        DefPoFVarPicked = [[1,1]]
    else:
        # DefPoFVarPicked = [[3,3],[1,3],[2,3],[4,3]]
        ## 1,3 not working for some reason
        DefPoFVarPicked = [[3,3],[2,3],[4,3]]
DefPoFTvarRef = DefPoFVarPicked[0]
PoFTvarPicked = ['PoF'+str(PoFShifts)+'to'+str(iPoF[0])+'dt'+str(iPoF[1]) for iPoF in DefPoFVarPicked]

if OnlySelVar:
    PoFTvarList = PoFTvarPicked
    DefPoFVarList = DefPoFVarPicked
else:
    DefPoFVarList = DeftodtList
    PoFTvarList = [itvar.replace('CM','PoF'+str(PoFShifts)) for itvar in DefTvarList]

    
# DefPoFVarList = [18,2]
# DefPoFVarList = [19,2]
# DefPoFVarList = [20,2]
# DefPoFVarList = [22,2]
AllPoFTSinkList = {str(kappa):[13]}
PoFTSinkList = AllPoFTSinkList[str(kappa)]
PoFTSinkStrList = ['tsink'+str(its) for its in PoFTSinkList]
##DEBUG##
# PoFTvarList = ['PoF'+PoFShifts+'Test']
# PoFReadTvarList = ['PoFTestnD'+str(PoFShifts)]
# ##uncomment below to restore##
## warning, below has been hard coded to fix the to16 is to17 problem i had
PoFDirTvarList = ['PoFto'+str(DefPoFVarList[0][0]-1)+'dt'+str(DefPoFVarList[0][1])]
PoFReadTvarList = ['PoFto'+str(DefPoFVarList[0][0]-1)+'dt'+str(DefPoFVarList[0][1])+'nD'+str(PoFShifts)]

PoFFlagList = [PickedStateStr+iPoF for iPoF in PoFTvarList]


# REvecPar26 = [[ 0.0004799, -0.0119381, 0.9999286 ],[0,0,0],[0,0,0]]
# REvecPar32 = [[ 0.0007567612, -0.0182265391, 0.9998335964],[0,0,0],[0,0,0]]

# if kappa == 12:
TSFFileFlags = ['CM','Tsink','test32','Small']
OSFFileFlags = ['CM','Tsink']
CfunCheckList = ['PoFRead']

MethodList = ['RF','Fits','SumMeth']+['TSF'+iTSF for iTSF in TSFFileFlags] + ['OSF'+iOSF for iOSF in OSFFileFlags]


def GetRenorm(thisstring):
    if 'giDi' in thisstring:
        return MomFracNorm
    elif re.search('g[1234]g5',thisstring):
        return PsVectorNorm
    elif re.search('g[1234]g[1234]',thisstring):
        return TensorNorm
    elif re.search('g[1234]',thisstring):
        return VectorNorm
    elif re.search('g5',thisstring):
        return PsScalarNorm
    elif re.search('I',thisstring):
        return ScalarNorm
    else:
        return 1.0


def mprint(*string):
    if isinstance(string,str):
        if not DoMulticore or Debug: print string
    else:
        if not DoMulticore or Debug: print ' '.join(map(str,list(string)))
        

TSinkDictList = {'REvec' : REvecTSinkList,
                 'REvecRead' : REvecTSinkList,
                 'CM' : CMTSinkList,
                 'cmRead' : CMTSinkList,
                 'Tsink' : AllTSinkList,
                 'tsinkRead' : AllTSinkListNoCM}


TSinkStrDictList = {'REvec' : REvecTSinkStrList,
                    'REvecRead' : REvecTSinkStrList,
                    'CM' : CMTSinkStrList,
                    'cmRead' : CMTSinkStrList,
                    'Tsink' : AllTSinkStrList,
                    'tsinkRead' : AllTSinkStrListNoCM}

if len(PoFTSinkList) > 0:
    TSinkDictList['PoF'] = PoFTSinkList+range(PoFTSinkList[-1]+1,PoFTSinkList[-1]+1+PoFShifts)
    TSinkDictList['PoFRead'] = PoFTSinkList+range(PoFTSinkList[-1]+1,PoFTSinkList[-1]+1+PoFShifts)
    TSinkStrDictList['PoF'] = PoFTSinkStrList+['tsink'+str(its) for its in range(PoFTSinkList[-1]+1,PoFTSinkList[-1]+1+PoFShifts)]
    TSinkStrDictList['PoFRead'] = PoFTSinkStrList+['tsink'+str(its) for its in range(PoFTSinkList[-1]+1,PoFTSinkList[-1]+1+PoFShifts)]


SmearDictList = {'PoF' : PoFFlagList,
                 'PoFRead' : DefSmearList,
                 'REvec' : REvecFlagList,
                 'REvecRead' : DefSmearList,
                 'CM' : DefSmList,
                 'cmRead' : DefSmearList,
                 'Tsink' : SingSmList,
                 'tsinkRead' : SingSmearList}

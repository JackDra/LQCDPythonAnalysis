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
                elif 'kappa' in thisread:
                    kappa = int(thisline)
                

        
PoFC2C3Dis = '665.'

myeps = np.finfo(0.0).eps
ScalarNorm = 0.6822 # normalisation for Scalar current
PsScalarNorm = 0.4948 # normalisation for Pseudo Scalar current
VectorNorm = 0.8574 # normalisation for Vector Current
PsVectorNorm = 0.8728 # normalisation for Pseudo Vector current
TensorNorm = 0.9945 # normalisation for Tensor Current
MultiCoreFitting = False # Multicore for Boot Fitting, not needed in current build
DoMulticore = True # Runs multicore wherever implemented
DoContentsCheck = False # True makes sure the xml file has the correct momenta first field, turn off for more performance
OnlySelVar = True # Selects "ThePickedSumVar" (see below) variable for all the method calculations instead of all
DoNorm = False # normalises the 2 point function (see CMSTech.py)

##DEBUG toggles (True/False):
ScaleByP4g4 = False ## scales out all operators by P4g4 instead of 2 point correlator at tsink for Ratio value (RF)
DebugCM = False # for debugging Correlation Matrix stuff
ShowConfNum = False # debugging, show number of configs during read
PrintRead = not DoMulticore # Screws up output if on and doing mulitcore reading
DoCmplx = True # reads complex opperator values as well as real values, should be on
DefWipe = False # Wipes sets before running RunMcorr, only doing if debugging, if working, should be False
    

VarMassCutoff = 0.4 # used in correlation matrix for cutting artifacts out of eigenmass sorting.

dirread = datadir+'/cfuns/k'+str(kappa)
outputdir = datadir+'results/'+ListOrSet+'k'+str(kappa)+'/'
logdir = scriptdir+'../logdir/k'+str(kappa)+'/'
momlistdir = datadir+'momdir/'
pickledir = datadir+"pickledir/"
REvecDir = scriptdir+'REvecSave/k'+str(kappa)+'/'
# RunMomList = qvecSet 
# For Debuggin, only use zero momenta
# RunMomList = [qvecSet[iqTOip(0)],qvecSet[qvecTOip([0,0,1])]]
RunMomList = [qvecSet[iqTOip(0)]]
mkdir_p(outputdir)
mkdir_p(pickledir)
mkdir_p(logdir)
mkdir_p(REvecDir)
mkdir_p(momlistdir)
nt = 64
nx = 32
ndim = [nx,nx,nx,nt]
latspace = 0.074 ## in fm
ns = 4


if ListOrSet == 'ReadList':
    nboot = 2
elif ListOrSet == 'ReadSet':
    nboot = 200
tsource = 16

# note: dim of StateSet < dim of SmearSet
GammaSet = ['I','g1','g2','g3','g4','g1g2','g1g3','g1g4','g2g3','g2g4','g3g4','g1g5','g2g5','g3g5','g4g5','g5']
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
if kappa == 12104:
    FileStruct = "b5p50kp"+str(kappa)+"0kp"+str(kappa)+".687."
    filelist = dirread+"/source1/@/"+FileStruct+"*"
elif kappa == 12090:
    ##CHECK THIS DEBUG##
    # FileStruct = "b5p50kp"+str(kappa)+"0kp"+str(kappa)+"0."
    if ListOrSet == 'ReadList':
        FileStruct = "b5p50kp"+str(kappa)+"0kp"+str(kappa)+"0."+PoFC2C3Dis
    elif ListOrSet == 'ReadSet':
        FileStruct = "b5p50kp"+str(kappa)+"0kp"+str(kappa)+"0."    
    filelist = dirread+"/source1/@/"+FileStruct+"*"
conflist = ['00105',
            '00115',
            '00125',
            '00265',
            '00275',
            '00285']
            

SourceList = ['source10',
              'source1',
              'source2',
              'source3',
              'source4',
              'source5',
              'source6',
              'source7',
              'source8',
              'source9']



DefProjGammaList = { 'GMA4' : AllGammaSet ,'GMA3' : AllGammaSet }
ReadProjDerList = {'GMA4' :['g4D4','g1D1','g2D2','g3D3']}
DefProjDerList = {'GMA4' :['giDi']}
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
DefCombGammaList,dump = UDIndex(DefGammaList)
DefRedGammaList = DefCombGammaList[len(DefGammaList):]

# DeftoList = [17,18,19,20,21,22,23]
DeftoList = [16,17,18,19,20]
DefdtList = [1,2,3,4]
# DeftodtPicked = (18,2)
DeftodtPicked = (18,2)
DefTvarPicked = 'CMto'+str(DeftodtPicked[0])+'dt'+str(DeftodtPicked[1])

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
    AnatodtList = [DeftodtPicked]
    AnaTvarList = [DefTvarPicked]
else:
    AnatodtList = DeftodtList
    AnaTvarList = DefTvarList
    
# DefSmearList = ['32','64','128']
# DefSmearList = ['8','16','32','64','128','256']
# DefSmearList = ['32','64','128']
DefSmearList = ['32','64','128']
SingSmearList = ['32']
# DefSmearList = ['32']
# DefSmearList = ['8','16','32']
StateSet = map(str,range(1,(PoFShifts+1)*len(DefSmearList)+1))
PickedState = 1
PickedStateStr = 'state'+str(PickedState)
DefInterpList = ['nucleon']
# DefInterpList = ['nucleon','nucleon2']
DefSmList = ['sm'+ism for ism in DefSmearList]
DefInterpSmearList = ElongateName(DefInterpList,DefSmList)
CMTSinkList = [29]
AllTSinkList = [26,29,32,35,38]
AllTSinkListVar = [26,27,29,32,35,38]
AllTSinkShift = [it-tsource for it in AllTSinkList]
AllTSinkStrList = ['tsink'+str(its) for its in AllTSinkList]
AllTSinkStrListVar = ['tsink'+str(its) for its in AllTSinkListVar]

# AllREvecTSinkList = {'12104':[29],'12090':[26,32]}
AllREvecTSinkList = {'12104':[29],'12090':[]}
REvecTSinkList = AllREvecTSinkList[str(kappa)]
REvecTSinkStrList = ['tsink'+str(its) for its in REvecTSinkList]
DefREvecVarList = [18,2]
REvecTvarList = ['REvecto'+str(DefREvecVarList[0])+'dt'+str(DefREvecVarList[1])]
# REvecTvarList = []
DefPoFVarList = [16,2]
AllPoFTSinkList = {'12104':[],'12090':[26,27]}
PoFTSinkList = AllPoFTSinkList[str(kappa)]
PoFTSinkStrList = ['tsink'+str(its) for its in PoFTSinkList]
##DEBUG##
# PoFTvarList = ['PoF'+PoFShifts+'Test']
# PoFReadTvarList = ['PoFTestnD'+str(PoFShifts)]
# ##uncomment below to restore##
## warning, below has been hard coded to fix the to16 is to17 problem i had
PoFTvarList = ['PoF'+str(PoFShifts)+'to'+str(DefPoFVarList[0])+'dt'+str(DefPoFVarList[1])]
PoFDirTvarList = ['PoFto'+str(DefPoFVarList[0]-1)+'dt'+str(DefPoFVarList[1])]
PoFReadTvarList = ['PoFto'+str(DefPoFVarList[0]-1)+'dt'+str(DefPoFVarList[1])+'nD'+str(PoFShifts)]

# REvecPar26 = [[ 0.0004799, -0.0119381, 0.9999286 ],[0,0,0],[0,0,0]]
# REvecPar32 = [[ 0.0007567612, -0.0182265391, 0.9998335964],[0,0,0],[0,0,0]]


if kappa == 12090:
    TSFFileFlags = ['CM','Tsink','test32','Small']
    OSFFileFlags = ['CM','Tsink']
elif kappa == 12104:
    TSFFileFlags = ['REvec']
    OSFFileFlags = ['REvec']
MethodList = ['RF','Fits','SumMeth']+['TSF'+iTSF for iTSF in TSFFileFlags] + ['OSF'+iOSF for iOSF in OSFFileFlags]

##DEAULT TAKEN FROM OSF from: /raid/jdragos/data/results/ReadSetk12090/cfuns/twopt/OSFCM/state1sum22to23dt4twoptm0.txt
##date/time calculated: 00:29 12 Oct 2015
DefMassVal = {'22-31' : {'Avg':0.4662535526,'Std': 0.0028592256 ,'Chi': 0.4851841332}}

def GetRenorm(thisstring):
    if re.search('g[1234]g5',thisstring):
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
        if not DoMulticore: print string
    else:
        if not DoMulticore: print ' '.join(map(str,list(string)))
        

    

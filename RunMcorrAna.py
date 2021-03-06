#!/usr/bin/env python

from array import array
import os
import numpy as np
import sys
from Params import *
from SetLists import *
from Try2ptPickCMS import CreateTwoPt
from Try3ptPickCMS import CreateRF
from MiscFuns import touch
from OppFuns import Wipe2pt,WipeSet
from MiscFuns import *
from MomParams import *
from FFParams import *
from FFFuns import *
from ReadTxt import Get3ptSetMoms,Get2ptSetMoms
from InputArgs import InputParams
import copy
from multiprocessing import Pool

def CreateSubOppList(igamma):
    GammaList = []
    for iProj in DefProjGammaList.keys():
        for iDS in DefDSList:
            for icmplx in ['cmplx','']:
                GammaList.append(iDS+'P'+iProj[-1]+igamma+icmplx)
    return GammaList

if kappa == 12090:
    # ReadColList = [['TSink',"26 32 35 38"],['CM',"29"],['REvec',"26 32"],['PoF',"26 27"]]
    ReadColList = [['TSink',"26 32 35 38"],['CM',"29"],['PoF',"26 27"]]
    # ReadColList = [['CM',"29"],['REvec',"'26 32'"]]
else:
    ReadColList = [['REvec',"29"]]


def CRFDWrap(RunType,itsinkList,thisSmearList,iPrefList,thisDPList):
    CreateRF(RunType,itsinkList,thisSmearList,iPrefList,
             [iqTOip(0)],thisPDList=thisDPList,giDi=True)
    

def CRFWrap(RunType,itsinkList,thisSmearList,iPrefList,thisMomList,iProj,igamma):
    DRZ = False
    if 'q = 0 0 0' not in thisMomList:
        DRZ = True
        thisMomList = ['q = 0 0 0']+thisMomList
    CreateRF(RunType,itsinkList,thisSmearList,iPrefList,
             DragpZ([qstrTOip(iq) for iq in thisMomList]),thisPGList={iProj:[igamma]},DontWriteZero=DRZ)


def RunOffCorrs(thisPool,Curr,RunType,RunTSinkList=None,WipeThisSet=False,feedin=None):

    # print "running " + Curr + ' ' + thisCol + ' tsinks: ' + ' '.join(RunTSinkList)
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__

    if 'CM' in RunType:
        thisTSinkList = RunTSinkList
        thisSmearList = DefSmearList
        wipeSL = thisSmearList
        thisPrefList = ['cm' for s in thisTSinkList]
        thisTvarList = AnaTvarList
        thisREvecTvarList = []
        thisPoFTvarList = []
        thisStateSet = CMStateSet
    elif 'REvec' in RunType:
        thisTSinkList = RunTSinkList
        thisSmearList = DefSmearList
        wipeSL = []
        thisPrefList = ['REvec' for s in thisTSinkList]
        thisTvarList = []
        thisREvecTvarList = REvecTvarList
        thisPoFTvarList = []
        thisStateSet = [PickedState]
    elif 'PoF' in RunType:
        thisTSinkList = RunTSinkList
        thisSmearList = DefSmearList
        wipeSL = []
        thisPrefList = ['PoF' for s in thisTSinkList]
        thisTvarList = []
        thisREvecTvarList = []
        thisPoFTvarList = PoFTvarList
        thisStateSet = [PickedState]
    elif 'TSink' in RunType:
        thisTSinkList = RunTSinkList
        thisSmearList = ['32']
        wipeSL = thisSmearList
        thisPrefList = []
        for itsink in thisTSinkList:
            if itsink == '29': thisPrefList.append('cm')
            else: thisPrefList.append('tsink')
        thisTvarList = []
        thisREvecTvarList = []
        thisPoFTvarList = []
        thisStateSet = [PickedState]
    elif 'TwoPt' in RunType:
        DoTwoPt = True
    else:
        raise IOError('Choose CM , REvec , TSink or TwoPt along with Tsinks')

    print ''
    print '----------------------------------------------------------------------------------'
    if RunType == 'TwoPt':
        print 'Two Point Analysis'
        Wipe2pt(outputdir,tvarlist=TwoPtDefTvarList,smlist=DefSmearList)
        thisMomList = Get2ptSetMoms(outputdir,RunMomList,tvarlist=TwoTotDefTvarList,smlist=DefSmearList,tsrclist=PoFtsourceList)
        CreateTwoPt(DragpZ([qstrTOip(iq) for iq in thisMomList]),DefSmearList,feedin=feedin)
        print 'Two Point Analysis Complete'
    else:
        print 'Three Point Analysis '+Curr + ' ' + RunType + ' tsinks: ' + ' '.join(RunTSinkList)
        for it,itsink in enumerate(thisTSinkList):
            if RunType == 'PoF':
                itsinkList = range(int(itsink),int(itsink)+PoFShifts+1)
                iPrefList = [thisPrefList[it]]*len(itsinkList)
            else:
                itsinkList,iPrefList = [int(itsink)],[thisPrefList[it]]
            thisSetList,this2ptSetList,dump = CreateSet(thisSmearL=wipeSL,thisCMTSinkL=[str(itsink)]
                                                        ,thisStateL=thisStateSet,thisSingSmearL=wipeSL,
                                                        thisTvarL=thisTvarList,thisTSinkL=[str(itsink)],
                                                        thisREvecTvarL=thisREvecTvarList,thisREvecTSinkL=[str(itsink)],
                                                        thisPoFTvarL=thisPoFTvarList,thisPoFTSinkL=[str(itsink)])
            if 'giDi' in Curr:
                if WipeThisSet:
                    WipeSet(outputdir,['doubP4giDi','singP4giDi'],thisSetList)
                    WipeSet(outputdir+'cfuns/',['doubP4giDi','singP4giDi'],thisSetList)
                if thisPool == False:
                    CRFDWrap(RunType,itsinkList,thisSmearList,iPrefList,DictCurrOpps[Curr])
                else:
                    thisPool.apply_async(CRFDWrap,(RunType,itsinkList,thisSmearList,iPrefList,DictCurrOpps[Curr]))
                # print 'Tsink = ' + str(itsink)+ ' giDi Run:'
            else:
                for iProj,iGammaList in DictCurrOpps[Curr].iteritems():
                    for icg,igamma in enumerate(iGammaList):
                        PiOpp = 'P'+iProj[-1]+igamma
                        gammalist = ['doub'+PiOpp,'sing'+PiOpp]
                        gammalistcmplx = ['doub'+PiOpp+'cmplx','sing'+PiOpp+'cmplx']
                        if WipeThisSet:
                            wipegammalist = gammalist + gammalistcmplx
                            WipeSet(outputdir,wipegammalist,thisSetList)
                            WipeSet(outputdir+'cfuns/',wipegammalist,thisSetList)
                        MomDone = Get3ptSetMoms(outputdir,gammalist,RunMomList,thisSetList) 
                        MomDoneCmplx = Get3ptSetMoms(outputdir,gammalistcmplx,RunMomList,thisSetList) 
                        runmomlist,runmomlistcmplx = [],[]
                        for iq in RunMomList:
                            iqvec = np.array(qstrTOqvec(iq))*qunit
                            if Curr != 'Test':
                                dump,rcheck,ccheck = CurrFFs[Curr](PiOpp,iqvec.tolist(),[0,0,0],1.0)
                            else:
                                rcheck,ccheck = True,False
                            if rcheck:
                                if iq in MomDone:
                                    print 'Adding   tsink='+str(itsink)+ ' ' +PiOpp+' '+iq
                                    runmomlist.append(iq)
                                # else:
                                #     print 'Skipping tsink='+str(itsink)+ ' ' +PiOpp+' '+iq
                            if ccheck:
                                if iq in MomDoneCmplx:
                                    print 'Adding   tsink='+str(itsink)+ ' ' +PiOpp+'cmplx '+iq
                                    runmomlistcmplx.append(iq)
                                # else:
                                #     print 'Skipping tsink='+str(itsink)+ ' ' +PiOpp+'cmplx '+iq
                        sys.stdout = sys.__stdout__
                        sys.stderr = sys.__stderr__
                        if len(runmomlist) > 0:
                            if thisPool == False:
                                CRFWrap(RunType,itsinkList,thisSmearList,iPrefList,copy.deepcopy(runmomlist),iProj,igamma)
                            else:
                                thisPool.apply_async(CRFWrap,(RunType,itsinkList,thisSmearList,iPrefList,copy.copy(runmomlist),iProj,igamma))
                        if len(runmomlistcmplx) > 0:
                            if thisPool == False:
                                CRFWrap(RunType,itsinkList,thisSmearList,iPrefList,copy.deepcopy(runmomlistcmplx),iProj,igamma+'cmplx')
                            else:
                                thisPool.apply_async(CRFWrap,(RunType,itsinkList,thisSmearList,iPrefList,copy.copy(runmomlistcmplx),iProj,igamma+'cmplx'))
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__
    if 'TwoPt' not in RunType: print 'Three Point Analysis '+Curr + ' ' + RunType + ' tsinks: ' + ' '.join(RunTSinkList) + ' Added to jobs'

if len(sys.argv) < 2: raise IOError("input current type as first argument")
print sys.argv[1]
CurrIn = sys.argv[1]

with open( logdir+'LogAll.log.start','a') as f:
    f.write('\n')
with open( logdir+'LogAll.log.end','a') as f:
    f.write('\n')
thisPool = False
if CurrIn == 'TwoPt':
    feedin = InputParams(sys.argv[2:])
    RunOffCorrs(False,CurrIn,CurrIn,WipeThisSet=DefWipe,feedin=feedin)
else:
    if len(sys.argv) < 3: raise IOError("input Collection of Data To compute as second argument (CM,TSink,REvec)")
    thisColIn = sys.argv[2]
    if thisColIn != 'All':
        if len(sys.argv) < 4: raise IOError("input tsinks in third input parameter (CM,TSink,REvec)")
        thisTSinkIn = sys.argv[3].split()
        feedin = InputParams(sys.argv[4:])
    else:
        feedin = InputParams(sys.argv[3:])
    if DoMulticore and feedin['anaproc'] > 1:
        thisPool = Pool(processes=feedin['anaproc'])
    else:
        thisPool = False
    if CurrIn == 'All':
        # RunOffCorrs(thisPool,'TwoPt','TwoPt')
        for iCurr in AllCurrTypes:
            # if iCurr != 'Tensor': continue
            if thisColIn == 'All':
                for iCol in ReadColList:
                    RunOffCorrs(thisPool,iCurr,iCol[0],RunTSinkList=iCol[1].split(),WipeThisSet=DefWipe)
            else:
                RunOffCorrs(thisPool,iCurr,thisColIn,RunTSinkList=thisTSinkIn,WipeThisSet=DefWipe)
    else:
        if thisColIn == 'All':
            for iCol in ReadColList:
                RunOffCorrs(thisPool,CurrIn,iCol[0],RunTSinkList=iCol[1].split(),WipeThisSet=DefWipe)
        else:
            RunOffCorrs(thisPool,CurrIn,thisColIn,RunTSinkList=thisTSinkIn,WipeThisSet=DefWipe)
        
if thisPool != False:
    thisPool.close()
    thisPool.join()
print 'All Complete'

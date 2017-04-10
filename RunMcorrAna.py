#!/usr/bin/env python

from array import array
import os
import numpy as np
import sys
from Params import *
from SetLists import *
from Try2ptPickCMS import CreateTwoPt,CreateTwoPtOverDet
from TryTopCharge import CreateTwoPtTop
from TryTop3pt import CreateRFTop
from Try3ptPickCMS import CreateRF
from MiscFuns import touch
from OppFuns import Wipe2pt,WipeSet
from MiscFuns import *
from MomParams import *
from FFParams import *
from FFFuns import *
from ReadTxt import Get3ptSetMoms,Get2ptSetMoms,GetTopSetMoms
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

# if kappa == 12090:
#     # ReadColList = [['TSink',"26 32 35 38"],['CM',"29"],['REvec',"26 32"],['PoF',"26 27"]]
#     ReadColList = [['TSink',"26 32 35 38"],['CM',"29"],['PoF',"26 27"]]
#     # ReadColList = [['CM',"29"],['REvec',"'26 32'"]]
# else:
#     ReadColList = [['REvec',"29"]]

ReadColList = [['TSink',"6"]]


def CRFDWrap(RunType,itsinkList,thisSmearList,iPrefList,thisDPList):
    CreateRF(RunType,itsinkList,thisSmearList,iPrefList,
             [iqTOip(0)],thisPDList=thisDPList,giDi=True)
    

def CRFWrap(RunType,itsinkList,thisiSmearList,thisjSmearList,iPrefList,thisMomList,iProj,igamma,DoTop):
    DRZ = False
    if 'q = 0 0 0' not in thisMomList:
        DRZ = True
        thisMomList = ['q = 0 0 0']+thisMomList
    if 'Wein' == DoTop:
        CreateRFTop(RunType,itsinkList,thisiSmearList,thisjSmearList,iPrefList,
                    DragpZ([qstrTOip(iq) for iq in thisMomList]),thisPGList={iProj:[igamma]},DontWriteZero=DRZ,Wein=True)
    elif DoTop == True:
        CreateRFTop(RunType,itsinkList,thisiSmearList,thisjSmearList,iPrefList,
                    DragpZ([qstrTOip(iq) for iq in thisMomList]),thisPGList={iProj:[igamma]},DontWriteZero=DRZ)
    else:
        CreateRF(RunType,itsinkList,thisiSmearList,thisjSmearList,iPrefList,
                 DragpZ([qstrTOip(iq) for iq in thisMomList]),thisPGList={iProj:[igamma]},DontWriteZero=DRZ)
        
def RunOffCorrs(thisPool,Curr,RunType,RunTSinkList=None,WipeThisSet=False,feedin=None,DoTop=False):

    # print "running " + Curr + ' ' + thisCol + ' tsinks: ' + ' '.join(RunTSinkList)
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__

    if 'CM' in RunType:
        thisTSinkList = RunTSinkList
        thisiSmearList = DefiSmearList
        thisjSmearList = DefjSmearList
        wipeiSL = thisiSmearList
        wipejSL = thisjSmearList
        wipeSSL = []
        thisPrefList = ['cm' for s in thisTSinkList]
        thisTvarList = AnaTvarList
        thisREvecTvarList = []
        thisPoFTvarList = []
        thisStateSet = CMStateSet
    elif 'REvec' in RunType:
        thisTSinkList = RunTSinkList
        thisiSmearList = DefiSmearList
        thisjSmearList = DefjSmearList
        wipeiSL = []
        wipejSL = []
        wipeSSL = []
        thisPrefList = ['REvec' for s in thisTSinkList]
        thisTvarList = []
        thisREvecTvarList = REvecTvarList
        thisPoFTvarList = []
        thisStateSet = [PickedState]
    elif 'PoF' in RunType:
        thisTSinkList = RunTSinkList
        thisiSmearList = DefiSmearList
        thisjSmearList = DefjSmearList
        wipeiSL = []
        wipejSL = []
        wipeSSL = []
        thisPrefList = ['PoF' for s in thisTSinkList]
        thisTvarList = []
        thisREvecTvarList = []
        thisPoFTvarList = PoFTvarList
        thisStateSet = [PickedState]
    elif 'TSink' in RunType:
        thisTSinkList = RunTSinkList
        thisSmearList = ['32']
        wipeiSL = []
        wipejSL = []
        wipeSSL = SingSmearList
        thisPrefList = []
        for itsink in thisTSinkList:
            if itsink == '29': thisPrefList.append('cm')
            else: thisPrefList.append('tsink')
        thisTvarList = []
        thisREvecTvarList = []
        thisPoFTvarList = []
        thisStateSet = [PickedState]
    elif 'TwoPt' in RunType or 'TopAlpha' in RunType or 'WeinAlpha' in RunType:
        DoTwo = True
    else:
        raise IOError('Choose CM , REvec , TSink or TwoPt along with Tsinks')

    print ''
    print '----------------------------------------------------------------------------------'
    if RunType == 'TwoPt':
        print 'Two Point Analysis'
        if DefWipe:
            thisMomList = RunAvgMomList
        else:
            thisMomList = Get2ptSetMoms(outputdir[0],RunAvgMomList,tvarlist=PoFTvarList,ismlist=DefiSmearList,jsmlist=DefjSmearList,tsrclist=PoFtsourceList)
        if len(thisMomList) == 0 :
            print 'Already Done, (if you want to overwrite, set DefWipe to True and rerun)'
        else:
            if OverDetRun:
                CreateTwoPtOverDet([qstrTOip(imom) for imom in DragpZstr(thisMomList)],DefiSmearList,DefjSmearList,feedin=feedin)
            else:
                CreateTwoPt([qstrTOip(imom) for imom in DragpZstr(thisMomList)],DefiSmearList,DefjSmearList,feedin=feedin)
            print 'Two Point Analysis Complete'
    elif RunType == 'TopAlpha':
        print 'Topological Charge Analysis'
        # Wipe2pt(outputdir[0],tvarlist=TwoPtDefTvarList,smlist=DefSmearList)
        # thisMomList = Get2ptSetMoms(outputdir[0],RunAvgMomList,tvarlist=TwoTotDefTvarList,smlist=DefSmearList,tsrclist=PoFtsourceList)
        if DefWipe:
            thisMomList = RunAvgMomList
        else:
            thisMomList = GetTopSetMoms(outputdir[0]+'/Top/',RunAvgMomList,tvarlist=PoFTvarList,ismlist=DefiSmearList,jsmlist=DefjSmearList,tsrclist=PoFtsourceList)
        if len(thisMomList) == 0 :
            print 'Already Done, (if you want to overwrite, set DefWipe to True and rerun)'
        else:
            CreateTwoPtTop([qstrTOip(imom) for imom in DragpZstr(thisMomList)],DefiSmearList,DefjSmearList,feedin=feedin)
            print 'Two Point Analysis Complete'
    elif RunType == 'WeinAlpha':
        print 'Weinberg Operator Analysis'
        # Wipe2pt(outputdir[0],tvarlist=TwoPtDefTvarList,smlist=DefSmearList)
        # thisMomList = Get2ptSetMoms(outputdir[0],RunAvgMomList,tvarlist=TwoTotDefTvarList,smlist=DefSmearList,tsrclist=PoFtsourceList)
        if DefWipe:
            thisMomList = RunAvgMomList
        else:
            thisMomList = GetTopSetMoms(outputdir[0]+'/Wein/',RunAvgMomList,tvarlist=PoFTvarList,ismlist=DefiSmearList,jsmlist=DefjSmearList,tsrclist=PoFtsourceList)
        if len(thisMomList) == 0 :
            print 'Already Done, (if you want to overwrite, set DefWipe to True and rerun)'
        else:
            CreateTwoPtTop([qstrTOip(imom) for imom in DragpZstr(thisMomList)],DefiSmearList,DefjSmearList,feedin=feedin,Wein=True)
            print 'Two Point Analysis Complete'
    else:
        print 'Three Point Analysis '+Curr + ' ' + RunType + ' tsinks: ' + ' '.join(RunTSinkList)
        for it,itsink in enumerate(thisTSinkList):
            if RunType == 'PoF':
                itsinkList = range(int(itsink),int(itsink)+PoFShifts+1)
                iPrefList = [thisPrefList[it]]*len(itsinkList)
            else:
                itsinkList,iPrefList = [int(itsink)],[thisPrefList[it]]
            thisSetList,this2ptSetList,dump = CreateSet(thisiSmearL=wipeiSL,thisjSmearL=wipejSL,thisSingSmearL=wipeSSL,thisCMTSinkL=[str(itsink)]
                                                        ,thisStateL=thisStateSet,
                                                        thisTvarL=thisTvarList,thisTSinkL=[str(itsink)],
                                                        thisREvecTvarL=thisREvecTvarList,thisREvecTSinkL=[str(itsink)],
                                                        thisPoFTvarL=thisPoFTvarList,thisPoFTSinkL=[str(itsink)])
            if Debug:
                print
                print 'SetList:'
                for iset in thisSetList:
                    print iset
                print
            if 'giDi' in Curr:
                if DoTop: print 'Warning, giDi not implemented with topological charge'
                if WipeThisSet:
                    WipeSet(outputdir[0],['doubP4giDi','singP4giDi'],thisSetList)
                    WipeSet(outputdir[0]+'cfun/',['doubP4giDi','singP4giDi'],thisSetList)
                if thisPool == False:
                    CRFDWrap(RunType,itsinkList,thisSmearList,iPrefList,DictCurrOpps[Curr])
                else:
                    thisPool.apply_async(CRFDWrap,(RunType,itsinkList,thisSmearList,iPrefList,DictCurrOpps[Curr]))
                # print 'Tsink = ' + str(itsink)+ ' giDi Run:'
            else:
                for iProj,iGammaList in DictCurrOpps[Curr].iteritems():
                    for icg,igamma in enumerate(iGammaList):
                        sys.stdout = sys.__stdout__
                        sys.stderr = sys.__stderr__
                        PiOpp = 'P'+iProj[-1]+igamma
                        if DoTop == 'Wein':
                            PiOpp += 'Wein'
                        elif DoTop :
                            PiOpp += 'Top'
                        gammalist = ['doub'+PiOpp,'sing'+PiOpp]
                        gammalistcmplx = ['doub'+PiOpp+'cmplx','sing'+PiOpp+'cmplx']
                        if WipeThisSet:
                            wipegammalist = gammalist + gammalistcmplx
                            WipeSet(outputdir[0],wipegammalist,thisSetList)
                            WipeSet(outputdir[0]+'cfun/',wipegammalist,thisSetList)    
                        MomDone = Get3ptSetMoms(outputdir[0],gammalist,RunMomList,thisSetList) 
                        MomDoneCmplx = Get3ptSetMoms(outputdir[0],gammalistcmplx,RunMomList,thisSetList) 
                        runmomlist,runmomlistcmplx = [],[]
                        for iq in RunMomList:
                            iqvec = np.array(qstrTOqvec(iq))*qunit
                            if Curr != 'Test':
                                checkPiOpp = PiOpp
                                if DoTop=='Wein': checkPiOpp = PiOpp+'Wein'
                                elif DoTop: checkPiOpp = PiOpp+'Top'
                                dump,rcheck,ccheck = CurrFFs[Curr](checkPiOpp,iqvec.tolist(),[0,0,0],1.0)
                            else:
                                rcheck,ccheck = True,False
                            if rcheck and ccheck:
                                print 'Warning, FormFactor result should be real or complex result, not both'
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
                                CRFWrap(RunType,itsinkList,thisiSmearList,thisjSmearList,iPrefList,copy.deepcopy(runmomlist),iProj,igamma,DoTop)
                            else:
                                thisPool.apply_async(CRFWrap,(RunType,itsinkList,thisiSmearList,thisjSmearList,iPrefList,copy.copy(runmomlist),iProj,igamma,DoTop))
                        if len(runmomlistcmplx) > 0:
                            if thisPool == False:
                                CRFWrap(RunType,itsinkList,thisiSmearList,thisjSmearList,iPrefList,copy.deepcopy(runmomlistcmplx),iProj,igamma+'cmplx',DoTop)
                            else:
                                thisPool.apply_async(CRFWrap,(RunType,itsinkList,thisiSmearList,thisjSmearList,iPrefList,copy.copy(runmomlistcmplx),iProj,igamma+'cmplx',DoTop))
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__
    if 'TwoPt' not in RunType and 'Alpha' not in RunType:
        print 'Three Point Analysis '+Curr + ' ' + RunType + ' tsinks: ' + ' '.join(RunTSinkList) + ' Added to jobs'

if len(sys.argv) < 2: raise IOError("input current type as first argument")
print sys.argv[1]
CurrIn = sys.argv[1]

with open( logdir+'LogAll.log.start','a') as f:
    f.write('\n')
with open( logdir+'LogAll.log.end','a') as f:
    f.write('\n')
thisPool = False
if CurrIn == 'TwoPt' or CurrIn == 'TopAlpha' or CurrIn == 'WeinAlpha':
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
            DoTop=False
            if 'Top' in iCurr: DoTop=True
            if 'Wein' in iCurr: DoTop='Wein'
            # if iCurr != 'Tensor': continue
            if thisColIn == 'All':
                for iCol in ReadColList:
                    RunOffCorrs(thisPool,iCurr,iCol[0],RunTSinkList=iCol[1].split(),WipeThisSet=DefWipe,DoTop=DoTop)
            else:
                RunOffCorrs(thisPool,iCurr,thisColIn,RunTSinkList=thisTSinkIn,WipeThisSet=DefWipe,DoTop=DoTop)
    else:
        DoTop=False
        if 'Top' in CurrIn: DoTop=True
        if 'Wein' in CurrIn: DoTop='Wein'
        if thisColIn == 'All':
            for iCol in ReadColList:
                RunOffCorrs(thisPool,CurrIn,iCol[0],RunTSinkList=iCol[1].split(),WipeThisSet=DefWipe,DoTop=DoTop)
        else:
            RunOffCorrs(thisPool,CurrIn,thisColIn,RunTSinkList=thisTSinkIn,WipeThisSet=DefWipe,DoTop=DoTop)
        
if thisPool != False:
    thisPool.close()
    thisPool.join()
print 'All Complete'

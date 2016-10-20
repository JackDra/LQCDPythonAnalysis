#!/usr/bin/env python

from Params import *
from FitParams import *
from collections import OrderedDict
from CombParams import CombList,CombFFList
from FFParams import NoFFList
import re

def CutDupSet(SetList):
    return list(OrderedDict.fromkeys([iSL[7:] for iSL in SetList]))



##data = [ igamma , ip , iset , it ]
#dataout = [ itsink , istate , igamma , ip , it]

def FlagList(AllSetList,*flag):
    SLOut = []
    for iset in AllSetList:
        if all([str(iflag) in iset for iflag in flag]):
            SLOut.append(iset)
    return SLOut

def PickCM(thisSetList):
    thisCM = None
    while thisCM not in DefTvarPicked:
        thisCM = raw_input("please choose a CM out of: " + ','.join(DefTvarPicked)+'  \n')
    SetListOut = []
    for iSet in thisSetList:
        if any([ism in iSet for ism in DefSmList]) or thisCM in iSet:
            SetListOut.append(iSet)
    return SetListOut
        
def SortSum(SLin):
    SLout = []
    SLpre = []
    for iSL in SLin:
        if 'SumMeth' not in iSL:
            SLpre.append(iSL)
    for ifitr in SumFitRListFlags:
        for icut in SumCutList:
            for iSLin in SLin:
                if 'SumMeth' in iSLin and icut in iSLin and ifitr in iSLin:
                    SLout.append(iSLin)
    return SLpre+SLout

def SortMySet(SLin,massset=False):
    # SLin.sort()
    SLout = []
    TSinkLout = []
    if massset:
        for ism in DefSmList:
            for itsrc in PoFtsourceList:
                for iSLin in SLin:
                    if (itsrc in iSLin and ism in iSLin) and iSLin not in SLout:
                        SLout.append(iSLin)
        for ism in TwoTotDefTvarList:
            for iSLin in SLin:
                if ism in iSLin and iSLin not in SLout:
                    SLout.append(iSLin)
    else:
        for ism in DefSmList+PoFTvarList+REvecTvarList+DefTvarList:
            for iSLin in SLin:
                if ism in iSLin and 'tsink' not in iSLin and iSLin not in SLout:
                    SLout.append(iSLin)
        for itsink,itint in zip(AllTSinkStrList,AllTSinkList):
            for iSLin in SLin:
                if itsink in iSLin and 'sm' not in iSLin and 'to' not in iSLin and iSLin not in SLout:
                    SLout.append(iSLin)
                    TSinkLout.append(itint)
        for itsink,itint in zip(map(str,range(20,40)),range(20,40)):
        # for itsink,itint in zip(AllTSinkStrList,AllTSinkList):
            # for ism in DefSmList + DefTvarList:
            for ism in DefSmList+DefTvarList+PoFTvarList+REvecTvarList:
                for iSLin in SLin:
                    if 'tsink'+itsink in iSLin and ism in iSLin and iSLin not in SLout:
                        TSinkLout.append(itint)
                        SLout.append(iSLin)
            # for ism in DefTvarList:
            #     for iSLin in SLin:
            #         if 'tsink'+itsink in iSLin and ism in iSLin and 'PoF' not in iSLin and 'REvec' not in iSLin and iSLin not in SLout:
            #             TSinkLout.append(itint)
            #             SLout.append(iSLin)
        # for itsink,itint in zip(REvecTSinkStrList,REvecTSinkList):
        #     for ism in REvecTvarList:
        #         for iSLin in SLin:
        #             if itsink in iSLin and ism in iSLin and VarPref not in iSLin:
        #                 TSinkLout.append(itint)
        #                 SLout.append(iSLin)
        # for itsink,itint in zip(PoFTSinkStrList,PoFTSinkList):
        #     for ism in PoFTvarList:
        #         for iSLin in SLin:
        #             if itsink in iSLin and ism in iSLin and 'PoF' in iSLin:
        #                 TSinkLout.append(itint)
        #                 SLout.append(iSLin)
                    
    if any(['SumMeth' in iSL for iSL in SLin]): SLout = SortSum(SLout)
    return SLout,TSinkLout

            
def SplitSet(data,SetList,TsinkList):
    dataout = []
    for itsink , thetsink in enumerate(TsinkList):
        dataout.append([])
        iset = -1
        for isetdatain,theset in enumerate(SetList):
            if 'tsink'+str(thetsink) in theset:
                iset += 1
                dataout[itsink].append([])
                for igamma,gammadata in enumerate(data):
                    dataout[itsink][iset].append([])
                    for pdata in gammadata:
                        dataout[itsink][iset][igamma].append(pdata[isetdatain])
    return dataout


def RemoveToDt(thisstring):
    for ivar in TwoPtDefTvarList:
        thisstring = thisstring.replace(ivar,'')
    return thisstring


def CreateStateSet(smL,stateL,tvarL):
    Lout = ['sm'+str(ism) for ism in smL]
    for istate in stateL:
        Lout += ['state'+str(istate)+itvar for itvar in tvarL]
    return Lout

def CreateTSinkStateSet(thisTsink,smL,stateL,tvarL):
    return ['tsink'+str(thisTsink)+istate for istate in CreateStateSet(smL,stateL,tvarL)]

def CreateStateTsinkSet(thisState,tsinkL):
    return ['tsink'+str(itsink)+thisState for itsink in tsinkL]


def CreateREvecSet(TSinkList,thisStateList,TvarList):
    SetGraph = []
    SetTsink = []
    for itsink in TSinkList:
        for istate in thisStateList:
            for itvar in TvarList:
                SetGraph.append('tsink'+str(itsink)+'state'+str(istate)+itvar)
                SetTsink.append(int(itsink))
    return SetGraph,SetTsink

def CreateGenericSet(thisTSinkList,thisSmearList,thisStateList,thisTvarList):
    output = []
    for itsink in thisTSinkList:
        output += CreateTSinkStateSet(itsink,thisSmearList,[],[])
        output += CreateTSinkStateSet(itsink,[],thisStateList,thisTvarList)
    return output

def CreateSet(thisSmearL=DefSmearList,thisSingSmearL=SingSmearList,
              thisStateL=[str(PickedState)],thisTvarL=AnaTvarList,
              thisTSinkL=AllTSinkList,thisCMTSinkL=CMTSinkList,
              thisREvecTvarL=REvecTvarList,thisREvecTSinkL=REvecTSinkList,
              thisPoFTvarL=PoFTvarList,thisPoFTSinkL=PoFTSinkList):
    SetGraph,SetMassGraph = [],[]
    SetTsink = []
    SetMassGraph += CreateMassSet(thisSmearL,thisStateL,thisTvarL)
    SetMassGraph += CreateMassSet(thisSingSmearL,[],[])
    SetGraph += CreateGenericSet(thisCMTSinkL,thisSmearL,[],[])
    SetGraph += CreateGenericSet(thisTSinkL,thisSingSmearL,[],[])
    SetGraph += CreateGenericSet(thisCMTSinkL,[],thisStateL,thisTvarL)        
    SetGraph += CreateGenericSet(thisREvecTSinkL,[],thisStateL,thisREvecTvarL)        
    SetGraph += CreateGenericSet(thisPoFTSinkL,[],thisStateL,thisPoFTvarL)
    SetGraph,SetTsink = SortMySet(SetGraph,massset=False)
    return [SetGraph,SortMySet(SetMassGraph,massset=True)[0],SetTsink]



def CreateMassSet(thisSmearL,thisStateList,thisTvarList,flipord=False,tsrc=False):
    SetGraph = []
    for ismear in thisSmearL:
        SetGraph.append('sm'+ismear)
    if flipord:
        for itvar in thisTvarList:
            for istate in thisStateList:
                SetGraph.append('state'+str(istate)+itvar)
    else:
        for istate in thisStateList:
            for itvar in thisTvarList:
                SetGraph.append('state'+str(istate)+itvar)
    return SortMySet(SetGraph,massset=True)[0]


## dataCM [ istate , todt ]
## data [ ism ]
## dataout = [ (istatetodt / ismism) ]


def CreateDataSet(dataCM,data,thisSmearList,thisStateList,thisTvarList,Interps=['nucleon']):
    dataCM,data = np.array(dataCM),np.array(data)
    dataCMflat = dataCM.reshape((dataCM.shape[0]*dataCM.shape[1],)+dataCM.shape[2:])
    return [np.array(data.tolist()+dataCMflat.tolist()),CreateMassSet(thisSmearList,thisStateList,thisTvarList)]

def CreateDataTsinkSet(dataCM,data,thisSmearList,thisStateList,thisTvarList,tsink):
    dataCMflat= dataCM.reshape((dataCM.shape[0]*dataCM.shape[1],)+dataCM.shape[2:])
    thisSetList = CreateMassSet(thisSmearList,thisStateList,thisTvarList)
    data = np.array(data)
    return [np.array(data.tolist()+dataCMflat.tolist()),['tsink'+str(tsink)+iS for iS in thisSetList]]


def SplitTSinkString(string): 
    tsinkstring = re.search('tsink..',string)
    if tsinkstring != None: tsinkstring = tsinkstring.group()
    smstring = None
    for ism in DefSmList:
        if ism in string:
            smstring = ism
    if smstring == None:
        if 'state' in string:
            for istate in map(str,StateSet):
                if 'state'+istate in string:
                    smstring = 'state'+istate
        for itvar in PoFTvarList+REvecTvarList+DefTvarList:
            if itvar in string:
                smstring += itvar
                break
    return tsinkstring,smstring

def PickTSinkSet(SetList,Tsink):
    ListOut = []
    for iset in SetList:
        if 'tsink'+str(Tsink) in iset:
            ListOut.append(iset)
    return ListOut

def PickSmSet(SetList,Sm):
    ListOut = []
    for iset in SetList:
        if Sm in iset:
            ListOut.append(Sm)
    return ListOut

def ReduceTsink(listin,NoCM=False,NoREvec=False,NoPoF=False):
    listout = []
    for ilist in listin:
        tsinkstr,smstr = SplitTSinkString(ilist)
        if smstr not in listout:
            listout.append(smstr)
    if NoPoF:
        for i,ilo in enumerate(listout):
            if 'PoF' in ilo:del listout[i]
    if NoCM:
        for i,ilo in enumerate(listout):
            if 'CM' in ilo:del listout[i]
    if NoREvec:
        for i,ilo in enumerate(listout):
            if 'REvec' in ilo:del listout[i]
    return listout


def GetTsinkSmLists(listin,NoREvec=False,Reduced=True):
    smlist,tsinklist = [],[]
    for istr in listin:
        tsinkstr,smstr = SplitTSinkString(istr)
        if 'REvec' in smstr and NoREvec: continue
        if tsinkstr not in tsinklist or not Reduced:
            tsinklist.append(tsinkstr)
        if smstr not in smlist or not Reduced:
            smlist.append(smstr)
    return tsinklist,smlist

def RemoveTSink(string):
    return SplitTSinkString(string)[1]

def GetintTSink(string):
    return int(SplitTSinkString(string)[0].replace('tsink',''))

def RemoveSet(string):
    return SplitTSinkString(string)[0]

def ReduceTooMassSet(thisSetList):
    notsinkList = ReduceTsink(thisSetList)
    setout = []
    for iset in notsinkList:
        if 'sm' in iset and 'state' not in iset:
            setout.append('tsrc'+tsource+iset)
        else:
            setout.append(iset)
            # for ic,ism in enumerate(notsinkList):
    #     if 'REvec' in ism: del notsinkList[ic]
    return setout


def PickSetForMethod(thismethod,thisSetList):
    outSetList = thisSetList
    if 'Tsink' in thismethod or 'Small' in thismethod or 'test32' in thismethod:
        outSetList = []
        for itsink in AllTSinkStrList:
            if itsink+SingSmList[0] in thisSetList:
                outSetList.append(itsink+SingSmList[0])
    elif 'CM' in thismethod:
        outSetList = []
        for iset in thisSetList:
            for ikey in ['CM','PoF']:
                for itsink in TSinkStrDictList[ikey]:
                    for ism in SmearDictList[ikey]:
                        if itsink+ism in iset:
                            outSetList.append(iset)
    elif 'PoF' in thismethod:
        outSetList = []
        for iset in thisSetList:
            ikey = 'PoF'
            for itsink in TSinkStrDictList[ikey]:
                for ism in SmearDictList[ikey]:
                    if itsink+ism in iset:
                        outSetList.append(iset)
    elif 'REvec' in thismethod:
        ikey = 'REvec'
        for itsink in TSinkStrDictList[ikey]:
            for ism in SmearDictList[ikey]:
                if itsink+ism in iset:
                    outSetList.append(iset)
    if 'TSF' in thismethod:
        outSetList = ReduceTsink(outSetList,NoPoF=True)
    elif 'SumMeth' in thismethod:
        outSetList = SingSmList
    return outSetList
                




if kappa == 12090:
    DefSetCol = CreateSet()
elif kappa == 12104:
    DefSetCol = CreateSet(thisSmearL=[],thisTSinkL=[],thisTvarL=[])
    DefSetCol[1] = CreateSet(thisTSinkL=[29],thisTvarL=[])[1]

AllCMSetList = CreateGenericSet(CMTSinkList,[],[PickedState],DefTvarList)
    
DefSetList,DefMassSetList,DefTSinkSetList = DefSetCol

def CreateOSFfitKey(smear):
    smearindex,deltashift = RemoveToDt(smear),0
    if 'PoF' in smear:
        deltashift = PoFShifts*2
        smearindex = PickedStateStr+'PoF'+str(PoFShifts)
    elif 'REvec' in smear:
        deltashift = 0
        smearindex = PickedStateStr+'REvec'
    elif 'CM' in smear:
        for itvar in DefTvarPicked:
            if itvar in smear:
                deltashift = 0
                smearindex = PickedStateStr+itvar
    return smearindex,deltashift


def SplitDSCurr(thisstr):
    for iDS in DefDSList+CombList:
        hold = thisstr.replace(iDS,'',1)
        if iDS in thisstr and hold in NoFFList.keys():
            return iDS,hold,''
        for iFFComb in CombFFList:
            hold = thisstr.replace(iDS,'',1).replace(iFFComb,'',1)
            if iFFComb in thisstr and iDS in thisstr and hold in NoFFList.keys():
                return iDS,hold,iFFComb.replace('/','')
    for iFFComb in CombFFList:
        hold = thisstr.replace(iFFComb,'',1)
        if iFFComb in thisstr and hold in NoFFList.keys():
            return '',hold,iFFComb.replace('/','')
    return '',thisstr,''
        
def CreateCurrCombFn(thisstr,spacing=''):
    DS,Curr,FFComb = SplitDSCurr(thisstr)
    return spacing.join([FFComb,DS,Curr])

def SplitKappa(thisstr):
    thiskappa = re.search('k.....',thisstr)
    if thiskappa == None:
        return None,thisstr
    else:
        thiskappa = thiskappa.group()
        return thiskappa,thisstr.replace(thiskappa,'')

 #!/usr/bin/env python

from Params import *
from FitParams import *
from collections import OrderedDict
import re

def CutDupSet(SetList):
    return list(OrderedDict.fromkeys([iSL[7:] for iSL in SetList]))



##data = [ igamma , ip , iset , it ]
#dataout = [ itsink , istate , igamma , ip , it]

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
        for ism in DefSmList+TwoTotDefTvarList:
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
            for ism in DefSmList+PoFTvarList+REvecTvarList:
                for iSLin in SLin:
                    if 'tsink'+itsink in iSLin and ism in iSLin and iSLin not in SLout:
                        TSinkLout.append(itint)
                        SLout.append(iSLin)
            for ism in DefTvarList:
                for iSLin in SLin:
                    if 'tsink'+itsink in iSLin and ism in iSLin and 'PoF' not in iSLin and 'REvec' not in iSLin and iSLin not in SLout:
                        TSinkLout.append(itint)
                        SLout.append(iSLin)
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
    SetGraph,SetMassGraph = [],set([])
    SetTsink = []
    SetMassGraph = CreateMassSet(thisSmearL,thisStateL,thisTvarL)
    SetGraph += CreateGenericSet(thisCMTSinkL,thisSmearL,thisStateL,[])
    SetGraph += CreateGenericSet(thisTSinkL,thisSingSmearL,thisStateL,[])
    SetGraph += CreateGenericSet(thisCMTSinkL,[],thisStateL,thisTvarL)        
    SetGraph += CreateGenericSet(thisREvecTSinkL,[],thisStateL,thisREvecTvarL)        
    SetGraph += CreateGenericSet(thisPoFTSinkL,[],thisStateL,thisPoFTvarL)
    SetGraph,SetTsink = SortMySet(SetGraph,massset=False)
    return [SetGraph,SortMySet(list(SetMassGraph))[0],SetTsink]



def CreateMassSet(thisSmearL,thisStateList,thisTvarList,flipord=False):
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

def ReduceTsink(listin,NoCM=False):
    listout = []
    for ilist in listin:
        tsinkstr,smstr = SplitTSinkString(ilist)
        if smstr not in listout:
            listout.append(smstr)
    if NoCM:
        for i,ilo in enumerate(listout):
            if 'state' in ilo:del listout[i]
    return listout


def GetTsinkSmLists(listin):
    smlist,tsinklist = [],[]
    for istr in listin:
        tsinkstr,smstr = SplitTSinkString(istr)
        if tsinkstr not in tsinklist:
            tsinklist.append(tsinkstr)
        if smstr not in smlist:
            smlist.append(smstr)
    return tsinklist,smlist

def RemoveTSink(string):
    return SplitTSinkString(string)[1]

def GetintTSink(string):
    return int(SplitTSinkString(string)[0].replace('tsink',''))

def RemoveSet(string):
    return SplitTSinkString(string)[0]

if kappa == 12090:
    DefSetCol = CreateSet()
elif kappa == 12104:
    DefSetCol = CreateSet(thisSmearL=[],thisTSinkL=[],thisTvarL=[])
    DefSetCol[1] = CreateSet(thisTSinkL=[29],thisTvarL=[])[1]
    
DefSetList,DefMassSetList,DefTSinkSetList = DefSetCol

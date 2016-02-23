#!/usr/bin/env python

from SetLists import *
import glob,os

def ReadSetDir(thisdir,thisgamma=''):
    thisSetList = glob.glob(thisdir+'/*.txt')
    thisSetList = [(iS.replace(thisdir,'')
                    .replace(thisgamma,'')
                    .replace('.txt','')
                    .replace('.boot','')) for iS in thisSetList]
    if OnlySelVar:
        SetOut = []
        for iSet in thisSetList:
            if (VarPref not in iSet) or (DefTvarPicked in iSet):
                SetOut.append(iSet)
    else:
        SetOut = thisSetOut        
    return SetOut


def ReadSetBootDir(thisdir,thisgamma='',ReadBoot=True):
    tSL = ReadSetDir(thisdir,thisgamma=thisgamma)
    if ReadBoot: 
        tSLBoot = ReadSetDir(thisdir+'boots/',thisgamma=thisgamma)
        return list(set(tSL) & set(tSLBoot))
    else:
        return tSL

def ReadAllDir(thisdir,thisgamma='',ReadBoot=True):
    thisSetList = []
    for (dirname,dirs,files) in os.walk(thisdir):
        if 'boots' in dirname: continue
        thisSetList += [dirname.replace(thisdir,'').replace('/','')+iset for iset in ReadSetBootDir(dirname,thisgamma=thisgamma,ReadBoot=ReadBoot)]
    return thisSetList


## thisCurrDict = [ iCurr : iSets ]
def GetCurrDict(thisCurrTypes,ReadBoot=False):
    thisCurrDict = {}
    thisSetList = None
    for iCurr in thisCurrTypes:
        print 'Reading ' , iCurr , ' Form Factors '
        thisCurrDict[iCurr] = ReadAllDir(outputdir+'FormFactors/'+iCurr+'/',thisgamma=iCurr,ReadBoot=ReadBoot)
        if thisSetList == None: 
            thisSetList = set(thisCurrDict[iCurr])
        else:
            thisSetList = thisSetList.intersection(thisCurrDict[iCurr])
        # print 'Sets Found:\n','\n'.join(thisCurrDict[iCurr])
        # print ''

    thisSetList = list(thisSetList)
    thisSetList.sort()
    print 'Sets Found:\n','\n'.join(thisSetList)
    return thisCurrDict

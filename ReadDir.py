#!/usr/bin/env python

from SetLists import *
import glob,os

def ReadSetDir(thisdir,thisgamma=''):
    thisSetList = glob.glob(thisdir+'/*.xml')
    thisSetList = [(iS.replace(thisdir,'')
                    .replace(thisgamma,'')
                    .replace('.xml','')) for iS in thisSetList]
    if OnlySelVar:
        SetOut = []
        for iSet in thisSetList:
            if (any([ism in iSet for ism in DefSmList]) or (DefTvarPicked in iSet)
                or any([iPoF in iSet for iPoF in PoFTvarList]) or any([iREvec in iSet for iREvec in REvecTvarList])):
                SetOut.append(iSet)
    else:
        SetOut = thisSetOut        
    return SetOut


# def ReadSetBootDir(thisdir,thisgamma='',ReadBoot=True):
#     tSL = ReadSetDir(thisdir,thisgamma=thisgamma)
#     if ReadBoot: 
#         tSLBoot = ReadSetDir(thisdir+'boots/',thisgamma=thisgamma)
#         return list(set(tSL) & set(tSLBoot))
#     else:
#         return tSL

def ReadAllDir(thisdir,thisgamma=''):
    thisSetList = []
    for (dirname,dirs,files) in os.walk(thisdir):
        if 'boots' in dirname: continue
        thisSetList += [dirname.replace(thisdir,'').replace('/','')+iset for iset in ReadSetDir(dirname,thisgamma=thisgamma)]
    return thisSetList


## thisCurrDict = [ iCurr : iSets ]
def GetCurrDict(thisCurrTypes):
    thisCurrDict = {}
    thisSetList = None
    for iCurr in thisCurrTypes:
        print 'Reading ' , iCurr , ' Form Factors '
        thisCurrDict[iCurr] = ReadAllDir(outputdir+'FormFactors/'+iCurr+'/',thisgamma=iCurr.replace('/',''))
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

def CheckCurrentSets(thisCurrDict):
    outCurrDict = {}
    for icurr,currset in thisCurrDict.iteritems():
        if 'doub' in icurr:
            singcurr = icurr.replace('doub','sing')
            if singcurr not in thisCurrDict.keys():
                raise LookupError(singcurr +' not found in current list')
            else:
                if sorted(currset) != sorted(thisCurrDict[singcurr]):
                    if Debug:
                        for icount in range(len(currset)):
                            print currset[icount], thisCurrDict[singcurr][icount],currset[icount] ==  thisCurrDict[singcurr][icount]
                    raise LookupError(singcurr +' has different set list as ' + icurr)
                else:
                    outCurrDict[icurr.replace('doub','')] = currset
        elif 'sing' in icurr:
            doubcurr = icurr.replace('sing','doub')
            if doubcurr not in thisCurrDict.keys():
                raise LookupError(doubcurr +' not found in current list')
            else:
                if sorted(currset) != sorted(thisCurrDict[doubcurr]):
                    if Debug:
                        for icount in range(len(currset)):
                            print currset[icount], thisCurrDict[doubcurr][icount],currset[icount] ==  thisCurrDict[doubcurr][icount]
                    raise LookupError(doubcurr +' has different set list as ' + icurr)
                else:
                    outCurrDict[icurr.replace('sing','')] = currset
        else:
            raise IOError('Depreciated file structure for current '+icurr+' , please check current list')
    return outCurrDict

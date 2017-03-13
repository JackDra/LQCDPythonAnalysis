#!/usr/bin/env python

from SetLists import *
import glob,os

def ReadSetDir(thisdir,thisgamma=''):
    thisSetList = glob.glob(thisdir+'/*.xml')
    if '/' in thisgamma:
        gl,gr = thisgamma.split('/')
    else:
        gl,gr = '',''
    thisSetList = [(iS.replace(thisdir,'')
                    .replace(thisgamma,'')
                    .replace(gl,'')
                    .replace(gr,'')
                    .replace('.xml','')) for iS in thisSetList]
    if OnlySelVar:
        SetOut = []
        for iSet in thisSetList:
            if (any([ism in iSet for ism in DefSmList]) or any([itvar in iSet for itvar in DefTvarPicked])
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
def GetCurrDict(thisCurrTypes,klist=False):
    thisCurrDict = OrderedDict()
    for iCurr in thisCurrTypes:
        thisSetList = None
        if klist == False:
            print 'Reading ' , iCurr , ' Form Factors '
            # thisCurrDict[iCurr] = ReadAllDir(outputdir[0]+'FormFactors/'+iCurr+'/',thisgamma=iCurr)
            thisCurrDict[iCurr] = ReadSetDir(outputdir[0]+'FormFactors/'+iCurr+'/',thisgamma=iCurr)
            if thisSetList == None: 
                thisSetList = set(thisCurrDict[iCurr])
            else:
                thisSetList = thisSetList.intersection(thisCurrDict[iCurr])
            # print 'Sets Found:\n','\n'.join(thisCurrDict[iCurr][ikappa])
            # print ''
        else:            
            thisCurrDict[iCurr] = {}
            for ikappa in klist:
                print 'Reading ' , iCurr , ' Form Factors '
                # thisCurrDict[iCurr] = ReadAllDir(outputdir[0]+'FormFactors/'+iCurr+'/',thisgamma=iCurr)
                thisCurrDict[iCurr][ikappa] = ReadSetDir(outputdir[0].replace(str(kappa),ikappa)+'FormFactors/'+iCurr+'/',thisgamma=iCurr)
                if thisSetList == None: 
                    thisSetList = set(thisCurrDict[iCurr][ikappa])
                else:
                    thisSetList = thisSetList.intersection(thisCurrDict[iCurr][ikappa])
                # print 'Sets Found:\n','\n'.join(thisCurrDict[iCurr][ikappa])
                # print ''

        thisSetList = list(thisSetList)
        thisSetList.sort()
        if Debug: print 'Sets Found:\n','\n'.join(thisSetList)
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
                        for icount in xrange(len(currset)):
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
                        for icount in xrange(len(currset)):
                            print currset[icount], thisCurrDict[doubcurr][icount],currset[icount] ==  thisCurrDict[doubcurr][icount]
                    raise LookupError(doubcurr +' has different set list as ' + icurr)
                else:
                    outCurrDict[icurr.replace('sing','')] = currset
        else:
            raise IOError('Depreciated file structure for current '+icurr+' , please check current list')
    return outCurrDict

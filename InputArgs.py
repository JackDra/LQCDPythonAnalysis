#!/usr/bin/env python

from Params import *
from SetLists import *
from FFParams import *
from MiscFuns import DelDubs
import sys
    
    
def DefWipeWarning():
    if DefWipe:
        thisinput = raw_input("Warning: DefWipe is true, Do you want to wipe existing data? (y/n)")
        if thisinput != 'y': sys.exit()

def ShowSetLists(thissetlist):
    print 
    print 'Set Lists:'
    print '\n'.join(thissetlist)
    print 

def ShowMethodList(thismethodlist):
    print 
    print 'Method Lists:'
    print '\n'.join(thismethodlist)
    print 

def ExpandSetList(thisSL):
    SLout = []
    for iset in thisSL:
        if 'SMSET' in iset:
            SLout += CreateGenericSet(CMTSinkList,DefSmearList,[],[])
        elif 'TSINKSET' in iset:
            SLout += CreateStateTsinkSet('sm32',AllTSinkList)
        elif 'PoFSET' in iset:
            SLout += CreateREvecSet(PoFTSinkList,[PickedState],PoFTvarList)[0]
        elif 'REvecSET' in iset:
            SLout += CreateREvecSet(REvecTSinkList,[PickedState],REvecTvarList)[0]
        elif 'CMSET' in iset:
            SLout += CreateREvecSet(CMTSinkList,[PickedState],AnaTvarList)[0]
        elif 'ALL' in iset:
            SLout += DefSetList
        else:
            SLout.append(iset)
    return SortMySet(SLout)

def ExpandMethodList(thisML):
    MLout = []
    for imethod in thisML:
        if imethod == 'OSF':
            MLout += ['OSF'+iOSF for iOSF in OSFFileFlags]
        elif imethod == 'TSF':
            MLout += ['TSF'+iTSF for iTSF in TSFFileFlags]
        elif imethod == 'Tsink':
            MLout += ['RF','Fits','SumMeth','OSFTsink','TSFTsink','TSFtest32','TSFSmall']
        elif imethod == 'CM' or imethod == 'PoF':
            MLout += ['RF','Fits','OSFCM','TSFCM']
        else:
            MLout.append(imethod)
    return DelDubs(MLout)

def InputParams(inputparams):
    feedout = {}
    feedout['anaproc'] = AnaProc
    feedout['gamma'] = ''
    feedout['set'] = DefSetList
    feedout['method'] = MethodList
    feedout['current'] = NoFFList.keys()
    feedout['mom'] = RunMomList
    SkipDefWipe = False
    for isys in inputparams:
        if isys[0] != '-':
            raise IOError("input arguments are specified with -, see -h for help")
        elif '-h' in isys:
            print 'commands are (with comma separated lists):'
            print '-g= specifies gamma matricies, choose from:'
            print 'Form is "projector""operator""real=blank, complex=cmplx"'
            print 'e.g. P3g1g2cmplx or P4g4'
            print 
            print '-s= specifies set list to use, choose from:\n' + '\n'.join(DefSetList)+'\n'
            print '-m= specifies Method used, choose from:\n' + '\n'.join(MethodList)+'\n'
            print '-c= specifies Current to look at, choose from:\n' + '\n'.join(CurrentDSList)+'\n'
            print "-p= specifies the momentium list to use, form is 'q = X Y Z', X,Y,Z = -3,-2,-1,0,1,2,3"
            print "-np= specifies the maximum number of processors used for this job"
            print "-noprompt does not display any prompts (for long runs)"
            print 
            exit()
        elif '-g' in isys:
            feedout['gamma'] = isys.replace('-g=','').split(',')
        elif '-s' in isys:
            feedout['set'] = ExpandSetList(isys.replace('-s=','').split(','))[0]
            for isl in feedout['set']:
                if isl not in DefSetList:
                    print 'Warning, ' + isl + ' not found in set list, skipping.'
                    feedout['set'].remove(isl)
            if len(feedout['set']) == 0:
                print 'Nothing found for set list, using default list'
                feedout['set'] = DefSetList
        elif '-p' in isys:
            momhold = isys.replace('-p=','').split(',')
            feedout['mom'] = []
            for imom in momhold:
                if 'zmom' in imom:
                    feedout['mom'] += ['q = 0 0 0']
                else:
                    feedout['mom'] += [' '.join(imom.replace('q','q=')).replace('- ','-')]
            for ipl in feedout['mom']:            
                if ipl not in qvecSet:                        
                    feedout['mom'].remove(ipl)
                    print 'Warning, ' + ipl + ' not found in qvecSet list, skipping.'
            if len(feedout['mom']) == 0:
                print 'Nothing found for mom list, using default list'
                feedout['set'] = RunMomList
        elif '-m' in isys:
            feedout['method'] = ExpandMethodList(isys.replace('-m=','').split(','))
            for iml in feedout['method']:
                if iml not in MethodList + ['CM','Tsink','OSF','TSF']:
                    print 'Warning, ' + iml + ' not found in method list, skipping.'
                    feedout['method'].remove(iml)
            if len(feedout['method']) == 0:
                print 'Nothing found for method list, using default list'
                feedout['method'] = MethodList
        elif '-c' in isys:
            feedout['current'] = isys.replace('-c=','').split(',')
            for icl in feedout['current']:
                if icl not in CurrentDSList:
                    print 'Warning, ' + icl + ' not found in current list, skipping.'
                    feedout['current'].remove(icl)
            if len(feedout['current']) == 0:
                print 'Nothing found for current list, using default list'
                feedout['current'] = NoFFList.keys()
        elif '-np' in isys:
            thisAnaProc = int(isys.replace('-np=',''))
            if AnaProc < thisAnaProc:
                print 'number of processors is larger than specified default in setup.cfg, using default'
            else:
                print 'number of processors = ',thisAnaProc
                feedout['anaproc'] = thisAnaProc
        elif '-noprompt' in isys:
            SkipDefWipe = True
    if not SkipDefWipe: DefWipeWarning()
    return feedout


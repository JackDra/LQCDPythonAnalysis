#!/usr/bin/env python

from Params import *
from SetLists import *
from FFParams import *
from MiscFuns import DelDubs
from CombParams import *
from FitParams import FitCutArgs
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

def ShowCombList(thiscomblist):
    print 
    print 'Combinging Lists:'
    print '\n'.join(thiscomblist)
    print 

def ExpandSetList(thisSL):
    SLout = []
    for iset in thisSL:
        if 'SMSET' in iset:
            SLout += CreateGenericSet(CMTSinkList,DefiSmearList,DefjSmearList,[],[])
        elif 'TSINKSET' in iset:
            SLout += CreateStateTsinkSet(SingSmList[0],AllTSinkList)
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
    feedout['cut'] = FitCutArgs
    feedout['method'] = MethodList
    feedout['current'] = NoFFList.keys()
    feedout['mom'] = RunMomList
    feedout['comb'] = DefDSList + CombList
    feedout['FFcomb'] = CombFFList
    feedout['FFcombNS'] = CombNSFFList
    feedout['klist'] = kappalist
    feedout['Wein'] = False
    feedout['DoCurr'] = True
    feedout['ffgraph'] = 'All'
    feedout['ForceTitle'] = False
    SkipDefWipe = False
    StillTitle = False
    for isys in inputparams:
        if isys[0] != '-':
            if StillTitle:
                feedout['ForceTitle'] += ' ' + isys
            else:
                raise IOError("input arguments are specified with -, see -h for help")
        else:
            StillTitle = False
            if '-h' in isys:
                print 'commands are (with comma separated lists):'
                print '-g= specifies gamma matricies, choose from:'
                print 'Form is "projector""operator""real=blank, complex=cmplx"'
                print 'e.g. P3g1g2cmplx or P4g4'
                print 
                print '-s= specifies set list to use, choose from:\n' + '\n'.join(DefSetList)+'\n'
                print '-m= specifies Method used, choose from:\n' + '\n'.join(MethodList)+'\n'
                print '-c= specifies Current to look at, choose from:\n' + '\n'.join(CurrentDSList)+'\n'
                print '-cut= specifies cut list, choose from:\n' + '\n'.join(FitCutArgs)+'\n'
                print '-DS= specifies how to combine DS for TryCombine.py, choose from:\n' + '\n'.join(DefDSList + CombList)+'\n'
                print '-FF= specifies how to combine Form Factors for TryFFCombine.py, choose from:\n' + '\n'.join(CombNSFFList)+'\n'
                print '-DoList= specifies a particular FF to plot in GraphFFs.py, choose ONE from:\n' + '\n'.join(DefGraphDoList+['All'])+'\n'
                print "-p= specifies the momentium list to use, form is 'q = X Y Z', X,Y,Z = -3,-2,-1,0,1,2,3"
                print "-np= specifies the maximum number of processors used for this job"
                print 
                print '-k= specifies kappas to compare when running GraphMKFFs.py from:\n' + '\n'.join(kappalist)+'\n'
                print 
                print "-noprompt does not display any prompts (for long runs)"
                print "-NoCurr ignores all Current calculations (vector, isovector etc..)"
                print "-FT= forces the title for all graphs"
                print 
                exit()
            elif '-NoCurr' in isys:
                feedout['DoCurr'] = False
            elif '-Wein' in isys:
                feedout['Wein'] = True
            elif '-g' in isys:
                feedout['gamma'] = isys.replace('-g=','').split(',')
            elif '-k' in isys:
                feedout['klist'] = isys.replace('-k=','').split(',')
            elif '-FT' in isys:
                feedout['ForceTitle'] = isys.replace('-FT=','')
                StillTitle = True
            elif '-cut' in isys:
                feedout['cut'] = isys.replace('-cut=','').split(',')
                for icl in feedout['cut']:
                    if icl not in FitCutArgs:
                        print 'Warning, ' + icl + ' not found in current list, skipping.'
                        feedout['cut'].remove(icl)
                if len(feedout['cut']) == 0:
                    print 'Nothing found for current list, using default list'
                    feedout['cut'] = FitCutArgs
            elif '-s' in isys:
                feedout['set'] = ExpandSetList(isys.replace('-s=','').split(','))[0]
                for isl in feedout['set']:
                    if isl not in DefSetList+AllCMSetList:
                        print 'Warning, ' + isl + ' not found in set list, skipping.'
                        feedout['set'].remove(isl)
                if len(feedout['set']) == 0:
                    print 'Nothing found for set list, using default list'
                    feedout['set'] = DefSetList
            elif '-p' in isys:
                if isys.replace('-p=','') == 'MassList':
                    feedout['mom'] = RunAvgMomList
                else:
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
            elif '-DS' in isys:
                feedout['comb'] = isys.replace('-DS=','').split(',')
                for iDS in feedout['comb']:
                    if iDS not in DefDSList+CombList:
                        print 'Warning, ' + iDS + ' not found in comb list, skipping.'
                        feedout['comb'].remove(iDS)
                if len(feedout['comb']) == 0:
                    print 'Nothing found for comb list, using default list'
                    feedout['comb'] = DefDSList + CombList            
            elif '-FF' in isys:
                feedout['FFcombNS'] = isys.replace('-FF=','').split(',')
                feedout['FFcomb'] = []
                for iFFc in feedout['FFcombNS']:
                    if len(iFFc) > 0:
                        feedout['FFcomb'].append('/'+iFFc)
                    else:
                        feedout['FFcomb'].append('')
                    if '/'+iFFc not in CombFFList:
                        if iFFc != '':
                            print 'Warning, ' + iFFc + ' not found in FFcomb list, skipping.'
                            feedout['FFcombNS'].remove(iFFc)
                            feedout['FFcomb'].remove('/'+iFFc)
                if len(feedout['FFcomb']) == 0:
                    print 'Nothing found for FFcomb list, using default list'
                    feedout['FFcombNS'] = CombNSFFList         
                    feedout['FFcomb'] = CombFFList         
            elif '-DoList' in isys:
                feedout['ffgraph'] = isys.replace('-DoList=','')
                if feedout['ffgraph'] not in DefGraphDoList:
                    print 'Warning, ' + iDS + ' not found in graphing list list, Doing all.'
                    feedout['ffgraph'] = 'All'
    if not SkipDefWipe: DefWipeWarning()
    if Debug: feedout['anaproc'] = 1
    return feedout


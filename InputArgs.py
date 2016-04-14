#!/usr/bin/env python

from Params import *
from SetLists import *
from FFParams import *


def ExpandSetList(thisSL):
    SLout = []
    for iset in thisSL:
        if 'SMSET' in iset:
            SLout += CreateTSinkStateSet(DefTSinkList[0],DefSmearList,[],[])
        elif 'TSINKSET' in iset:
            SLout += CreateStateTsinkSet('sm32',AllTSinkList)
        elif 'PoFSET' in iset:
            SLout += CreateREvecSet(PoFTSinkList,[1],PoFTvarList)[0]
        elif 'REvecSET' in iset:
            SLout += CreateREvecSet(REvecTSinkList,[1],REvecTvarList)[0]
    return SortMySet(SLout)[0]


def InputParams(inputparams):
    feedout = {}
    feedout['gamma'] = ''
    feedout['set'] = DefSetList
    feedout['method'] = MethodList
    feedout['current'] = CurrOpps.keys()
    feedout['mom'] = RunMomList
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
            print '-c= specifies Current to look at, choose from:\n' + '\n'.join(CurrOpps.keys())+'\n'
            print "-p= specifies the momentium list to use, form is 'q = X Y Z', X,Y,Z = -3,-2,-1,0,1,2,3"
            print 
            exit()
        elif '-g' in isys:
            feedout['gamma'] = isys.replace('-g=','').split(',')
        elif '-s' in isys:
            feedout['set'] = ExpandSetList(isys.replace('-s=','').split(','))
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
                feedout['mom'] += [' '.join(imom).replace('- ','-')]
            for ipl in feedout['mom']:            
                if ipl not in qvecSet:                        
                    feedout['mom'].remove(ipl)
                    if 'zmom' in ipl:
                        feedout['mom'] += ['q = 0 0 0']
                    else:
                        print 'Warning, ' + ipl + ' not found in qvecSet list, skipping.'
            if len(feedout['mom']) == 0:
                print 'Nothing found for set list, using default list'
                feedout['set'] = RunMomList
        elif '-m' in isys:
            feedout['method'] = isys.replace('-m=','').split(',')
            for iml in feedout['method']:
                if iml not in MethodList:
                    print 'Warning, ' + iml + ' not found in method list, skipping.'
                    feedout['method'].remove(iml)
            if len(feedout['method']) == 0:
                print 'Nothing found for method list, using default list'
                feedout['method'] = MethodList
        elif '-c' in isys:
            feedout['current'] = isys.replace('-c=','').split(',')
            for icl in feedout['current']:
                if icl not in CurrOpps.keys():
                    print 'Warning, ' + icl + ' not found in current list, skipping.'
                    feedout['current'].remove(icl)
            if len(feedout['current']) == 0:
                print 'Nothing found for current list, using default list'
                feedout['current'] = CurrOpps.keys()
    return feedout


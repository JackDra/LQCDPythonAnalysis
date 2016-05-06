#!/usr/bin/env python
import numpy as np
import sys
from Params import *
from FitParams import *
from Fitting import OneStateSetFit
from ReadTxt import ReadCfunsnp
from SetLists import *
from OutputData import PrintOSFMassToFile,PrintOSFSetToFile
from CreateCombs import MakeUmD
from OppFuns import CreateGammaList,PrintOpps,WipeSFSet
from FFParams import *
# from multiprocessing import Pool
# from multiprocessing.map import ThreadPool
# from MultiWrap import *
import time
import datetime
import cPickle as pickle
# import pickle
from InputArgs import *

# sys.stdout = open(logfile,'a',0)
# sys.stderr = sys.stdout

print '----------------------------------------------------------------------------------'

if len(sys.argv) < 2: raise IOError('Input CM, Tsink, JustPoF or REvec as first arg')
outfile = sys.argv[1]

feedin = InputParams(sys.argv[2:])

print 'Gamma Input (For re-running): -g=' , feedin['gamma']
ReadGammaList = CreateGammaList(feedin['gamma'],twopt=True)

OSFColList = ['Tsink','CM','JustPoF','REvec']

if outfile == 'Tsink':
    ReadSmearList = ['32']
    ReadTSinkList = AllTSinkList
    CaptString = ['TSINKSET']
    ReadTvarList = []
    ReadREvecTSinkList = []
    ReadREvecTvarList = []
elif outfile == 'CM':
    ReadSmearList = DefSmearList
    ReadTSinkList = [29]
    # CaptString = ['SMSET','CMSET','PoFSET','REvecSET']
    CaptString = ['SMSET','CMSET','PoFSET']
    ReadTvarList = AnaTvarList
    ReadREvecTSinkList = []
    ReadREvecTvarList = []
elif outfile == 'JustPoF':
    outfile = 'CM'
    ReadSmearList = []
    ReadTSinkList = [29]
    ReadTvarList = AnaTvarList
    CaptString = ['PoFSET']
    ReadREvecTSinkList = []
    ReadREvecTvarList = []
elif outfile == 'REvec':
    ReadSmearList = []
    ReadTSinkList = []
    ReadTvarList = []
    CaptString = ['REvecSET']
    ReadREvecTSinkList = REvecTSinkList
    ReadREvecTvarList = REvecTvarList
elif outfile == 'All':
    for iCol in OSFColList:
        os.system(scriptdir+"TryOneStateFit.py " + iCol+' '+' '.join(sys.argv[2:]) )
    exit()
else:
    raise SyntaxError('Input CM or Tsink')

thisFitOSFR = CreateFitList(OSF2ptMinStart,OSF2ptMinEnd,OSF2ptMaxStart,OSF2ptMaxEnd,OSF3ptCutMin,OSF3ptCutMax)

print 'Creating SetList'
[ReadSetList,SetTsink] = ExpandSetList(CaptString)
print ''
print 'nboot = ' + str(nboot)
print 'All Sets:\n' + '\n'.join(ReadSetList)+'\n'


def DoOSF(thisSetList,thisGammaList,OSF2ptarray,twoptGammaMomList):
    print 'Reading Data'
    [data3pt,dump,thisGammaMomList,BorA] = ReadCfunsnp(thisGammaList,thisSetList,thisMomList=feedin['mom'])
    thisGammaMomList['twopt'] = twoptGammaMomList['twopt']
    thisGammaList = thisGammaMomList.keys()
    PrintOpps(thisGammaList)
    print 'Data Read is: ' + BorA
    thisGammaList.remove('twopt')
    ## data2pt = [ ip , iset2pt , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
    ## OSF2ptarray = [ OneFit2pt,OneFit2ptAvg,OneFitChi ]
    #OneFit2pt    = [ ifit2pt , ip , ism  , params ] bs1
    #OneFit2ptAvg = [ ifit2pt , ip , ism  , params ]
    #OneFit2ptChi = [ ifit2pt , ip , ism  ]
    ## data3pt = [ igamma , ip , iset , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)


    start = time.time()
    # thisFitOSFR = [thisFitOSFR[0]]
    for icf,ifit2pt in enumerate(thisFitOSFR):
        thispicklefile = pickledir+'tempOSF'+outfile+'fit'+'to'.join(map(str,ifit2pt))+'.p'
        if not os.path.isfile(thispicklefile):
            perdone = (icf+1)/float(len(thisFitOSFR))
            thisOSF2ptarray = [OSF2ptarray[0][ifir2pt],OSF2ptarray[1][ifir2pt],OSF2ptarray[2][ifir2pt]]
            tempout = OneStateSetFit(thisOSF2ptarray,,data3pt,OSF3ptCutList,thisSetList,thisGammaMomList,[ifit2pt,int(perdone*100)])
            pfile = open( thispicklefile, "wb" )
            pickle.dump( tempout, pfile )
            pfile.close()
            timeleft = (time.time()-start)*((1-perdone)/perdone)
            print 'Current Fit Time: ' , str(datetime.timedelta(seconds=(time.time()-start))) ,' h:m:s  Time Remaining: ' , str(datetime.timedelta(seconds=timeleft)) , ' h:m:s'
            print ''
    print 'Fitting took: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s'

    del data3pt
    del data2pt

    start = time.time()
    OneFit2pt,OneFit2ptAvg,OneFit2ptChi = OSF2ptarray
    OneFit3pt = []
    OneFit3ptAvg = []
    OneFit3ptChi = []
    for icf,ifit2pt in enumerate(thisFitOSFR):
        thispicklefile = pickledir+'tempOSF'+outfile+'fit'+'to'.join(map(str,ifit2pt))+'.p'
        print 'Reading Picked file: ' , thispicklefile , '                                \r',
        if os.path.isfile(thispicklefile):
            pfile = open( thispicklefile, "rb" )
            [OSF2pt,OSF3pt,OSF2ptAvg,OSF3ptAvg,OSF2ptChi,OSF3ptChi] = pickle.load( pfile )
            pfile.close()
            OneFit3pt.append(OSF3pt)
            OneFit3ptAvg.append(OSF3ptAvg)
            OneFit3ptChi.append(OSF3ptChi)
        else:
            raise IOError('Pickled file not found: ' + thispicklefile)

    print 'Reading Pickled files took: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s'

    #OneFit3pt    = [ ifit2pt , igamma , ip , iset , ifit3pt , params ] bs1
    #OneFit3ptAvg = [ ifit2pt , igamma , ip , iset , ifit3pt , params ]
    #OneFit3ptChi = [ ifit2pt , igamma , ip , iset , ifit3pt ]

    start = time.time()
    print 'Printing OSF Results to file: \r',
    WipeSFSet(outputdir,thisGammaList+['twopt'],'OSF'+outfile,'One',setlist=thisSetList)
    PrintOSFMassToFile(OneFit2pt,OneFit2ptChi,thisSetList,thisFitOSFR,outfile,thisGammaMomList['twopt'])
    PrintOSFSetToFile(OneFit3pt,OneFit3ptChi,thisGammaMomList,thisSetList,thisFitOSFR,outfile)

    for icf,ifit2pt in enumerate(thisFitOSFR):
        thispicklefile = pickledir+'tempOSF'+outfile+'fit'+'to'.join(map(str,ifit2pt))+'.p'
        print 'Removing Picked file: ' , thispicklefile , '                                \r',
        os.remove(thispicklefile)

    print 'Printting OSF Results to file  took: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s'

print 'Reading and fitting 2 point correlator data'
[dump,data2pt,twoptGammaMomList,dump3] = ReadCfunsnp(['twopt'],ReadSetList,thisMomList=feedin['mom'])
for icf,ifit2pt in enumerate(thisFitOSFR):
    OSF2ptarray = OneStateSet2pt(data2pt,thisSetList,twoptGammaMomList,ifit2pt)
print 'Reading and fitting 2 point correlators finished'
print ''

for igamma in ReadGammaList:
    if 'doub' not in igamma and 'sing' not in igamma:
        print 'Running ' + igamma
        DoOSF(ReadSetList,[igamma,'doub'+igamma,'sing'+igamma],OSF2ptarray,twoptGammaMomList)
    elif igamma.replace('doub','').replace('sing','') not in ReadGammaList:
        print 'Running ' + igamma
        DoOSF(ReadSetList,[igamma],OSF2ptarray,twoptGammaMomList)

#!/usr/bin/env python
import numpy as np
import sys
from Params import *
from FitParams import *
from Fitting import OneStateSetFit,OneStateSet2pt
from ReadTxt import ReadCfunsnp
from SetLists import *
from OutputData import PrintOSFMassToFile,PrintOSFSetToFile
from CreateCombs import MakeUmD
from OppFuns import CreateGammaList,PrintOpps,WipeSFSet
from FFParams import *
from multiprocessing import Pool
# from multiprocessing.map import ThreadPool
from MultiWrap import *
import time
import datetime
import cPickle as pickle
# import pickle
from InputArgs import *
from CheckXml import *
from MiscFuns import GetTimeStr

# sys.stdout = open(logfile,'a',0)
# sys.stderr = sys.stdout

print '----------------------------------------------------------------------------------'

if len(sys.argv) < 2: raise IOError('Input CM, Tsink, JustPoF or REvec as first arg')
outfile = sys.argv[1]

feedin = InputParams(sys.argv[2:])
DefWipeWarning()

print 'Gamma Input (For re-running): -g=' , feedin['gamma']
ReadGammaList = CreateGammaList(feedin['gamma'],twopt=True)

# OSFColList = ['Tsink','CM','JustPoF','REvec']
OSFColList = ['Tsink','CM','JustPoF']

twoptoutfile = outfile
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
elif outfile == 'JustCM':
    outfile = 'CM'
    ReadSmearList = DefSmearList
    ReadTSinkList = [29]
    # CaptString = ['SMSET','CMSET','PoFSET','REvecSET']
    CaptString = ['CMSET']
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


def DoOSF(thisSetList,thisGammaList,OSF2ptarray,twoptGammaMomList,thisMomList):
    print 'Running ' + thisGammaList[0] + ' ' + thisMomList[0]
    totstart = time.time()
    mprint( 'Reading Data')
    [data3pt,dump,thisGammaMomList,BorA,infolistRF,infolist2pt] = ReadCfunsnp(thisGammaList,thisSetList,thisMomList=thisMomList)
    thisMom = qstrTOqcond(thisMomList[0])
    thisGammaMomList['twopt'] = twoptGammaMomList['twopt']
    thisGammaList = thisGammaMomList.keys()
    # PrintOpps(thisGammaList)
    mprint( 'Data Read is: ' + BorA)
    thisGammaList.remove('twopt')
    ## data3pt = [ igamma , ip , iset , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)


    start = time.time()
    # thisFitOSFR = [thisFitOSFR[0]]
    thispicklelist = []
    for icf,ifit2pt in enumerate(thisFitOSFR):
        thispickledir = pickledir+thisGammaList[0]+'/'+thisMom
        mkdir_p(thispickledir)
        thispicklelist.append(thispickledir+'/tempOSF'+outfile+'fit'+'to'.join(map(str,ifit2pt))+thisGammaList[0]+thisMom+'.p')
        if not os.path.isfile(thispicklelist[icf]):
            perdone = (icf+1)/float(len(thisFitOSFR))
            tempout = OneStateSetFit(OSF2ptarray[icf],data3pt,OSF3ptCutList,thisSetList,thisGammaMomList,[ifit2pt,int(perdone*100)])
            pfile = open( thispicklelist[icf], "wb" )
            pickle.dump( tempout, pfile )
            pfile.close()
            timeleft = (time.time()-start)*((1-perdone)/perdone)
            mprint( 'Current Fit Time: ' , str(datetime.timedelta(seconds=(time.time()-start))) ,' h:m:s  Time Remaining: ' , str(datetime.timedelta(seconds=timeleft)) , ' h:m:s')
            mprint( '')
    mprint( 'Fitting took: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s')

    del data3pt

    start = time.time()
    OneFit3pt = []
    OneFit3ptAvg = []
    OneFit3ptChi = []
    for icf,(ifit2pt,thispicklefile) in enumerate(zip(thisFitOSFR,thispicklelist)):
        mprint( 'Reading Picked file: ' , thispicklefile , '                                \r',)
        if os.path.isfile(thispicklefile):
            pfile = open( thispicklefile, "rb" )
            [OSF3pt,OSF3ptAvg,OSF3ptChi] = pickle.load( pfile )
            pfile.close()
            OneFit3pt.append(OSF3pt)
            OneFit3ptAvg.append(OSF3ptAvg)
            OneFit3ptChi.append(OSF3ptChi)
        else:
            raise IOError('Pickled file not found: ' + thispicklefile)

    mprint( 'Reading Pickled files took: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s')

    #OneFit3pt    = [ ifit2pt , igamma , ip , iset , ifit3pt , params ] bs1
    #OneFit3ptAvg = [ ifit2pt , igamma , ip , iset , ifit3pt , params ]
    #OneFit3ptChi = [ ifit2pt , igamma , ip , iset , ifit3pt ]

    start = time.time()
    mprint( 'Printing OSF Results to file: \r',)
    # WipeSFSet(outputdir,thisGammaList+['twopt'],'OSF'+outfile,'One',setlist=thisSetList)
    print 'Debugging: OSF ' + thisGammaList[0]+' ' + thisMom + ' took ' , str(datetime.timedelta(seconds=time.time()-totstart)) , ' h:m:s'
    PrintOSFSetToFile(OneFit3pt,OneFit3ptChi,thisGammaMomList,thisSetList,thisFitOSFR,outfile,infolistRF)

    for icf,(ifit2pt,thispicklefile) in enumerate(zip(thisFitOSFR,thispicklelist)):
        mprint( 'Removing Picked file: ' , thispicklefile , '                                \r',)
        os.remove(thispicklefile)

    print 'OSF ' + thisGammaList[0]+' ' + thisMom + ' took ' , str(datetime.timedelta(seconds=time.time()-totstart)) , ' h:m:s'




picklefile2pt = pickledir+'tempOSF'+twoptoutfile+'fittwopt.p'
if os.path.isfile(picklefile2pt):
    print '2 point picked file found, reading in'
    with open( picklefile2pt, "rb" ) as pfile:
        OSF2ptarray,twoptGammaMomList = pickle.load( pfile )
    print '2 point picked file read in'        
else:
    print 'Reading and fitting 2 point correlator data'
    [dump,data2pt,twoptGammaMomList,dump3,dump4,infolist2pt] = ReadCfunsnp(['twopt'],ReadSetList,thisMomList=feedin['mom'])
    ## data2pt = [ ip , iset2pt , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
    OSF2ptarray = []
    OneFit2pt = []
    OneFit2ptChi = []
    start = time.time()
    for icf,ifit2pt in enumerate(thisFitOSFR):
        perdone = (icf+1)/float(len(thisFitOSFR))
        OSF2ptarray.append(OneStateSet2pt(data2pt,ReadSetList,twoptGammaMomList,[ifit2pt,int(perdone*100)]))
        OneFit2pt.append(OSF2ptarray[-1][0])
        OneFit2ptChi.append(OSF2ptarray[-1][2])
        timeleft = (time.time()-start)*((1-perdone)/perdone)
        print 'Current Fit Time: ' , str(datetime.timedelta(seconds=(time.time()-start))) ,' h:m:s  Time Remaining: ' , str(datetime.timedelta(seconds=timeleft)) , ' h:m:s'
        ## OSF2ptarray = [ ifit2pt , OneFit2pt/OneFit2ptAvg/OneFitChi ]
        #OneFit2pt    = [ ifit2pt , ip , ism  , params ] bs1
        #OneFit2ptAvg = [ ifit2pt , ip , ism  , params ]
        #OneFit2ptChi = [ ifit2pt , ip , ism  ]
    del data2pt
    print 'Pickling 2 point correlators'
    with open( picklefile2pt, "wb" ) as pfile:
        pickle.dump( [OSF2ptarray,twoptGammaMomList], pfile )
    print 'Printing 2 point correlators to file'
    PrintOSFMassToFile(OneFit2pt,OneFit2ptChi,ReadSetList,thisFitOSFR,twoptoutfile,twoptGammaMomList['twopt'],infolist2pt)


inputparams = []
for igamma in ReadGammaList:
    if 'twopt' in igamma: continue
    if 'doub' not in igamma and 'sing' not in igamma:
        if DefWipe:
            QueMomList = feedin['mom']
        else:
            QueMomList = Check3ptAllSets([igamma,'doub'+igamma,'sing'+igamma],ReadSetList,thisMomList=feedin['mom'],CheckType='OSF'+outfile)
            QueMomList = QueMomList[igamma]
        for imom in QueMomList:
            # print 'adding to que: ' , igamma , imom
            inputparams.append((ReadSetList,[igamma,'doub'+igamma,'sing'+igamma],OSF2ptarray,twoptGammaMomList,[imom]))
    elif igamma.replace('doub','').replace('sing','') not in ReadGammaList:
        if DefWipe:
            QueMomList = feedin['mom']
        else:
            QueMomList = Check3ptAllSets([igamma],ReadSetList,thisMomList=feedin['mom'],CheckType='OSF'+outfile)
            QueMomList = QueMomList[igamma]
        for imom in QueMomList:
            # print 'adding to que: ' , igamma , imom
            inputparams.append((ReadSetList,[igamma],OSF2ptarray,twoptGammaMomList,[imom]))


RunStart = time.time()
if len(inputparams) > 0:
    if DoMulticore:
        print 'Running Multicore'
        makeContextFunctions(DoOSF)
        thisPool = Pool(min(len(inputparams),feedin['anaproc']))
        thisPool.map(DoOSF.mapper,inputparams)
        thisPool.close()
        thisPool.join()
    else:
        print 'Running Single Core'
        for iin in inputparams: DoOSF(*iin)
else:
    print 'nothing to calculate'
# print 'removing pickled 2pt file'
# if os.path.isfile(picklefile2pt): os.remove(picklefile2pt)
print 'all finished'
print 'total time: ', GetTimeStr(time.time() - RunStart)

#!/usr/bin/env python
import numpy as np
import sys
from Params import *
from FitParams import *
from Fitting import MomTSSetFit,MomTSSetFit2pt
from ReadTxt import ReadCfunsnp
from SetLists import *
from OutputData import PrintTSFMassToFile,PrintTSFSetToFile
from OppFuns import CreateGammaList,PrintOpps,WipeSF
from FFParams import *
# from multiprocessing.pool import ThreadPool
from multiprocessing import Pool
from MultiWrap import *
import time
import datetime
import os
import cPickle as pickle
# import pickle
import warnings
from InputArgs import *
from CheckXml import *
from MiscFuns import GetTimeStr

print '----------------------------------------------------------------------------------'

# def fxn():
#     warnings.warn("", RuntimeWarning)
    
# with warnings.catch_warnings():
#     warnings.simplefilter("ignore")
#     fxn()


if len(sys.argv) < 2: raise IOError('Input CM, Tsink or Sm as first arg')
outfile = sys.argv[1]
feedin = InputParams(sys.argv[2:])

print 'Gamma Input (For re-running): -g=' , feedin['gamma']
ReadGammaList = CreateGammaList(feedin['gamma'],twopt=True)

TSFColList = ['Tsink','test32','Small','CM']

if outfile == 'Tsink':
    print 'Tsink run'
    CaptString = ['TSINKSET']
    ReadSmearList = ['32']
    ReadTSinkList = AllTSinkList
    ReadStateList = []
    ReadTvarList = []
    ReadREvecTSinkList = []
    ReadREvecTvarList = []
elif outfile == 'test32':
    print 'test32 run'
    ReadSmearList = ['32']
    ReadTSinkList = [32,35,38]
    CaptString = ['tsink'+str(its)+'sm32' for its in ReadTSinkList]
    ReadStateList = []
    ReadTvarList = []
    ReadREvecTSinkList = []
    ReadREvecTvarList = []
elif outfile == 'Small':
    print 'Small run'
    ReadSmearList = ['32']
    ReadTSinkList = [26,29,32]
    CaptString = ['tsink'+str(its)+'sm32' for its in ReadTSinkList]
    ReadStateList = []
    ReadTvarList = []
    ReadREvecTSinkList = []
    ReadREvecTvarList = []
elif outfile == 'CM':
    print 'CM run'
    ReadSmearList = DefSmearList
    ReadTSinkList = [29]
    CaptString = ['SMSET','CMSET']
    ReadStateList = []
    ReadTvarList = []
    ReadREvecTSinkList = []
    ReadREvecTvarList = []
elif outfile == 'JustCM':
    print 'JustCM run'
    outfile = 'CM'
    ReadSmearList = DefSmearList
    ReadTSinkList = [29]
    CaptString = ['CMSET']
    ReadStateList = []
    ReadTvarList = []
    ReadREvecTSinkList = []
    ReadREvecTvarList = []
elif outfile == 'REvec':
    print 'REvec run'
    CaptString = ['REvecSET']
    ReadSmearList = []
    ReadTSinkList = []
    ReadTvarList = []
    ReadREvecTSinkList = REvecTSinkList
    ReadREvecTvarList = REvecTvarList
elif outfile == 'PoF'+str(PoFShifts):
    thisPoF = 'PoF'+str(PoFShifts)
    print thisPoF+' run'
    CaptString = ['PoFSET']
    ReadSmearList = []
    ReadTSinkList = []
    ReadTvarList = []
    ReadREvecTSinkList = PoFTSinkList
    ReadREvecTvarList = PoFTvarList
elif outfile == 'All':
    for iCol in TSFColList:
        os.system(scriptdir+"TryTwoStateFit.py " + iCol+' ' +' '.join(sys.argv[2:]))
    exit()
else:
    raise SyntaxError('Input CM, Tsink or Sm for first arg')



def DoTSF(thisSetList,thisGammaList,TSF2ptarray,twoptGammaMomList,thisMomList):
    print 'Running ' + thisGammaList[0] + ' ' + thisMomList[0]
    totstart = time.time()
    mprint( 'Reading Data')

    [data3pt,dump,thisGammaMomList,BorA,infoRF,info2pt] = ReadCfunsnp(thisGammaList,thisSetList,thisMomList=thisMomList)
    thisGammaMomList['twopt'] = twoptGammaMomList['twopt']
    thisGammaList = thisGammaMomList.keys()
    thisMom = qstrTOqcond(thisMomList[0])
    mprint( 'Data Read is: ' + BorA)

    ## data2pt = [ ip , iset2pt , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
    ## data3pt = [ igamma , ip , iset , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)

    thisGammaList.remove('twopt')

    # thisFitTSFR = [thisFitTSFR[0]]
    start = time.time()
    thispicklelist = []
    for icf,ifit2pt in enumerate(thisFitTSFR):
        thispickledir = pickledir+thisGammaList[0] +'/'+ thisMom
        mkdir_p(thispickledir)
        thispicklelist.append(thispickledir+'/tempTSF'+outfile+'fit'+'to'.join(map(str,ifit2pt))+thisGammaList[0] + thisMom+'.p')
        if not os.path.isfile(thispicklelist[icf]):
            perdone = (icf+1)/float(len(thisFitTSFR))
            tempout = MomTSSetFit(TSF2ptarray[icf],data3pt,TSF3ptCutList,thisSetList,thisGammaMomList,[ifit2pt,int(perdone*100)])
            pfile = open( thispicklelist[icf], "wb" )
            pickle.dump( tempout, pfile )
            pfile.close()
            timeleft = (time.time()-start)*((1-perdone)/perdone)
            mprint( 'Current Fit Time: ' , str(datetime.timedelta(seconds=(time.time()-start))) ,' h:m:s  Time Remaining: ' , str(datetime.timedelta(seconds=timeleft)) , ' h:m:s')
            mprint( '')
    # mprint( 'Fitting took: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s'

    del data3pt

    start = time.time()
    TwoFit3pt = []
    TwoFit3ptAvg = []
    TwoFit3ptChi = []
    for icf,(ifit2pt,thispicklefile) in enumerate(zip(thisFitTSFR,thispicklelist)):
        mprint( 'Reading Picked File: ' , thispicklefile , '                         \r',)
        if os.path.isfile(thispicklefile):
            pfile = open( thispicklefile, "rb" )
            [TSF3pt,TSF3ptAvg,TSF3ptChi] = pickle.load( pfile )
            pfile.close()
            # [TSF2pt,TSF3pt,TSF2ptAvg,TSF3ptAvg,TSF2ptChi,TSF3ptChi] = TwoStateSetFit(np.rollaxis(np.array(data2pt),1),data3pt,ifit2pt,TSF3ptCutList,ReadTSinkList)
            TwoFit3pt.append(TSF3pt)
            TwoFit3ptAvg.append(TSF3ptAvg)
            TwoFit3ptChi.append(TSF3ptChi)
        else:
            raise IOError('Pickled file not found: ' + thispicklefile)

    mprint( 'Reading Pickled files took: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s          ')



    #TwoFit2pt    = [ ifit2pt , ip , params ] bs1
    #TwoFit2ptAvg = [ ifit2pt , ip , params ]
    #TwoFit2ptChi = [ ifit2pt , ip ]
    #TwoFit3pt    = [ ifit2pt , ip , igamma , istate , ifit3pt , params ] bs1
    #TwoFit3ptAvg = [ ifit2pt , ip , igamma , istate , ifit3pt , params ]
    #TwoFit3ptChi = [ ifit2pt , ip , igamma , istate , ifit3pt ]

    start = time.time()
    mprint( 'Printing TSF Results to file: ')
    # WipeSF(outputdir,thisGammaList+['twopt'],'TSF'+outfile,'Two',statelist=ReadStateList,todtlist=ReadTvarList,smlist=ReadSmearList)
    PrintTSFSetToFile(TwoFit3pt,TwoFit3ptChi,thisGammaMomList,thisSetList,thisFitTSFR,outfile,infoRF)

    for icf,(ifit2pt,thispicklefile) in enumerate(zip(thisFitTSFR,thispicklelist)):
        mprint( 'Removing Picked File: ' , thispicklefile , '                         \r',)
        os.remove(thispicklefile)

    print 'TSF ' + thisGammaList[0]+' ' + thisMom +' took ' , str(datetime.timedelta(seconds=time.time()-totstart)) , ' h:m:s'
    


thisFitTSFR = CreateFitList(TSF2ptMinStart,TSF2ptMinEnd,TSF2ptMaxStart,TSF2ptMaxEnd,TSF3ptCutMin,TSF3ptCutMax)

print 'Creating SetList'
[ReadSetList,SetTsink] = ExpandSetList(CaptString)
if outfile == 'CM': ReadSetList = PickCM(ReadSetList)
print ''
print 'nboot = ' + str(nboot)
print 'All Sets:\n' + '\n'.join(ReadSetList)+'\n'



    
picklefile2pt = pickledir+'tempTSF'+outfile+'fittwopt.p'
print 'looking for ',picklefile2pt
if os.path.isfile(picklefile2pt):
    print '2 point picked file found, reading in'
    with open( picklefile2pt, "rb" ) as pfile:
        TSF2ptarray,twoptGammaMomList = pickle.load( pfile )
    print '2 point picked file read in'        
else:
    print 'Reading and fitting 2 point correlator data'
    [dump,data2pt,twoptGammaMomList,dump3,infoRF,info2pt] = ReadCfunsnp(['twopt'],ReadSetList,thisMomList=feedin['mom'])
    ## data2pt = [ ip , iset2pt , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
    TSF2ptarray = []
    TwoFit2pt = []
    TwoFit2ptChi = []
    for icf,ifit2pt in enumerate(thisFitTSFR):
        perdone = (icf+1)/float(len(thisFitTSFR))
        TSF2ptarray.append(MomTSSetFit2pt(data2pt,ReadSetList,twoptGammaMomList,[ifit2pt,int(100*perdone)]))
        TwoFit2pt.append(TSF2ptarray[-1][0])
        TwoFit2ptChi.append(TSF2ptarray[-1][2])
        ## TSF2ptarray = [ ifit2pt , TwoFit2pt/TwoFit2ptAvg/TwoFitChi ]
        #TwoFit2pt    = [ ifit2pt , ip , ism  , params ] bs1
        #TwoFit2ptAvg = [ ifit2pt , ip , ism  , params ]
        #TwoFit2ptChi = [ ifit2pt , ip , ism  ]
    del data2pt
    print 'Pickling 2 point correlators'
    with open( picklefile2pt, "wb" ) as pfile:
        pickle.dump( [TSF2ptarray,twoptGammaMomList], pfile )
    print 'Printing 2 point correlators to file'
    PrintTSFMassToFile(TwoFit2pt,TwoFit2ptChi,ReadSetList,thisFitTSFR,outfile,twoptGammaMomList['twopt'],info2pt)


inputparams = []
for igamma in ReadGammaList:
    if 'twopt' in igamma: continue
    if 'doub' not in igamma and 'sing' not in igamma:
        if DefWipe:
            QueMomList = GetMomFromGamma(igamma,thisMomList=feedin['mom'])
        else:
            QueMomList = Check3ptAllSets(['doub'+igamma,'sing'+igamma],ReadSetList,thisMomList=feedin['mom'],CheckType='TSF'+outfile)
            QueMomList = QueMomList[igamma]
        for imom in QueMomList:
            # print 'adding to que: ' , igamma , imom
            inputparams.append((ReadSetList,['doub'+igamma,'sing'+igamma],TSF2ptarray,twoptGammaMomList,[imom]))
    elif igamma.replace('doub','').replace('sing','') not in ReadGammaList:
        if DefWipe:
            QueMomList = GetMomFromGamma(igamma,thisMomList=feedin['mom'])
        else:
            QueMomList = Check3ptAllSets([igamma],ReadSetList,thisMomList=feedin['mom'],CheckType='TSF'+outfile)
            QueMomList = QueMomList[igamma.replace('doub','').replace('sing','')]
        for imom in QueMomList:
            # print 'adding to que: ' , igamma , imom
            inputparams.append((ReadSetList,[igamma],TSF2ptarray,twoptGammaMomList,[imom]))

RunStart = time.time()
if len(inputparams) > 0:
    if DoMulticore and feedin['anaproc'] > 1 and len(inputparams) > 1:
        print 'Running Multicore'
        makeContextFunctions(DoTSF)
        thisPool = Pool(min(len(inputparams),feedin['anaproc']))
        thisPool.map(DoTSF.mapper,inputparams)
        thisPool.close()
        thisPool.join()
    else:
        print 'Running Single Core'
        for iin in inputparams: DoTSF(*iin)
else:
    print 'nothing to do'        
# print 'removing pickled 2pt file'
# if os.path.isfile(picklefile2pt): os.remove(picklefile2pt)
print 'all finished'
print 'total time: ', GetTimeStr(time.time()- RunStart)

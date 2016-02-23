#!/usr/bin/env python
import numpy as np
import sys
from Params import *
from FitParams import *
from Fitting import MomTSSetFit
from ReadTxt import ReadCfunsnp
from SetLists import *
from OutputData import PrintTSFMassToFile,PrintTSFToFile
from CreateCombs import MakeUmD
from OppFuns import CreateGammaList,PrintOpps,WipeSF
from FFParams import *
# from multiprocessing.pool import ThreadPool
# from multiprocessing import Pool
# from MultiWrap import *
import time
import datetime
import os
import cPickle as pickle
# import pickle
import warnings

print '----------------------------------------------------------------------------------'

def fxn():
    warnings.warn("", RuntimeWarning)
    
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    fxn()


if len(sys.argv) < 2: raise IOError('Input CM, Tsink or Sm as first arg')
thisReadMomList = qvecSet
# thisReadMomList = ['q = 0 0 0']
outfile = sys.argv[1]
print 'Gamma Input (For re-running): ' , sys.argv[2:]
thisGammaList = CreateGammaList(sys.argv[2:],twopt=True)

TSFColList = ['Tsink','test32','Small','CM']
thisFitTSFR = CreateFitList(TSF2ptMinStart,TSF2ptMinEnd,TSF2ptMaxStart,TSF2ptMaxEnd,TSF3ptCutMin,TSF3ptCutMax)

if outfile == 'Tsink':
    print 'Tsink run'
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
    ReadStateList = []
    ReadTvarList = []
    ReadREvecTSinkList = []
    ReadREvecTvarList = []
elif outfile == 'Small':
    print 'Small run'
    ReadSmearList = ['32']
    ReadTSinkList = [26,29,32]
    ReadStateList = []
    ReadTvarList = []
    ReadREvecTSinkList = []
    ReadREvecTvarList = []
elif outfile == 'CM':
    print 'CM run'
    ReadSmearList = DefSmearList
    ReadTSinkList = [29]
    ReadStateList = []
    ReadTvarList = []
    ReadREvecTSinkList = []
    ReadREvecTvarList = []
elif outfile == 'REvec':
    print 'REvec run'
    ReadSmearList = []
    ReadTSinkList = []
    ReadTvarList = []
    ReadREvecTSinkList = REvecTSinkList
    ReadREvecTvarList = REvecTvarList
elif outfile == 'All':
    for iCol in TSFColList:
        os.system(scriptdir+"TryTwoStateFit.py " + iCol+' ' +' '.join(sys.argv[2:]))
    exit()
else:
    raise SyntaxError('Input CM, Tsink or Sm for first arg')


print 'Creating SetList'
[ReadSetList,ReadSet2pt,SetTsink] = CreateSet(thisSmearL=ReadSmearList,thisTvarL=ReadTvarList,thisTSinkL=ReadTSinkList,
                                              thisREvecTvarL=ReadREvecTvarList,thisREvecTSinkL=ReadREvecTSinkList)

print ''
print 'nboot = ' + str(nboot)
print 'All T Sinks: '+ ', '.join(map(str,ReadTSinkList))
# print 'All Operators: \n'+'\n'.join(thisGammaList)+'\n'
print 'All Sets:\n' + '\n'.join(ReadSetList)+'\n'
# print 'All 2pt Sets:\n' + '\n'.join(map(str,ReadSet2pt))
# print 'All Fit Ranges:\n' + '\n'.join(map(str,thisFitTSFR))+'\n'

print 'Reading Data'

[data3pt,data2pt,thisGammaMomList,BorA] = ReadCfunsnp(thisGammaList,ReadSetList,thisMomList=thisReadMomList)
thisGammaList = thisGammaMomList.keys()
PrintOpps(thisGammaList)
print 'Data Read is: ' + BorA

## data2pt = [ ip , iset2pt , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)
## data3pt = [ igamma , ip , iset , it ] = bootstrap1 class (.Avg, .Std, .values, .nboot)

# data3pt,thisGammaList = MakeUmD(data3pt,thisGammaList)
thisGammaList.remove('twopt')
 
# thisFitTSFR = [thisFitTSFR[0]]
start = time.time()
for icf,ifit2pt in enumerate(thisFitTSFR):
    thispicklefile = pickledir+'tempTSF'+outfile+'fit'+'to'.join(map(str,ifit2pt))+'.p'
    if not os.path.isfile(thispicklefile):
        perdone = (icf+1)/float(len(thisFitTSFR))
        tempout = MomTSSetFit(data2pt,data3pt,TSF3ptCutList,ReadSetList,thisGammaMomList,[ifit2pt,int(perdone*100)])
        pfile = open( thispicklefile, "wb" )
        pickle.dump( tempout, pfile )
        pfile.close()
        timeleft = (time.time()-start)*((1-perdone)/perdone)
        print 'Current Fit Time: ' , str(datetime.timedelta(seconds=(time.time()-start))) ,' h:m:s  Time Remaining: ' , str(datetime.timedelta(seconds=timeleft)) , ' h:m:s'
        print ''
# print 'Fitting took: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s'

del data3pt
del data2pt

start = time.time()
TwoFit2pt = []
TwoFit3pt = []
TwoFit2ptAvg = []
TwoFit3ptAvg = []
TwoFit2ptChi = []
TwoFit3ptChi = []
for icf,ifit2pt in enumerate(thisFitTSFR):
    thispicklefile = pickledir+'tempTSF'+outfile+'fit'+'to'.join(map(str,ifit2pt))+'.p'
    print 'Reading Picked File: ' , thispicklefile , '                         \r',
    if os.path.isfile(thispicklefile):
        pfile = open( thispicklefile, "rb" )
        [TSF2pt,TSF3pt,TSF2ptAvg,TSF3ptAvg,TSF2ptChi,TSF3ptChi] = pickle.load( pfile )
        pfile.close()
        # [TSF2pt,TSF3pt,TSF2ptAvg,TSF3ptAvg,TSF2ptChi,TSF3ptChi] = TwoStateSetFit(np.rollaxis(np.array(data2pt),1),data3pt,ifit2pt,TSF3ptCutList,ReadTSinkList)
        TwoFit2pt.append(TSF2pt)
        TwoFit3pt.append(TSF3pt)
        TwoFit2ptAvg.append(TSF2ptAvg)
        TwoFit3ptAvg.append(TSF3ptAvg)
        TwoFit2ptChi.append(TSF2ptChi)
        TwoFit3ptChi.append(TSF3ptChi)
    else:
        raise IOError('Pickled file not found: ' + thispicklefile)

print 'Reading Pickled files took: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s          '
 

 
#TwoFit2pt    = [ ifit2pt , ip , params ] bs1
#TwoFit2ptAvg = [ ifit2pt , ip , params ]
#TwoFit2ptChi = [ ifit2pt , ip ]
#TwoFit3pt    = [ ifit2pt , ip , igamma , istate , ifit3pt , params ] bs1
#TwoFit3ptAvg = [ ifit2pt , ip , igamma , istate , ifit3pt , params ]
#TwoFit3ptChi = [ ifit2pt , ip , igamma , istate , ifit3pt ]

start = time.time()
print 'Printing TSF Results to file: '
WipeSF(outputdir,thisGammaList+['twopt'],'TSF'+outfile,'Two',statelist=ReadStateList,todtlist=ReadTvarList,smlist=ReadSmearList)
PrintTSFMassToFile(TwoFit2pt,TwoFit2ptChi,ReadSetList,thisFitTSFR,outfile,thisGammaMomList['twopt'])
PrintTSFToFile(TwoFit3pt,TwoFit3ptChi,thisGammaMomList,ReadSetList,thisFitTSFR,outfile)

for icf,ifit2pt in enumerate(thisFitTSFR):
    thispicklefile = pickledir+'tempTSF'+outfile+'fit'+'to'.join(map(str,ifit2pt))+'.p'
    print 'Removing Picked File: ' , thispicklefile , '                         \r',
    os.remove(thispicklefile)

print 'Printting TSF Results to file took: ' , str(datetime.timedelta(seconds=time.time()-start)) , ' h:m:s                  '






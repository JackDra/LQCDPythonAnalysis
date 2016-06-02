#!/usr/bin/env python

from Params import *
from FitParams import *
import numpy as np
from BootTest import BootStrap1
from MiscFuns import *
from SetLists import GetTsinkSmLists
from OppFuns import *
from FormFactors import NoFFPars
from OutputXmlData import *



##C3set [ igamma , iset , it ] bs1

def PrintCfunToFile(C3set,thisSetList,thisMomList, thisGammaList,AddDict={}):
    cfundir = outputdir + 'cfuns/'
    for thegamma,gammadata in zip(thisGammaList,C3set):
        gammadir = cfundir+CreateOppDir(thegamma)+'/'
        for iset,setdata in zip(thisSetList,gammadata):
            print 'Printing : ' , thegamma , iset , '                \r',
            PrintToFile(np.array(setdata),gammadir,iset+thegamma,range(64),thisMomList,AddDict=AddDict,frmtflag='e')

##dataset [ igamma , iset , ip , it ] bs1

def PrintSetToFile(dataset,thisSetList,thisMomList, thisGammaList,tsink,AddDict={}):
    for thegamma,gammadata in zip(thisGammaList,dataset):
        gammadir = outputdir+CreateOppDir(thegamma)+'/'
        for iset,setdata in zip(thisSetList,gammadata):
            print 'Printing : ' , thegamma , iset , '                \r',
            if thegamma == 'Mass':
                calcflag = 'Mass'
                setdata = cfunTOmass(setdata)
                tlist = range(64)
            else:
                calcflag = 'Ratio_Factor'
                tlist = range(tsource,int(tsink)+1)
            PrintToFile(setdata,gammadir,iset+thegamma,tlist,thisMomList,AddDict=AddDict)


##sumdata [ igamma , ip , icut , itsink ] bs1
##sumfits [ igamma , ip , icut , fitvals , par ] bs1
##sumfitsChi [ igamma , ip , icut , fitvals ]

def PrintSumSetToFile(sumdata,sumfits,sumfitschi,thisFitList,thissm, thisGammaMomList,thisTSinkList,thisCutList,infoRF):
    for (thegamma,thisMomList),gammadata,gammafitdata,gfdchi,gfitlist in zip(thisGammaMomList.iteritems(),sumdata,sumfits,sumfitschi,thisFitList):
        print 'Printing : ' , thegamma , '                \r',
        gammadir = outputdir+CreateOppDir(thegamma)+'/SumMeth/'
        filename = thissm+thegamma
        infosetRF = [ip[igamma] for ip in infosetRF]
        PrintSumToFile(gammadata,gammafitdata,gfdchi,gammadir,filename,gfitlist,thisMomList,thisTSinkList,thisCutList,infosetRF)



#FitData = [ igamma , ip , icut , iset ]

def PrintFitSetToFile(dataset,datasetChi,thisGammaMomList,thisSetList,thisCutList,infosetRF):
    for igamma,(thisgamma,thismomlist) in enumerate(thisGammaMomList.iteritems()):
        gammadir = outputdir+CreateOppDir(thisgamma)+'/Fits/'
        for iset,thisset in enumerate(thisSetList):
            infosetRF = [ip[igamma][iset] for ip in infosetRF]
            print 'Printing : ' , thisgamma , thisset , '                \r',
            filename = thisset+thisgamma
            PrintFitToFile(dataset[igamma],datasetChi[igamma],iset,gammadir,filename,thismomlist,thisCutList,mominfoRF)


#dataset    = [ cuts , ip , istate ] bs1
#datasetChi = [ cuts , ip , istate ]
#dataset (after roll)   = [ istate , ip , cuts ] bs1
#datasetChi (after roll)= [ istate , ip , cuts ]
#statedata = [ ip , icut ]

def PrintFitMassSetToFile(dataset,datasetChi,thisMomList,thisStateList,thisFitR,AddDict={}):
    dataset = np.rollaxis(np.rollaxis(dataset,1),2)
    datasetChi = np.rollaxis(np.rollaxis(datasetChi,1),2)
    gammadir = outputdir+'Mass/fits/'
    mkdir_p(gammadir)
    for thisstate,statedata,statedataChi in zip(thisStateList,dataset,datasetChi):
        filename = (thisstate+'Mass')
        PrintFitMassToFile(statedata,statedataChi,thisstate,gammadir,filename,thisMomList,thisFitR,AddDict=AddDict)


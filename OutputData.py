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

def PrintCfunToFile(C3set,thisSetList,thisMomList, thisGammaList):
    cfundir = outputdir + 'cfuns/'
    for thegamma,gammadata in zip(thisGammaList,C3set):
        gammadir = cfundir+CreateOppDir(thegamma)+'/'
        bootgammadir = gammadir + 'boots/'
        mkdir_p(bootgammadir)
        for iset,setdata in zip(thisSetList,gammadata):
            print 'Printing : ' , thegamma , iset , '                \r',
            PrintToFile(np.array(setdata),gammadir,iset+thegamma,range(64),thisMomList,frmtflag='e')
            # filename = (bootgammadir + iset+thegamma)
            # PrintBootToFile(np.array(setdata),filename,range(64),thisMomList)

##dataset [ igamma , iset , ip , it ] bs1

def PrintSetToFile(dataset,thisSetList,thisMomList, thisGammaList,tsink):
    for thegamma,gammadata in zip(thisGammaList,dataset):
        gammadir = outputdir+CreateOppDir(thegamma)+'/'
        bootgammadir = gammadir + 'boots/'
        mkdir_p(bootgammadir)
        for iset,setdata in zip(thisSetList,gammadata):
            print 'Printing : ' , thegamma , iset , '                \r',
            if thegamma == 'Mass':
                calcflag = 'Mass'
                setdata = cfunTOmass(setdata)
                tlist = range(64)
            else:
                calcflag = 'Ratio_Factor'
                tlist = range(tsource,int(tsink)+1)
            PrintToFile(setdata,gammadir,iset+thegamma,tlist,thisMomList)
            # filename = (bootgammadir +iset+thegamma)
            # PrintBootToFile(setdata,filename,tlist,thisMomList)


##sumdata [ igamma , ip , icut , itsink ] bs1
##sumfits [ igamma , ip , icut , fitvals , par ] bs1
##sumfitsChi [ igamma , ip , icut , fitvals ]

def PrintSumSetToFile(sumdata,sumfits,sumfitschi,thisFitList,thissm, thisGammaMomList,thisTSinkList,thisCutList):
    for (thegamma,thisMomList),gammadata,gammafitdata,gfdchi,gfitlist in zip(thisGammaMomList.iteritems(),sumdata,sumfits,sumfitschi,thisFitList):
        print 'Printing : ' , thegamma , '                \r',
        gammadir = outputdir+CreateOppDir(thegamma)+'/SumMeth/'
        mkdir_p(bootgammadir)
        filename = thissm+thegamma
        PrintSumToFile(gammadata,gammafitdata,gfdchi,gammadir,filename,gfitlist,thisMomList,thisTSinkList,thisCutList)



#FitData = [ igamma , ip , icut , iset ]

def PrintFitSetToFile(dataset,datasetChi,thisGammaMomList,thisSetList,thisCutList):
    for igamma,(thisgamma,thismomlist) in enumerate(thisGammaMomList.iteritems()):
        gammadir = outputdir+CreateOppDir(thisgamma)+'/Fits/'
        for iset,thisset in enumerate(thisSetList):
            print 'Printing : ' , thisgamma , thisset , '                \r',
            filename = thisset+thisgamma
            PrintFitToFile(dataset[igamma],datasetChi[igamma],iset,gammadir,filename,thismomlist,thisCutList)
            # PrintFitBootToFile(dataset[igamma],bootfilename,iset,thismomlist,thisCutList)


#dataset    = [ cuts , ip , istate ] bs1
#datasetChi = [ cuts , ip , istate ]
#dataset (after roll)   = [ istate , ip , cuts ] bs1
#datasetChi (after roll)= [ istate , ip , cuts ]
#statedata = [ ip , icut ]

def PrintFitMassSetToFile(dataset,datasetChi,thisMomList,thisStateList,thisFitR):
    dataset = np.rollaxis(np.rollaxis(dataset,1),2)
    datasetChi = np.rollaxis(np.rollaxis(datasetChi,1),2)
    gammadir = outputdir+'Mass/fits/'
    mkdir_p(gammadir)
    for thisstate,statedata,statedataChi in zip(thisStateList,dataset,datasetChi):
        filename = (thisstate+'Mass')
        PrintFitMassToFile(statedata,statedataChi,thisstate,gammadir,filename,thisMomList,thisFitR)


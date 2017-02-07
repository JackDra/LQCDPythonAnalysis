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

## IN DEV
##C3setTop [ igamma , iset , iflow , ip , it ] bs1

def PrintTopCfunToFile(C3setTop,thisSetList,thisMomList, thisGammaList,thisTopList,AddDict={},Wein=False):
    if Wein:
        cfundir = outputdir[0] + 'Wein/cfun/'
    else:
        cfundir = outputdir[0] + 'Top/cfun/'
        
    for thegamma,gammadata in zip(thisGammaList,C3setTop):
        gammadir = cfundir+CreateOppDir(thegamma)+'/'
        if 'cmplx' in thegamma:
            thegamma = thegamma.replace('cmplx','Topcmplx')
        else:
            thegamma = thegamma+'Top'
        for iset,setdata in zip(thisSetList,gammadata):
            print 'Printing cfuns Top : ' , thegamma , iset , '                \r',
            Print3ptTopToFile(np.rollaxis(np.array(setdata),1),gammadir,iset+thegamma,thisTopList,range(nt),thisMomList,AddDict=AddDict,frmtflag='e')


##dataset [ igamma , iset , iflow , ip , it ] bs1

def PrintTopSetToFile(C3setTop,thisSetList,thisMomList, thisGammaList,tsink,thisTopList,AddDict={},Wein=False):
    if Wein:
        TopOrWein = 'Wein'
    else:
        TopOrWein = 'Top'
    for thegamma,gammadata in zip(thisGammaList,C3setTop):
        if Wein:
            gammadir = outputdir[0]+'Wein/'+CreateOppDir(thegamma)+'/'
        else:
            gammadir = outputdir[0]+'Top/'+CreateOppDir(thegamma)+'/'
        if 'cmplx' in thegamma:
            thegamma = thegamma.replace('cmplx',TopOrWein+'cmplx')
        else:
            thegamma = thegamma+TopOrWein
        for iset,setdata in zip(thisSetList,gammadata):
            print 'Printing ',TopOrWein,' : ' , thegamma , iset , '                \r',
            calcflag = 'Ratio_Factor_'+TopOrWein
            tlist = range(tsource,int(tsink)+1)
            Print3ptTopToFile(np.rollaxis(np.array(setdata),1),gammadir,iset+thegamma,thisTopList,tlist,thisMomList,AddDict=AddDict)








            

##topdataset [ iset , ttop, ip ,it ] bs1
##dataset [iset , ip , it ] bs1

def PrintAlphaSetToFile(topdataset,dataset,thisSetList,thisMomList, thisTopList, AddDict={},Wein=False):
    if Wein:
        topdir = outputdir[0] + 'Wein/Alpha/'
    else:        
        topdir = outputdir[0] + 'Top/Alpha/'
    mkdir_p(topdir)
    mkdir_p(topdir.replace('Alpha/','cfun/NNQ/'))
    for iset,setdata,topsetdata in zip(thisSetList,dataset,topdataset):
        print 'Printing : ' , iset ,'                \r',
        PrintTopToFile(topsetdata,setdata,topdir,iset,thisTopList,range(nt),thisMomList,AddDict=AddDict)



##C3set [ igamma , iset , ip , it ] bs1

def PrintCfunToFile(C3set,thisSetList,thisMomList, thisGammaList,AddDict={},Top=False):
    if Top == 'Wein':
        cfundir = outputdir[0] + 'Wein/cfun/'        
    elif Top:
        cfundir = outputdir[0] + 'Top/cfun/'
    else:
        cfundir = outputdir[0] + 'cfun/'
    for thegamma,gammadata in zip(thisGammaList,C3set):
        gammadir = cfundir+CreateOppDir(thegamma)+'/'
        for iset,setdata in zip(thisSetList,gammadata):
            print 'Printing : ' , thegamma , iset , '                \r',
            PrintToFile(np.array(setdata),gammadir,iset+thegamma,range(nt),thisMomList,AddDict=AddDict,frmtflag='e')

##dataset [ igamma , iset , ip , it ] bs1

def PrintSetToFile(dataset,thisSetList,thisMomList, thisGammaList,tsink,AddDict={}):
    for thegamma,gammadata in zip(thisGammaList,dataset):
        gammadir = outputdir[0]+CreateOppDir(thegamma)+'/'
        for iset,setdata in zip(thisSetList,gammadata):
            print 'Printing : ' , thegamma , iset , '                \r',
            if thegamma == 'Mass':
                calcflag = 'Mass'
                setdata = cfunTOmass(setdata)
                tlist = range(nt)
            else:
                calcflag = 'Ratio_Factor'
                tlist = range(tsource,int(tsink)+1)
            PrintToFile(setdata,gammadir,iset+thegamma,tlist,thisMomList,AddDict=AddDict)


##sumdata [ igamma , ip , icut , itsink ] bs1
##sumfits [ igamma , ip , icut , fitvals , par ] bs1
##sumfitsChi [ igamma , ip , icut , fitvals ]

def PrintSumSetToFile(sumdata,sumfits,sumfitschi,thisFitList,thissm, thisGammaMomList,thisTSinkList,thisCutList,infoRF):
    for igamma,((thegamma,thisMomList),gammadata,gammafitdata,gfdchi,gfitlist) in enumerate(zip(thisGammaMomList.iteritems(),sumdata,sumfits,sumfitschi,thisFitList)):
        mprint('Printing : ' , thegamma , '                \r',)
        gammadir = outputdir[0]+CreateOppDir(thegamma)+'/SumMeth/'
        filename = thissm+thegamma
        infosetRF = [ip for ip in infoRF[igamma]]
        PrintSumToFile(gammadata,gammafitdata,gfdchi,gammadir,filename,gfitlist,thisMomList,thisTSinkList,thisCutList,infosetRF)



#FitData = [ igamma , ip , icut , iset ]

def PrintFitSetToFile(dataset,datasetChi,thisGammaMomList,thisSetList,thisCutList,infosetRF,flowlist):
    for igamma,(thisgamma,thismomlist) in enumerate(thisGammaMomList.iteritems()):
        gammadir = outputdir[0]+CreateOppDir(thisgamma)+'/Fits/'
        for iset,thisset in enumerate(thisSetList):
            mominfoRF = [ip[iset] for ip in infosetRF[igamma]]
            print 'Printing : ' , thisgamma , thisset , '                \r',
            filename = thisset+thisgamma
            PrintFitToFile(dataset[igamma],datasetChi[igamma],iset,gammadir,filename,thismomlist,thisCutList,mominfoRF,flowlist)


#dataset    = [ cuts , ip , istate ] bs1
#datasetChi = [ cuts , ip , istate ]
#dataset (after roll)   = [ istate , ip , cuts ] bs1
#datasetChi (after roll)= [ istate , ip , cuts ]
#statedata = [ ip , icut ]

def PrintFitMassSetToFile(dataset,datasetChi,thisMomList,thisStateList,thisFitR,AddDict={}):
    dataset = np.rollaxis(np.rollaxis(dataset,1),2)
    datasetChi = np.rollaxis(np.rollaxis(datasetChi,1),2)
    gammadir = outputdir[0]+'Mass/fits/'
    mkdir_p(gammadir)
    for thisstate,statedata,statedataChi in zip(thisStateList,dataset,datasetChi):
        filename = (thisstate+'Mass')
        PrintFitMassToFile(statedata,statedataChi,thisstate,gammadir,filename,thisMomList,thisFitR,AddDict=AddDict)


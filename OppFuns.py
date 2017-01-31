#!/usr/bin/env python

from MiscFuns import *
from FitParams import StateParList
from FFParams import *
from Params import *
from SetLists import GetTsinkSmLists
from CombParams import CombList
from copy import copy

def DoubSingList(listin):
    return ['doub'+il for il in listin] + ['sing'+il for il in listin]

def DoubSingCmplxList(listin):
    return ['doub'+il+'cmplx' for il in listin] + ['sing'+il+'cmplx' for il in listin]

                
def CreateOppDir(Opp):
    if Opp in ['Mass','twopt']: return Opp+'/'
    thisSplitOpp = SplitOpp(Opp)
    contents = thisSplitOpp.keys()
    iscmplx = ''
    DS = ''
    prefolder = ''
    Proj = ''
    if 'Gamma' in contents:
        igamma = thisSplitOpp['Gamma']
        if 'Der' in contents:
            igamma = igamma + thisSplitOpp['Gamma']
    else:
        raise  IOError('Invalid Opp string '+Opp)            
    if 'DS' in contents:
        DS = thisSplitOpp['DS']
    if 'Proj' in contents:
        Proj = thisSplitOpp['Proj']
    else:
        raise  IOError('Invalid Opp string '+Opp)
    if 'Run' in contents: iscmplx = 'cmplx'
    if 'Top' in contents: prefolder = 'Top/'
    return prefolder+'/'+igamma + '/' +DS+Proj + iscmplx + '/'

def SetOpps(AllList):
    [Extra,thisOppList,thisProjList,thisDSList,RunList] = set([]),set([]),set([]),set([]),set([])
    iscmplx = ''
    TopRun = False
    for iOpp in AllList:        
        if iOpp in ['Mass','twopt','P4giDi','doubP4giDi','singP4giDi']: 
            Extra = Extra.union([iOpp])
            continue
        SplitOpps = SplitOpp(iOpp)
        contents = SplitOpps.keys()
        if 'Gamma' in contents:
            if 'Der' in contents:
                thisOppList = thisOppList.union([SplitOpps['Gamma']+SplitOpps['Der']])
            else:
                thisOppList = thisOppList.union([SplitOpps['Gamma']])            
        if 'DS' in contents:
            thisDSList = thisDSList.union([SplitOpps['DS']])
        else:
            thisDSList = thisDSList.union(['UmD'])
        if 'Proj' in contents:
            thisProjList = thisProjList.union([SplitOpps['Proj']])
        if 'Run' in contents: iscmplx = 'cmplx'
        if 'Top' in contents: TopRun=True
    return sorted(Extra),sorted(thisOppList),sorted(thisDSList),sorted(thisProjList),'real, '+iscmplx,TopRun

def Wipe2pt(thisoutputdir,tvarlist=[],ismlist=[],jsmlist=[],thisMomList=RunMomList,tsrclist = PoFtsourceList):
    thistvarlist = ['PoF'+str(PoFShifts)+itvar for itvar in tvarlist]
    thistvarlist += ['CM'+itvar for itvar in tvarlist]
    xmlMomList = map(qstrTOqcond,thisMomList)
    for iflag in ['cfun/twopt','Mass']:
        for ip in xmlMomList:
            thisdir = thisoutputdir+iflag+MakeMomDir(ip)
            for itvar in thistvarlist:
                ifile = thisdir+itvar+'LREM'+ip+'.xml'
                if os.path.isfile(ifile): os.remove(ifile)
                ifile = thisdir+'boots/'+itvar+'LREM'+ip+'.boot.p'
                if os.path.isfile(ifile): os.remove(ifile)
                # for istate in GetStateSet(itvar):
                for istate in StateSet:
                    ifile = thisdir+'state'+istate+itvar+iflag.replace('cfun/','')+ip+'.xml'
                    if os.path.isfile(ifile): os.remove(ifile)
                    ifile = thisdir+'boots/state'+istate+itvar+iflag.replace('cfun/','')+ip+'.boot.p'
                    if os.path.isfile(ifile): os.remove(ifile)
            for itsrc in tsrclist:
                for ism in ismlist:
                    for jsm in jsmlist:
                        ifile = thisdir+'tsrc'+itsrc+'ism'+ism+'jsm'+jsm+iflag.replace('cfun/','')+ip+'.xml'
                        if os.path.isfile(ifile): os.remove(ifile)
                        ifile = thisdir+'boots/tsrc'+itsrc+'ism'+ism+'jsm'+jsm+iflag.replace('cfun/','')+ip+'.boot.p'
                        if os.path.isfile(ifile): os.remove(ifile)
    

def WipeSet(thisoutputdir,thisGammaList,setlist,thisMomList=RunMomList,filepref=''):
    xmlMomList = map(qstrTOqcond,thisMomList)
    for igamma in thisGammaList:
        thisdir = thisoutputdir+CreateOppDir(igamma)
        for iset in setlist:
            for ip in xmlMomList:
                ifile = thisdir+MakeMomDir(ip)+filepref+iset+igamma+ip+'.xml'
                if os.path.isfile(ifile): os.remove(ifile)
                ifile = thisdir+MakeMomDir(ip)+'boots/'+filepref+iset+igamma+ip+'.boot.p'
                if os.path.isfile(ifile): os.remove(ifile)

            
# def WipeSet(outputdir[0],thisGammaList,tlist=[],treveclist=[],statelist=[],revectodtlist=[],todtlist=[],smlist=[],filepref=''):
#     for igamma in thisGammaList:
#         thisdir = outputdir[0]+CreateOppDir(igamma) + filepref        
#         if len(tlist) > 0 or len(treveclist) > 0:
#             for it in treveclist:
#                 for istate in statelist:
#                     for itodt in revectodtlist:
#                         ifile = thisdir+'tsink'+it+'state'+istate+itodt+igamma+'.txt'
#                         ifb = thisdir+'boots/tsink'+it+'state'+istate+itodt+igamma+'.boot.txt'
#                         # print ifile
#                         # print ifb
#                         if os.path.isfile(ifile): os.remove(ifile)
#                         if os.path.isfile(ifb): os.remove(ifb)
#             for it in tlist:
#                 for istate in statelist:
#                     for itodt in todtlist:
#                         ifile = thisdir+'tsink'+it+'state'+istate+itodt+igamma+'.txt'
#                         ifb = thisdir+'boots/tsink'+it+'state'+istate+itodt+igamma+'.boot.txt'
#                         # print ifile
#                         # print ifb
#                         if os.path.isfile(ifile): os.remove(ifile)
#                         if os.path.isfile(ifb): os.remove(ifb)
#                 for ism in smlist:
#                     ifile = thisdir+'tsink'+it+'sm'+ism+igamma+'.txt'
#                     ifb = thisdir+'boots/tsink'+it+'sm'+ism+igamma+'.boot.txt'
#                     # print ifile
#                     # print ifb
#                     if os.path.isfile(ifile): os.remove(ifile)
#                     if os.path.isfile(ifb): os.remove(ifb)
#         else:
#             for istate in statelist:
#                 for itodt in todtlist:
#                     ifile = thisdir+'state'+istate+itodt+igamma+'.txt'
#                     ifb = thisdir+'boots/state'+istate+itodt+igamma+'.boot.txt'
#                     # print ifile
#                     # print ifb
#                     if os.path.isfile(ifile): os.remove(ifile)
#                     if os.path.isfile(ifb): os.remove(ifb)
#                 for itodt in revectodtlist:
#                     ifile = thisdir+'state'+istate+itodt+igamma+'.txt'
#                     ifb = thisdir+'boots/state'+istate+itodt+igamma+'.boot.txt'
#                     # print ifile
#                     # print ifb
#                     if os.path.isfile(ifile): os.remove(ifile)
#                     if os.path.isfile(ifb): os.remove(ifb)
#             for ism in smlist:
#                 ifile = thisdir+'sm'+ism+igamma+'.txt'
#                 ifb = thisdir+'boots/sm'+ism+igamma+'.boot.txt'
#                 # print ifile
#                 # print ifb
#                 if os.path.isfile(ifile): os.remove(ifile)
#                 if os.path.isfile(ifb): os.remove(ifb)

def WipeSF(thisoutputdir,thisGammaList,RunName,OoT,statelist=[],todtlist=[],smlist=[],tsinklist=['']):
    for igamma in thisGammaList:
        if igamma == 'twopt':
            thisdir = thisoutputdir+'cfun/'+CreateOppDir(igamma)+RunName+'/'
            thisParList = StateParList[OoT]['C2']
            thistsinklist = ['']
        else:
            thisdir = thisoutputdir+CreateOppDir(igamma)+RunName+'/'
            thisParList = StateParList[OoT]['C3']
            thistsinklist = tsinklist
        for ip in thisParList:
            if igamma == 'twopt' and 'Two' in OoT:
                ifile = thisdir+igamma+ip+'.txt'
                ifb = thisdir+'boots/'+igamma+ip+'.boot.txt'
                if os.path.isfile(ifile): os.remove(ifile)
                if os.path.isfile(ifb): os.remove(ifb)
            for its in thistsinklist:
                for istate in statelist:
                    for itodt in todtlist:
                        ifile = thisdir+its+'state'+istate+itodt+igamma+ip+'.txt'
                        ifb = thisdir+'boots/'+its+'state'+istate+itodt+igamma+ip+'.boot.txt'
                        if os.path.isfile(ifile): os.remove(ifile)
                        if os.path.isfile(ifb): os.remove(ifb)
                for ism in smlist:
                    ifile = thisdir+its+'sm'+ism+igamma+ip+'.txt'
                    ifb = thisdir+'boots/'+its+'sm'+ism+igamma+ip+'.boot.txt'
                    if os.path.isfile(ifile): os.remove(ifile)
                    if os.path.isfile(ifb): os.remove(ifb)

def WipeSFSet(thisoutputdir,thisGammaList,RunName,OoT,setlist=[]):
    for igamma in thisGammaList:
        if igamma == 'twopt':
            thisdir = thisoutputdir+'cfun/'+CreateOppDir(igamma)+RunName+'/'
            thisParList = StateParList[OoT]['C2']
            thissetlist = GetTsinkSmLists(setlist)[1]
        else:
            thisdir = thisoutputdir+CreateOppDir(igamma)+RunName+'/'
            thisParList = StateParList[OoT]['C3']
            thissetlist = setlist
        for ip in thisParList:
            if igamma == 'twopt' and 'Two' in OoT:
                ifile = thisdir+igamma+ip+'.txt'
                ifb = thisdir+'boots/'+igamma+ip+'.boot.txt'
                if os.path.isfile(ifile): os.remove(ifile)
                if os.path.isfile(ifb): os.remove(ifb)
            else:
                for iset in thissetlist:
                    ifile = thisdir+iset+igamma+ip+'.txt'
                    ifb = thisdir+'boots/'+iset+igamma+ip+'.boot.txt'
                    if os.path.isfile(ifile): os.remove(ifile)
                    if os.path.isfile(ifb): os.remove(ifb)

def SplitOpp(All):
    outputDict = OrderedDict()
    Split,OrdSplit,contents = [],[],[]
    if any([igamma in All for igamma in GammaSet + ['gi']]):
        contents.append('Gamma')
        gammalen,thisgamma = 0,''
        for igamma in GammaSet+['gi']:
            if igamma in All:
                if len(igamma) > gammalen:
                    thisgamma = igamma
                gammalen = len(igamma)
        # Split.append(thisgamma)
        outputDict['Gamma'] = thisgamma
    if 'D' in All:
        for iDer in DerSet+['Di']:
            if iDer in All: outputdict['Der'] = iDer
    if any([iDS in All for iDS in DefDSList]) or any([icomb in All for icomb in CombList]):
        for iDS in DefDSList+CombList:
            if iDS in All and 'Iso'+iDS not in All:
                outputDict['DS'] = iDS
    if 'P4' in All or 'P3' in All:
        if 'P4' in All: outputDict['Proj'] = 'P4'
        if 'P3' in All: outputDict['Proj'] = 'P3'
    if 'cmplx' in All:
        outputDict['Run'] = 'cmplx'
    if 'Top' in All:
        outputDict['Top'] = 'Top'
    return outputDict

def PrintOpps(AllList):
    Extra,thisGS,thisDSS,thisProjS,RunRS,TopRun = SetOpps(AllList)
    print 'All Opperators: \n'+'\n'.join(thisGS)
    print ''
    print 'Projectors: '+', '.join(thisProjS)
    # print 'DS: '+', '.join(thisDSS)
    if TopRun:
        print 'Run: ' +RunRS
    else:
        print 'Run: ' +RunRS + ', TopCharge '        
    print 'Extras: ' + ', '.join(Extra)
    print ''

#input arguments can be:
# specific operators (e.g. 'doubP4g4')
# no doub/sing operators give up-down results (e.g. 'P4g4')
# current types (e.g. 'Vector')
# 'giDi' for momenta fraction
# 'twopt' or 'mass' for mass analysis
# 'SmallSet' for small set defined below


def CreateGammaList(thislist,twopt=False):
    if len(thislist) == 0:
        print 'No Gamma Inputted, using whole set (see Params.py DefGammaList)'
        GLout = DefCombGammaList
    else:
        GLout = []
        for ig in thislist:
            if 'NoTensor' in ig :
                for icurr in CurrTypes:
                    if icurr in 'Tensor': continue
                    GLout += DoubSingList(CurrOpps[icurr])
                    GLout += DoubSingCmplxList(CurrOpps[icurr])      
            elif ig in CurrTypes:
                GLout += DoubSingList(CurrOpps[ig])
                GLout += DoubSingCmplxList(CurrOpps[ig])      
            elif ig in [icurr + 'OnlyGamma' for icurr in CurrTypes]:
                GLout += CurrOpps[ig.replace('OnlyGamma','')]
                GLout += [icurr+'cmplx' for icurr in CurrOpps[ig.replace('OnlyGamma','')]]
            elif 'SmallSet' in ig :
                GLout += DoubSingList(['P4g4','P3g3g5','P4I','P3g1g2','P4giDi'])
            elif 'OnlyDS' in ig :
                GLout += DefGammaList
            elif 'OnlyGamma' in ig :
                GLout += DefNoDSGammaList
            elif ig in DerCurrTypes:
                GLout += DoubSingList(['P4'+ig])
            elif ig in DefCombGammaList:
                GLout += [ig]
            elif ig in [(igamma+'Top').replace('cmplxTop','Topcmplx') for igamma in DefCombGammaList]:
                GLout += [ig]
            elif ig in ['twopt','Mass']:
                GLout += ['twopt']
            else:
                print 'Warning, opperator not found: ' , ig
    
    if twopt: GLout += ['twopt']
    GLout = DelDubs(GLout)
    PrintOpps(GLout)
    return GLout

#!/usr/bin/env python

from MiscFuns import *
from FitParams import StateParList
from FFParams import *
from Params import DefGammaList
from SetLists import GetTsinkSmLists
from copy import copy

def DoubSingList(listin):
    return ['doub'+il for il in listin] + ['sing'+il for il in listin]

def DoubSingCmplxList(listin):
    return ['doub'+il+'cmplx' for il in listin] + ['sing'+il+'cmplx' for il in listin]

                
def CreateOppDir(Opp):
    if Opp in ['Mass','twopt']: return Opp+'/'
    thisSplitOpp,contents = SplitOpp(Opp)
    print thisSplitOpp,contents
    count,iscmplx = 0,''
    DS = ''
    if 'Gamma' in contents:
        igamma = thisSplitOpp[count]
        count += 1
        if 'Der' in contents:
            igamma = igamma + thisSplitOpp[count]
            count += 1
    if 'DS' in contents:
        DS = thisSplitOpp[count]
        count += 1
    if 'Proj' in contents:
        Proj = thisSplitOpp[count]
    if 'Run' in contents: iscmplx = 'cmplx'
    try:
        return igamma + '/' +DS+Proj + iscmplx + '/'
    except:
        raise  IOError('Invalid Opp string '+Opp)

def SetOpps(AllList):
    [Extra,thisOppList,thisProjList,thisDSList,RunList] = set([]),set([]),set([]),set([]),set([])
    iscmplx = ''
    for iOpp in AllList:        
        if iOpp in ['Mass','twopt','P4giDi','doubP4giDi','singP4giDi']: 
            Extra = Extra.union([iOpp])
            continue
        SplitOpps,contents = SplitOpp(iOpp)
        count = 0
        if 'Gamma' in contents:
            if 'Der' in contents:
                thisOppList = thisOppList.union([SplitOpps[count]+SplitOpps[count+1]])
                count += 1
            else:
                thisOppList = thisOppList.union([SplitOpps[count]])
            count += 1
            
        if 'DS' in contents:
            thisDSList = thisDSList.union([SplitOpps[count]])
            count += 1
        else:
            thisDSList = thisDSList.union(['UmD'])
        if 'Proj' in contents:
            thisProjList = thisProjList.union([SplitOpps[count]])
        if 'Run' in contents: iscmplx = 'cmplx'
    return sorted(Extra),sorted(thisOppList),sorted(thisDSList),sorted(thisProjList),'real, '+iscmplx

def Wipe2pt(outputdir,statelist=[],todtlist=[],smlist=[],thisMomList=RunMomList):
    thistodtlist = ['PoF'+str(PoFShifts)+itodt for itodt in todtlist]
    thistodtlist += ['CM'+itodt for itodt in todtlist]
    xmlMomList = map(qstrTOqcond,thisMomList)
    for iflag in ['cfuns/twopt','Mass']:
        for ip in xmlMomList:
            thisdir = outputdir+iflag+MakeMomDir(ip)
            for itodt in thistodtlist:
                ifile = thisdir+itodt+'LREM.xml'
                if os.path.isfile(ifile): os.remove(ifile)
                for istate in statelist:
                    ifile = thisdir+'state'+istate+itodt+iflag.replace('cfuns/','')+ip+'.xml'
                    if os.path.isfile(ifile): os.remove(ifile)
            for ism in smlist:
                ifile = thisdir+'sm'+ism+iflag.replace('cfuns/','')+ip+'.xml'
                if os.path.isfile(ifile): os.remove(ifile)
    

def WipeSet(outputdir,thisGammaList,setlist,thisMomList=RunMomList,filepref=''):
    xmlMomList = map(qstrTOqcond,thisMomList)
    for igamma in thisGammaList:
        thisdir = outputdir+CreateOppDir(igamma)
        for iset in setlist:
            for ip in xmlMomList:
                ifile = thisdir+MakeMomDir(ip)+filepref+iset+igamma+ip+'.xml'
                if os.path.isfile(ifile): os.remove(ifile)

            
# def WipeSet(outputdir,thisGammaList,tlist=[],treveclist=[],statelist=[],revectodtlist=[],todtlist=[],smlist=[],filepref=''):
#     for igamma in thisGammaList:
#         thisdir = outputdir+CreateOppDir(igamma) + filepref        
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

def WipeSF(outputdir,thisGammaList,RunName,OoT,statelist=[],todtlist=[],smlist=[],tsinklist=['']):
    for igamma in thisGammaList:
        if igamma == 'twopt':
            thisdir = outputdir+'cfuns/'+CreateOppDir(igamma)+RunName+'/'
            thisParList = StateParList[OoT]['C2']
            thistsinklist = ['']
        else:
            thisdir = outputdir+CreateOppDir(igamma)+RunName+'/'
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

def WipeSFSet(outputdir,thisGammaList,RunName,OoT,setlist=[]):
    for igamma in thisGammaList:
        if igamma == 'twopt':
            thisdir = outputdir+'cfuns/'+CreateOppDir(igamma)+RunName+'/'
            thisParList = StateParList[OoT]['C2']
            thissetlist = GetTsinkSmLists(setlist)[1]
        else:
            thisdir = outputdir+CreateOppDir(igamma)+RunName+'/'
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
    Split,OrdSplit,contents = [],[],[]
    SearchFlags = [['g','I'],['D'],['d','s'],['P'],['c']]
    gdone = False
    if any([iDS in All for iDS in DefDSList]) or any([icomb in All for icomb in DefCombList]):
        contents.append('DS')
        for iDS in DefDSList+DefCombList:
            if iDS in All:
                Split.append(iDS)
    if 'cmplx' in All:
        contents.append('Run')
        Split.append('cmplx')
    if 'P4' in All or 'P3' in All:
        contents.append('Proj')
        if 'P4' in All: Split.append('P4')
        if 'P3' in All: Split.append('P3')
    if 'D' in All:
        contents.append('Der')
        for iDer in DerSet+['Di']:
            if iDer in All: Split.append(iDer)
    if any([igamma in All for igamma in GammaSet]):
        contents.append('Gamma')
        gammalen,thisgamma = 0,''
        for igamma in GammaSet:
            if igamma in All:
                if len(igamma) > gammalen:
                    thisgamma = igamma
                gammalen = len(igamma)
        Split.append(thisgamma)
        
    # for ichar,char in enumerate(All):
    #     if char in ['d','s']:
    #         Split.append(All[ichar:ichar+4])
    #         contents.append('DS')
    #     if char == 'c':            
    #         Split.append(All[ichar:ichar+5])
    #         contents.append('Run')
    #     elif char == 'P':
    #         Split.append(All[ichar:ichar+2])
    #         contents.append('Proj')
    #     elif char == 'D':
    #         Split.append(All[ichar:ichar+2])
    #         contents.append('Der')
    #     elif char == 'g' and not gdone:
    #         if All[ichar-1] == 'n': continue
    #         contents.append('Gamma')
    #         cutlen = 2
    #         gdone = True
    #         if len(All[ichar+2:]) > 1:
    #             for char2 in All[ichar+2:]:
    #                 if char2 == 'g':
    #                     cutlen += 2
    #                 else: 
    #                     break
    #         Split.append(All[ichar:ichar+cutlen])
    #     elif char == 'I' and not gdone:
    #         Split.append(All[ichar])
    #         gdone == True
    #         contents.append('Gamma')
    for iflag in SearchFlags:
        for ichar,(char,icont) in enumerate(zip(Split,contents)):
            if char[0] in iflag:
                OrdSplit.append(char)
                break
    return OrdSplit,contents

def PrintOpps(AllList):
    Extra,thisGS,thisDSS,thisProjS,RunRS = SetOpps(AllList)
    print 'All Opperators: \n'+'\n'.join(thisGS)
    print ''
    print 'Projectors: '+', '.join(thisProjS)
    print 'DS: '+', '.join(thisDSS)
    print 'Run: ' +RunRS
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
        GLout = DefGammaList
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
            elif 'SmallSet' in ig :
                GLout += DoubSingList(['P4g4','P3g3g5','P4I','P3g1g2','P4giDi'])
            elif ig in DerCurrTypes:
                GLout += DoubSingList(['P4'+ig])
            elif ig in DefGammaList:
                GLout += [ig]
            elif ig in ['twopt','Mass']:
                GLout += ['twopt']
            else:
                print 'Warning, opperator not found: ' , ig
    
    if twopt: GLout += ['twopt']
    GLout = DelDubs(GLout)
    PrintOpps(GLout)
    return GLout

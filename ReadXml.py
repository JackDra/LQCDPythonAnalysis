#!/usr/bin/env python

import xmltodict
from XmlFuns import *
from XmlFormatting import *
import os
from BootTest import BootStrap1
from Params import *
from FitParams import *
import cPickle as pickle
import warnings


def ReadPickleBoot(filein):
    try:
        with open(filein , "rb") as pfile:
            dataout = pickle.load( pfile )
    except:
        print 'Reading pickle file fail: ' + filein
        dataout = {}
    return dataout


def ReadXmlDict(filein,Boot=True):
    try:
        with open(filein,'r') as f:
            xmldata = RecFTDAS(xmltodict.parse(f.read()))
        if Boot:
            bootfile = xmldata[xmldata.keys()[0]]['Boots']
        else:
            bootfile = None
    except:
        if Debug: print 'Reading xml file fail: ' + filein
        xmldata = {}
        bootfile = 'NoBootDir'
    return xmldata,bootfile

def ReadXmlAndPickle(filein):
    if Debug: print 'reading: ' + filein
    xmldata,bootfile = ReadXmlDict(filein)
    firstkey = xmldata.keys()
    if len(firstkey) > 0:
        xmldata[firstkey[0]]['Boots'] = ReadPickleBoot(bootfile)
        if 'Info' not in xmldata[firstkey[0]].keys():
            xmldata[firstkey[0]]['Info'] = CreateDefaultInfodict()
        # else:
            # if not ('nconfig' in xmldata[firstkey[0]]['Info'].keys() or any(['qsqrd' in ixm for ixm in xmldata[firstkey[0]]['Info'].keys()])):
            #     xmldata[firstkey[0]]['Info'] = CreateDefaultInfodict()
    return xmldata,bootfile

def CheckNconfFile(filein):
    if Debug: print 'reading: ' + filein
    if not os.path.isfile(filein): return 'File Missing'
    mprint('reading: ' + filein)
    Nconf = 10e17
    with open(filein,'r') as thisfile:
        for line in thisfile:
            strpline = line.strip()
            if '</Info>' in strpline:
                return Nconf
            elif '<nMeasure>' in strpline and '</nMeasure>' in strpline:
                thisNconf = int(strpline.replace('<nMeasure>','').replace('</nMeasure>',''))
                if thisNconf > 0: Nconf = min(Nconf,thisNconf)
    if Nconf > 10e10 or Nconf < 0:
        return 'Dep'
    else:
        return Nconf


def CheckMomFile(filein,nconftest = False):
    mprint(filein)
    if nconftest != False:
        filenconf = CheckNconfFile(filein)
        mprint(filein)
        mprint('has nconf = ' + str(filenconf))
        if Debug:
            print filein
            print 'has nconf = ' + str(filenconf)
        if nconftest == filenconf:
            return True
        else:
            return False
    else:
        if not os.path.isfile(filein): return False
        if not DoContentsCheck: return True
        xmldata,bootfile = ReadXmlDict(filein)
        if GetqcondFromFilename(filein) in xmldata.keys()[0]:
            return True
        else:
            return False


                    
def GetDPFitValue(SearchSet,iFF,thisCurr,thiskappa=str(kappa)):
    if thiskappa == None: thiskappa = kappa 
    FFdir = outputdir[0].replace('k'+str(kappa),'k'+str(thiskappa)) +'/FormFactors/DPfits/'
    filename = FFdir +thisCurr
    if Debug:
        print filename+'.xml'
        if not os.path.isfile(filename+'.xml'):
            print 'Warrning, file not found' , filename+'.xml'
    if not os.path.isfile(filename+'.xml'): return [],[]
    data = ReadXmlDict(filename+'.xml')[0]
    thisDPpAvg = []
    thisDPpStd = []
    if Debug: print data.keys(),'DP_Fits'
    if Debug: print data[data.keys()[0]].keys(),'Values'
    if Debug: print iFF, SearchSet
    if iFF in data['DP_Fits']['Values'].keys():
        # if any([SearchSet in iset  for iset in data['DP_Fits']['Values'][iFF]['Fzero'].keys()]):
        if SearchSet in data['DP_Fits']['Values'][iFF]['Fzero'].keys() and SearchSet in data['DP_Fits']['Values'][iFF]['mEM'].keys() :
            thisDPpAvg.append(data['DP_Fits']['Values'][iFF]['Fzero'][SearchSet]['Avg'])
            thisDPpAvg.append(data['DP_Fits']['Values'][iFF]['mEM'][SearchSet]['Avg'])
            thisDPpAvg.append(data['DP_Fits']['Values'][iFF]['Radius'][SearchSet]['Avg'])
            thisDPpStd.append(data['DP_Fits']['Values'][iFF]['Fzero'][SearchSet]['Std'])
            thisDPpStd.append(data['DP_Fits']['Values'][iFF]['mEM'][SearchSet]['Std'])
            thisDPpStd.append(data['DP_Fits']['Values'][iFF]['Radius'][SearchSet]['Std'])
    return thisDPpAvg,thisDPpStd


def ReadTopFileKappas(thiskappalist,filedirlist,iset,thisMomList=RunMomList,Wein=False):
    dictout = OrderedDict()
    for ikappa,ifile in zip(thiskappalist,filedirlist):
        dictout[ikappa] = ReadTopFile(ifile,iset,thisMomList=thisMomList,Wein=Wein)
    return dictout
        
def ReadTopFile(filedir,iset,thisMomList=RunMomList,Wein=False):
    dictout = OrderedDict()
    if Wein:TopOrWein ='Wein'
    else: TopOrWein = 'Top'
    dictout['RF'] = OrderedDict()
    dictout['NNQ'] = OrderedDict()
    dictout['cfun'] = OrderedDict()
    for thismom in thisMomList:
        ip = qstrTOqcond(thismom)
        readfile = filedir+TopOrWein+'/Alpha/'+MakeMomDir(ip)+iset+ip
        if Debug: print 'Reading TopCharge :' ,readfile
        if os.path.isfile(readfile+'.xml'):
            dictout['RF'][thismom] = ReadXmlAndPickle(readfile+'.xml')[0][ip]
        thisreadfile = filedir+TopOrWein+'/cfun/twopt/'+MakeMomDir(ip)+iset+'twopt'+ip
        if os.path.isfile(thisreadfile+'.xml'):
            dictout['cfun'][thismom] = ReadXmlAndPickle(thisreadfile+'.xml')[0][ip]
        thisreadfile = readfile.replace(TopOrWein+'/Alpha/',TopOrWein+'/cfun/NNQ/')
        # readfile = filedir+'Top/NNQ/'+MakeMomDir(ip)+iset+ip
        if os.path.isfile(thisreadfile+'.xml'):
            dictout['NNQ'][thismom] = ReadXmlAndPickle(thisreadfile+'.xml')[0][ip]
    return dictout


## dictout = { 'RF' , imom , Info/ Boots: itflow , tsink }  
def ReadAlphaFile(filedir,iset,thisMomList=RunMomList,Wein=False):
    dictout = OrderedDict()
    if Wein:TopOrWein ='Wein'
    else: TopOrWein = 'Top'
    dictout['RF'] = OrderedDict()
    for thismom in thisMomList:
        ip = qstrTOqcond(thismom)
        readfile = filedir+TopOrWein+'/Alpha/'+MakeMomDir(ip)+iset+ip
        if Debug: print 'Reading TopCharge :' ,readfile+'.xml'
        if os.path.isfile(readfile+'.xml'): dictout['RF'][thismom] = ReadXmlAndPickle(readfile+'.xml')[0][ip]
    return dictout


def ReadAlphaFitFileKappas(thiskappalist,filedirlist,iset,thisMomList=RunMomList,Wein=False):
    dictout = OrderedDict()
    for ikappa,ifile in zip(thiskappalist,filedirlist):
        dictout[ikappa] = ReadAlphaFitFile(ifile,iset,thisMomList=thisMomList,Wein=Wein)
    return dictout
        
## dictout = { imom , Info/ Boots: itflow , fitr }  
def ReadAlphaFitFile(filedir,iset,thisMomList=RunMomList,givefile=False,Wein=False):
    dictout = OrderedDict()
    if Wein:TopOrWein ='Wein'
    else: TopOrWein = 'Top'
    for thismom in thisMomList:
        ip = qstrTOqcond(thismom)
        readfile = filedir+TopOrWein+'/Alpha/Fits/'+MakeMomDir(ip)+iset+ip
        if Debug: print 'Reading TopCharge Fits :' ,readfile+'.xml'
        if os.path.isfile(readfile+'.xml'): dictout[thismom] = ReadXmlAndPickle(readfile+'.xml')[0][ip]
    if givefile:
        return dictout,readfile
    else:
        return dictout

##Also works for cfuns##
##xmlinput = { Ratio_Factor , Boots/Values , thismomlist , tlist } 
##outputdict = { thismom , [tVals] / [Vals] / [Valserr] / [Boot] bs }
def ReadRFFile(filedir,filename,thisMomList=RunMomList):
    renorm = GetRenorm(filename)
    dictout = {}
    TopFile = 'Top' in filename.replace('Topk'+str(kappa),'') or 'Wein' in filename.replace('Weink'+str(kappa),'')
    for thismom in thisMomList:
        ip = qstrTOqcond(thismom)
        readfile = filedir+MakeMomDir(ip)+filename.replace('.xml',ip+'.xml')
        # if os.path.isfile(readfile):
        data = ReadXmlAndPickle(readfile)[0]
        if len(data.keys()) > 0:
            data = data[data.keys()[0]]
            if 'Boots' in data.keys():
                bootdata = data['Boots']
                dictout[thismom] = OrderedDict()
                dictout[thismom]['Info'] = data['Info']
                if TopFile:
                    for iflow,flowdata in bootdata.iteritems():
                        dictout[thismom][iflow] = OrderedDict()
                        dictout[thismom][iflow]['tVals'] = map(untstr,flowdata.keys())
                        dictout[thismom][iflow]['Boot'] = []
                        dictout[thismom][iflow]['Vals'] = []
                        dictout[thismom][iflow]['Valserr'] = []
                        for bootlist in flowdata.itervalues():
                            dictout[thismom][iflow]['Boot'].append(BootStrap1(nboot,0))
                            if any(np.isnan(bootlist)): warnings.warn('Warning: NaN found in a file, please dont run this set: '+ readfile)
                            dictout[thismom][iflow]['Boot'][-1].values = np.array([iboot*renorm for iboot in bootlist])                            
                            dictout[thismom][iflow]['Boot'][-1].Stats()
                            dictout[thismom][iflow]['Vals'].append(dictout[thismom][iflow]['Boot'][-1].Avg)
                            dictout[thismom][iflow]['Valserr'].append(dictout[thismom][iflow]['Boot'][-1].Std)

                else:
                    dictout[thismom]['tVals'] = map(untstr,bootdata.keys())
                    dictout[thismom]['Boot'] = []
                    dictout[thismom]['Vals'] = []
                    dictout[thismom]['Valserr'] = []
                    for bootlist in bootdata.itervalues():
                        dictout[thismom]['Boot'].append(BootStrap1(nboot,0))
                        if any(np.isnan(bootlist)): warnings.warn('Warning: NaN found in a file, please dont run this set: '+ readfile)
                        dictout[thismom]['Boot'][-1].values = np.array([iboot*renorm for iboot in bootlist])
                        dictout[thismom]['Boot'][-1].Stats()
                        dictout[thismom]['Vals'].append(dictout[thismom]['Boot'][-1].Avg)
                        dictout[thismom]['Valserr'].append(dictout[thismom]['Boot'][-1].Std)
                        
            elif 'Avg' in data.keys():
                bootdata = data['Values']
                dictout[thismom] = {}
                dictout[thismom]['Info'] = data['Info']
                if TopFile:
                    for iflow,flowdata in bootdata.iteritems():
                        dictout[thismom][iflow] = OrderedDict()
                        dictout[thismom][iflow]['tVals'] = map(untstr,flowdata.keys())
                        dictout[thismom][iflow]['Vals'] = []
                        dictout[thismom][iflow]['Valserr'] = []
                        for tdata in flowdata.itervalues():
                            dictout[thismom][iflow]['Vals'].append(tdata['Avg']*renorm)
                            dictout[thismom][iflow]['Valserr'].append(tdata['Std'])
                else:
                    dictout[thismom]['tVals'] = map(untstr,bootdata.keys())
                    dictout[thismom]['Vals'] = []
                    dictout[thismom]['Valserr'] = []
                    for tdata in bootdata.itervalues():
                        dictout[thismom]['Vals'].append(tdata['Avg']*renorm)
                        dictout[thismom]['Valserr'].append(tdata['Std'])
                    
    return dictout
        
                    

def ReadFitFile(filedir,filename,thisMomList=RunMomList):
    dictout = {}
    for thismom in thisMomList:
        ip = qstrTOqcond(thismom)
        readfile = filedir+MakeMomDir(ip)+filename.replace('.xml',ip+'.xml')
        # if os.path.isfile(readfile):
        dictout[thismom] = {}
        data = ReadXmlAndPickle(readfile)[0]
        if len(data.keys()) > 0:
            data = data[data.keys()[0]]
            if 'Boots' in data.keys():
                bootdata = data['Boots']
                if 'Top' in readfile.replace('Top'+kappalist[0],'') or 'Wein' in readfile.replace('Wein'+kappalist[0],''):
                    dictout[thismom] = OrderedDict()
                    for iflow,flowdata in bootdata.iteritems():
                        dictout[thismom][iflow] = OrderedDict()
                        for icut,cutdata in flowdata.iteritems():
                            dictout[thismom][iflow][icut] = OrderedDict()
                            dictout[thismom][iflow][icut]['Boot'] = BootStrap1(nboot,0)
                            dictout[thismom][iflow][icut]['Boot'].values = np.array(cutdata)
                            dictout[thismom][iflow][icut]['Boot'].Stats()
                            dictout[thismom][iflow][icut]['Avg'] = dictout[thismom][iflow][icut]['Boot'].Avg
                            dictout[thismom][iflow][icut]['Std'] = dictout[thismom][iflow][icut]['Boot'].Std
                            dictout[thismom][iflow][icut]['Chi'] = data['Values'][iflow][icut]['Chi']
                else:
                    for icut,cutdata in bootdata.iteritems():
                        dictout[thismom][icut] = OrderedDict()
                        dictout[thismom][icut]['Boot'] = BootStrap1(nboot,0)
                        dictout[thismom][icut]['Boot'].values = np.array(cutdata)
                        dictout[thismom][icut]['Boot'].Stats()
                        dictout[thismom][icut]['Avg'] = dictout[thismom][icut]['Boot'].Avg
                        dictout[thismom][icut]['Std'] = dictout[thismom][icut]['Boot'].Std
                        dictout[thismom][icut]['Chi'] = data['Values'][icut]['Chi']
            else:
                dictout[thismom] = data
            dictout[thismom]['Info'] = data['Info']
    return dictout

##outputdict = { thismom , cutpar , tsinkrpar/tsinkval , Avg / Std / Chi / Boot (bs) }
def ReadSumFile(filedir,filename,thisMomList=RunMomList):
    dictout = {}
    for thismom in thisMomList:
        ip = qstrTOqcond(thismom)
        readfile = filedir+MakeMomDir(ip)+filename.replace('.xml',ip+'.xml')
        # if os.path.isfile(readfile):
        data = ReadXmlAndPickle(readfile)[0]
        if len(data.keys()) > 0:
            data = data[data.keys()[0]]
            if 'Boots' in data.keys():
                bootdata = data['Boots']
                dictout[thismom] = {}
                dictout[thismom]['Info'] = data['Info']
                for icut,cutdata in bootdata.iteritems():
                    dictout[thismom][icut] = {}
                    for it,tdata in cutdata.iteritems():
                        if 'tsink' in it:
                            dictout[thismom][icut][it] = {}
                            dictout[thismom][icut][it]['Boot'] = BootStrap1(nboot,0)
                            dictout[thismom][icut][it]['Boot'].values = np.array(tdata)
                            dictout[thismom][icut][it]['Boot'].Stats()
                            dictout[thismom][icut][it]['Avg'] = dictout[thismom][icut][it]['Boot'].Avg
                            dictout[thismom][icut][it]['Std'] = dictout[thismom][icut][it]['Boot'].Std
                        else:
                            for ir,rdata in tdata.iteritems():
                                thisfitr = FitFlagXmlToOld(it,ir)
                                dictout[thismom][icut][thisfitr] = {}
                                dictout[thismom][icut][thisfitr]['Boot'] = BootStrap1(nboot,0)
                                dictout[thismom][icut][thisfitr]['Boot'].values = np.array(rdata)
                                dictout[thismom][icut][thisfitr]['Boot'].Stats()
                                dictout[thismom][icut][thisfitr]['Avg'] = dictout[thismom][icut][thisfitr]['Boot'].Avg
                                dictout[thismom][icut][thisfitr]['Std'] = dictout[thismom][icut][thisfitr]['Boot'].Std
                                dictout[thismom][icut][thisfitr]['Chi'] = data['Values'][icut][it][ir]['Chi']

            else:
                valdata = data['Values']
                dictout[thismom] = {}
                dictout[thismom]['Info'] = data['Info']
                for icut,cutdata in bootdata.iteritems():
                    dictout[thismom][icut] = {}
                    for it,tdata in cutdata.iteritems():
                        if 'tsink' in it:
                            dictout[thismom][icut][it] = tdata
                        else:
                            for ir,rdata in tdata.iteritems():
                                thisfitr = FitFlagXmlToOld(it,ir)
                                dictout[thismom][icut][thisfitr] = rdata
    return dictout

## dataout = { Mass:Set/Avg/Std/Chi/Boot , FF#:qsqrd:Avg/Std/Boot , Chi:qsqrd}
def ReadFFFile(filename):
    dataout = False
    if '.txt' in filename: filename = filename.replace('.txt','.xml')
    # if os.path.isfile(filename):
    dataout = {}
    data = ReadXmlAndPickle(filename)[0]
    if len(data.keys()) > 0:
        data = data[data.keys()[0]]
        dataout = OrderedDict()
        dataout['Info'] = data['Info']
        if 'Mass' in data['Values'].keys():
            dataout['Mass'] = data['Values']['Mass']
            thismass = dataout['Mass']['Avg']
        else:
            thismass = 1.0
        dataout['Chi'] = OrderedDict()
        if 'Boots' in data.keys():
            for iq,qdata in data['Boots'].iteritems():
                dataout['Chi'][iq] = data['Values'][iq]['Chi']
                for iff,ffdata in qdata.iteritems():
                    if 'Chi' in iff: continue
                    if iff not in dataout.keys():
                        dataout[iff] = OrderedDict()
                    dataout[iff][iq] = OrderedDict()
                    dataout[iff][iq]['Boot'] = BootStrap1(nboot,0)
                    if F3Div2M and 'VectorTop' in filename and '3' in iff:
                        dataout[iff][iq]['Boot'].values = np.array(ffdata*hbarc/(2*thismass*hbarcdivlat)) ## in fm
                    else:
                        dataout[iff][iq]['Boot'].values = np.array(ffdata)
                    dataout[iff][iq]['Boot'].Stats()
                    dataout[iff][iq]['Avg'] = dataout[iff][iq]['Boot'].Avg
                    dataout[iff][iq]['Std'] =dataout[iff][iq]['Boot'].Std

        else:
            for iq,qdata in data['Values'].iteritems():
                for iff,ffdata in qdata.iteritems():
                    if iff not in dataout.keys():
                        dataout[iff] = OrderedDict()
                    dataout[iff][iq] = ffdata
    return dataout


## dataout = { Info:qsqrd#:nconf , Mass:Set/Avg/Std/Chi/Boot , qsqrd#:Avg/Std/Boot/Chi}
def ReadFFCombFile(filename):
    dataout = False
    if '.txt' in filename: filename = filename.replace('.txt','.xml')
    if os.path.isfile(filename):
        dataout = {}
        data = ReadXmlAndPickle(filename)[0]
        data = data[data.keys()[0]]
        dataout = OrderedDict()
        dataout['Info'] = data['Info']
        if 'Mass' in data['Values'].keys(): dataout['Mass'] = data['Values']['Mass']
        if 'Boots' in data.keys():
            for iq,qdata in data['Boots'].iteritems():
                dataout[iq] = OrderedDict()
                dataout[iq]['Chi'] = data['Values'][iq]['Chi']
                dataout[iq]['Boot'] = BootStrap1(nboot,0)
                dataout[iq]['Boot'].values = np.array(qdata.values)
                dataout[iq]['Boot'].Stats()
                dataout[iq]['Avg'] = dataout[iq]['Boot'].Avg
                dataout[iq]['Std'] = dataout[iq]['Boot'].Std
        else:
            dataout = data['Values']     
        return MakeFFCombLikeFF(dataout)
    return dataout

## dataout = { Mass:Set/Avg/Std/Chi/Boot , FF#:qsqrd:Avg/Std/Boot , Chi:qsqrd}
def MakeFFCombLikeFF(data):
    dataout = {}
    dataout['Info'] = data['Info']
    if 'Mass' in data.keys(): dataout['Mass'] = data['Mass']
    dataout['FF1'] = OrderedDict()
    dataout['Chi'] = OrderedDict()
    for iq in data.iterkeys():
        if 'qsqrd' not in iq: continue
        dataout['FF1'][iq] = {}
        dataout['FF1'][iq]['Avg'] = data[iq]['Avg']
        dataout['FF1'][iq]['Std'] = data[iq]['Std']
        dataout['Chi'][iq] = data[iq]['Chi']
    return dataout
        
##outputdict = { thismom , fitpar , 2corfitr , 3corcutr , Avg / Std / Chi / Boot (bs) }
## put ## as parameter
def ReadSFFile(filedir,filename,OneOrTwo='Two',thisMomList=RunMomList):
    dictout = {}
    if '.txt' in filename: filename = filename.replace('.txt','.xml')
    if 'twopt' in filename:
        twoptread = True
        thisCorr = 'C2'
    else:
        twoptread = False
        thisCorr = 'C3'
    for ipar in StateParList[OneOrTwo][thisCorr]:
        if ipar in ['m0','Dm'] and OneOrTwo == 'Two':
            thisfilename = re.sub('sm.*twopt','twopt',filename)
            thisfilename = re.sub('state.*twopt','twopt',thisfilename).replace('##',ipar)
        else:
            thisfilename = filename.replace('##',ipar)
        for thismom in thisMomList:
            ip = qstrTOqcond(thismom)
            readfile = filedir+MakeMomDir(ip)+thisfilename.replace('.xml',ip+'.xml')
            if 'Two' in OneOrTwo: readfile = readfile.replace('tsrc'+str(tsource),'')
            # if os.path.isfile(readfile):
            data = ReadXmlAndPickle(readfile)[0]
            if len(data.keys()) > 0:
                data = data[data.keys()[0]]
                if 'Boots' in data.keys():
                    bootdata = data['Boots']
                    if thismom not in dictout.keys():dictout[thismom] = {}
                    dictout[thismom][ipar] = {}
                    dictout[thismom][ipar]['Info'] = data['Info']
                    for ifit,fitrdata in bootdata.iteritems():
                        thisfit = FitFlagXmlToOldSF(ifit)
                        dictout[thismom][ipar][thisfit] = {}
                        if twoptread:
                            dictout[thismom][ipar][thisfit]['Boot'] = BootStrap1(nboot,0)
                            dictout[thismom][ipar][thisfit]['Boot'].values = np.array(fitrdata)
                            dictout[thismom][ipar][thisfit]['Boot'].Stats()
                            dictout[thismom][ipar][thisfit]['Avg'] = dictout[thismom][ipar][thisfit]['Boot'].Avg
                            dictout[thismom][ipar][thisfit]['Std'] = dictout[thismom][ipar][thisfit]['Boot'].Std
                            dictout[thismom][ipar][thisfit]['Chi'] = data['Values'][ifit]['Chi']
                        else:
                            for icut,cutdata in fitrdata.iteritems():
                                dictout[thismom][ipar][thisfit][icut] = {}
                                dictout[thismom][ipar][thisfit][icut]['Boot'] = BootStrap1(nboot,0)
                                dictout[thismom][ipar][thisfit][icut]['Boot'].values = np.array(cutdata)
                                dictout[thismom][ipar][thisfit][icut]['Boot'].Stats()
                                dictout[thismom][ipar][thisfit][icut]['Avg'] = dictout[thismom][ipar][thisfit][icut]['Boot'].Avg
                                dictout[thismom][ipar][thisfit][icut]['Std'] = dictout[thismom][ipar][thisfit][icut]['Boot'].Std
                                dictout[thismom][ipar][thisfit][icut]['Chi'] = data['Values'][ifit][icut]['Chi']

                else:
                    bootdata = data['Values']
                    if thismom not in dictout.keys():dictout[thismom] = {}
                    dictout[thismom][ipar] = {}
                    dictout[thismom]['Info'] = data['Info']
                    for ifit,fitrdata in bootdata.iteritems():
                        thisfit = FitFlagXmlToOldSF(ifit)
                        dictout[thismom][ipar][thisfit] = dictout[thismom][ifit]
    return dictout
            
            

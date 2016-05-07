#!/usr/bin/env python

import xmltodict
from XmlFuns import *
from XmlFormatting import *
import os
from BootTest import BootStrap1
from Params import *
from FitParams import *
import cPickle as pickle


def ReadPickleBoot(filein):
    try:
        with open(filein , "rb") as pfile:
            dataout = pickle.load( pfile )
    except:
        print 'Reading pickle file fail: ' + filein
        dataout = {}
    return dataout


def ReadXmlDict(filein):
    try:
        with open(filein,'r') as f:
            xmldata = RecFTDAS(xmltodict.parse(f.read()))
        bootfile = xmldata[xmldata.keys()[0]]['Boots']
    except:
        # print 'Reading xml file fail: ' + filein
        xmldata = {'Null':{'Values':{},'Boots':{}}}
        bootfile = 'NoBootDir'
    return xmldata,bootfile

def ReadXmlAndPickle(filein):
    xmldata = ReadXmlDict(filein)
    firstkey,bootfile = xmldata.keys()[0]
    if firstkey != 'Null':
        xmldata[firstkey]['Boots'] = ReadPickleBoot(bootfile)
    return xmldata,bootfile

    
def CheckMomFile(filein):
    ip = GetqcondFromFilename(filein)
    xmldata,bootfile = ReadXmlDict(filein)
    firstkey = xmldata.keys()[0]
    print filein
    if os.path.isfile(bootfile) and ip in firstkey:
        print 'True'
        return True
    print 'False'
    return False


##Also works for cfuns##
##xmlinput = { Ratio_Factor , Boots/Values , thismomlist , tlist } 
##outputdict = { thismom , [tVals] / [Vals] / [Valserr] / [Boot] bs }
def ReadRFFile(filedir,filename,thisMomList=RunMomList):
    renorm = GetRenorm(filename)
    dictout = {}
    for thismom in thisMomList:
        ip = qstrTOqcond(thismom)
        readfile = filedir+MakeMomDir(ip)+filename.replace('.xml',ip+'.xml')
        if os.path.isfile(readfile):
            data = ReadXmlAndPickle(readfile)[0]
            data = data[data.keys()[0]]
            if 'Boots' in data.keys():
                bootdata = data['Boots']
                dictout[thismom] = {}
                dictout[thismom]['tVals'] = map(untstr,bootdata.keys())
                dictout[thismom]['Boot'] = []
                dictout[thismom]['Vals'] = []
                dictout[thismom]['Valserr'] = []
                for bootlist in bootdata.itervalues():
                    dictout[thismom]['Boot'].append(BootStrap1(nboot,0))
                    dictout[thismom]['Boot'][-1].values = [iboot*renorm for iboot in bootlist]
                    dictout[thismom]['Boot'][-1].Stats()
                    dictout[thismom]['Vals'].append(dictout[thismom]['Boot'][-1].Avg)
                    dictout[thismom]['Valserr'].append(dictout[thismom]['Boot'][-1].Std)
            else:
                bootdata = data['Values']
                dictout[thismom] = {}
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
        if os.path.isfile(readfile):
            data = ReadXmlAndPickle(readfile)[0]
            data = data[data.keys()[0]]
            if 'Boots' in data.keys():
                bootdata = data['Boots']
                for icut,cutdata in bootdata.iteritems():
                    dictout[thismom][icut] = {}
                    dictout[thismom][icut]['Boot'] = BootStrap1(nboot,0)
                    dictout[thismom][icut]['Boot'].values = cutdata
                    dictout[thismom][icut]['Boot'].Stats()
                    dictout[thismom][icut]['Avg'] = dictout[thismom][icut]['Boot'].Avg
                    dictout[thismom][icut]['Std'] = dictout[thismom][icut]['Boot'].Std
                    dictout[thismom][icut]['Chi'] = data['Values'][icut]['Chi']
            else:
                dictout[thismom] = data
    return dictout

##outputdict = { thismom , cutpar , tsinkrpar/tsinkval , Avg / Std / Chi / Boot (bs) }
def ReadSumFile(filedir,filename,thisMomList=RunMomList):
    dictout = {}
    for thismom in thisMomList:
        ip = qstrTOqcond(thismom)
        readfile = filedir+MakeMomDir(ip)+filename.replace('.xml',ip+'.xml')
        if os.path.isfile(readfile):
            data = ReadXmlAndPickle(readfile)[0]
            data = data[data.keys()[0]]
            if 'Boots' in data.keys():
                bootdata = data['Boots']
                dictout[thismom] = {}
                for icut,cutdata in bootdata.iteritems():
                    dictout[thismom][icut] = {}
                    for it,tdata in cutdata.iteritems():
                        if 'tsink' in it:
                            dictout[thismom][icut][it] = {}
                            dictout[thismom][icut][it]['Boot'] = BootStrap1(nboot,0)
                            dictout[thismom][icut][it]['Boot'].values = tdata
                            dictout[thismom][icut][it]['Boot'].Stats()
                            dictout[thismom][icut][it]['Avg'] = dictout[thismom][icut][it]['Boot'].Avg
                            dictout[thismom][icut][it]['Std'] = dictout[thismom][icut][it]['Boot'].Std
                        else:
                            for ir,rdata in tdata.iteritems():
                                thisfitr = FitFlagXmlToOld(it,ir)
                                dictout[thismom][icut][thisfitr] = {}
                                dictout[thismom][icut][thisfitr]['Boot'] = BootStrap1(nboot,0)
                                dictout[thismom][icut][thisfitr]['Boot'].values = rdata
                                dictout[thismom][icut][thisfitr]['Boot'].Stats()
                                dictout[thismom][icut][thisfitr]['Avg'] = dictout[thismom][icut][thisfitr]['Boot'].Avg
                                dictout[thismom][icut][thisfitr]['Std'] = dictout[thismom][icut][thisfitr]['Boot'].Std
                                dictout[thismom][icut][thisfitr]['Chi'] = data['Values'][icut][it][ir]['Chi']

            else:
                valdata = data['Values']
                dictout[thismom] = {}
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
    dataout = {}
    if '.txt' in filename: filename = filename.replace('.txt','.xml')
    if not os.path.isfile(filename):
        mprint(filename + ' not found')
    else:
        data = ReadXmlAndPickle(filename)[0]
        data = data[data.keys()[0]]
        dataout = OrderedDict()
        dataout['Mass'] = data['Values']['Mass']
        dataout['Chi'] = OrderedDict()
        if 'Boots' in data.keys():
            for iq,qdata in data['Boots'].iteritems():
                dataout['Chi'][iq] = data['Values'][iq]['Chi']
                for iff,ffdata in qdata.iteritems():
                    if iff not in dataout.keys():
                        dataout[iff] = OrderedDict()
                    dataout[iff][iq] = OrderedDict()
                    dataout[iff][iq]['Boot'] = BootStrap1(nboot,0)
                    dataout[iff][iq]['Boot'].values = ffdata
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
            if os.path.isfile(readfile):
                data = ReadXmlAndPickle(readfile)[0]
                data = data[data.keys()[0]]
                if 'Boots' in data.keys():
                    bootdata = data['Boots']
                    if thismom not in dictout.keys():dictout[thismom] = {}
                    dictout[thismom][ipar] = {}
                    for ifit,fitrdata in bootdata.iteritems():
                        thisfit = FitFlagXmlToOldSF(ifit)
                        dictout[thismom][ipar][thisfit] = {}
                        if twoptread:
                            dictout[thismom][ipar][thisfit]['Boot'] = BootStrap1(nboot,0)
                            dictout[thismom][ipar][thisfit]['Boot'].values = fitrdata
                            dictout[thismom][ipar][thisfit]['Boot'].Stats()
                            dictout[thismom][ipar][thisfit]['Avg'] = dictout[thismom][ipar][thisfit]['Boot'].Avg
                            dictout[thismom][ipar][thisfit]['Std'] = dictout[thismom][ipar][thisfit]['Boot'].Std
                            dictout[thismom][ipar][thisfit]['Chi'] = data['Values'][ifit]['Chi']
                        else:
                            for icut,cutdata in fitrdata.iteritems():
                                dictout[thismom][ipar][thisfit][icut] = {}
                                dictout[thismom][ipar][thisfit][icut]['Boot'] = BootStrap1(nboot,0)
                                dictout[thismom][ipar][thisfit][icut]['Boot'].values = cutdata
                                dictout[thismom][ipar][thisfit][icut]['Boot'].Stats()
                                dictout[thismom][ipar][thisfit][icut]['Avg'] = dictout[thismom][ipar][thisfit][icut]['Boot'].Avg
                                dictout[thismom][ipar][thisfit][icut]['Std'] = dictout[thismom][ipar][thisfit][icut]['Boot'].Std
                                dictout[thismom][ipar][thisfit][icut]['Chi'] = data['Values'][ifit][icut]['Chi']

                else:
                    bootdata = data['Values']
                    if thismom not in dictout.keys():dictout[thismom] = {}
                    dictout[thismom][ipar] = {}
                    for ifit,fitrdata in bootdata.iteritems():
                        thisfit = FitFlagXmlToOldSF(ifit)
                        dictout[thismom][ipar][thisfit] = dictout[thismom][ifit]
    return dictout
            
            

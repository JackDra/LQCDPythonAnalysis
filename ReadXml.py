#!/usr/bin/env python

import xmltodict
from XmlFuns import *
from XmlFormatting import *
import os
from BootTest import BootStrap1
from Params import *


def ReadXmlDict(filein):
    with open(filein,'r') as f:
        data = RecFTDAS(xmltodict.parse(f.read()))
    return data


def ReadXmlArrayed(filein):
    data = ReadXmlDict(filein)
    

##Also works for cfuns##
##xmlinput = { Ratio_Factor , Boots/Values , thismomlist , tlist } 
##outputdict = { thismom , [tVals] / [Vals] / [Valserr] / [Boot] bs }
def ReadRFFile(filename,bootfn='',thisMomList=[]):
    renorm = GetRenorm(filename)
    dictout = {}
    if '.txt' in filename: filename = filename.replace('.txt','.xml')
    if not os.path.isfile(filename):
        mprint(filename + ' not found')
    else:
        data = ReadXmlDict(filename)
        data = data[data.keys()[0]]
        if 'Boots' in data.keys():
            bootdata = data['Boots']
            for imom,momdata in bootdata.iteritems():
                thismom = qcondTOqstr(imom)
                dictout[thismom] = {}
                dictout[thismom]['tVals'] = map(untstr,momdata.keys())
                dictout[thismom]['Boot'] = []
                dictout[thismom]['Vals'] = []
                dictout[thismom]['Valserr'] = []
                for bootlist in momdata.itervalues():
                    dictout[thismom]['Boot'].append(BootStrap1(nboot,0))
                    dictout[thismom]['Boot'][-1].values = [iboot*renorm for iboot in bootlist]
                    dictout[thismom]['Boot'][-1].Stats()
                    dictout[thismom]['Vals'].append(dictout[thismom]['Boot'][-1].Avg)
                    dictout[thismom]['Valserr'].append(dictout[thismom]['Boot'][-1].Std)
        else:
            bootdata = data['Values']
            for imom,momdata in bootdata.iteritems():
                thismom = qcondTOqstr(imom)
                dictout[thismom] = {}
                dictout[thismom]['tVals'] = map(untstr,momdata.keys())
                dictout[thismom]['Vals'] = []
                dictout[thismom]['Valserr'] = []
                for tdata in momdata.itervalues():
                    dictout[thismom]['Vals'].append(tdata['Avg']*renorm)
                    dictout[thismom]['Valserr'].append(tdata['Std'])
    return dictout
        
                    

def ReadFitFile(filename,bootfn='',thisMomList=[]):
    dictout = {}
    if '.txt' in filename: filename = filename.replace('.txt','.xml')
    if not os.path.isfile(filename):
        mprint(filename + ' not found')
    else:
        data = ReadXmlDict(filename)
        data = data[data.keys()[0]]
        if 'Boots' in data.keys():
            bootdata = data['Boots']
            for imom,momdata in bootdata.iteritems():
                thismom = qcondTOqstr(imom)
                dictout[thismom] = {}
                for icut,cutdata in momdata.iteritems():
                    dictout[thismom][icut] = {}
                    dictout[thismom][icut]['Boot'] = BootStrap1(nboot,0)
                    dictout[thismom][icut]['Boot'].values = cutdata
                    dictout[thismom][icut]['Boot'].Stats()
                    dictout[thismom][icut]['Avg'] = dictout[thismom][icut]['Boot'].Avg
                    dictout[thismom][icut]['Std'] = dictout[thismom][icut]['Boot'].Std
                    dictout[thismom][icut]['Chi'] = data['Values'][imom][icut]['Chi']
        else:
            bootdata = data['Values']
            for imom,momdata in bootdata.iteritems():
                thismom = qcondTOqstr(imom)
                dictout[thismom] = momdata
    return dictout

##outputdict = { thismom , cutpar , tsinkrpar/tsinkval , Avg / Std / Chi / Boot (bs) }
def ReadSumFile(filename,bootfn='',thisMomList=[]):
    dictout = {}
    if '.txt' in filename: filename = filename.replace('.txt','.xml')
    if not os.path.isfile(filename):
        mprint(filename + ' not found')
    else:
        data = ReadXmlDict(filename)
        data = data[data.keys()[0]]
        if 'Boots' in data.keys():
            bootdata = data['Boots']
            for imom,momdata in bootdata.iteritems():
                thismom = qcondTOqstr(imom)
                dictout[thismom] = {}
                for icut,cutdata in momdata.iteritems():
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
                                dictout[thismom][icut][thisfitr]['Chi'] = data['Values'][imom][icut][it][ir]['Chi']

        else:
            valdata = data['Values']
            for imom,momdata in valdatadata.iteritems():
                thismom = qcondTOqstr(imom)
                dictout[thismom] = {}
                for icut,cutdata in momdata.iteritems():
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
def ReadFFFile(filename,bootfn=''):
    outputdict = OrderedDict()
    if os.path.isfile(filename):
        with open(filename,'r') as f:
            for line in f:
                rdata = line.strip()
                if len(rdata) == 0: continue
                rdata = rdata.split()
                if 'Mass' in rdata[0]:
                    outputdict['Mass'] = {'Set':rdata[1],
                                          'Avg':float(rdata[2]),
                                          'Std':float(rdata[3]),
                                          'Chi':float(rdata[4])}
                elif 'q**2' in rdata[0]:
                    for idata in rdata:
                        if 'FF' in idata and 'Err' not in idata :
                            outputdict[idata] = OrderedDict()
                        elif 'Chi' in idata:
                            outputdict['Chi'] = OrderedDict()
                elif 'qsqrd' in rdata[0]:
                    thisikey = 0
                    for ikey,dictkey in enumerate(outputdict.keys()):
                        if 'Mass' in dictkey: continue
                        thisikey += 1
                        if 'FF' in dictkey:
                            outputdict[dictkey][rdata[0]] = {'Avg':float(rdata[thisikey]),
                                                             'Std':float(rdata[thisikey+1])}
                            thisikey +=1
                        else:
                            outputdict[dictkey][rdata[0]] = float(rdata[thisikey])                            
        if os.path.isfile(bootfn):
            RFF,qsqrdread = None,None
            with open(bootfn,'r') as f:
                for line in f:
                    rdata = line.strip()
                    if len(rdata) == 0: continue
                    rdata = rdata.split()
                    if 'nboot' in rdata[0]:
                        nboot = int(rdata[2])
                    elif 'Mass' in rdata[0] and 'Mass' in outputdict.keys():
                        outputdict['Mass']['Boot'] = BootStrap1(nboot,1)
                        RFF = 'Mass'
                    elif 'FF' in rdata[0]:
                        RFF = rdata[0]
                    elif 'qsqrd' in rdata[0] and RFF in outputdict.keys():
                        qsqrdread = rdata[0]
                        outputdict[RFF][qsqrdread]['Boot'] = BootStrap1(nboot,1)
                    else:
                        if 'Mass' in RFF:
                            outputdict['Mass']['Boot'].values[int(rdata[0])] = float(rdata[1])
                        else:
                            outputdict[RFF][qsqrdread]['Boot'].values[int(rdata[0])] = float(rdata[1])
            BootNdimDict(outputdict)
    return outputdict

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
    dictout = {}
    if '.txt' in filename: filename = filename.replace('.txt','.xml')
    if os.path.isfile(filename):
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
                    dictout[thismom]['Boot'][-1].values = bootlist
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
                    dictout[thismom]['Vals'].append(tdata['Avg'])
                    dictout[thismom]['Valserr'].append(tdata['Std'])
    return dictout
        
                    

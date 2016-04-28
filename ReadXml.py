#!/usr/bin/env python

import xmltodict
from XmlFuns import *

def ReadMyXml(filein):
    with open(filein+'.xml','r') as f:
        data = RecursiveFTDASxmltodict.parse(f.read()))
    return data


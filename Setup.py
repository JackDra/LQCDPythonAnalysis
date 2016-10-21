#!/usr/bin/env python

import os
import sys


thisscriptdir = os.getcwd()

if 'LQCDPythonAnalysis' not in thisscriptdir:
    raise EnvironmentError('Please run from the script directory .../LQCDPythonAnalysis ')

print 'setting scriptdir to ',thisscriptdir

if 'default' in sys.argv[-1]:
    thisdatadir = thisscriptdir
    thisAnaProc = 1
    thisListOrSet = 'ReadSet'
    thiskappa = 12090
    thisNShifts = 1
    with open(thisscriptdir+'/setup.cfg','w') as f:
        f.write('\nscriptdir:\n')
        f.write(thisscriptdir+'/\n')
        f.write('\ndatadir:\n')
        f.write(thisdatadir+'\n')
        f.write('\nAnaProc:\n')
        f.write(str(thisAnaProc)+'\n')
        f.write('\nListOrSet:\n')
        f.write(thisListOrSet+'\n')
        f.write('\nkappa:\n')
        f.write(str(thiskappa)+'\n')
        f.write('\nPoFShifts:\n')
        f.write(str(thisNShifts)+'\n')
    sys.exit()
    
print ''
print 'Please Specify data directory:'
print 'This should include correlators in .../cfuns/k######/source#/.../....3cf'
print '(see README.md for more information)'
thisdatadir = raw_input("")

print ''
print "Please Specify number of processors for multiprocessing:"
thisAnaProc = int(raw_input("(note, don't kill your machine)\n"))

# print ''
# thisListOrSet = raw_input("Read All Correlators in directorys? (y/n):\n")

# if thisListOrSet == 'y':
#     thisListOrSet = 'ReadSet'
# elif thisListOrSet == 'n':
#     thisListOrSet = 'ReadList'
# else:
#     raise IOError('Input y or n')

print ''
thisListOrSet = raw_input("ReadList or ReadSet (add previx/suffix for different runs if you like)\n")

if 'ReadList' in thisListOrSet:
    print 'Using ReadList'
elif 'ReadSet' in thisListOrSet:
    print 'Using ReadSet'
else:
    raise IOError('Input string that contains ReadList or ReadSet')


print ''
thiskappa = int(raw_input("Kappa Value: (0.XXXXX)\n"))

print ''
thisNShifts = int(raw_input("How Many Shifts?\n"))

with open(thisscriptdir+'/setup.cfg','w') as f:
    f.write('\nscriptdir:\n')
    f.write(thisscriptdir+'/\n')
    f.write('\ndatadir:\n')
    f.write(thisdatadir+'\n')
    f.write('\nAnaProc:\n')
    f.write(str(thisAnaProc)+'\n')
    f.write('\nListOrSet:\n')
    f.write(thisListOrSet+'\n')
    f.write('\nkappa:\n')
    f.write(str(thiskappa)+'\n')
    f.write('\nPoFShifts:\n')
    f.write(str(thisNShifts)+'\n')


from FFParams import DumpAllMomLists
DumpAllMomLists()

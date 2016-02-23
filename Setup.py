#!/usr/bin/env python

import os

thisscriptdir = os.getcwd()

if 'LQCDPythonAnalysis' not in thisscriptdir:
    raise EnvironmentError('Please run from the script directory .../LQCDPythonAnalysis ')

print 'setting scriptdir to ',thisscriptdir


print ''
print 'Please Specify data directory:'
print 'This should include correlators in .../cfuns/k######/source#/.../....3cf'
print '(see README.md for more information)'
thisdatadir = raw_input("")

print ''
print "Please Specify number of processors for multiprocessing:"
thisAnaProc = int(raw_input("(note, don't kill your machine)\n"))

print ''
thisListOrSet = raw_input("Read All Correlators in directorys? (y/n):\n")

if thisListOrSet == 'y':
    thisListOrSet = 'ReadSet'
elif thisListOrSet == 'n':
    thisListOrSet = 'ReadList'
else:
    raise IOError('Input y or n')
    
print ''
thiskappa = int(raw_input("Kappa Value:\n"))


print ''
thisBoolPoF = raw_input("Do Pencil of Function method (y/n)?\n")

if thisBoolPoF == 'y':
    print ''
    thisNShifts = int(raw_input("How Many Shifts?\n"))
else:
    thisNShifts = 0



f = open(thisscriptdir+'/setup.cfg','w')
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
if thisBoolPoF == 'y':
    f.write('\nPoFOrSum:\n')
    f.write('PoF\n')
    f.write('\nPoFShifts:\n')
    f.write(str(thisNShifts)+'\n')
    f.write('\nVarPref:\n')
    f.write('PoF'+str(thisNShifts)+'\n')
else:
    f.write('\nPoFOrSum:\n')
    f.write('sum\n')
    f.write('\nPoFShifts:\n')
    f.write(str(thisNShifts)+'\n')
    f.write('\nVarPref:\n')
    f.write('sum\n')

f.close()



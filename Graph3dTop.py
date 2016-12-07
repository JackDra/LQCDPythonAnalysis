#!/usr/bin/env python

import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

import cPickle as pickle
import os.path
import sys

##datadict = { gamma } { mom } { method } { set }
##thisdatadict = { method } { set }
##thisGammaList format: D/S P4/3 Opp (e.g. doubP4g4)
##thisMomList in qstr format (e.g. 'q = 0 0 0')
##NB q = 0 0 0 MUST BE INCLUDED FOR TSF ETC

PickleFile = '/home/jackdra/PHD/CHROMA/TestVar/scratch/pickledir/TopChargePlot.p'

trange = 3,21
valrange = -10,10

pfile = open( PickleFile, "rb" )
thiscollist,thistlist,thistflowlist,thisplotAvg,thisplotUp,thisplotDown =pickle.load( pfile )
pfile.close()
fig = pl.figure()
ax = fig.gca(projection='3d')
for icol,it,itflow,iavg,iup,idown in zip(thiscollist,thistlist,thistflowlist,thisplotAvg,thisplotUp,thisplotDown):
    for jcol,jt,jtflow,javg,jup,jdown in zip(icol,it,itflow,iavg,iup,idown):
        # print
        # for pt,ptflow,pavg in zip(jt,jtflow,javg):
        #     print pt,ptflow,pavg
        jtp,jtflowp,javgp,jupp,jdownp = [],[],[],[],[]
        for tplot,tflowplot,avgplot,upplot,downplot in zip(jt,jtflow,javg,jup,jdown):
            if tplot < trange[1] and tplot > trange[0]:
                if avgplot < valrange[1] and avgplot > valrange[0]:
                    jtp.append(tplot)
                    jtflowp.append(tflowplot)
                    javgp.append(avgplot)
                    jupp.append(upplot)
                    jdownp.append(downplot)
        ax.scatter(jtp,jtflowp,javgp,color=jcol) 
        for i in np.arange(0, len(javgp)):
            ax.plot([jtp[i], jtp[i]], [jtflowp[i], jtflowp[i]], [jupp[i], jdownp[i]],color=jcol, marker="_")
ax.set_xlabel(r'$t_{sink}$')
ax.set_ylabel(r'$t_{flow}$')
ax.set_zlabel(r'$\frac{\langle NNQ \rangle>}{\langle NN\rangle}$')
ax.set_xlim3d(*trange)
ax.set_zlim3d(*valrange)
pl.show()


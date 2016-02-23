#!/usr/bin/env python

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
from SetLists import *
from ReadTxt import *
from OppFuns import *
import FitFunctions as ff
import cPickle as pickle
import os.path

colourset8 = [ '#000080','#B22222','#00B800','#8B008B', '#0000FF','#FF0000','#000000','#32CD32','#FF0066']
colcyc = itertools.cycle(colourset8)
thisGammaList = CreateGammaList(sys.argv[1:],twopt=True)
thisMomList = ['q = 0 0 0']
thisMethodList = ['RF','TSFTsink','TSFtest32','TSFSmall']
# thisMethodList = ['RF','TSFTsink']
loopML = thisMethodList[1:]
fixedTSFset = 'sm32'
thisTSFcut = 'cut3'
intcut = int(thisTSFcut.replace('cut',''))
thisres = 1
thisalpha = .4
thisSetList = CreateSet(thisSmearL=['32'],thisStateL=[])[0]



if 'erwin' in THISMACHINE:
    datadict = ReadSetFitRFDict(outputdir,thisSetList,thisGammaList,thisMethodList,thisMomList=thisMomList)
    for imom in thisMomList:
        for igamma in thisGammaList:
            if any([idst in igamma for idst in ['doub','sing','twopt']]): continue
            PickleRFFile = '/home/accounts/jdragos/scripts/Plot3D/RF'+igamma+imom.replace('=','').replace(' ','')+'.p'
            print igamma, ' ', imom
            pxdata = []
            pydata = []
            pzdata = []
            pzedata = []
            for thisset,setdata in datadict[igamma][imom]['RFTSFTsink'].iteritems():
                thists = setdata['tVals'][0]
                xval = np.array(setdata['tVals'])-thists-(GetintTSink(thisset)-thists)/2.0
                yval = [GetintTSink(thisset)-thists]*len(setdata['tVals'])
                zval = setdata['Vals']
                zerr = setdata['Valserr']
                pxdata.append(xval[1:-1])
                pydata.append(yval[1:-1])
                pzdata.append(zval[1:-1])
                pzedata.append(zerr[1:-1])
            pfile = open(PickleRFFile , "wb" )
            pickle.dump( [pxdata,pydata,pzdata,pzedata], pfile )
            pfile.close()
            for iTSF in loopML:
                PickleTSFFile = '/home/accounts/jdragos/scripts/Plot3D/'+iTSF+igamma+imom.replace('=','').replace(' ','')+'.p'
                data3pt = datadict[igamma][imom][iTSF][fixedTSFset]
                data2pt = datadict['twopt'][imom][iTSF][fixedTSFset]
                pars2pt,pars3pt = [],[]
                for ipar in StateParList['Two']['C2']:
                    pars2pt.append(data2pt[ipar][TSFfitr]['Avg'])
                for ipar in StateParList['Two']['C3']:
                    pars3pt.append(data3pt[ipar][TSFfitr][thisTSFcut]['Avg'])
                def RFFun(thistsink,thistau):
                    return ff.C3TSFLineFun(pars2pt+pars3pt,thistau,thistsink)/ff.C2TSFLineFun(thistsink,pars2pt)
                mintsink,maxtsink = 10.0,26.0
                thistsink = np.arange(mintsink,maxtsink+thisres,thisres)
                tsinkmesh,taumesh,RFplot = [],[],[]
                for its in thistsink:
                    for itau in np.arange(-its/2.0+intcut,its/2.0+thisres-intcut,thisres):
                        tsinkmesh.append(its)
                        taumesh.append(itau)
                        RFplot.append(RFFun(its,itau+its/2.0))
                pfile = open(PickleTSFFile , "wb" )
                pickle.dump( [taumesh,tsinkmesh,RFplot], pfile )
                pfile.close()

else:
    for imom in thisMomList:
        for igamma in thisGammaList:
            if any([idst in igamma for idst in ['doub','sing','twopt']]): continue
            print igamma, ' ', imom
            PickleRFFile = '/home/jackdra/PHD/DataAnalysis/Plot3D/RF'+igamma+imom.replace('=','').replace(' ','')+'.p'
            pfile = open( PickleRFFile, "rb" )
            xvallist,yvallist,zvallist,zerrlist = pickle.load( pfile )
            pfile.close()
            for iTSF in loopML:
                fig = plt.figure()
                ax = fig.gca(projection='3d')
                colcyc = itertools.cycle(colourset8)
                for xval,yval,zval,zerr in zip(xvallist,yvallist,zvallist,zerrlist):
                    thiscol = next(colcyc)
                    ax.scatter(xval,yval,zval,color=thiscol)
                    for i in np.arange(0, len(xval)):
                        ax.plot([xval[i], xval[i]], [yval[i], yval[i]], [zval[i]+zerr[i], zval[i]-zerr[i]],color=thiscol, marker="_")
                PickleTSFFile = '/home/jackdra/PHD/DataAnalysis/Plot3D/'+iTSF+igamma+imom.replace('=','').replace(' ','')+'.p'
                if os.path.isfile(PickleTSFFile):
                    pfile = open( PickleTSFFile, "rb" )
                    taumesh,tsinkmesh,RFplot = pickle.load( pfile )
                    ax.plot_trisurf(taumesh,tsinkmesh,RFplot,linewidth=0,alpha=thisalpha)
                    # ax.plot_surface(taumesh,tsinkmesh,RFplot,linewidth=0,alpha=thisalpha)
                ax.set_xlabel(r'$\tau$')
                ax.set_ylabel(r'$t$')
                ax.set_zlabel(r'$RF$')
                ax.set_title(iTSF+' '+igamma + ' ' + imom)
                plt.show()
    
        

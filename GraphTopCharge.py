#!/usr/bin/env python

# from Params import *
# import numpy as np
from ReadXml import *
# from MiscFuns import *
# from FitParams import *
# from StringFix import *
# from SetLists import *
# from OppFuns import CreateGammaList
# import time,datetime
# from MultiWrap import *
# from multiprocessing import Pool
from InputArgs import *
# from CreateCombs import CreateDictOldCombs
# from CombParams import *
import cPickle as pickle


##datadict = { gamma } { mom } { method } { set }
##thisdatadict = { method } { set }
##thisGammaList format: D/S P4/3 Opp (e.g. doubP4g4)
##thisMomList in qstr format (e.g. 'q = 0 0 0')
##NB q = 0 0 0 MUST BE INCLUDED FOR TSF ETC


def PlotTopCharge(data,iSet,iMom):
    tflowlist,tlist,plotAvg,plotUp,plotDown = [],[],[],[],[]
    momdata = data['RF'][iMom]['Boots']
    for itflow,flowdata in momdata.iteritems():        
        # tflowlist.append([])
        # tlist.append([])
        # plotAvg.append([])
        # plotUp.append([])
        # plotDown.append([])
        for it,tdata in flowdata.iteritems():
            # tflowlist[-1].append(untflowstr(itflow))
            # tlist[-1].append(untstr(it))
            tavg = np.mean(tdata)
            tstd = np.std(tdata)
            # plotAvg[-1].append(tdata.Avg)
            # plotUp[-1].append(tdata.Avg+tdata.Std)
            # plotDown[-1].append(tdata.Avg-tdata.Std)
            # plotAvg[-1].append(tavg)
            # plotUp[-1].append(tavg+tstd)
            # plotDown[-1].append(tavg-tstd)
            tflowlist.append(untflowstr(itflow))
            tlist.append(untstr(it))
            plotAvg.append(tavg)
            plotUp.append(tavg+tstd)
            plotDown.append(tavg-tstd)

    if len(plotAvg) == 0: return
    return tlist,tflowlist,plotAvg,plotUp,plotDown
    # ax.errorbar(np.array(plotlist),Pullflag(plotdata,'Avg'),Pullflag(plotdata,'Std'),color=thiscolor,label=LegLab(iSet+'\ '+qstrTOqcond(iMom)))
    
def PlotTopSetCharge(data,thisSetList,thisMomList,FT):
    # global ax
    ForceTitle = FT
    thiscolcyc = itertools.cycle([ '#000080','#B22222','#00B800','#8B008B', '#0000FF','#FF0000','#000000','#32CD32','#FF0066'])
    # fig = pl.figure()
    # ax = fig.gca(projection='3d')
    thiscollist = []
    thistlist,thistflowlist,thisplotAvg,thisplotUp,thisplotDown = [],[],[],[],[]
    for iset,setdata in zip(thisSetList,data):
        thiscollist.append([])
        thistlist.append([])
        thistflowlist.append([])
        thisplotAvg.append([])
        thisplotUp.append([])
        thisplotDown.append([])
        for imom in thisMomList:
            if CheckDict(setdata,'RF',imom,'Boots'):
                print 'plotting ', iset, imom
                thiscollist[-1].append(thiscolcyc.next())
                tlist,tflowlist,plotAvg,plotUp,plotDown = PlotTopCharge(setdata,iset,imom)
                thistlist[-1].append(tlist)
                thistflowlist[-1].append(tflowlist)
                thisplotAvg[-1].append(plotAvg)
                thisplotUp[-1].append(plotUp)
                thisplotDown[-1].append(plotDown)
                
                # for it,itflow,iAvg in zip(tlist,tflowlist,plotAvg):
                #     print
                #     for jt,jtflow,jAvg in zip(it,itflow,iAvg):
                #         print jt,jtflow,jAvg
                # ax.scatter(tlist,tflowlist,plotAvg,color=thiscolcyc.next())
    # SetTopAxies('Ratio NNQdivNN ')
    # ax.savefig(CreateFile('','twopt','q = 0 0 0',,subdir='Top')+'.pdf')
    # fig.show()
    return thiscollist,thistlist,thistflowlist,thisplotAvg,thisplotUp,thisplotDown




def ReadAndPlotDis(thisSetList,thisMomList):
    # thisAllSetList = thisSmearList+thisSetList
    # for isetlist,dump in thisSetPoFLists[:len(thisSetPoFLists)/2]:
    #     thisAllSetList += isetlist
    iterSetList = SortMySet(ReduceTooMassSet(thisSetList))[0]
    thiscollist = []
    thistlist,thistflowlist,thisplotAvg,thisplotUp,thisplotDown = [],[],[],[],[]
    for iset in iterSetList:
        datadict = ReadTopFile(outputdir[0],iset,thisMomList=thisMomList)
        setstart = time.time()
        thiscol,tlist,tflowlist,plotAvg,plotUp,plotDown = PlotTopSetCharge([datadict],[iset],thisMomList,feedin['ForceTitle'])
        print 'Getting ' , iset, 'Took: ' , str(datetime.timedelta(seconds=(time.time()-setstart))) ,' h:m:s                      '
        thiscollist.append(thiscol[0])
        thistlist.append(tlist[0])
        thistflowlist.append(tflowlist[0])
        thisplotAvg.append(plotAvg[0])
        thisplotUp.append(plotUp[0])
        thisplotDown.append(plotDown[0])
    return thiscollist,thistlist,thistflowlist,thisplotAvg,thisplotUp,thisplotDown


# if all('-m' not in iin for iin in sys.argv[1:]):
#     feedin = InputParams(sys.argv[1:] + ['-noprompt'] + ['-m=Fits,TSFTsink,TSFtest32,OSF'])
# else:
feedin = InputParams(sys.argv[1:] + ['-noprompt'])
    


feedin['set'] = ReduceTooMassSet(feedin['set'])
ShowSetLists(feedin['set'])
feedin['mom'] = GetAvgMomList(feedin['mom'])

if DoMulticore and len(feedin['set']) > 1  and feedin['anaproc'] > 1:
    # inputparams = [([iset],feedin['mom'],feedin['method']) for iset in feedin['set']]
    inputparams = []
    for iset in feedin['set']:
        for imom in feedin['mom']:
            inputparams.append([iset],[imom])
    makeContextFunctions(ReadAndPlotDis)
    thisPool = Pool(min(len(inputparams),feedin['anaproc']))
    output = thisPool.map(ReadAndPlotDis.mapper,inputparams)
    thisPool.close()
    thisPool.join()
    counter = 0
    thiscollist = []
    thistlist,thistflowlist,thisplotAvg,thisplotUp,thisplotDown = [],[],[],[],[]
    for iset in feedin['set']:
        thiscollist.append([])
        thistlist.append([])
        thistflowlist.append([])
        thisplotAvg.append([])
        thisplotUp.append([])
        thisplotDown.append([])
        for imom in feedin['mom']:
            collist,tlist,tflowlist,plotAvg,plotUp,plotDown = output[counter]
            thiscollist[-1].append(collist[0])
            thistlist[-1].append(tlist[0])
            thistflowlist[-1].append(tflowlist[0])
            thisplotAvg[-1].append(plotAvg[0])
            thisplotUp[-1].append(plotUp[0])
            thisplotDown[-1].append(plotDown[0])
            counter += 1
else:
    thiscollist,thistlist,thistflowlist,thisplotAvg,thisplotUp,thisplotDown = ReadAndPlotDis(feedin['set'],feedin['mom'])


print pickledir+'TopChargePlot.p'    
pfile = open(pickledir+'TopChargePlot.p' , "wb" )
pickle.dump( [thiscollist,thistlist,thistflowlist,thisplotAvg,thisplotUp,thisplotDown], pfile )
pfile.close()

# fig = pl.figure()
# ax = fig.gca(projection='3d')
# for icol,it,itflow,iavg,iup,idown in zip(thiscollist,thistlist,thistflowlist,thisplotAvg,thisplotUp,thisplotDown):
#     for jcol,jt,jtflow,javg,jup,jdown in zip(icol,it,itflow,iavg,iup,idown):
#         ax.scatter(jt,jtflow,javg,color=jcol) 
# ax.set_xlabel(r'$t_{sink}$')
# ax.set_ylabel(r'$t_{flow}$')
# ax.set_zlabel(r'$\frac{\langle NNQ \rangle>}{\langle NN\rangle}$')
# fig.show()
         
# print 'Graphing all complete'
    

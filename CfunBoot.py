#!/usr/bin/env python
from array import array
import numpy as np
import BootTest as bt
from Params import *
from MiscFuns import *
from ReadBinaryCfuns import *
from FitFunctions import CreateOORFF,OORNFFDer
from LLSBoot import LSFit
from StringFix import GetCfgNumb
import matplotlib.pyplot as pl

##NB CreateBoot take array [ iconf , it ]

#data = [ iconf , igamma , ip , it ]
#dataout = [igamma , ip , it ] bs
def BootSet3pt(data,thisMomList,thisGammaList,nboot,printstr='',randlist=[]):
    dataout = []
    for ig,igamma in enumerate(thisGammaList):
        dataout.append([])
        for ip,imom in enumerate(thisMomList):
            # print 'Booting '+printstr+igamma+' {}%               \r'.format(int((ip*100)/float(len(thisMomList)))),
            # for it in range(16):
            #     print 'pie3pt'
            #     for idata in np.array(data)[:,ig,ip,:]:
            #         print it,idata[it]
            bootdata,randlist = bt.CreateBoot(np.array(data)[:,ig,ip,:],nboot,0,randlist=randlist)
            dataout[ig].append(bootdata)
    # print '                             \r',
    return dataout


#data = [ iconf , ip , it]
#dataout = [ ip , it ]. bs
def BootSet2pt(data,thisMomList,nboot,randlist=[]):
    dataout = []
    randlist = []
    for ip,imom in enumerate(thisMomList):
        # print 'Booting {}%  \r'.format(int((ip*100)/float(len(thisMomList)))),
        # for icfg,cfgdata in enumerate(np.array(data)[:,ip,:]):
        #     print ''
        #     print 'icfg=',icfg
        #     for it,tdata in enumerate(cfgdata):
        #         print tdata
        # print 'pie2pt'
        # for iboot,idata in enumerate(np.array(data)[:,ip,12]):
        #     print iboot,idata
        bootdata,randlist = bt.CreateBoot(data[:,ip,:],nboot,0,randlist=randlist)
        dataout.append(bootdata)
    # print '                              \r',
    return dataout,randlist

#data = [ iconf , iflow , igamma , ip , it ]
#dataout = [igamma , iflow,  ip , it ] bs
def BootSet3ptTC(data,thisMomList,thisGammaList,nboot,tflowlist,printstr='',randlist=[]):
    dataout = []
    for ig,igamma in enumerate(thisGammaList):
        dataout.append([])
        for icf,iflow in enumerate(tflowlist):
            dataout[ig].append([])
            for ip,imom in enumerate(thisMomList):
                # print 'Booting '+printstr+igamma+' {}%               \r'.format(int((ip*100)/float(len(thisMomList)))),
                # for it in range(16):
                #     print 'pie3pt'
                #     for idata in np.array(data)[:,ig,ip,:]:
                #         print it,idata[it]
                bootdata,randlist = bt.CreateBoot(np.array(data)[:,icf,ig,ip,:],nboot,0,randlist=randlist)
                dataout[ig][icf].append(bootdata)
    # print '                             \r',
    return dataout


#data = [ iconf , iflow , ip , it]
#dataout = [ iflow, ip , it ]. bs
def BootSet2ptTC(data,thisMomList,nboot,tflowlist,randlist=[]):
    dataout = []
    randlist = []
    for icf,iflow in enumerate(tflowlist):
        dataout.append([])
        for ip,imom in enumerate(thisMomList):
        # print 'Booting {}%  \r'.format(int((ip*100)/float(len(thisMomList)))),
        # for icfg,cfgdata in enumerate(np.array(data)[:,ip,:]):
        #     print ''
        #     print 'icfg=',icfg
        #     for it,tdata in enumerate(cfgdata):
        #         print tdata
        # print 'pie2pt'
        # for iboot,idata in enumerate(np.array(data)[:,ip,12]):
        #     print iboot,idata
            bootdata,randlist = bt.CreateBoot(data[:,icf,ip,:],nboot,0)
            # if icf == 40 and ip == 0:
            #     print 'BootValue'
            #     print data[:,icf,ip,7]
            #     myseed=1234*len(data)/nboot
            #     np.random.seed(myseed)
            #     locrandint=np.random.random_integers
            #     thislist = locrandint(0,len(data)-1,len(data))
            #     print thislist
            #     print np.mean(np.array(data)[thislist,icf,ip,7])
            #     print bootdata[7].values[0]
            #     # print bootdata[7].Avg, bootdata[7].Std
            #     print
            dataout[-1].append(bootdata)
    # print '                              \r',
    return dataout

#dataout = [ ip , it ]. bs
def ReadAndBoot2pt(readfilelist,thisMomList,thisnboot,randlist=[]):
    tempdata = []
    shiftlist = []
    for ifilepref,ifileList in readfilelist.iteritems():
        # print 'Reading {}%  \r'.format(int((iconf*100)/float(len(readfilelist)))),
        try:
            if CHROMA:
                if NoXAvg:
                    for ifile in ifileList:
                        data = R2CChromaXMLFileList([ifile],thisMomList,Dog5=False)
                        tempdata.append(data.data)
                else:
                    data = R2CChromaXMLFileList(ifileList,thisMomList,Dog5=False)
                    tempdata.append(data.data)
                shiftlist.append(data.tshiftlist)
        except NaNCfunError as e:
            print 
            print 'Deleting file ' + ifile
            print e
            print 'MUST RE-RUN AFTER THIS TO EXCLUDE BAD CONFIGS'
            print
            # os.remove(ifile)
    return BootSet2pt(np.array(tempdata),thisMomList,thisnboot,randlist=randlist),shiftlist


#dataout [ igamma , ip , it ]. bs
def ReadAndBoot3pt(readfilelist,thisMomList,thisGammaList,thisDerGammaList,thisnboot,printstr='',randlist=[]):
    tempdata = []
    for ifileList in readfilelist.itervalues():
        # print 'Reading '+printstr+' {}%            \r'.format(int((iconf*100)/float(len(readfilelist)))),
        try:
            if CHROMA:
                if NoXAvg:
                    for ifile in ifileList:
                        if len(thisGammaList) > 0:
                            data = ReadFSCfunPickCHROMA([ifile],thisMomList,thisGammaList)
                            tempdata.append(data.data)
                        if len(thisDerGammaList) > 0:
                            raise IOError('Chroma version does not do derivatives, make DerList in Params.py be empty')
                else:
                    if len(thisGammaList) > 0:
                        data = ReadFSCfunPickCHROMA(ifileList,thisMomList,thisGammaList)
                        tempdata.append(data.data)
                    if len(thisDerGammaList) > 0:
                        raise IOError('Chroma version does not do derivatives, make DerList in Params.py be empty')
            else:
                raise IOError('Top Charge not implemented for non chroma results')
        except NaNCfunError as e:
            print 
            print 'Deleting file ' + ifile
            print e
            print 'MUST RE-RUN AFTER THIS TO EXCLUDE BAD CONFIGS'
            print
            # os.remove(ifile)
    # print 'Values'
    # for idata in tempdata:
    #     print idata[0][0][5]
    # print 
    if len(thisGammaList) > 0:
        return BootSet3pt(tempdata,thisMomList,thisGammaList,thisnboot,printstr='',randlist=randlist)
    elif len(thisDerGammaList) > 0:
        return BootSet3pt(tempdata,thisMomList,thisDerGammaList,thisnboot,printstr='',randlist=randlist)


    
def ReadAndBoot2ptTop(readfilelist,thisMomList,thisnboot,chargedata,chargecfglist,flowlist,randlist=[]):
    tempdata = []
    tempdataTop = []
    shiftlist = []
    for ifilepref,ifileList in readfilelist.iteritems():
        # print 'Reading {}%  \r'.format(int((iconf*100)/float(len(readfilelist)))),
        try:
            if CHROMA:
                chargeindex = FileToChargeCfg(ifilepref,chargecfglist)
                if NoXAvg:
                    for ifile in ifileList:
                        data = R2CChromaXMLFileList([ifile],thisMomList,Dog5=True)
                        tempdataTop.append([])
                        tempdata.append(data.data)
                        # tempdata.append(np.array(data.datag5)*chargedata[chargeindex][40])
                        for iflowdata in chargedata[chargeindex]:
                            # tempdataTop[-1].append(np.array(data.data)*iflowdata)
                            # tempdataTop[-1].append((np.array(data.datag5)*iflowdata).tolist())
                            # tempdataTop[-1].append(np.array(data.datag5)*iflowdata)
                            tempdataTop[-1].append(np.array(data.datag5)*iflowdata)
                            # tempdataTop[-1].append(np.abs(np.array(data.datag5)*iflowdata))
                        # print ifile
                        # print data.datag5[0][7],chargedata[chargeindex][40]
                        # print 
                else:
                    data = R2CChromaXMLFileList(ifileList,thisMomList,Dog5=True)
                    tempdataTop.append([])
                    tempdata.append(data.data)
                    for iflowdata in chargedata[chargeindex]:
                        tempdataTop[-1].append(np.array(data.datag5)*iflowdata)
                        # tempdataTop[-1].append(np.array(data.datag5))
                shiftlist.append(data.tshiftlist)
                
            else:
                raise IOError('Top Charge not implemented for non chroma results')
        except NaNCfunError as e:
            print 
            print 'Deleting file ' + ifile
            print e
            print 'MUST RE-RUN AFTER THIS TO EXCLUDE BAD CONFIGS'
            print
            # os.remove(ifile)
    # pl.hist(setlist,bins=BinList,color=collist,label=leglist,stacked=Stacked,histtype=HistType,normed=Normed)
    # pl.hist(np.array(tempdataTop)[:,40,0,7])
    TCBdata,(Bdata,rlist),shiftlist = (BootSet2ptTC(np.array(tempdataTop),thisMomList,thisnboot,flowlist,randlist=randlist),
                                       BootSet2pt(np.array(tempdata),thisMomList,thisnboot,randlist=randlist),shiftlist)
    if PlotMonte:
        mkdir_p('./montegraphs')
        xlist,yavg,yerr,yavgNNQ,yerrNNQ = [],[],[],[],[]
        # print 'values'
        # print np.mean(np.array(tempdataTop)[:,MonteFlow,0,MonteTime-1]),np.std(np.array(tempdataTop)[:,MonteFlow,0,MonteTime-1])
        # print
        bplotdata = np.array(Bdata)[0,MonteTime-1]
        bplotdataNNQ = np.array(TCBdata)[MonteFlow,0,MonteTime-1]
        plotdata = np.array(tempdata)[:,0,MonteTime-1]
        plotdataNNQ = np.array(tempdataTop)[:,MonteFlow,0,MonteTime-1]
        for ic,(icfg,iread) in enumerate(readfilelist.iteritems()):
            xlist.append(icfg)
            # xlist += iread
            if NoXAvg:
                yavgNNQ.append(np.mean(plotdataNNQ[ic*len(iread):(ic+1)*len(iread)]))
                yerrNNQ.append(np.std(plotdataNNQ[ic*len(iread):(ic+1)*len(iread)]))
                yavg.append(np.mean(plotdata[ic*len(iread):(ic+1)*len(iread)]))
                yerr.append(np.std(plotdata[ic*len(iread):(ic+1)*len(iread)]))
            else:
                yavgNNQ.append(plotdataNNQ[ic])
                yerrNNQ.append(0.0)
                yavg.append(plotdata[ic])
                yerr.append(0.0)
        # pl.scatter(map(GetCfgNumb,xlist),plotdataNNQ)
        pl.fill_between([2500,4500],[bplotdataNNQ.Avg+bplotdataNNQ.Std,bplotdataNNQ.Avg+bplotdataNNQ.Std],
                        [bplotdataNNQ.Avg-bplotdataNNQ.Std,bplotdataNNQ.Avg-bplotdataNNQ.Std],alpha=0.7,color='green',edgecolor='none')
        pl.fill_between([2500,4500],np.mean(plotdataNNQ)-np.std(plotdataNNQ),np.mean(plotdataNNQ)+np.std(plotdataNNQ),alpha=0.5,color='red',edgecolor='none')
        pl.errorbar(map(GetCfgNumb,xlist),yavgNNQ,yerrNNQ,fmt='o')
        # pl.ylim(np.min(np.array(yavgNNQ)-np.array(yerrNNQ)),np.max(np.array(yavgNNQ)+np.array(yerrNNQ)))
        pl.ylim(-2.5*10**-10,2.5*10**-10)
        # pl.ylim(np.min(plotdataNNQ),np.max(plotdataNNQ))
        pl.ylabel('C2')
        pl.xlabel('icfg')
        pl.title('Monte Carlo time dependence of NNQ')
        pl.savefig('./montegraphs/MonteNNQflow'+str(MonteFlow)+'ts'+str(MonteTime)+'INg5'+INg5+'.pdf')
        pl.clf()

        # pl.scatter(map(GetCfgNumb,xlist),plotdata)
        val = bplotdata.Avg
        err = bplotdata.Std
        up = [val+err,val+err]
        down = [val-err,val-err]
        pl.fill_between([2500,4500],up,down,alpha=0.7,color='green',edgecolor='none')
        pl.fill_between([2500,4500],np.mean(plotdata)-np.std(plotdata),np.mean(plotdata)+np.std(plotdata),alpha=0.5,color='red',edgecolor='none')
        pl.errorbar(map(GetCfgNumb,xlist),yavg,yerr,fmt='o')
        # pl.ylim(np.min(np.array(yavg)-np.array(yerr)),np.max(np.array(yavg)+np.array(yerr)))
        pl.ylim(0*10**-10,3*10**-10)
        # pl.ylim(np.min(plotdata),np.max(plotdata))
        pl.ylabel('C2')
        pl.xlabel('icfg')
        pl.title('Monte Carlo time dependence of NN')
        pl.savefig('./montegraphs/MonteNNts'+str(MonteTime)+'.pdf')
        pl.clf()
    if PlotXSrcDep:
        mkdir_p('./montegraphs')
        # print 'values'
        # print np.mean(np.array(tempdataTop)[:,MonteFlow,0,MonteTime-1]),np.std(np.array(tempdataTop)[:,MonteFlow,0,MonteTime-1])
        # print
        bplotdata = np.array(Bdata)[0,MonteTime-1]
        bplotdataNNQ = np.array(TCBdata)[MonteFlow,0,MonteTime-1]
        plotdata = np.array(tempdata)[:,0,MonteTime-1]
        plotdataNNQ = np.array(tempdataTop)[:,MonteFlow,0,MonteTime-1]
        plotxsrcNNQ,plotxsrcNNQerr = [],[]
        plotxsrc,plotxsrcerr = [],[]
        bootxsrcNNQ,bootxsrcNNQerr = [],[]
        bootxsrc,bootxsrcerr = [],[]
        xlist = []
        for iXSrc in range(1,XSrcLen+1):
            yavg,yerr,yavgNNQ,yerrNNQ = [],[],[],[]
            xlist.append(iXSrc)
            for ic,(icfg,iread) in enumerate(readfilelist.iteritems()):
                # xlist += iread
                if NoXAvg:
                    yavgNNQ += plotdataNNQ[ic*len(iread):ic*len(iread)+iXSrc].tolist()
                    yavg += plotdata[ic*len(iread):ic*len(iread)+iXSrc].tolist()
                else:
                    raise IOError('NoXAvg must be True if doing PlotXSrcDep (check Params.py)')
            plotxsrcNNQ.append(np.mean(yavgNNQ))
            plotxsrcNNQerr.append(np.std(yavgNNQ))
            plotxsrc.append(np.mean(yavg))
            plotxsrcerr.append(np.std(yavg))
            bootdata,randlist = bt.CreateBoot(np.array(yavg).reshape(len(yavg),1),nboot,0)
            bootdataNNQ,randlist = bt.CreateBoot(np.array(yavgNNQ).reshape(len(yavg),1),nboot,0)
            bootxsrcNNQ.append(bootdataNNQ[0].Avg)
            bootxsrcNNQerr.append(bootdataNNQ[0].Std)
            bootxsrc.append(bootdata[0].Avg)
            bootxsrcerr.append(bootdata[0].Std)
        ncfg = len(readfilelist.keys())
        thisfitfun = CreateOORFF(bootxsrcNNQerr[0])
        Bcoeff,covar,chisqpdf = LSFit(1,np.array([xlist]),[1]*len(bootxsrcNNQerr),thisfitfun,bootxsrcNNQerr,derfun=OORNFFDer,iGuess=[bootxsrcNNQerr[0]])
        fitx = np.arange(xlist[0],xlist[-1],0.1)
        fitline = thisfitfun([fitx],Bcoeff)
        # for ix,ifit in zip(fitx,fitline):
        #     print ix,ifit
        pl.plot(fitx,fitline,color='blue',label=r'$f(n)=A+\frac{B}{n^{1/2}}\ ,\ \chi^{2}_{pdf}=$'+'{:.2e}'.format(chisqpdf))
        # pl.errorbar(xlist,[0]*len(xlist),plotxsrcNNQerr,fmt='o',label='Reg')
        pl.errorbar(np.array(xlist),[0]*len(xlist),bootxsrcNNQerr,fmt='o',label='boot='+str(nboot),color='blue')
        # pl.ylim(0,np.max(plotxsrcNNQerr+bootxsrcNNQerr))
        pl.ylim(0,np.max(bootxsrcNNQerr))
        pl.ylabel('Error')
        pl.xlabel('iXSrc#')
        pl.legend()
        pl.title('XSrc number Error of NNQ, ncfg = '+str(ncfg))
        pl.savefig('./montegraphs/XSrcErrNNQflow'+str(MonteFlow)+'ts'+str(MonteTime)+'INg5'+INg5+'.pdf')
        pl.clf()

        thisfitfun = CreateOORFF(bootxsrcerr[0])
        Bcoeff,covar,chisqpdf = LSFit(1,np.array(xlist),[1]*len(bootxsrcerr),thisfitfun,bootxsrcerr,derfun=OORNFFDer,iGuess = [bootxsrcerr[0]])
        fitx = np.arange(xlist[0],xlist[-1],0.1)
        fitline = thisfitfun([fitx],Bcoeff)
        pl.plot(fitx,fitline,color='blue',label=r'$f(n)=A+\frac{B}{n^{1/2}}\ ,\ \chi^{2}_{pdf}=$'+'{:.2e}'.format(chisqpdf))
        # pl.errorbar(xlist,[0]*len(xlist),plotxsrcerr,fmt='o',label='Reg')
        pl.errorbar(np.array(xlist),[0]*len(xlist),bootxsrcerr,fmt='o',label='boot='+str(nboot),color='blue')
        # pl.ylim(0,np.max(plotxsrcerr+bootxsrcerr))
        pl.ylim(0,np.max(bootxsrcerr))
        pl.ylabel('Error')
        pl.xlabel('iXSrc#')
        pl.title('XSrc number Error of NN, ncfg = '+str(ncfg))
        pl.legend()
        pl.savefig('./montegraphs/XSrcErrNNts'+str(MonteTime)+'.pdf')
        pl.clf()
        
        # pl.errorbar(xlist,plotxsrcNNQ,plotxsrcNNQerr,fmt='o',label='Reg')
        pl.errorbar(np.array(xlist)+0.25,bootxsrcNNQ,bootxsrcNNQerr,fmt='o',label='boot='+str(nboot))
        # pl.ylim(0,np.max(plotxsrcNNQerr+bootxsrcNNQerr))
        pl.ylabel('C2')
        pl.xlabel('iXSrc#')
        pl.legend()
        pl.title('XSrc number dependance of NNQ, ncfg = '+str(ncfg))
        pl.savefig('./montegraphs/XSrcNNQflow'+str(MonteFlow)+'ts'+str(MonteTime)+'INg5'+INg5+'.pdf')
        pl.clf()

        # pl.errorbar(xlist,plotxsrc,plotxsrcerr,fmt='o',label='Reg')
        pl.errorbar(np.array(xlist)+0.25,bootxsrc,bootxsrcerr,fmt='o',label='boot='+str(nboot))
        # pl.ylim(0,np.max(plotxsrcerr+bootxsrcerr))
        pl.ylabel('C2')
        pl.xlabel('iXSrc#')
        pl.title('XSrc number dependance of NN, ncfg = '+str(ncfg))
        pl.legend()
        pl.savefig('./montegraphs/XSrcNNts'+str(MonteTime)+'.pdf')
        pl.clf()
    if DoPlotAuto:
        plotdata = np.array(tempdata)[:,0,MonteTime-1]
        plotdataNNQ = np.array(tempdataTop)[:,MonteFlow,0,MonteTime-1]
        PlotAutoCorrDetailed(plotdata,plotdataNNQ)

        bplotdata = np.array(Bdata)[0,:]
        bplotdataNNQ = np.array(TCBdata)[MonteFlow,0,:]
        plotdata = np.array(tempdata)[:,0,:]
        plotdataNNQ = np.array(tempdataTop)[:,MonteFlow,0,:]
        PlotAutoCorr(plotdata,plotdataNNQ,'t',bplotdata,bplotdataNNQ)

        bplotdata = np.array(Bdata)[0,MonteTime-1]
        bplotdataNNQ = np.array(TCBdata)[:,0,MonteTime-1]
        plotdata = np.array(tempdata)[:,0,MonteTime-1]
        plotdataNNQ = np.array(tempdataTop)[:,:,0,MonteTime-1]
        PlotAutoCorr(plotdata,plotdataNNQ,'flow',bplotdata,bplotdataNNQ)

    return TCBdata,(Bdata,rlist),shiftlist



##NNQdata [ icfg ] 
##NNdata [ icfg  ] 
## auto_gamma = [ W ] 
def PlotAutoCorrDetailed(NNdata,NNQdata):
    mkdir_p('./montegraphs')
    auto_gamma,Cw,Gfun,Wpick,auto_error = GammaAlpha_estimate(NNQdata,NNdata,Norm=True)
    if Wpick == -1:
        print
        print 'Optimal W not found'
    pl.plot(range(len(Gfun[:3*Wpick+1])),Gfun[:3*Wpick+1],'.-')
    pl.axvline(Wpick, color='k', linestyle='-')
    pl.axhline(0.0, color='k', linestyle='--')
    pl.ylabel(r'$ \Gamma$')
    pl.xlabel('W')
    pl.title('Autocorrelation of $ \\alpha $ for nconf=' + str(len(NNdata)))
    pl.savefig('./montegraphs/AutoCorrflow'+str(MonteFlow)+'ts'+str(MonteTime)+'.pdf')
    pl.clf()
    if Debug:
        for iW,(itau,itauerr) in enumerate(zip(auto_gamma,auto_error)):
            print iW,itau,itauerr
    pl.errorbar(range(len(auto_gamma[:3*Wpick+1])),auto_gamma[:3*Wpick+1],auto_error[:3*Wpick+1],label='Error={:.2g}'.format(Cw[Wpick]))
    pl.axvline(Wpick, color='k', linestyle='-')
    pl.axhline(0.5, color='k', linestyle='--')
    pl.ylabel(r'$ \tau_{int}$')
    pl.xlabel('W')
    pl.legend()
    pl.title('Integrated Autocorrelation function for nconf=' + str(len(NNdata)))
    pl.savefig('./montegraphs/IntAutoCorrflow'+str(MonteFlow)+'ts'+str(MonteTime)+'.pdf')
    pl.clf()


##NNQdata [ icfg, t/flow ] 
##NNdata [ icfg, t/flow  ] 
##NNQboot [ icfg, t/flow ] 
##NNboot [ icfg, t/flow  ] 
def PlotAutoCorr(NNdata,NNQdata,TorFlow,NNboot,NNQboot):
    mkdir_p('./montegraphs')
    meanlist,taulist,tauerrlist,alphaerr = [],[],[],[]
    thisshift = 0.1
    if 'flow' in TorFlow:
        iNN = NNdata
        for iNNQ in np.rollaxis(np.array(NNQdata),1):
            auto_gamma,Cw,Gfun,Wpick,auto_error = GammaAlpha_estimate(iNNQ,iNN,Norm=True)
            taulist.append( auto_gamma[Wpick])
            tauerrlist.append(auto_error[Wpick])
            alphaerr.append(Cw[Wpick])
            meanlist.append(np.mean(iNNQ)/np.mean(iNN))
    else:
        for iNN, iNNQ in zip(np.rollaxis(np.array(NNdata),1),np.rollaxis(np.array(NNQdata),1)):
            auto_gamma,Cw,Gfun,Wpick,auto_error = GammaAlpha_estimate(iNNQ,iNN,Norm=True)
            taulist.append( auto_gamma[Wpick])
            tauerrlist.append(auto_error[Wpick])
            alphaerr.append(Cw[Wpick])
            meanlist.append(np.mean(iNNQ)/np.mean(iNN))

        
    pl.errorbar(range(len(taulist)),taulist,tauerrlist)
    pl.axhline(0.5, color='k', linestyle='--')
    pl.ylabel(r'$ \tau_{int}$')
    
    pl.title('$\\tau (\\alpha)$ for nconf=' + str(len(NNdata)))
    if 'flow' in TorFlow:
        mkdir_p('./montegraphs/OverFlow/')
        pl.xlabel(r'$t_{flow}$')
        pl.savefig('./montegraphs/OverFlow/IntAutoCorrts'+str(MonteTime)+'.pdf')
    else:
        mkdir_p('./montegraphs/Overts/')
        pl.xlabel(r'$t_{sep}$')
        pl.savefig('./montegraphs/Overts/IntAutoCorrflow'+str(MonteFlow)+'.pdf')
    pl.clf()
    pl.errorbar(range(len(meanlist)),meanlist,alphaerr,label='Autocorr')
    pl.errorbar(np.arange(len(NNboot))+thisshift,Pullflag(NNboot,'Avg'),Pullflag(NNboot,'Std'),label='Bootstrap')
    pl.ylabel(r'$ \alpha$')
    pl.legend()
    pl.title('$\\alpha$ for nconf=' + str(len(NNdata)))
    if 'flow' in TorFlow:
        pl.xlabel(r'$t_{flow}$')
        pl.savefig('./montegraphs/OverFlow/AutoAlphats'+str(MonteTime)+'.pdf')
    else:
        pl.xlabel(r'$t_{sep}$')
        pl.savefig('./montegraphs/Overts/AutoAlphaflow'+str(MonteFlow)+'.pdf')
    pl.clf()


    

#dataout [ iflow, igamma , ip , it ]. bs
def ReadAndBoot3ptTop(readfilelist,thisMomList,thisGammaList,thisDerGammaList,thisnboot,chargedata,chargecfglist,flowlist,printstr='',randlist=[]):
    tempdata = []
    tempdataTop = []
    for ifilepref,ifileList in readfilelist.iteritems():
        # print 'Reading '+printstr+' {}%            \r'.format(int((iconf*100)/float(len(readfilelist)))),
        try:
            if CHROMA:
                chargeindex = FileToChargeCfg(ifilepref,chargecfglist)
                if NoXAvg:
                    for ifile in ifileList:
                        if len(thisGammaList) > 0:
                            data = ReadFSCfunPickCHROMA([ifile],thisMomList,thisGammaList)
                            tempdata.append(data.data)
                            tempdataTop.append([])
                            for iflowdata in chargedata[chargeindex]:
                                tempdataTop[-1].append(np.array(data.data)*iflowdata)
                        if len(thisDerGammaList) > 0:
                            raise IOError('Chroma version does not do derivatives, make DerList in Params.py be empty')
                else:
                    if len(thisGammaList) > 0:
                        data = ReadFSCfunPickCHROMA(ifileList,thisMomList,thisGammaList)
                        tempdata.append(data.data)
                        tempdataTop.append([])
                        for iflowdata in chargedata[chargeindex]:
                            tempdataTop[-1].append(np.array(data.data)*iflowdata)
                    if len(thisDerGammaList) > 0:
                        raise IOError('Chroma version does not do derivatives, make DerList in Params.py be empty')
            else:
                raise IOError('Top Charge not implemented for non chroma results')
        except NaNCfunError as e:
            print 
            print 'Deleting file ' + ifile
            print e
            print 'MUST RE-RUN AFTER THIS TO EXCLUDE BAD CONFIGS'
            print
            # os.remove(ifile)
    # print 'Values'
    # for idata in tempdata:
    #     print idata[0][0][5]
    # print 
    if len(thisGammaList) > 0:
        # print 
        # print 'TC',np.array(BootSet3ptTC(tempdataTop,thisMomList,thisGammaList,thisnboot,flowlist,printstr='',randlist=randlist)).shape
        # print '3pt',np.array(BootSet3pt(tempdata,thisMomList,thisGammaList,thisnboot,printstr='',randlist=randlist)).shape
        return (BootSet3ptTC(tempdataTop,thisMomList,thisGammaList,thisnboot,flowlist,printstr='',randlist=randlist),
                BootSet3pt(tempdata,thisMomList,thisGammaList,thisnboot,printstr='',randlist=randlist))
    else:
        raise IOError('No Gammas?')
    # elif len(thisDerGammaList) > 0:
    #     return BootSet3pt(tempdata,thisMomList,thisDerGammaList,thisnboot,printstr='',randlist=randlist)


def FileToChargeCfg(ifile,chargecfglist):
    for ic,icharge in enumerate(chargecfglist):
        if icharge in ifile:
            return ic
    print 'charge cfg not found'
    return None



# #dataout [ igamma , ip , it ]. bs
# def ReadAndBoot3pt(readfilelist,thisMomList,thisGammaList,thisDerGammaList,thisnboot,shiftlist,printstr='',randlist=[]):
#     tempdata = []
#     counter = -1
#     for iconf,ifile in enumerate(readfilelist):
#         # print 'Reading '+printstr+' {}%            \r'.format(int((iconf*100)/float(len(readfilelist)))),
#         try:
#             if CHROMA:
#                 if xsrcList[0] in ifile or not XAvg:
#                     counter += 1
#                     if len(thisGammaList) > 0:
#                         tempdata.append(ReadFSCfunPickCHROMA(ifile,thisMomList,thisGammaList,srcshift=shiftlist[counter]).data)
#                     if len(thisDerGammaList) > 0:
#                         raise IOError('Chroma version does not do derivatives, make DerList in Params.py be empty')
#             else:
#                 if len(thisGammaList) > 0:
#                     tempdata.append(ReadFSCfunPick(ifile,thisMomList,thisGammaList).data)
#                 if len(thisDerGammaList) > 0:
#                     tempdata.append(ReadFSDerCfunPick(ifile,thisMomList,thisDerGammaList).data)
#         except NaNCfunError as e:
#             print 
#             print 'Deleting file ' + ifile
#             print e
#             print 'MUST RE-RUN AFTER THIS TO EXCLUDE BAD CONFIGS'
#             print
#             # os.remove(ifile)
#     if len(thisGammaList) > 0:
#         return BootSet3pt(tempdata,thisMomList,thisGammaList,thisnboot,printstr='',randlist=randlist)
#     elif len(thisDerGammaList) > 0:
#         return BootSet3pt(tempdata,thisMomList,thisDerGammaList,thisnboot,printstr='',randlist=randlist)


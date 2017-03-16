#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl
from copy import deepcopy

def MyCorrelate(x,y,Norm=True,MinAvg=True):
    if MinAvg:
        x = x-np.mean(x)
        y = y-np.mean(y)
    listout = []
    for it in xrange(len(x)):
        listout.append(0.0)
        for index in xrange(len(x)):
            if it+index < len(x):
                listout[-1] = listout[-1] + x[index]*y[index+it]
        if Norm: listout[-1] = listout[-1]/(len(x)-it)
    return np.array(listout)


### autocorrelation work taken from https://arxiv.org/pdf/hep-lat/0306017.pdf
def autocorr(x,y):
   """
   http://stackoverflow.com/q/14297012/190597
   http://en.wikipedia.org/wiki/Autocorrelation#Estimation
   """
   x = np.array(x)
   y = np.array(y)
   
   n = len(x)
   # variance = x.var()
   x = x-x.mean()
   y = y-y.mean()
   r = np.correlate(x, y, mode = 'full')[-n:]
   # assert np.allclose(r, np.array([(x[:n-k]*x[-(n-k):]).sum() for k in xrange(n)]))
   result = r/np.arange(n, 0, -1)
   return result


def gW(tauW,thisN):
    ## using auto fitting window method used from (52)
    for iW,it in enumerate(tauW):
        if iW == 0: continue
        val = np.exp(-iW/it)-it/np.sqrt(iW*thisN)
        if val < 0.0:
            return iW
    return -1
        
def VarTau(tau,N):
    ## Using aproximate formula (42) from paper
    return [4/N * (iW + 0.5 - itau) * itau**2 for iW,itau in enumerate(tau)]

def BiasCorrect(CfW,N):
    ## Bias corrections using (49)
    W = np.arange(len(CfW))
    return CfW*(1+((2*W+1)/float(N)))
    # ## Testing
    # return CfW
    

#data = [ variable , monte time ]
#fun(variables) output is value
#funder(variables) output is partial derivatives (w.r.t variables)
## AllOut = False
## average, error[Woptimal], tau[Woptimal], tauerr[Woptimal]
## AllOut = True
## average, error(W), tau(W), tauerr(W), GFt(W), Wopt
def uWerrMine(data,fun,funder,Sparam=1.5,AllOut=False,plot=False):
    data = np.array(data)
    glen = len(data[0])
    avgdata = [idata.mean() for idata in data]

    ## (31) matrix of autocorrelations w.r.t ab= variables
    G_ab_t = []
    for adat in data:
        G_ab_t.append([ autocorr(adat,bdat) for bdat in data])
        # G_ab_t.append([ MyCorrelate(adat,bdat) for bdat in data])

        
    ## (33) alpha function derivates (w.r.t variables)
    f_a = np.array(funder(*avgdata))
    f_ab = []
    for if_a in f_a:
        f_ab.append([if_b*if_a for if_b in f_a])
    f_ab = np.array(f_ab)

    ## (33)
    GFt = [np.sum(f_ab * G_ab) for G_ab in np.rollaxis(np.array(G_ab_t),-1) ]

    ## (35)
    CFW = np.array([GFt[0]]+[GFt[0] + 2*np.sum(GFt[1:W+1]) for W in xrange(1,len(GFt))])

    # ## Bias corrections (49)
    # CFW = BiasCorrect(CFW,glen)


    ## (34)
    nuF = CFW[0]

    
    ## equation (41)
    tauint = CFW / (2*nuF)


    
    ## From Paper
    tauintpass = deepcopy(tauint)
    for it,itauint in enumerate(tauint):
        if itauint <= 0.5:
            tauintpass[it] = 0.00000001
    ## (51) 
    tau = Sparam/np.log((2*tauintpass+1)/(2*tauintpass-1))

    # ## From Matlab
    # tau = []
    # for iGFt in GFt:
    #     if iGFt <= 0.0:
    #         tau.append(0.0000001)
    #     else:
    #         tau.append(Sparam/np.log((iGFt+1)/iGFt))

    ## (52)
    Wopt = gW(np.array(tau),glen)

    ## (42)
    dtauint = VarTau(tauint,glen)

    if plot != False:
        xmax = int(Wopt*3)
        step = int(np.ceil(Wopt/20)) or 1
        thisshift = xmax*0.002
        fig = pl.figure(1)
        Gplt= fig.add_subplot(211)
        Gplt.set_ylabel(r'$\Gamma$')
        Gplt.set_xlabel('$W$')
        GFtplot = GFt[:xmax:step]/GFt[0]
        # GFtplot = GFt[:xmax:step]
        pl.errorbar(range(len(GFtplot)), GFtplot,fmt="o", color='b')
        pl.axvline(Wopt+thisshift, color='r')
        tplt = fig.add_subplot(212)
        tplt.set_ylabel(r'$\tau_{\mathrm{int}}$')
        tplt.set_xlabel(r'$W$')
        # tauintplot = tauint[:xmax:step]
        tauintplot = CFW[:xmax:step]
        pl.errorbar(range(len(tauintplot)), tauintplot,
                     dtauint[:xmax:step], fmt="o", color='b')
        pl.axvline(Wopt+thisshift, color='r')
        if plot == True:
            pl.show()
        else:
            pl.savefig(plot+'.pdf')
            pl.clf()
        
    if AllOut:
        return fun(*avgdata),np.sqrt(np.abs(CFW)/float(glen)),tauint,dtauint,GFt,Wopt
    else:
        return fun(*avgdata),np.sqrt(np.abs(CFW[Wopt])/float(glen)),tauint[Wopt],dtauint[Wopt]
    
    
# def GammaAlpha_estimate(gQ,gN):
#    gQ = np.array(gQ)
#    gN = np.array(gN)
#    glen = len(gQ)
#    gQAvg = gQ.mean()
#    gNAvg = gN.mean()
   
#    GQQt=autocorr(gQ,gQ)
#    GNQt=autocorr(gN,gQ)
#    GQNt=autocorr(gQ,gN)
#    GNNt=autocorr(gN,gN)
   
#    ##alpha function derivates wrt NNQ and NN for NNQ/NN
#    fQ = gNAvg**(-1)
#    fN = -gQAvg/(gNAvg**2)
   
#    ## equation (33)
#    Gat = fQ**2 * GQQt + (fQ*fN * (GQNt + GNQt)) + fN**2 * GNNt
   
#    # if Norm: Gat = np.array(Gat)/Gat[0]

#    ## equation (35)
#    CaW = [Gat[0] + 2*np.sum(Gat[1:W]) for W in xrange(1,len(Gat))]
#    ## equation (41)
#    tau = np.array(CaW) / (2*Gat[0])
#    Wopp = gW(tau)
#    return np.array(tau),np.sqrt(np.abs(CaW[Wopp])/float(glen)),Gat,Wopp, np.array(VarTau(tau))


# def Gamma1D_est(data):
#    data = np.array(data)
#    glen = len(data)
#    dataAvg = data.mean()
   
#    GQQt=autocorr(data,data)
   
#    ##alpha function derivates wrt NNQ and NN for NNQ/NN
#    fQ = 1
   
#    ## equation (33)
#    Gat = GQQt   
#    # if Norm: Gat = np.array(Gat)/Gat[0]

#    ## equation (35)
#    CaW = [Gat[0] + 2*np.sum(Gat[1:W]) for W in xrange(1,len(Gat))]
#    ## equation (41)
#    tau = np.array(CaW) / (2*Gat[0])
#    Wopp = gW(tau)
#    return np.array(tau),np.sqrt(np.abs(CaW[Wopp])/float(glen)),Gat,Wopp, np.array(VarTau(tau))

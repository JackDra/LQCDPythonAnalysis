#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl

# def MyCorrelate(x,y,Norm=True,MinAvg=True):
#     if MinAvg:
#         x = x-np.mean(x)
#         y = y-np.mean(y)
#     listout = []
#     for it in xrange(len(x)):
#         listout.append(0.0)
#         for index in xrange(len(x)):
#             if it+index < len(x):
#                 listout[-1] = listout[-1] + x[index]*y[index+it]
#         if Norm: listout[-1] = listout[-1]/(len(x)-it)
#     return np.array(listout)


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


def gW(tauW):
    thisN = len(tauW)
    ## using auto fitting window method used from (52)
    for iW,it in enumerate(tauW):
        if iW == 0: continue
        val = np.exp(-iW/it)-it/np.sqrt(iW*thisN)
        if val < 0.0:
            return iW
    return -1
        
def VarTau(tau):
    ## Using aproximate formula (42) from paper
    ## starting from W = 1, (iW = W -1, need to add 1 to start from iW = 1)
    N = float(len(tau))
    return [0] + [4/N * (iW + 1.5 - itau) * itau**2 for iW,itau in enumerate(tau[1:])]

def BiasCorrect(CfW,W,N):
    return CfW*(1+(2*W+1)/N)

#data = [ variable , monte time ]
#fun(variables) output is value
#funder(variables) output is partial derivatives (w.r.t variables)
## AllOut = False
## average, error[Woptimal], tau[Woptimal], tauerr[Woptimal]
## AllOut = True
## average, error(W), tau(W), tauerr(W), Gat(W), Wopt
def uWerrMine(data,fun,funder,Sparam=1.5,AllOut=False,plot=False):
    data = np.array(data)
    glen = len(data[0])
    avgdata = [idata.mean() for idata in data]
    
    G_ab = []
    for adat in data:
        G_ab.append([ autocorr(adat,bdat) for bdat in data])

    ##alpha function derivates (w.r.t variables)
    f_a = np.array(funder(*avgdata))
    f_ab = []
    for if_a in f_a:
        f_ab.append([if_b*if_a for if_b in f_a])
    f_ab = np.array(f_ab)

    ## equation (33)
    Gat = [np.sum(f_ab * G_ab_t) for G_ab_t in np.rollaxis(np.array(G_ab),-1) ]
    
    ## equation (35)
    CaW = BiasCorrect(np.array([Gat[0] + 2*np.sum(Gat[1:W]) for W in xrange(1,len(Gat))]),np.arange(1,len(Gat)),glen)
    ## equation (41)
    tauint = np.array(CaW) / (2*Gat[0])
    tauintpass = tauint
    for it,itauint in enumerate(tauint):
        if itauint <= 0.5:
            tauintpass[it] = 0.000001
    tau = Sparam/np.log((2*tauint+1)/(2*tauint-1))
    Wopt = gW(tau)
    dtauint = VarTau(tauint)
    if plot != False:
        xmax = int(Wopt*2)
        step = int(np.ceil(Wopt/20)) or 1
        thisshift = xmax*0.002
        fig = pl.figure(1)
        Gplt= fig.add_subplot(211)
        Gplt.set_ylabel(r'$\Gamma$')
        Gplt.set_xlabel('$W$')
        Gatplot = Gat[:xmax:step]/Gat[0]
        pl.errorbar(range(Gatplot), Gatplot,
                     fmt="o", color='b')
        pl.axvline(Wopt+thisshift, color='r')
        tplt = fig.add_subplot(212)
        tplt.set_ylabel(r'$\tau_{\mathrm{int}}$')
        tplt.set_xlabel(r'$W$')
        tauintplot = tauint[:xmax:step]
        pl.errorbar(range(tauintplot), tauintplot,
                     dtauint[:xmax:step], fmt="o", color='b')
        pl.axvline(Wopt+thisshift, color='r')
        if plot == True:
            pl.show()
        else:
            pl.savefig(plot+'.pdf')
            pl.clf()
        
    if AllOut:
        return fun(*avgdata),np.sqrt(np.abs(CaW)/float(glen)),tauint,dtauint,Gat,Wopt
    else:
        return fun(*avgdata),np.sqrt(np.abs(CaW[Wopt])/float(glen)),tauint[Wopt],dtauint[Wopt]
    
    
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

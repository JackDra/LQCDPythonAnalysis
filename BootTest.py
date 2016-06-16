#!/usr/bin/env python
#from pylab import *
# from Defs import *
#from Functions import *
from numpy import zeros, size, sort, exp,average,split,random,append,multiply,sum,sqrt,log
#from Numeric import zeros
#from ROOT import TF1, TGraphErrors, TVirtualFitter
#import random
#from scipy import exp
from numpy import array as narray
from array import array

startseed=1234

# This could be useful
# v = numpy.random.normal(mu,sigma,10000)

def UpdateAvgStd(data):
    data['FitBoot'].Stats()
    data['FitAvg'],data['FitStd'] = data['FitBoot'].Avg,data['FitBoot'].Std
    return data

def CreateBoot(RawData,nboot,confidence,ifold=0,bin=1,weights=[],randlist= []):
    nconf=len(RawData)
    if len(weights)==nconf:
        print "Found Weights"
    nt=len(RawData[0])
    bootlist=[]
    bappend=bootlist.append
    for it in xrange(nt):
        tmpboot=BootStrap1(nboot,confidence)
        if it==0 or ifold==0:
            randlist = tmpboot.Import([RawData[iconf][it] for iconf in xrange(nconf)],bin,weights,randlist=randlist)
        elif it!=0 and ifold==1:
            randlist = tmpboot.Import([0.5*(RawData[iconf][it]+RawData[iconf][nt-it]) for iconf in xrange(nconf)],bin,randlist=randlist)
        elif it!=0 and ifold==2:
            randlist = tmpboot.Import([0.5*(RawData[iconf][it]-RawData[iconf][nt-it]) for iconf in xrange(nconf)],bin,randlist=randlist)
        bappend(tmpboot)
        bootlist[-1].Stats()
    return bootlist,randlist

def CreateCplxBoot(RawData,nboot,confidence,ifold=0,bin=1):
    nsp=len(RawData)/2
    nconf=len(RawData[0])
    nt=len(RawData[0][0])
    bootlist=[]
    bootlist1=[]
    bootlist2=[]
    bappend1=bootlist1.append
    bappend2=bootlist2.append
    for it in xrange(nt):
        sumbootr=BootStrap1(nboot,confidence)
        sumbooti=BootStrap1(nboot,confidence)
        sumbootr.Constant(0.0)
        sumbooti.Constant(0.0)
        for isp in xrange(nsp):
            tmpbootr=BootStrap1(nboot,confidence)
            tmpbooti=BootStrap1(nboot,confidence)
            if it==0 or ifold==0:
               tmpbootr.Import([RawData[2*isp][iconf][it].real for iconf in xrange(nconf)],bin)
               tmpbooti.Import([RawData[2*isp][iconf][it].imag for iconf in xrange(nconf)],bin)
            elif it!=0 and ifold==1:
                tmpbootr.Import([0.5*(RawData[2*isp][iconf][it].real+RawData[2*isp][iconf][nt-it].real) for iconf in xrange(nconf)],bin)
                tmpbooti.Import([0.5*(RawData[2*isp][iconf][it].imag+RawData[2*isp][iconf][nt-it].imag) for iconf in xrange(nconf)],bin)
            elif it!=0 and ifold==2:
                assert 1==0,["Not supported"]
            sumbootr=sumbootr+(1./nsp)*tmpbootr
            sumbooti=sumbooti+(1./nsp)*tmpbooti
        bappend1([sumbootr,sumbooti])
        sumbootr=BootStrap1(nboot,confidence)
        sumbooti=BootStrap1(nboot,confidence)
        sumbootr.Constant(0.0)
        sumbooti.Constant(0.0)
        for isp in xrange(nsp):
            tmpbootr=BootStrap1(nboot,confidence)
            tmpbooti=BootStrap1(nboot,confidence)
            if it==0 or ifold==0:
                tmpbootr.Import([RawData[2*isp+1][iconf][it].real for iconf in xrange(nconf)],bin)
                tmpbooti.Import([RawData[2*isp+1][iconf][it].imag for iconf in xrange(nconf)],bin)
            elif it!=0 and ifold==1:
                tmpbootr.Import([0.5*(RawData[2*isp+1][iconf][it].real+RawData[2*isp+1][iconf][nt-it].real) for iconf in xrange(nconf)],bin)
                tmpbooti.Import([0.5*(RawData[2*isp+1][iconf][it].imag+RawData[2*isp+1][iconf][nt-it].imag) for iconf in xrange(nconf)],bin)
            elif it!=0 and ifold==2:
                assert 1==0,["Not supported"]
            sumbootr=sumbootr+(1./nsp)*tmpbootr
            sumbooti=sumbooti+(1./nsp)*tmpbooti
        bappend2([sumbootr,sumbooti])
        bootlist1[-1][0].Stats()
        bootlist1[-1][1].Stats()
        bootlist2[-1][0].Stats()
        bootlist2[-1][1].Stats()
    bootlist.append(bootlist1)
    bootlist.append(bootlist2)
    return bootlist

def CreateDoubleBoot(RawData1,RawData2,nboot,confidence,ifold=0,bin=1):
    nconf=len(RawData1)
    assert nconf==len(RawData2)
    nt=len(RawData1[0])
    assert nt==len(RawData2[0])
    bootlist=[]
    for it in xrange(nt):
        tmpboot=BootStrap1(nboot,confidence)
        if it==0 or ifold==0:
            if ifold==3:
                tmpboot.Import([RawData2[iconf][it] for iconf in range(nconf)])
            else:
                tmpboot.Import([RawData1[iconf][it] for iconf in range(nconf)])
        elif it!=0 and ifold==1:
            tmpboot.Import([0.5*(RawData1[iconf][it]+RawData2[iconf][nt-it]) for iconf in range(nconf)])
        elif it!=0 and ifold==2:
            tmpboot.Import([0.5*(RawData1[iconf][it]-RawData2[iconf][nt-it]) for iconf in range(nconf)])
        elif it!=0 and ifold==3:
            tmpboot.Import([RawData2[iconf][nt-it] for iconf in range(nconf)])
        bootlist.append(tmpboot)
        bootlist[-1].Stats()
    return bootlist

class BootStrap1:
    def __init__(self,Nboot,Confidence):
        self.nboot=Nboot
        self.confidence=Confidence
        self.Avg=0.0
        self.Std=0.0
        self.values=zeros((self.nboot))
        self.Raw=zeros(0)
    def Stats(self):
        self.Avg=average(self.values)
        tmp=(self.values-self.Avg)**2
        self.Std=sqrt(average(tmp))
        return
    def OldStats(self):
        tmp=sort(self.values,0)
        omit=(100-self.confidence)/2
        ilo=(omit*self.nboot)/100
        ihi=(self.nboot-1-(omit*self.nboot)/100)
        self.Std=(tmp[ihi]-tmp[ilo])/2
        return
    def Import(self,data,bin=1,weights=[],randlist = []):
        if len(weights)>0 and bin!=1:
            print "Weights only currently supported for binsize=1"
            assert 1==0
        self.nconf=len(data)/bin
        if bin>1:
            tmpdata=narray([average(ien) for ien in split(narray(data[0:bin*(len(data)/bin)]),self.nconf)])
        else:
            tmpdata=narray(data)
#        self.Raw=append(self.Raw,tmpdata)
#        self.nconf=len(self.Raw)
#        rint=zeros( (self.nboot,self.nconf))
        # if randlist == []:
        myseed=startseed*self.nconf/self.nboot
        random.seed(myseed)
        #        locranint=random.randint
        locranint=random.random_integers
        # else:
        #     locranint = randlist
        if len(weights)>0:
            tmpweight=narray(weights)
            self.Avg=average(multiply(tmpdata,tmpweight))/sum(tmpweight)
        else:
            self.Avg=average(tmpdata)
        for iboot in range(self.nboot):
            rint=locranint(0,self.nconf-1,self.nconf)
#            tmp=0.0
#            for iconf in range(self.nconf):
#                rint=locranint(0,self.nconf-1)
#                tmp+=tmpdata[rint]/self.nconf
#            self.values[iboot]=tmp
            if len(weights)>0:
                tw2=tmpweight[rint]                
                td2=multiply(tmpdata[rint],tw2)
                self.values[iboot]=average(td2)/sum(tw2)
#average(multiply(tmpdata[rint],tmpweight[rint]))/sum(tmpweight[rint])
            else:
                self.values[iboot]=average(tmpdata[rint])
        return locranint
    def Create(self,mean,stdev):
        self.Avg=mean
        self.Std=stdev
#        for iboot in range(self.nboot):
#            self.values[iboot]=random.gauss(mean,stdev)
        self.values=random.normal(mean,stdev,self.nboot)
    def Constant(self,const):
        self.Avg=const
        self.values.fill(const)
        return
    def write(self,file):
        iclose=False
        try:
            # If file is not open, open it as a fileobject
            fo=open(file,'wb')
            iclose=True
        except:
            # If already a fileobject, we will append to it
            fo=file
        narray(self.nboot).tofile(fo)
        narray(self.Avg).tofile(fo)
        narray(self.Std).tofile(fo)
        self.values.tofile(fo)
        
        # If file was a new file, close it
        if iclose:
            fo.close()
    def read(self,file):
        iclose=False
        try:
            # If file is not open, open it as a fileobject
            fo=open(file,'rb')
            iclose=True
        except:
            # If already a fileobject, we will append to it
            fo=file
        ti=array('i')
        ti.read(fo,1)
        self.nboot=ti[0]
        tf=array('d')
        tf.read(fo,2)
        self.Avg=tf[0]
        self.Std=tf[1]
        tf=array('d')
        tf.read(fo,self.nboot)
        self.values=narray(tf)
        
        # If file was a new file, close it
        if iclose:
            fo.close()
    def sqrt(self):
        self.Avg = sqrt(self.Avg)
        self.values = sqrt(self.values)
    def __mul__(self,fac):
        result=BootStrap1(self.nboot, self.confidence)
        try:
            real=float(fac)
            result.Avg=self.Avg*real
            result.values=self.values*real            
        except:
            try:
                tnboot=fac.nboot
                result.Avg=self.Avg*fac.Avg
                result.values=self.values*fac.values
            except:
                print "ERROR: UNknown boot multiply"
                print fac.nboot
                assert 1==0
        return result
    def __rmul__(self,real=float(1)):
        result=BootStrap1(self.nboot, self.confidence)
        result.Avg=self.Avg*real
        result.values=self.values*real
        return result
    def __imul__(self,real=float(1)):
        result=BootStrap1(self.nboot, self.confidence)
        result.Avg=self.Avg*real
        result.values=self.values*real
        return result
    def __div__(self,fac):
        result=BootStrap1(self.nboot, self.confidence)
        try:
            real=float(fac)
            result.Avg=self.Avg/real
            result.values=self.values/real            
        except:
            try:
                print fac.nboot
                tnboot=fac.nboot
                result.Avg=self.Avg/fac.Avg
                result.values = []
                for ival,ival2 in zip(self.values,fac.values):
                    result.values=ival/ival2
            except:
                print "ERROR: UNknown boot divide"
                print tnboot,fac.nboot
                assert 1==0
        return result
    def __add__(self,fac):
        result=BootStrap1(self.nboot, self.confidence)
        try:
            real=float(fac)
            result.Avg=self.Avg+real
            for iboot in range(self.nboot):
                result.values[iboot]=self.values[iboot]+real            
        except:
            try:
                real=complex(fac)
                result.Avg=self.Avg+real
                for iboot in range(self.nboot):
                    result.values[iboot]=self.values[iboot]+real            
            except:
                try:
                    tnboot=fac.nboot
                    result.Avg=self.Avg+fac.Avg
                    for iboot in range(self.nboot):
                        result.values[iboot]=self.values[iboot]+fac.values[iboot]
                except:
                    assert 1==0
                    print "ERROR: UNknown boot addition"
        return result
    def __radd__(self,real=float(1)):
        result=BootStrap1(self.nboot, self.confidence)
        result.Avg=self.Avg+real
        for iboot in range(self.nboot):
            result.values[iboot]=self.values[iboot]+real
        return result
    def __iadd__(self,real=float(1)):
        result=BootStrap1(self.nboot, self.confidence)
        result.Avg=self.Avg+real
        for iboot in range(self.nboot):
            result.values[iboot]=self.values[iboot]+real
        return result
    def __sub__(self,fac):
        result=BootStrap1(self.nboot, self.confidence)
        try:
            real=float(fac)
            result.Avg=self.Avg-real
            for iboot in range(self.nboot):
                result.values[iboot]=self.values[iboot]-real            
        except:
            try:
                tnboot=fac.nboot
                result.Avg=self.Avg-fac.Avg
                for iboot in range(self.nboot):
                    result.values[iboot]=self.values[iboot]-fac.values[iboot]
            except:
                assert 1==0
                print "ERROR: UNknown boot multiply"
        return result
    def __rsub__(self,real=float(1)):
        result=BootStrap1(self.nboot, self.confidence)
        result.Avg=real-self.Avg
        for iboot in range(self.nboot):
            result.values[iboot]=real-self.values[iboot]
        return result
    def __isub__(self,real=float(1)):
        result=BootStrap1(self.nboot, self.confidence)
        result.Avg=self.Avg-real
        for iboot in range(self.nboot):
            result.values[iboot]=self.values[iboot]-real
        return result
    def __pow__(self,real):
        result=BootStrap1(self.nboot, self.confidence)
        result.Avg=self.Avg**real
#        for iboot in range(self.nboot):
        result.values=self.values**real
        return result
    def exp(self,real):
        result=BootStrap1(self.nboot, self.confidence)
        result.Avg=exp(self.Avg*real)
        for iboot in range(self.nboot):
            result.values[iboot]=exp(self.values[iboot]*real)
        return result
    def log(self):
        result=BootStrap1(self.nboot, self.confidence)
        result.Avg=log(self.Avg)
        for iboot in range(self.nboot):
            result.values[iboot]=log(self.values[iboot])
        return result
    def __abs__(self):
        result=BootStrap1(self.nboot, self.confidence)
        result.Avg=abs(self.Avg)
        result.values=abs(self.values)
        return result
    def __ge__(self,real):
        tav=False
        tboot=False
        if self.Avg>=real:
            tav=True
        if (self.values>=real).all():
            tboot=True
        result=tav*tboot
        return result

class BootStrap2:
    def __init__(self,Nboot,Ndim,Confidence):
        self.nboot=Nboot
        self.ndim=Ndim
        self.confidence=Confidence
        self.Avg=zeros((self.ndim))
        self.Std=zeros((self.ndim))
        self.values=zeros( (self.nboot,self.ndim))
    def Import(self,data):
        self.nconf=size(data,0)
        rint=zeros( [self.nboot,self.nconf])
        myseed=startseed*self.nconf/self.nboot
        random.seed(myseed)
        locranint=random.random_integers
#
# Note: These changes haven't been tested!
#
#        for iconf in range(self.nconf):
#            self.Avg+=data[iconf]/self.nconf
        self.Avg=average(data,axis=0)
        for iboot in range(self.nboot):
            rint=locranint(0,self.nconf-1,self.nconf)
#            tmp=zeros((self.ndim))
#            for iconf in range(self.nconf):
#                rint=random.randint(0,self.nconf-1)
#                tmp+=data[rint]/self.nconf
#            self.values[iboot]=tmp
            self.values[iboot]=average(tmpdata[rint],axis=0)
    def Stats(self):
        tmp=sort(self.values,0)
        omit=(100-self.confidence)/2
        ilo=(omit*self.nboot)/100
        ihi=(self.nboot-1-(omit*self.nboot)/100)
        for it in range(self.ndim):
            self.Std[it]=(tmp[ihi][it]-tmp[ilo][it])/2
        return
    def write(self,file):
        iclose=False
        try:
            # If file is not open, open it as a fileobject
            fo=open(file,'wb')
            iclose=True
        except:
            # If already a fileobject, we will append to it
            fo=file
        narray(self.nboot).tofile(fo)
        narray(self.ndim).tofile(fo)
        self.Avg.tofile(fo)
        self.Std.tofile(fo)
        self.values.tofile(fo)
        
        # If file was a new file, close it
        if iclose:
            fo.close()
    def read(self,file):
        iclose=False
        try:
            # If file is not open, open it as a fileobject
            fo=open(file,'rb')
            iclose=True
        except:
            # If already a fileobject, we will append to it
            fo=file
        ti=array('i')
        ti.read(fo,2)
        self.nboot=ti[0]
        self.ndim=ti[1]
        self.values.resize(self.nboot,self.ndim)
        tf=array('d')
        tf.read(fo,self.ndim)
        self.Avg=narray(tf)
        tf=array('d')
        tf.read(fo,self.ndim)
        self.Std=narray(tf)
        for iboot in xrange(self.nboot):
            tf=array('d')
            tf.read(fo,self.ndim)
            self.values[iboot]=narray(tf)
        
        # If file was a new file, close it
        if iclose:
            fo.close()
    def FitInit(self,xin,xinerr,xinmin,xinmax,fitflg):
        self.__x=array('f')
        self.__xerr=array('f')
        for ix in range(xinmin,xinmax+1):
            self.__x.append(xin[ix])
            self.__xerr.append(xinerr[ix])
        self.ndat=len(self.__x)
        self.__xmin=min(self.__x)
        self.__xmax=max(self.__x)
        print self.ndat
        if fitflg==POLE:
            self.fitfnc=monopole()
        elif fitflg==DIPOLE:
            self.fitfnc=dipole()
        elif fitflg==TRIPOLE:
            self.fitfnc=tripole()
        elif fitflg==CNST:
            self.fitfnc=constant()
        elif fitflg==LINR:
            self.fitfnc=linear()
        elif fitflg==QUAD:
            self.fitfnc=quad()
        elif fitflg==CUBIC:
            self.fitfnc=cubic()
        return
    def Fit(self):
        fitparams=[]
        fitt = TF1( "fitt", self.fitfnc.fnc, self.__xmin, self.__xmax,self.fitfnc.npar)
        for ipar in range(self.fitfnc.npar):
            fitt.SetParameter( ipar,self.fitfnc.initpar[ipar])
            tmpboot=BootStrap1(self.nboot,self.confidence)
            fitparams.append(tmpboot)

# Initialise the graph
        gra = TGraphErrors(self.ndat, self.__x, self.Avg, self.__xerr, self.Std)
# Fit the graph with the predefined "fitt" function
        gra.Fit(fitt,"V","R")
# Access the fit resuts
        results=fitt.GetParameters()
        errors=fitt.GetParErrors()        
# Access elements of correlation matrix
        self.vij=[]
        fitter = TVirtualFitter.GetFitter()
        for ires in range(self.fitfnc.npar):
            fitparams[ires].Avg=results[ires]
            mres=ires+1
            while mres<self.fitfnc.npar:
                self.vij.append(fitter.GetCovarianceMatrixElement(ires,mres)/sqrt(fitter.GetCovarianceMatrixElement(ires,ires)*fitter.GetCovarianceMatrixElement(mres,mres)))
                mres=mres+1
            
        self.chi2ndf=fitt.GetChisquare()/fitt.GetNDF()        
# Fit for each bootstrap ensemble
        for iboot in range(self.nboot):
            tmpy=array('f')
            for idim in range(self.ndim):
                tmpy.append(self.values[iboot][idim])
# Initialise the graph
            gra = TGraphErrors(self.ndat, self.__x, tmpy, self.__xerr, self.Std)
# Fit the graph with the predefined "fitt" function
            gra.Fit(fitt,"R")
# Access the fit resuts
            results=fitt.GetParameters()
            for ires in range(self.fitfnc.npar):
                fitparams[ires].values[iboot]=results[ires]
        for ires in range(self.fitfnc.npar):
            fitparams[ires].Stats()
        return fitparams

#!/usr/bin/env python
from array import array
import numpy as np
from Params import *
from FFSympy import GammaTOChroma

SeekIncSize = 8
ChromaSIS = 16
bar3ptChromaByteSwap = True
# complextype = np.complex64
complextype = np.complex64
CheckMagic = False

def GetBar3ptLoc(ig,ip,tsinklen,momlistlen):
    intlen = np.dtype(np.int32).itemsize
    strlen = np.dtype('S13').itemsize
    cmplxlen = np.dtype(complextype).itemsize
    cmplxblocklen = cmplxlen*tsinklen+intlen
    # cmplxblocklen = intlen*16+intlen
    headerlen = intlen*27+strlen
    magicNumshift = intlen*6
    magicNumshiftCons = intlen*5
    momblocklen = cmplxblocklen + magicNumshift
    momblocklenCons = cmplxblocklen*2 + magicNumshiftCons
    momblockList = [momblocklen,momblocklenCons,momblocklenCons,momblocklen,momblocklenCons,
                    momblocklen,momblocklen,momblocklenCons,momblocklenCons,momblocklen,
                    momblocklen,momblocklenCons,momblocklen,momblocklenCons,momblocklenCons,
                    momblocklen,momblocklen,momblocklen,momblocklen,momblocklen]
    
    gammablocklen = momblocklen*momlistlen + intlen*2
    gammablocklenCons = momblocklenCons*momlistlen + intlen*2
    gammablockList = [gammablocklen,gammablocklenCons,gammablocklenCons,gammablocklen,gammablocklenCons,
                      gammablocklen,gammablocklen,gammablocklenCons,gammablocklenCons,gammablocklen,
                      gammablocklen,gammablocklenCons,gammablocklen,gammablocklenCons,gammablocklenCons,
                      gammablocklen,gammablocklen,gammablocklen,gammablocklen,gammablocklen]
    outloc = headerlen + sum(gammablockList[:ig]) + momblockList[ig]*ip
    return outloc , outloc + cmplxblocklen,magicNumshift
        
    
## noDer[igamma , ip , it]
## forcent is total length, not array index!
class ReadFSCfunPickCHROMA:
    def __init__(self,thisfile,thisMomList,thisGammaList,forcent=nt,srcshift=None):
        self.data = []
        if XAvg: thisxsrcList = xsrcList
        else: thisxsrcList = [xsrcList[0]]
        if srcshift == None: srcshift = [0]*len(thisxsrcList)
        for xsrc,isshift in zip(thisxsrcList,srcshift):
            # if thistsink == 'FileName':
            #     thistsink = int(re.findall('tsink.*?p',thisfile)[0].replace('tsink','').replace('p',''))
            datahold = []
            corr_shift = isshift
                
            f = open(thisfile.replace(xsrcList[0],xsrc),'rb')
            for igamma,thisgamma in enumerate(thisGammaList):
                igammaloc = GammaTOChroma(thisgamma)
                datahold.append([])
                for ip,iploc in enumerate(thisMomList):      
                    loc,loccons,magmin = GetBar3ptLoc(igammaloc,iploc,forcent,len(qvecSet))
                    magicloc = loc-magmin
                    if 'Cons' in thisgamma or (RepWithCons and ('g1' in thisgamma or 'g2' in thisgamma or 'g3' in thisgamma or 'g4' in thisgamma )):
                        loc = loccons
                    if CheckMagic:
                        f.seek(magicloc)
                        MagicList =  np.fromfile(f,dtype=np.int32,count=1).byteswap()
                        if MagicList[0] != 20301:
                            f.seek(magicloc)
                            MagicList =  np.fromfile(f,dtype=np.int32,count=100).byteswap()
                            print MagicList
                            raise IOError('Magic Number not found for ' +thisgamma+' '+ipTOqstr(iploc))                
                    f.seek(loc)
                    tmpdata = np.fromfile(f,dtype=complextype,count=forcent)
                    if bar3ptChromaByteSwap:tmpdata = tmpdata.byteswap()
                    if 'cmplx' in thisgamma:
                        datahold[igamma].append(tmpdata.imag)
                    else:                    
                        datahold[igamma].append(tmpdata.real)
                    # if isshift> nt-min(AllTSinkList)-1:
                    datahold[igamma][-1][-corr_shift-1] = -datahold[igamma][-1][-corr_shift-1]
                    # print
                    # print thisfile.replace(xsrcList[0],xsrc)
                    # print tmpdata.real
                    # print corr_shift
                    # print datahold[igamma][-1]
                    # print
                    # datahold[igamma][-1] = np.abs(datahold[igamma][-1])
                    if any(np.isnan(datahold[igamma][ip])) and DeleteNanCfgs:
                        raise NaNCfunError('NaN Values: '+thisgamma+' ' +qvecSet[iploc]  )
            f.close()
            if len(self.data) == 0:                    
                self.data = np.array(datahold)
            else:
                self.data += np.array(datahold)
        self.data = self.data/len(thisxsrcList)

## Just a reference function for the extra parameters in the 3 point function file.
def holder(thisfile):
    f = open(thisfile,'rb')    
    d1 = np.fromfile(f,dtype=np.int32,count=10).byteswap()
    output_version = d1[0] 
    mom2max = d1[1]
    j_decay = d1[2] 
    nrow = d1[4:8]
    wilson_version = d1[8]
    Num_seqsrc = d1[9]
    
    nuc_type = np.fromfile(f,dtype='S13',count=1)[0]
    d2 = np.fromfile(f,dtype=np.int32,count=17).byteswap()
    tsource = d2[0] 
    tsink = d2[1] 
    sinkmom = d2[3:6]
    gammanumber = d2[6]
    FFoutput_version = d2[7]
    numFF = d2[8]
    gammaopp = d2[9]
    numMom = d2[10]
    MagicNumber = d2[11]
    qmom = d2[13:16]
    tcurrlen = d2[16]
    values = np.fromfile(f,dtype=complextype,count=tsink-tsource)
    
    f.close()
    
        

## Der[igamma , ip , it]
class ReadFSDerCfunPick:
    def __init__(self,thisfile,thisMomList,thisDGList):
        self.data = []
        f = open(thisfile,'rb')
        for igd,thisgamma in enumerate(thisDGList):
            cmplxflag = False
            if 'cmplx' in thisgamma: 
                cmplxflag = True
                thisgamma = thisgamma.replace('cmplx','')
            igdloc = DGSet.index(thisgamma)
            self.data.append([])
            for ip,iploc in enumerate(thisMomList):         
                self.data[igd].append([])
                for it in range(nt):
                    loc = 2*SeekIncSize*(igdloc + iploc*DGBSize + it*DGBSize*nmom)
                    if cmplxflag: loc += SeekIncSize
                    f.seek(loc)
                    tmpdata = array('d')
                    tmpdata.read(f,1)
                    tmpdata.byteswap()
                    self.data[igd][ip].append(tmpdata[0])
                    if np.isnan(self.data[igd][ip][it]) and DeleteNanCfgs:
                        raise NaNCfunError('NaN Values: '+thisgamma+' ' +qvecSet[iploc] + ' it='+str(it+1) )
        f.close()

## noDer[igamma , ip , it]
class ReadFSCfunPick:
    def __init__(self,thisfile,thisMomList,thisGammaList):
        self.data = []
        f = open(thisfile,'rb')
        for igamma,thisgamma in enumerate(thisGammaList):
            cmplxflag = False
            if 'cmplx' in thisgamma: 
                cmplxflag = True
                thisgamma = thisgamma.replace('cmplx','')
            igammaloc = GammaSet.index(thisgamma)
            self.data.append([])
            for ip,iploc in enumerate(thisMomList):      
                self.data[igamma].append([])
                for it in range(nt):
                    loc = 2*SeekIncSize*(igammaloc + iploc*GBSize + it*GBSize*nmom)
                    if cmplxflag: loc += SeekIncSize
                    f.seek(loc)
                    tmpdata = array('d')
                    tmpdata.read(f,1)
                    tmpdata.byteswap()
                    self.data[igamma][ip].append(tmpdata[0])
                    if np.isnan(self.data[igamma][ip][it]) and DeleteNanCfgs:
                        if Debug: self.data[igamma][ip][it] = 0.0
                        raise NaNCfunError('NaN Values: '+thisgamma+ ' ' +qvecSet[iploc] + ' it='+str(it+1) )
        f.close()



class Read2ptCfunPick:
    def __init__(self,thisfile,thisMomList):
        self.data = []
        f = open(thisfile,'rb')
        if XAvg: thisxsrcList = xsrcList
        else: thisxsrcList = [xsrcList[0]]
        self.tshiftlist = [0]*len(thsixsrclist)
        for ip,iploc in enumerate(thisMomList):
            self.data.append([])
            for it in range(nt):
        ##11 + 22 projected spin component ns = 4 works
                loc = SeekIncSize*2*(ns**2)*(it + iploc*nt)
                f.seek(loc)
                tmpdata = array('d')
                tmpdata.read(f,1)
                tmpdata.byteswap()
                val = tmpdata[0]
                f.seek(loc+10*SeekIncSize)
                tmpdata = array('d')
                tmpdata.read(f,1)
                tmpdata.byteswap()
                val += tmpdata[0]
                self.data[ip].append(val)
                if np.isnan(self.data[ip][it]) and DeleteNanCfgs:
                    raise NaNCfunError('NaN Values: ' +qvecSet[iploc] + ' it='+str(it+1) )
                # if 'cmplx' in RI:
                #     self.datai[ipread].append(tmpdata[loc+1]+tmpdata[loc+11])
            if self.data[ip][tsource-1] < 0.0:
                self.data[ip] = np.negative(self.data[ip])
        f.close()


class Read2ptCfunChroma:
    def __init__(self,thisfile,thisMomList):
        ##CURRENT: LIME_UNPACK THE FILE AND READ THE BINARY INSIDE##
        self.data = []
        barnum = 0
        # barnum = 21
        if XAvg: thisxsrcList = xsrcList
        else: thisxsrcList = [xsrcList[0]]
        self.tshiftlist = [0]*len(thsixsrclist)
        for ip,iploc in enumerate(thisMomList):
            self.data.append(np.memmap(thisfile,dtype=np.complex128,mode='r',offset=nt*ChromaSIS*ip,shape=(nt,)).byteswap())

class Read2ptCfunChromaXML:
    def __init__(self,thisfile,thisMomList):
        self.data = []
        if XAvg: thisxsrcList = xsrcList
        else: thisxsrcList = [xsrcList[0]]
        self.tshiftlist = []
        for xsrc in thisxsrcList:
            self.OutMomList = []
            TSRC_read = False
            datahold = []
            # print 'Reading ' ,thisfile.replace(xsrcList[0],xsrc)
            with open(thisfile.replace(xsrcList[0],xsrc),'r') as f:
                BarPart,InterpPart,ReadMom = False,False,False                
                for line in f:
                    strline = line.strip()
                    if 't_srce' in strline and not TSRC_read:
                        TSRC_read = True
                        thissrclist = strline.replace('<t_srce>','').replace('</t_srce>','').split()
                        # print thisfile.replace(xsrcList[0],xsrc)
                        # print int(thissrclist[-1])
                        self.tshiftlist.append(int(thissrclist[-1]))
                    elif strline == '<Shell_Shell_Wilson_'+MesOrBar+'s>':
                    # if strline == '<Shell_Shell_Wilson_Baryons>':
                        BarPart = True
                    elif InterpFlag in strline:
                        if strline.replace('<'+InterpFlag+'>','').replace('</'+InterpFlag+'>','') == InterpNumb:
                            InterpPart = True
                    elif BarPart and InterpPart:
                        if '<sink_mom_num>' in strline:
                            thismom = int(strline.replace('<sink_mom_num>','').replace('</sink_mom_num>',''))
                            if thismom in thisMomList:
                                datahold.append([])
                                self.OutMomList.append(thismom)
                                ReadMom = True
                            else:
                                ReadMom = False
                        elif '<re>' in strline and ReadMom:
                            datahold[-1].append(float(strline.replace('<re>','').replace('</re>','')))
                            if np.isnan(datahold[-1][-1]) and DeleteNanCfgs:
                                raise NaNCfunError('NaN Values: '+thisfile+'  ' +qvecSet[int(self.OutMomList[-1])]  )
                    if strline == '</momenta>' and InterpPart:
                        if len(datahold) > 0: break
            # print xsrc, ' data '
            # print np.array(datahold)
            # print self.data
            # print 
            if len(self.data) == 0:                    
                self.data = np.array(datahold)
            else:
                self.data += np.array(datahold)
        self.data = self.data/len(thisxsrcList)
        indicies =  np.searchsorted(self.OutMomList,thisMomList)
        # if Debug:
        #     print 
        #     print thisMomList
        #     print thisfile
        #     print indicies
        #     print self.data
        self.data = np.array(self.data)[indicies].tolist()
                    
                    
                
        
class NaNCfunError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)












# ## Der[igamma , ip , it]
# class ReadFSDerCfunPick:
#     def __init__(self,file,thisMomList,thisDGList):
#         self.data = []
#         f = open(file,'rb')
#         tmpdata = array('d')
#         tmpdata.read(f,DGBSize*nmom*nt*2)
#         tmpdata.byteswap()
#         for igd,thisgamma in enumerate(thisDGList):
#             cmplxflag = False
#             if 'cmplx' in thisgamma: 
#                 cmplxflag = True
#                 thisgamma = thisgamma.replace('cmplx','')
#             igdloc = DGSet.index(thisgamma)
#             self.data.append([])
#             for ip,iploc in enumerate(thisMomList):         
#                 self.data[igd].append([])
#                 for it in range(nt):
#                     loc = igdloc*2 + iploc*DGBSize*2 + it*DGBSize*nmom*2
#                     if cmplxflag: loc += 1
#                     self.data[igd][ip].append(tmpdata[loc])
#                     if np.isnan(self.data[igd][ip][it]):
#                         raise NaNCfunError('NaN Values: '+thisgamma+' ' +qvecSet[iploc] + ' it='+str(it+1) )

#         f.close()

# ## noDer[igamma , ip , it]
# class ReadFSCfunPick:
#     def __init__(self,file,thisMomList,thisGammaList):
#         self.data = []
#         f = open(file,'rb')
#         tmpdata = array('d')
#         try:
#             tmpdata.read(f,GBSize*nmom*nt*2)
#         except:
#             print file
#         tmpdata.byteswap()
#         for igamma,thisgamma in enumerate(thisGammaList):
#             cmplxflag = False
#             if 'cmplx' in thisgamma: 
#                 cmplxflag = True
#                 thisgamma = thisgamma.replace('cmplx','')
#             igammaloc = GammaSet.index(thisgamma)
#             self.data.append([])
#             for ip,iploc in enumerate(thisMomList):      
#                 self.data[igamma].append([])
#                 for it in range(nt):
#                     loc = igammaloc*2 + iploc*GBSize*2 + it*GBSize*nmom*2
#                     if cmplxflag: loc += 1
#                     self.data[igamma][ip].append(tmpdata[loc])
#                     if np.isnan(self.data[igamma][ip][it]):
#                         raise NaNCfunError('NaN Values: '+thisgamma+ ' ' +qvecSet[iploc] + ' it='+str(it+1) )

#         f.close()

# ##[ip,it]
# class Read2ptCfunPickOld:
#     def __init__(self,file,thisMomList):
#         self.data = []
#         f = open(file,'rb')
#         tmpdata = array('d')
#         tmpdata.read(f,nmom*ns*ns*nt*2)
#         tmpdata.byteswap()
#         for ip,iploc in enumerate(thisMomList):
#             self.data.append([])
#             for it in range(nt):
#         ##11 + 22 projected spin component ns = 4 works
#                 loc = it*2*ns**2 + iploc*nt*2*ns**2
#                 self.data[ip].append(tmpdata[loc]+tmpdata[loc+10])
#                 if np.isnan(self.data[ip][it]):
#                     raise NaNCfunError('NaN Values: ' +qvecSet[iploc] + ' it='+str(it+1) )
#                 # if 'cmplx' in RI:
#                 #     self.datai[ipread].append(tmpdata[loc+1]+tmpdata[loc+11])
#             if self.data[ip][tsource-1] < 0.0:
#                 self.data[ip] = np.negative(self.data[ip])
#         f.close()

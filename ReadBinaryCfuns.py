#!/usr/bin/env python
from array import array
import numpy as np
from Params import *
from FFSympy import GammaTOChroma

SeekIncSize = 8
ChromaSIS = 16
bar3ptChromaByteSwap = True
# complextype = np.complex64
complextype = np.complex128
CheckMagic = False

def GetBar3ptLoc(ig,ip,tsinklen,momlistlen,thiscomplextype,PolProj=False):
    intlen = np.dtype(np.int32).itemsize
    if PolProj:
        strlen = np.dtype('S11').itemsize
    else:
        strlen = np.dtype('S13').itemsize
    cmplxlen = np.dtype(thiscomplextype).itemsize
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
    def __init__(self,thisfileList,thisMomList,thisGammaList,forcent=nt):
        self.data = []
        for thisfile in thisfileList:
            # if thistsink == 'FileName':
            #     thistsink = int(re.findall('tsink.*?p',thisfile)[0].replace('tsink','').replace('p',''))
            datahold = []            
            if os.path.getsize(thisfile) == 826009: thiscomplextype = np.complex128
            elif os.path.getsize(thisfile) < 600000:
                # raise IOError(thisfile+ ' Is old 64 bit, please remove all smaller sized files')
                thiscomplextype = np.complex64
            else: thiscomplextype = complextype
            f = open(thisfile,'rb')
            for igamma,thisgamma in enumerate(thisGammaList):
                igammaloc = GammaTOChroma(thisgamma)
                datahold.append([])
                for ip,iploc in enumerate(thisMomList):      
                    loc,loccons,magmin = GetBar3ptLoc(igammaloc,iploc,forcent,len(qvecSet),thiscomplextype,PolProj=('GMA3' in thisfile))
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
                    tmpdata = np.fromfile(f,dtype=thiscomplextype,count=forcent)
                    if bar3ptChromaByteSwap:tmpdata = tmpdata.byteswap()
                    if 'cmplx' in thisgamma:
                        datahold[igamma].append(tmpdata.imag)
                    else:                    
                        datahold[igamma].append(tmpdata.real)
                    # if isshift> nt-min(AllTSinkList)-1:
                    # datahold[igamma][-1][-corr_shift-1] = -datahold[igamma][-1][-corr_shift-1]
                    # print
                    # print tmpdata.real
                    # print
                    # totsign = np.sign(np.sum(np.sign(datahold[igamma][-1][:CMTSinkList[0]])))
                    # datahold[igamma][-1] = totsign*np.abs(datahold[igamma][-1])
                    datahold[igamma][-1] = np.sign(datahold[igamma][-1][0])*np.abs(datahold[igamma][-1])
                    # print np.abs(datahold[igamma][-1][13]-datahold[igamma][-1][12])/datahold[igamma][-1][12]
                    # if np.abs(datahold[igamma][-1][13]-datahold[igamma][-1][12])/datahold[igamma][-1][12] > 1.0:
                    #     print
                    #     print thisfile.replace(xsrcList[0],xsrc)
                    #     print corr_shift
                    #     for it,tdata in enumerate(datahold[igamma][-1]):
                    #         print it, tdata
                    # print np.roll(datahold[igamma][-1],isshift)
                    # print np.roll(datahold[igamma][-1],-isshift)
                    if any(np.isnan(datahold[igamma][ip])):
                        print 
                        print thisfile
                        print 'thisgamma ',thisgamma,' iploc ', iploc
                        print 'data:'
                        for it,idata in enumerate(datahold[igamma][ip]):
                            print it,idata
                        if DeleteNanCfgs:
                            raise NaNCfunError('NaN Values: '+thisgamma+' ' +qvecSet[iploc]  )
            f.close()
            if len(self.data) == 0:                    
                self.data = np.array(datahold)
                datalen = 1
            else:
                self.data += np.array(datahold)
                datalen += 1
        self.data = self.data/float(datalen)

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
    
    nuc_type = np.fromfile(f,dtype='S11',count=1)[0]
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

    print 'output_version',output_version
    print 'mom2max',mom2max 
    print 'j_decay',j_decay
    print 'nrow',nrow
    print 'wilson_version',wilson_version
    print 'Num_seqsrc',Num_seqsrc
    
    print 'nuc_type',nuc_type
    print 'tsource',tsource
    print 'tsink ',tsink 
    print 'sinkmom ',sinkmom 
    print 'gammanumber ',gammanumber 
    print 'FFoutput_version ',FFoutput_version 
    print 'numFF ',numFF 
    print 'gammaopp ',gammaopp 
    print 'numMom ',numMom 
    print 'MagicNumber',MagicNumber
    print 'qmom ',qmom 
    print 'tcurrlen ',tcurrlen 
    print 'values',values

    
        

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
        self.tshiftlist = [0]*len(thisxsrcList)
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
        self.tshiftlist = [0]*len(thisxsrcList)
        for ip,iploc in enumerate(thisMomList):
            self.data.append(np.memmap(thisfile,dtype=np.complex128,mode='r',offset=nt*ChromaSIS*ip,shape=(nt,)).byteswap())

class Read2ptCfunChromaXML:
    def __init__(self,thisfile,thisMomList,Dog5=False):
        self.data = []
        self.datag5 = []
        if XAvg: thisxsrcList = xsrcList
        else: thisxsrcList = [xsrcList[0]]
        self.tshiftlist = []
        datalen,datag5len = 0,0
        for xsrc in thisxsrcList:
            self.OutMomList = []
            TSRC_read = False
            datahold = []
            datag5hold = []
            # print 'Reading ' ,thisfile.replace(xsrcList[0],xsrc)
            with open(thisfile.replace(xsrcList[0],xsrc),'r') as f:
                BarPart,InterpPart,InterpPartg5,ReadMom = False,False,False,False
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
                        if int(strline.replace('<'+InterpFlag+'>','').replace('</'+InterpFlag+'>','')) == int(InterpNumb):
                            InterpPartg5 = False
                            InterpPart = True
                        elif int(strline.replace('<'+InterpFlag+'>','').replace('</'+InterpFlag+'>','')) == int(INg5):
                            InterpPart = False
                            InterpPartg5 = True
                        else:
                            InterpPart = False
                            InterpPartg5 = False
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
                    elif BarPart and InterpPartg5:
                        if '<sink_mom_num>' in strline:
                            thismom = int(strline.replace('<sink_mom_num>','').replace('</sink_mom_num>',''))
                            if thismom in thisMomList:
                                datag5hold.append([])
                                ReadMom = True
                            else:
                                ReadMom = False
                        elif '<re>' in strline and ReadMom:
                            datag5hold[-1].append(float(strline.replace('<re>','').replace('</re>','')))
                            if np.isnan(datag5hold[-1][-1]) and DeleteNanCfgs:
                                raise NaNCfunError('NaN Values: '+thisfile+'  ' +qvecSet[int(self.OutMomList[-1])]  )
                    if strline == '</momenta>' and InterpPart:
                        InterpPart = False
                        if len(datahold) > 0 and ((len(datag5hold) > 0) or (not Dog5 )): break
            # print xsrc, ' data '
            # print np.array(datahold)
            # print self.data
            # print 
            if len(self.data) == 0:                    
                self.data = np.rollaxis(np.array(datahold),0,1)
                datalen = 1
            else:
                self.data += np.rollaxis(np.array(datahold),0,1)
                datalen += 1
            if Dog5:
                if len(self.datag5) == 0:                    
                    datag5len = 1
                    self.datag5 = np.rollaxis(np.array(datag5hold),0,1)
                else:
                    self.datag5 += np.rollaxis(np.array(datag5hold),0,1)
                    datag5len += 1
        self.data = self.data/float(datalen)
        if Dog5: self.datag5 = self.datag5/float(datag5len)
        indicies =  np.searchsorted(self.OutMomList,thisMomList)
        # if Debug:
        #     print 
        #     print thisMomList
        #     print thisfile
        #     print indicies
        #     print self.data
        self.data = np.array(self.data)[indicies].tolist()
        if Dog5: self.datag5 = np.array(self.datag5)[indicies].tolist()


class R2CChromaXMLFileList:
    def __init__(self,thisfileList,thisMomList,Dog5=False):
        self.data = []
        self.datag5 = []
        self.tshiftlist = []
        datalen,datag5len = 0,0
        for thisfile in thisfileList:
            self.OutMomList = []
            TSRC_read = False
            datahold = []
            datag5hold = []
            # print 'Reading ' ,thisfile.replace(xsrcList[0],xsrc)
            with open(thisfile,'r') as f:
                BarPart,InterpPart,InterpPartg5,ReadMom = False,False,False,False
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
                        if int(strline.replace('<'+InterpFlag+'>','').replace('</'+InterpFlag+'>','')) == int(InterpNumb):
                            InterpPartg5 = False
                            InterpPart = True
                        elif int(strline.replace('<'+InterpFlag+'>','').replace('</'+InterpFlag+'>','')) == int(INg5):
                            InterpPart = False
                            InterpPartg5 = True
                        else:
                            InterpPart = False
                            InterpPartg5 = False
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
                    elif BarPart and InterpPartg5:
                        if '<sink_mom_num>' in strline:
                            thismom = int(strline.replace('<sink_mom_num>','').replace('</sink_mom_num>',''))
                            if thismom in thisMomList:
                                datag5hold.append([])
                                ReadMom = True
                            else:
                                ReadMom = False
                        elif '<re>' in strline and ReadMom:
                            datag5hold[-1].append(float(strline.replace('<re>','').replace('</re>','')))
                            if np.isnan(datag5hold[-1][-1]) and DeleteNanCfgs:
                                raise NaNCfunError('NaN Values: '+thisfile+'  ' +qvecSet[int(self.OutMomList[-1])]  )
                    if strline == '</momenta>' and InterpPart:
                        if len(datahold) > 0 and ((len(datag5hold) > 0) or (not Dog5 )): break
            # print xsrc, ' data '
            # print np.array(datahold)
            # print self.data
            # print 
            if len(self.data) == 0:                    
                self.data = np.rollaxis(np.array(datahold),0,1)
                datalen = 1
            else:
                self.data += np.rollaxis(np.array(datahold),0,1)
                datalen += 1
            if Dog5:
                if len(self.datag5) == 0:                    
                    datag5len = 1
                    self.datag5 = np.rollaxis(np.array(datag5hold),0,1)
                else:
                    self.datag5 += np.rollaxis(np.array(datag5hold),0,1)
                    datag5len += 1
        if datalen != datag5len and Dog5: print 'Warning!, datag5len not equal to datalen'
        self.data = self.data/float(datalen)
        if Dog5: self.datag5 = self.datag5/float(datag5len)
        indicies =  np.searchsorted(self.OutMomList,thisMomList)
        # if Debug:
        #     print 
        #     print thisfile
        #     print thisMomList
        #     print self.OutMomList
        #     print indicies
        #     print self.data
        #     print self.datag5
        self.data = np.array(self.data)[indicies].tolist()
        if Dog5: self.datag5 = np.array(self.datag5)[indicies].tolist()

        
                    
                
        
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

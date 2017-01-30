#!/usr/bin/env python

import numpy as np
from Params import *
from MomParams import *
from MiscFuns import *
from FFFuns import *
from FitFunctions import *
from collections import OrderedDict
from BootTest import BootStrap1
from FFParams import *
from LLSBoot import *
import time

## data { { gamma } { mom } { Fit(Boot/Avg/Std/Chi) } }
## dataout { { momsqrd } { Boot/Avg/Chi }
##REMEBER deal with cmplx signals
def CreateFF(data,mass,iCurr,gammaflag='',Rfac=True,alphalist = [1.0]):
    thisdataout = {}
    DoTop = 'Top' in iCurr
    baseCurr = iCurr.replace(gammaflag,'')
    Opps = CurrOpps[baseCurr]
    thisdataout = OrderedDict()
    infodict = {}
    infodict['alpha'] = OrderedDict()
    infodict['alpha']['Avg'] = alphalist[0]
    infodict['alpha']['Std'] = np.std(alphalist[1:])                        
    for iqsqrd in MomSqrdSet:        
        momstart = time.time()
        iqs = 'qsqrd'+str(iqsqrd)
        thisdataout['qsqrd'+str(iqsqrd)] = {}
        datavals,FFcoeff = [],[]
        for ica in range(len(alphalist)):
            FFcoeff.append([])
            for iff in range(NoFFPars[baseCurr]):
                FFcoeff[ica].append([])
        opplist = []
        for iopp in Opps:                    
            flagopp = gammaflag+iopp
            RealVal,CmplxVal = False,False
            if flagopp in data.keys(): RealVal = True
            if flagopp+'cmplx' in data.keys(): CmplxVal = True
            if not RealVal and not CmplxVal: continue
            start = time.time()
            for iq in qvecINqsqrd(int(iqsqrd)):
                if RealVal: 
                    if iq not in data[flagopp].keys(): continue 
                if CmplxVal: 
                    if iq not in data[flagopp+'cmplx'].keys(): continue 
                for ica,ialpha in enumerate(alphalist):
                    FFcoeffhold,rcheck,ccheck = CurrFFs[baseCurr](iopp,np.array(qstrTOqvec(iq))*qunit,[0,0,0],mass,Rfac=Rfac,alpha=ialpha)
                    if CmplxVal and ccheck:
                        for iFF,iFFcof in enumerate(FFcoeffhold):
                            FFcoeff[ica][iFF].append(iFFcof.imag)
                    if RealVal and rcheck:
                        for iFF,iFFcof in enumerate(FFcoeffhold):
                            FFcoeff[ica][iFF].append(iFFcof.real)
                    if ica == 0:
                        if CmplxVal and ccheck:
                            datavals.append(data[flagopp+'cmplx'][iq]['Boot'])
                            infodict[iqs] = data[flagopp+'cmplx'][iq]['Info']
                            opplist.append(flagopp+'cmplx '+ iq)
                        if RealVal and rcheck:
                            datavals.append(data[flagopp][iq]['Boot'])
                            infodict[iqs] = data[flagopp][iq]['Info']
                            opplist.append(flagopp +' '+ iq)
            # print 'PullOutLHS iopp',iopp,', time taken:' , GetTimeStr(time.time()-start)
        if len(datavals) == 0: continue
        zboot,zvec = [BootStrap1(nboot,0)],[0.0]        
        if not DoTop:
            FFcoeff = FFcoeff[0]
        if 'Scalar' in baseCurr :
            if Debug:
                print 'Printing Form Factors debug:'
                for FF1,res in zip(FFcoeff[0],datavals):
                    print iqsqrd, '   ' , FF1,'FF1 = ',res.Avg
                print ''
            FFBoothold,FFAvghold,FFChihold = FitBoots(datavals,FFcoeff,FFFitFuns[baseCurr],tBooted=DoTop)
            thisdataout[iqs]['Boot'] = FFBoothold
            thisdataout[iqs]['Avg'] = FFAvghold
        elif DoTop:
            FFcoeff = np.array(FFcoeff)
            if Debug:
                print 'Printing Form Factors debug:'
                for FF1,FF2,FF3,res,iopp in zip(FFcoeff[0,0],FFcoeff[0,1],FFcoeff[0,2],datavals,opplist):
                    print iqsqrd, '   ' , iopp, ' ',FF1,'FF1 + ',FF2,'FF2 + ',FF3,'FF3 = ',res.Avg, res.Std
                print ''
            if len(datavals) == 1:
                if sum(ia == [0.0] for ia in FFcoeff[0]) != 2: continue
                if FFcoeff[0,0] != [0.0]:
                    FFBoothold,FFAvghold,FFChihold = FitBoots(datavals,FFcoeff[:,[0],:],FFFitFuns['Scalar'],tBooted=DoTop)
                    thisdataout[iqs]['Boot'] = FFBoothold+zboot+zboot
                    thisdataout[iqs]['Avg'] = FFAvghold+zvec+zvec
                elif FFcoeff[0,1] != [0.0]:
                    FFBoothold,FFAvghold,FFChihold = FitBoots(datavals,FFcoeff[:,[1],:],FFFitFuns['Scalar'],tBooted=DoTop)
                    thisdataout[iqs]['Boot'] = zboot+FFBoothold+zboot
                    thisdataout[iqs]['Avg'] = zvec+FFAvghold+zvec
                elif FFcoeff[0,2] != [0.0]:
                    FFBoothold,FFAvghold,FFChihold = FitBoots(datavals,FFcoeff[:,[2],:],FFFitFuns['Scalar'],tBooted=DoTop)
                    thisdataout[iqs]['Boot'] = zboot+zboot+FFBoothold
                    thisdataout[iqs]['Avg'] = zvec+zvec+FFAvghold
            elif len(datavals) == 2:
                if [0.0] not in FFcoeff[0]: continue
                if FFcoeff[0,0] == [0.0]:
                    FFBoothold,FFAvghold,FFChihold = FitBoots(datavals,FFcoeff[:,[1,2],:],FFFitFuns['Vector'],tBooted=DoTop)
                    thisdataout[iqs]['Boot'] = zboot+FFBoothold
                    thisdataout[iqs]['Avg'] = zvec+FFAvghold
                elif FFcoeff[0,1] == [0.0]:
                    FFBoothold,FFAvghold,FFChihold = FitBoots(datavals,FFcoeff[:,[0,2],:],FFFitFuns['Vector'],tBooted=DoTop)
                    thisdataout[iqs]['Boot'] = [FFAvghold[0]]+zvec+[FFAvghold[1]]
                    thisdataout[iqs]['Avg'] = [FFBoothold[0]]+zboot+[FFBoothold[1]]
                elif FFcoeff[0,2] == [0.0]:
                    FFBoothold,FFAvghold,FFChihold = FitBoots(datavals,FFcoeff[:,[0,1],:],FFFitFuns['Vector'],tBooted=DoTop)
                    thisdataout[iqs]['Boot'] = FFBoothold+zboot
                    thisdataout[iqs]['Avg'] = FFAvghold+zvec
            else:
                FFBoothold,FFAvghold,FFChihold = FitBoots(datavals,FFcoeff,FFFitFuns[baseCurr],tBooted=DoTop)
                thisdataout[iqs]['Boot'] = FFBoothold
                thisdataout[iqs]['Avg'] = FFAvghold
        elif 'Tensor' in baseCurr:
            if len(datavals) == 1:
                if sum(ia == [0.0] for ia in FFcoeff) != 2: continue
                if FFcoeff[0] != [0.0]:
                    FFBoothold,FFAvghold,FFChihold = FitBoots(datavals,FFcoeff[0],FFFitFuns['Scalar'],tBooted=DoTop)
                    thisdataout[iqs]['Boot'] = FFBoothold+zboot+zboot
                    thisdataout[iqs]['Avg'] = FFAvghold+zvec+zvec
                elif FFcoeff[1] != [0.0]:
                    FFBoothold,FFAvghold,FFChihold = FitBoots(datavals,FFcoeff[1],FFFitFuns['Scalar'],tBooted=DoTop)
                    thisdataout[iqs]['Boot'] = zboot+FFBoothold+zboot
                    thisdataout[iqs]['Avg'] = zvec+FFAvghold+zvec
                elif FFcoeff[2] != [0.0]:
                    FFBoothold,FFAvghold,FFChihold = FitBoots(datavals,FFcoeff[2],FFFitFuns['Scalar'],tBooted=DoTop)
                    thisdataout[iqs]['Boot'] = zboot+zboot+FFBoothold
                    thisdataout[iqs]['Avg'] = zvec+zvec+FFAvghold
            elif len(datavals) == 2:
                if [0.0] not in FFcoeff: continue
                if FFcoeff[0] == [0.0]:
                    FFBoothold,FFAvghold,FFChihold = FitBoots(datavals,[FFcoeff[1],FFcoeff[2]],FFFitFuns['Vector'],tBooted=DoTop)
                    thisdataout[iqs]['Boot'] = zboot+FFBoothold
                    thisdataout[iqs]['Avg'] = zvec+FFAvghold
                elif FFcoeff[1] == [0.0]:
                    FFBoothold,FFAvghold,FFChihold = FitBoots(datavals,[FFcoeff[0],FFcoeff[2]],FFFitFuns['Vector'],tBooted=DoTop)
                    thisdataout[iqs]['Boot'] = [FFAvghold[0]]+zvec+[FFAvghold[1]]
                    thisdataout[iqs]['Avg'] = [FFBoothold[0]]+zboot+[FFBoothold[1]]
                elif FFcoeff[2] == [0.0]:
                    FFBoothold,FFAvghold,FFChihold = FitBoots(datavals,[FFcoeff[0],FFcoeff[1]],FFFitFuns['Vector'],tBooted=DoTop)
                    thisdataout[iqs]['Boot'] = FFBoothold+zboot
                    thisdataout[iqs]['Avg'] = FFAvghold+zvec
            else:
                FFBoothold,FFAvghold,FFChihold = FitBoots(datavals,FFcoeff,FFFitFuns[baseCurr],tBooted=DoTop)
                thisdataout[iqs]['Boot'] = FFBoothold
                thisdataout[iqs]['Avg'] = FFAvghold
        elif 'Vector' in baseCurr:
            ## DEBUG ##
            if Debug:
                print 'Printing Form Factors debug:'
                for FF1,FF2,res in zip(FFcoeff[0],FFcoeff[1],datavals):
                    print iqsqrd, '   ' , FF1,'FF1 + ',FF2,'FF2 = ',res.Avg
                print ''
            if [0.0] not in FFcoeff and len(datavals) == 1: continue
            if FFcoeff[0] == [0.0]:
                FFBoothold,FFAvghold,FFChihold = FitBoots(datavals,FFcoeff[1],FFFitFuns['Scalar'],tBooted=DoTop)
                thisdataout[iqs]['Boot'] = zboot+FFBoothold
                thisdataout[iqs]['Avg'] = zvec+FFAvghold
            elif FFcoeff[1]== [0.0]:
                FFBoothold,FFAvghold,FFChihold = FitBoots(datavals,FFcoeff[0],FFFitFuns['Scalar'],tBooted=DoTop)
                thisdataout[iqs]['Boot'] = FFBoothold+zboot
                thisdataout[iqs]['Avg'] = FFAvghold+zvec
            else:
                FFBoothold,FFAvghold,FFChihold = FitBoots(datavals,FFcoeff,FFFitFuns[baseCurr],tBooted=DoTop)
                thisdataout[iqs]['Boot'] = FFBoothold
                thisdataout[iqs]['Avg'] = FFAvghold
        thisdataout[iqs]['Chi'] = FFChihold[0]
        # print 'All Done',iqs,', time taken:' , GetTimeStr(time.time()-momstart)
    return thisdataout,infodict


    
            # print 'qsqrd = ',iqsqrd            
            # for iFF,ival in enumerate(datavals):
            #     if len(FFcoeff) == 1:
            #         print FFcoeff[0][iFF] , ival.Avg 
            #     if len(FFcoeff) == 2:
            #         print FFcoeff[0][iFF] ,FFcoeff[1][iFF] , ival.Avg 
            #     if len(FFcoeff) == 3:
            #         print FFcoeff[0][iFF] ,FFcoeff[1][iFF] ,FFcoeff[2][iFF] , ival.Avg 
            #         if 'Tensor' in CT and iqsqrd == '0':
            #             print iopp,FFcoeffhold, data[iopp][iq]['FitBoot'].Avg
            # if 'Tensor' in CT and iqsqrd == '0':
            #     print datavals , NoFFPars[CT]



                # print 'Fit:'
                # for iboot in thisdataout[CT]['qsqrd'+str(iqsqrd)]['Boot']:
                #     print iboot.Avg,iboot.Std,thisdataout[CT]['qsqrd'+str(iqsqrd)]['Chi']
                # print ''

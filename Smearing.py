#!/usr/bin/env python

from Params import *
import numpy as np
from BootTest import BootStrap1
import operator as op
from collections import OrderedDict
import time,datetime


def SweepsToRMS(alpha,Nsmears,LatDim=ndim[:-1]):
   LatDim = np.array(LatDim)
   sDim = LatDim/2
   Smat = np.zeros(LatDim,dtype=np.float64)
   Smat[sDim[0],sDim[1],sDim[2]] = np.float64(1.0)
   cofone = np.float64(1.0-alpha)
   coftwo = np.float64(alpha/6.0)
   for ism in range(1,Nsmears+1):
      # print Smat[sDim[0],sDim[1],sDim[2]]
      # print Smat[sDim[0]+1,sDim[1],sDim[2]]
      # print Smat[sDim[0],sDim[1]+1,sDim[2]]
      # print Smat[sDim[0],sDim[1],sDim[2]+1]
      # print Smat[sDim[0]-1,sDim[1],sDim[2]]
      # print Smat[sDim[0],sDim[1]-1,sDim[2]]
      # print Smat[sDim[0],sDim[1],sDim[2]-1]
      # print ''
      Smat = (cofone * Smat + 
              coftwo * (np.roll(Smat,1,0) +
                        np.roll(Smat,1,1) +
                        np.roll(Smat,1,2) +
                        np.roll(Smat,-1,0) +
                        np.roll(Smat,-1,1) +
                        np.roll(Smat,-1,2)))
   SqrRad = np.float64(0.0)
   for ix in range(1,LatDim[0]):
      for iy in range(1,LatDim[1]):
         for iz in range(1,LatDim[2]):
            SqrRad += Smat[ix,iy,iz]**2 * ((ix-sDim[0])**2+(iy-sDim[1])**2+(iz-sDim[2])**2)
   RMS = np.sqrt(SqrRad / np.sum(Smat**2))

   print 'RMS lattice units: ' , RMS
   print 'RMS physical units: ' , '{0:5.3f} fm'.format(RMS*latspace)

   print 'S summed: ' , np.sum(Smat)
   


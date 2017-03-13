#!/usr/bin/env python

from OvDetEigen import OverdetEigen

from scipy.linalg import eig
import numpy as np
import numpy.random as npr

Scale = 0.1
A = np.arange(16).reshape(4,4)
# B = np.eye(*A.shape)
B = np.roll(np.arange(16),5).reshape(4,4)[::-1]

randlist1,randlist2 = npr.rand(*A.shape)*Scale,npr.rand(*A.shape)*Scale
randlist3,randlist4 = npr.rand(*B.shape)*Scale,npr.rand(*B.shape)*Scale
A = np.bmat([[A+randlist1],[A+randlist2]])
# A = np.bmat([[A+randlist1],[A+randlist1]])
B = np.bmat([[B+randlist3],[B+randlist4]])

Asq = A[:A.shape[1],:]
Bsq = B[:B.shape[1],:]
Iter = 1000



print 'input:'
print 'A=',A
print 'B=',B
print
print 'Asq=',Asq
print 'Bsq=',Bsq
print
print 'Niter=',Iter

Evals,REvec,Dump = OverdetEigen(A,B,Iter)
Evals3,REvec3,Dump = OverdetEigen(Asq,Bsq,Iter)
Evals2,REvec2 = eig(Asq,Bsq)
print
print
print 'Overdet'
print 'Evals'
print Evals
print 'Evec'
print REvec
print
print 'Overdet Square'
print 'Evals'
print Evals3
print 'Evec'
print REvec3
print
print
print 'reg'
print 'Evals'
print Evals2
print 'Evec'
print REvec2
                                                                                                

Asq = np.matrix(Asq)
Bsq = np.matrix(Bsq)
REvec = np.matrix(REvec)
REvec2 = np.matrix(REvec2)
REvec3 = np.matrix(REvec3)

# print
# print 'shapes'
# print Asq.shape
# print Bsq.shape
# print REvec.shape
# print REvec2.shape

print
print 'test Overdet'
# print Asq*REvec[0,:].T - Bsq*REvec[0,:].T*Evals[0]
print A*REvec[:,0] - B*REvec[:,0]*Evals[0]


print
print 'test Overdet'
# print Asq*REvec[0,:].T - Bsq*REvec[0,:].T*Evals[0]
print Asq*REvec3[:,0] - Bsq*REvec3[:,0]*Evals3[0]

print
print 'test eig'
print Asq*REvec2[:,0] - Bsq*REvec2[:,0]*Evals2[0]


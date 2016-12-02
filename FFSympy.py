#!/usr/bin/env python 

import os, errno
from sympy import *

NoLoad = False
try:
        from sympy import *
except:
        # print 'Sympy not present on machine, attempting to bypass'
        NoLoad = True

## HOW TO USE ##
#0. Open Python and execute "from FFSympy import *"
#1. Set Gamma Basis (CreateGamma('chiral' or 'sakurai')
#2. Find Form Factor (look at end of file)
#3. use vectorsubs or object.subs to imput your parameters

def GammaTOChroma(Opp):
    for i in [1,2,3,4]:
        for j in [1,2,3,4]:
            if 'g'+str(i)+'g'+str(j) in Opp:
                return 2**(i-1) + 2**(j-1)
    for i in [1,2,3,4]:
        if 'g'+str(i)+'g5' in Opp:
            return 15-2**(i-1)
    for i in [1,2,3,4]:
        if 'g'+str(i) in Opp:
            return 2**(i-1)
    if 'g5' in Opp:
        return 15
    elif 'I' in Opp:
        return 0
    else:
        return -1


if not NoLoad:
    GB = "UNSET"

    g0 = MatrixSymbol('g0', 4,4)
    g1 = MatrixSymbol('g1', 4,4)
    g2 = MatrixSymbol('g2', 4,4)
    g3 = MatrixSymbol('g3', 4,4)
    g4 = MatrixSymbol('g4', 4,4)
    g5 = MatrixSymbol('g5', 4,4)

    g_mu = [g0,g1,g2,g3,g4]

    sigma_mu_nu = [[zeros(4,4),-1j*g0*g1,-1j*g0*g2,-1j*g0*g3,zeros(4,4)],
                   [-1j*g1*g0,zeros(4,4),-1j*g1*g2,-1j*g1*g3,-1j*g1*g4],
                   [-1j*g2*g0,-1j*g2*g1,zeros(4,4),-1j*g2*g3,-1j*g2*g4],
                   [-1j*g3*g0,-1j*g3*g1,-1j*g3*g2,zeros(4,4),-1j*g3*g4],
                   [zeros(4,4),-1j*g4*g1,-1j*g4*g2,-1j*g4*g3,zeros(4,4)]]
    sigma_5_nu = [g5*g1,g5*g2,g5*g3,g5*g4]
    sigma_nu_5 = [g1*g5,g2*g5,g3*g5,g4*g5]

    G_unpol = (eye(4) + g4) / 2.
    G_pol = G_unpol* g5 * g3


    Ep, px, py, pz = symbols('Ep px py pz')
    Epp, pxp, pyp, pzp = symbols('Ep\' px\' py\' pz\'')
    # Eq, qx, qy, qz = symbols('Eq qx qy qz')

    Eq = Epp - Ep
    qx = pxp - px
    qy = pyp - py
    qz = pzp - pz


    m = symbols('m')
    F1, F2 ,F3 = symbols('F1 F2 F3')

    p = Matrix([1j*Ep,px,py,pz])
    pp = Matrix([1j*Epp,pxp,pyp,pzp])
    q = Matrix([1j*Eq,qx,qy,qz])

    paulie1 = ImmutableMatrix([[  0,  1],
                               [  1,  0]])

    paulie2 = ImmutableMatrix([[  0,-1j],
                               [ 1j,  0]])

    paulie3 = ImmutableMatrix([[  1,  0],
                               [  0, -1]])

    IMat2 = Identity(2)
    ZeroMat2 = ZeroMatrix(2,2)
    IMat4 = Identity(4)
    ZeroMat4 = ZeroMatrix(4,4)

    def GetRealCmplxMat(thisMat):
        return (thisMat+ImutableMatrix(np.conj(thisMat)))/2.,(thisMat-ImutableMatrix(np.conj(thisMat)))/2.

    def CreateGamma(GammaBasis):
        global GB
        global g0 
        global g1
        global g2 
        global g3
        global g4 
        global g5
        global g_mu 
        global sigma_mu_nu
        global sigma_5_nu
        global G_unpol
        global G_pol

        if GammaBasis == "sakurai":
            GB = GammaBasis
            g0 = diag(1,1,-1,-1)
            g1 = ImmutableMatrix([[  0,  0,  0,-1j],
                                  [  0,  0,-1j,  0],
                                  [  0, 1j,  0,  0],
                                  [ 1j,  0,  0,  0]])

            g2 = ImmutableMatrix([[  0,  0,  0, -1],
                                  [  0,  0,  1,  0],
                                  [  0,  1,  0,  0],
                                  [ -1,  0,  0,  0]])

            g3 = ImmutableMatrix([[  0,  0,-1j,  0],
                                  [  0,  0,  0, 1j],
                                  [ 1j,  0,  0,  0],
                                  [  0,-1j,  0,  0]])

            g4 = g0
            g5 = ImmutableMatrix([[  0,  0, -1,  0],
                                  [  0,  0,  0, -1],
                                  [ -1,  0,  0,  0],
                                  [  0, -1,  0,  0]])
        elif GammaBasis == "chiral":
            GB = GammaBasis
            g0 = ImmutableMatrix([[  0,  0,  1,  0],
                                  [  0,  0,  0,  1],
                                  [  1,  0,  0,  0],
                                  [  0,  1,  0,  0]])

            g1 = ImmutableMatrix([[  0,  0,-1j,  0],
                                  [  0,  0,  0,-1j],
                                  [ 1j,  0,  0,  0],
                                  [  0, 1j,  0,  0]])

            g2 = ImmutableMatrix([[  0,  0,  0, -1],
                                  [  0,  0,  1,  0],
                                  [  0,  1,  0,  0],
                                  [ -1,  0,  0,  0]])

            g3 = ImmutableMatrix([[  0,  0,-1j,  0],
                                  [  0,  0,  0, 1j],
                                  [ 1j,  0,  0,  0],
                                  [  0,-1j,  0,  0]])

            g4 = g0
            g5 = ImmutableMatrix([[  0,  0, -1,  0],
                                  [  0,  0,  0, -1],
                                  [ -1,  0,  0,  0],
                                  [  0, -1,  0,  0]])        
        else:
            # print "atempting symbolic gamma matricies"
            raise Exception("GammaBasis not recognised")
        g_mu = [g0,g1,g2,g3,g4]
        sigma_mu_nu = [[zeros(4,4),-1j*g0*g1,-1j*g0*g2,-1j*g0*g3,zeros(4,4)],
                       [-1j*g1*g0,zeros(4,4),-1j*g1*g2,-1j*g1*g3,-1j*g1*g4],
                       [-1j*g2*g0,-1j*g2*g1,zeros(4,4),-1j*g2*g3,-1j*g2*g4],
                       [-1j*g3*g0,-1j*g3*g1,-1j*g3*g2,zeros(4,4),-1j*g3*g4],
                       [zeros(4,4),-1j*g4*g1,-1j*g4*g2,-1j*g4*g3,zeros(4,4)]]
        sigma_5_nu = [g5*g1,g5*g2,g5*g3,g5*g4]
        sigma_nu_5 = [g1*g5,g2*g5,g3*g5,g4*g5]

        G_unpol = (eye(4) + g4) / 2.
        G_pol = G_unpol* g5 * g3


    def GetProj(thisProj):
        if 'P4' in thisProj:
            return G_unpol
        elif 'P3' in thisProj:
            return G_pol
        else:
            return -1

    def GetGamma(Opp):
        for i in [1,2,3,4]:
            for j in [1,2,3,4]:
                if 'g'+str(i)+'g'+str(j) in Opp:
                    return sigma_mu_nu[i][j]
        for i in [1,2,3,4]:
            if 'g'+str(i)+'g5' in Opp:
                return g_mu[i]*g5
        for i in [1,2,3,4]:
            if 'g'+str(i) in Opp:
                return g_mu[i]
        if 'g5' in Opp:
            return g5
        elif 'I' in Opp:
            return eye(4)
        else:
            return -1

    def GetProjGamma(Opp):
        return GetGamma(Opp),GetProj(Opp)

    ## Solves Equation without a projector (output spin matrix result) 
    ## Opp = opperator (e.g. 'g3g5')

    def FFunOpp(Opp):
        pplusm = (g4 - (1j/Ep) * (p[1]*g1 + p[2]*g2 + p[3]*g3) + (m / Ep) * eye(4))
        pprimeplusm = (g4 - (1j/Epp) * (pp[1]*g1 + pp[2]*g2 + pp[3]*g3) + (m / Epp )* eye(4))
        return simplify(pplusm * Opp * pprimeplusm)*sqrt(Epp*Ep/((Epp+m)*(Ep+m))) * 1/4.

    ## Solves Equation with projector and tracing
    ## Opp = 'P(4/3)Opp'

    def FFun(Opp):
        thisOpp,thisProj = GetProjGamma(Opp)
        return simplify(Trace(thisProj * FFunOpp(thisOpp)).doit())


    def FFunOppCheck(Opp,pmu,ppmu,mass,Rfac=True):
        p,pp = pmu,ppmu
        m = mass
        Ep, Epp = -1.0j*p[0],-1.0j*pp[0]
        pplusm = (g4 - (1j/Ep) * (p[1]*g1 + p[2]*g2 + p[3]*g3) + (m / Ep) * eye(4))
        pprimeplusm = (g4 - (1j/Epp) * (pp[1]*g1 + pp[2]*g2 + pp[3]*g3) + (m / Epp )* eye(4))
        if Rfac:
            return (pplusm * Opp * pprimeplusm)*sqrt(Epp*Ep/((Epp+m)*(Ep+m))) * 1/4.
        else:
            return (pplusm * Opp * pprimeplusm) * 1/4.

    def FFunCheck(Opp,pmu,ppmu,mass,Rfac=True):
        thisOpp,thisProj = GetProjGamma(Opp)
        return complex(Trace(thisProj * FFunOppCheck(thisOpp,pmu,ppmu,mass,Rfac=Rfac)).doit())


    # def FormFactor(Gproj,index):
    #     sigterm1 = 0
    #     for i in [1,2,3,4]:
    #         sigterm1 = sigterm1 + FFun(Gproj,sigma_mu_nu[index][i])*q[i]
    #     return FFun(Gproj,g_mu[index]) * F1 + (1/(2.*m)) * sigterm1 * F2


    def subsvector(vecen,vecx,vecy,vecz,vecennew,vecxnew,vecynew,vecznew, theobject):
        return theobject.subs(vecx,vecxnew).subs(vecy,vecynew).subs(vecz,vecznew).subs(vecen,vecennew)

    def subsZeroSinkMom(theobject):
        return subsvector(Epp,pxp,pyp,pzp,m,0,0,0,theobject)

    def subsZeroSourceMom(theobject):
        return subsvector(Ep,px,py,pz,m,0,0,0,theobject)

    def subsZeroMom(theobject):
        return subsZeroSinkMom(subsZeroSourceMom(theobject))

    #solves FFunOpp and substitutes in zero source and sink momenta
    def ZeroMomFFunOpp(Opp):
        return subsZeroMom(FFunOpp(Opp))


    def Slashed(vec):
        outvec = ZeroMat4
        for ivec,igma in zip(vec,g_mu):
            outvec += ivec*igma
        return outvec

    def IndexToGamma(*Vals):
        iout = []
        for iVal in Vals:
            if isinstance(iVal,int):
                iout.append(g_mu[iVal])
            else:
                iout.append(Slashed(iVals))
        return iout

    def CapGamma(*Vals):
        # outVals = IndexToGamma(Vals)
        outVals = Vals
        if len(Vals) == 0:
            return IMat4
        if len(Vals) == 1:
            return outVals[0]
        if len(Vals) == 2:
            return simplify((outVals[0]*outVals[1] - outVals[1]*outVals[0])/(2.))
        if len(Vals) == 3:
            output = ZeroMat4
            for ii in range(3):
                for ij in range(3):
                    for ik in range(3):
                        output += LeviCivita(ii,ij,ik)*outVals[ii]*outVals[ij]*outVals[ik]
            return simplify(output/(3.*2.))
        if len(Vals) == 4:
            output = ZeroMat4
            for ii in range(4):
                for ij in range(4):
                    for ik in range(4):
                        for il in range(4):
                            output += LeviCivita(ii,ij,ik,il)*outVals[ii]*outVals[ij]*outVals[ik]*outVals[il]
            return simplify(output/(4.*3.*2.))

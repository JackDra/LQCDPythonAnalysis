import numpy as np
alpha = -0.1587859292
mass = 0.5723293466
import MomParams as mp
E = np.sqrt(mass**2 +mp.qunit**2)
E
mp.qunit
qi = mp.qunit
F1 = 0.8169657582
F2 = 1.0741372468
F3 = -0.0052768705
F2coeff = (E+3*mass)/(2*mass)
F2coeff
F12coeff = -4*alpha*mass*qi
F12coeff
F12Term = F12coeff*(F1 + F2coeff*F2)
F12Term
F3coeff = -2*(E+mass)*qi
F3coeff
F1
F1 = 1.6077267731
F2 = 0.9821494063
F3 = 0.1029939487
F12Term = F12coeff*(F1 + F2coeff*F2)
F12Term
Rfac = 0.5
F12Term
F3Pred = (Rfac - F12Term)/F3coeff
F3Pred
F3
F12Term
F3coeff
Rfac
N =0.25*(E*mass*(E+mass)*(mass**2))**(-0.5)
N
F3Pred = (Rfac/N - F12Term)/F3coeff
F3Pred
N
alpha
F3Pred/(2*mass)
F3Pred = (Rfac/N + F12Term)/F3coeff
F3Pred
F3Pred = (Rfac/N - F12Term)/F3coeff
F3Pred
Rfac
Rfac = -Rfac
F3Pred = (Rfac/N - F12Term)/F3coeff
F3Pred
F3Pred = (Rfac/N + F12Term)/F3coeff
F3Pred
F3Pred = (Rfac/N - F12Term)/F3coeff
F3Pred
Rfac
Rfac = -0.3
F3Pred = (Rfac/N - F12Term)/F3coeff
F3Pred
F12Term
EqF1coeff = -0.0230956772517
EqF2coeff = -0.0468520317716
EqF1 = 1.6056013319
EqF2 = 1.0054554750
EqRfac = -0.138133521498
EqF3coeff = 0.169198825359
EqF3 = (EqRfac - EqF1coeff*EqF1 - EqF2coeff*EqF2)/EqF3coeff
N = EqF1coeff/(qi*mass*-4*alpha)
EqF3 = (EqRfac/N - EqF1coeff*EqF1 - EqF2coeff*EqF2)/EqF3coeff
EqF3
N
N*4
E+mass
(E+mass)/(2*mass)
(E+mass)**2/(2*mass)
import FFFuns as ff
opp = 'P3g4Top'
thisqvec = [0,0,1]
thisppvec = [0,0,0]
thismass = mass
thisalpha = alpha
[F1coeff,F2coeff,F3coeff],dump,dump2 = ff.VectorFFTop(opp,thisqvec,thisppvec,thismass,Rfac=True,alpha=thisalpha)
F1coeff
F2coeff
F3coeff
thismass
Erat = (E*mass*(E+mass)*mass**2)
Erat = (E*mass*(E+mass)*mass**2)**(-0.5)
Erat
Erat/4
Erat/4.
F1coeff/(-4*alpha*mass*qi)
Erat = ((E*mass)/(E+mass)*mass**2)**(0.5)
Erat
Erat = ((E*mass)/(E+mass)*2*mass)**(0.5)
Erat
Erat = ((E*mass)/((E+mass)*2*mass))**(0.5)
Erat
Erat = ((E*mass)/((E+mass)*2*mass))**(0.5)
F1coeff/(-4*alpha*mass*qi)
F1coeff/(-4*alpha*mass*qi*Erat)
Erat
Erat = ((E*mass)/((E+mass)*2*mass))**(0.5)/4.
Erat = ((E*mass)/((E+mass)*2*mass))**(0.5)
EqF1coeff/(-4*alpha*mass*qi*Erat)
N = (1/(4*E*mass)) * ((E*mass)/((E+mass)*mass**2))**(-0.5)
N
2*alpha*mass*qi*N
[F1coeff,F2coeff,F3coeff],dump,dump2 = ff.VectorFFTop(opp,thisqvec,thisppvec,thismass,Rfac=True,alpha=thisalpha)
[F1coeff,F2coeff,F3coeff],dump,dump2 = ff.VectorFFTop(opp,thisqvec,thisppvec,thismass,Rfac=False,alpha=thisalpha)
[F1coeffV,F2coeffV],dump,dump2 = ff.VectorFF('P4g4',thisqvec,thisppvec,thismass,Rfac=False)
F1coeffV
[F1coeffV,F2coeffV],dump,dump2 = ff.VectorFF('P4g4',thisqvec,thisppvec,thismass,Rfac=True)
F1coeffV
np.sqrt((E+mass)/(2*E))
[F1coeffV,F2coeffV],dump,dump2 = ff.VectorFF('P4g4',thisqvec,thisppvec,mass,Rfac=True)
F1coeffV
np.sqrt((E+mass)/(2*E))
E
reload(ff)
[F1coeffV,F2coeffV],dump,dump2 = ff.VectorFF('P4g4',thisqvec,thisppvec,mass,Rfac=True)
thispvec = [ipp-iq for ipp,iq in zip(thisppvec,thisqvec)]
E = ff.CreateEs(thispvec,thismass)
[F1coeffV,F2coeffV],dump,dump2 = ff.VectorFF('P4g4',thisqvec,thisppvec,mass,Rfac=True)
np.sqrt((E+mass)/(2*E))
F1coeff
F1coeffV
F1coeff
[F1coeffV,F2coeffV],dump,dump2 = ff.VectorFF('P4g4',thisqvec,thisppvec,mass,Rfac=True)
F1coeffV
[F1coeff,F2coeff,F3coeff],dump,dump2 = ff.VectorFFTop(opp,thisqvec,thisppvec,thismass,Rfac=False,alpha=thisalpha)
F1coeff
-4*alpha*mass*qi
-alpha*mass*qi/np.sqrt((E+mass)*2*mass**2*E)
F1coeff
[F1coeff,F2coeff,F3coeff],dump,dump2 = ff.VectorFFTop(opp,thisqvec,thisppvec,thismass,Rfac=True,alpha=thisalpha)
F1coeff
-alpha*mass*qi/np.sqrt((E+mass)*2*mass**2*E)
0.07965230114931679/0.13781128793927286
F1coeff
ExpF1coeff = -alpha*mass*qi/np.sqrt((E+mass)*2*mass**2*E)
F1coeff/ExpF1coeff
ExpF1coeff = -4*alpha*mass*qi/(4*mass*np.sqrt((E+mass)*2*E))
ExpF1coeff
F1coeff
F1coeff/ExpF1coeff
qi
ExpF1coeff = -4*alpha*mass/(4*mass*np.sqrt((E+mass)*2*E))
ExpF1coeff
F1coeff/ExpF1coeff
[F1coeff,F2coeff,F3coeff],dump,dump2 = ff.VectorFFTop(opp,np.array(thisqvec)*ff.qunit,thisppvec,thismass,Rfac=True,alpha=thisalpha)
F1coeff
ExpF1coeff = -4*alpha*mass*ff.qunit/(4*mass*np.sqrt((E+mass)*2*E))
ExpF1coeff
F1coeff
F1coeff/ExpF1coeff
thisqvec = np.array(thisqvec)*ff.qunit
thispvec = [ipp-iq for ipp,iq in zip(thisppvec,thisqvec)]
E = ff.CreateEs(thispvec,thismass)
ExpF1coeff = -4*alpha*mass*ff.qunit/(4*mass*np.sqrt((E+mass)*2*E))
F1coeff/ExpF1coeff
ExpF2coeff = ExpF1coeff*(E+3*mass)/(2*mass)
F2coeff
ExpF2coeff
F2coeff/ExpF2coeff
N = 1/(4*mass*np.sqrt((E+mass)*2*E))
ExpF3coeff = N*(-2*(E+mass)*ff.qunit)
ExpF3coeff
F3coeff
F1coeff
F2coeff
F3coeff
Rfac = -0.157461095832
F1
F3exp = ((Rfac/N) -F1coeff*F1 - F2coeff*F2)/(F3coeff)
F3exp
F3coeff
F2coeff
F1coeff
F1
F2
F3
F1 = 1.6056013319
F2 = 1.0054554750
F3 = 0.2780207625
F1coeff
F2coeff
F3coeff
F1coeff*F1 + F2coeff*F2 + F3coeff*F3
import readline
readline.write_history_file('./Testingff.int.py')

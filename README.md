(n.b. read files for more information)



THIS LIBARY USES xmltodict.py WHICH WAS PRODUCED BY:

https://github.com/martinblech/xmltodict.git


### Setup ###

./Setup.py

datadir is where your data is to be stored (n.b. your correlators should be stored in "datadir"/cfuns/k(YOUR KAPPA VALUE)/)
AnaProc is your max number of processors you want to use for analysis (run nproc to find your number of cpus, recommended using about half)

in "datadir"/results/ReadSetk....../ is where all the results get put.
All readable files are all in xml formmat, boots are stored in cPicked formatting. also 

Ratio Functions formatted as:
/Opperator/(doub/sing)Projector/qsqrd#/q###/files
/Operator/(doub/sing)Projector/qsqrd#/q###//boot/bootfiles

for example:
g4/doubP4/qsqrd0/q000/tsink26sm32doubP4g4q000.xml
g4/doubP4/qsqrd0/q000/boots/tsink26sm32doubP4g4q000.boot.p


Fitting formatted as:
/Opperator/(doub/sing)Projector/YOURFITTINGFLAG/qsqrd#/q###/files
/Operator/(doub/sing)Projector/YOURFITTINGFLAG/qsqrd#/q###//boot/bootfiles

for example:
g4/P4/OSFTsink/qsqrd0/q000/tsink26sm32doubP4g4B00q000.xml
g4/P4/OSFTsink/qsqrd0/q000/boots/tsink26sm32doubP4g4B00q000.boot.p




graphs are in "datadir"/results/ReadSetk....../graphs/ in similar file structure.


### Parameter Files ###

Params.py: Holds most parameters.

FitParams.py: Mainly fit ranges for different methods

MomParams.py: Momenta parameters and function

FFParams.py: Parameters for form factor calculations (mainly operators matricies in particular form factors)


### Bootstrapping ###


RunMcorrAna.py: Runs off correlators from binary correlators to bootstrapped results.
		Outputs ratio factors in results folder, and correlators in "results/cfuns" folder.

	./RunMcorrAna.py TwoPt
	Runs two point correlators off, creates correlation matrix results as well.

	./RunMcorrAna.py All All
	Runs off all three point constructions for all operators and sets of correlators

	./RunMcorrAna.py thisFF All
	Runs off all sets of correlators for opperators in "thisFF" (substitute thisFF for Vector/Scalar/PsScalar/PsVector/Tensor)
	
	./RunMcorrAna.py All CM '29'
	Runs off all operators for the correlation matrix set of data. 29 refers to the sink time for the set.
	  
	./RunMcorrAna.py All TSink '26 32 35 38'
	Runs off all operators for the sink time varying data. '26 32 35 38' refers to the sink times for the set.

	./RunMcorrAna.py All REvec '26 32'
	Runs off all operators for the Right eigenvector pre-inversion data. '26 32' refers to the sink times for the set.

	./RunMcorrAna.py All PoF '26 27'
	Runs off all operators for the Pencil of function correlation matrix set of data. '26 27' refers to the sink time for the set.



### Fitting ###

----- All scripts below follow a common input flag system defined in InputArgs.py. -----
----- Run Any script with -h for more information. -----


TryFits.py: Does constant fit to all data

	./TryFits.py
	Fits all operators

TrySummation.py: Does summation method analysis

	./TrySummation.py
	Summation analysis on all operators		 

TryOneStateFit.py: Does One State Fit

	./TryOneStateFit All
	Does all one state fit analysis (all operators and all collections)

	./TryOneStateFit CM 
	Does one state fit for correlation matrix results individually 

	./TryOneStateFit Tsink 
	Does one state fit for tsink time varying results individually 


TryTwoStateFit.py: Does Two State Fit

	./TryTwoStateFit All
	Does all two state fit analysis (all operators and all collections)

	./TryTwoStateFit CM 
	Does two state fit for correlation matrix results as combine fit 

	./TryTwoStateFit Tsink 
	Does two state fit for all tsink time varying results as combine fit 

	./TryTwoStateFit test32 
	Does two state fit for tsink time varying results 32,35,38 as combine fit 

	./TryTwoStateFit SmallSet 
	Does two state fit for tsink time varying results 26,29 as combine fit 


TryFF.py: Does Form Factor Analysis

	./TryFF.py
	Does Form Factor analysis for what ever is present in the results file (run this with no flags if you are lazy)


### Plotting ###

GraphRapNew.py: Graphs opperators (file comments to plot just zero momenta for all methods or lots of momenta for just RF)

	./GraphRapNew.py 

GraphRapSummary.py: Does the same as above, but just plots the extracted values instead of "R" function results

	./GraphRapSummary.py

GraphFFs.py: plots form factor results two ways, over q^2 and for each q^2 compared over methods

	./GraphFFs.py 





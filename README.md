(n.b. read files for more information)





### Parameter Files ###

Params.py: Holds most parameters.

FitParams.py: Mainly fit ranges for different methods

MomParams.py: Momenta parameters and function

FFParams.py: Parameters for form factor calculations (mainly operators matricies in particular form factors)







### Analysis ###

RunMcorrAna.py: Runs off correlators from binary correlators to bootstrapped results.

	./RunMcorrAna.py TwoPt
	Runs two point correlators off, creates correlation matrix results as well.

	./RunMcorrAna.py All All
	Runs off all three point constructions for all operators and sets of correlators

	./RunMcorrAna.py thisFF All
	Runs off all sets of correlators for opperators in "thisFF" (substitute thisFF for Vector/Scalar/PsScalar/PsVector/Tensor)
	
	./RunMcorrAna.py All CM 29
	Runs off all operators for the correlation matrix set of data. 29 refers to the sink time for the set.
	  
	./RunMcorrAna.py All TSink '26 32 35 38'
	Runs off all operators for the sink time varying data. '26 32 35 38' refers to the sink times for the set.

	./RunMcorrAna.py All REvec '26 32'
	Runs off all operators for the Right eigenvector pre-inversion data. '26 32' refers to the sink times for the set.


TryFits.py: Does constant fit to all data

	./TryFits.py
	Fits all operators

	./TryFits.py "THISOP"
	Fits all operators associated with "THISOP" (see CreateGammaList in OppFuns.py for more details)


TrySummation.py: Does summation method analysis

	./TrySummation.py
	Summation analysis on all operators

	./TrySummation.py "THISOP"
	Summation analysis on all operators associated with "THISOP" (see CreateGammaList in OppFuns.py for more details)
		 

TryOneStateFit.py: Does One State Fit

	./TryOneStateFit All
	Does all one state fit analysis (all operators and all collections)

	./TryOneStateFit CM "THISOP"
	Does one state fit for correlation matrix results individually ("THISOP" above)

	./TryOneStateFit Tsink "THISOP"
	Does one state fit for tsink time varying results individually ("THISOP" above)


TryTwoStateFit.py: Does Two State Fit

	./TryTwoStateFit All
	Does all two state fit analysis (all operators and all collections)

	./TryTwoStateFit CM "THISOP"
	Does two state fit for correlation matrix results as combine fit ("THISOP" above)

	./TryTwoStateFit Tsink "THISOP"
	Does two state fit for all tsink time varying results as combine fit ("THISOP" above)

	./TryTwoStateFit test32 "THISOP"
	Does two state fit for tsink time varying results 32,35,38 as combine fit ("THISOP" above)

	./TryTwoStateFit SmallSet "THISOP"
	Does two state fit for tsink time varying results 26,29 as combine fit ("THISOP" above)


TryFF.py: Does Form Factor Analysis

	./TryFF.py
	Does Form Factor analysis for what ever is present in the results file (run this if you are lazy)

	./TryFF.py "THISFF" "THISMETHOD"
	Does Form Factor analysis for type of form factor "THISFF" (Scalar Vector etc) and method types "THISMETHOD" (Fits/SumMeth/OSF../TSF..)


### Plotting ###

GraphRapNew.py: Graphs opperators (file comments to plot just zero momenta for all methods or lots of momenta for just RF)

	./GraphRapNew.py "THISOP"
	Graphs "THISOP" (see above)

GraphRapSummary.py: Does the same as above, but just plots the extracted values instead of "R" function results


GraphFFs.py: plots form factor results two ways, over q^2 and for each q^2 compared over methods

	./GraphFFs.py "THISFF"
	plots form factor "THISFF" (Scalar Vector etc)


### File Structure ###

in Params.py, add your machine to the "in THISMASHINE" caluses with your directory structures:
scriptdir is where scripts are
datadir is where your data is to be stored (n.b. your correlators should be stored in "datadir"/cfuns/k...../)
AnaProc is your max number of processors you want to use for analysis

in "datadir"/results/ReadSetk....../ is where all the results get put.
formatted as:
/Opperator/(doub/sing)Projector/files
/Opperator/(doub/sing)Projector/boot/bootfiles

files are all sets and each file has all the momenta results. also 

graphs are in "datadir"/results/ReadSetk....../graphs/ in similar formatting.



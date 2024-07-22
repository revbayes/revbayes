## name
mcmcmc
## title
Metropolis-Coupled MCMC analysis object
## description
The Mcmcmc analysis object keeps a model and the associated moves and monitors. The object is used to run Metropolis Couped Markov chain Monte Carlo (Mcmcmc) simulation on the model, using the provided moves, to obtain a sample of the posterior probability distribution. During the analysis, the monitors are responsible for sampling model parameters of interest.
## details
The Mcmcmc analysis object produced by a call to this function keeps copies of the model and the associated moves and monitors. The Mcmcmc analysis object is used to run Markov chain Monte Carlo (Mcmcmc) simulation on the model, using the provided moves, to obtain a sample of the posterior probability distribution. During the analysis, the monitors are responsible for sampling model parameters of interest.

The `mcmcmc.run()` method begins or continues an MCMCMC analysis. The run will continue for `generations`, or until a stopping rule is triggered: perhaps once the run has attained convergence, or after a certain amount of time has passed.  The run will be terminated once *all* convergence rules ([`srGelmanRubin()`], [`srGeweke()`], [`srMinESS()`], [`srStationarity()`]) have been fulfilled; or once *any* threshold rules ([`srMaxTime()`], [`srMaxIteration()`]) are met.

The parameters `checkpointFile` and `checkpointInterval` generate snapshots of the current state of the MCMCMC run from which the run can be continued if interrupted using the `mcmc.initializeFromCheckpoint()` method. An example is given on the documentation page for [`mcmc()`].

In node-dating analyses that employ calibration distributions, it is often of interest to run the analysis without character data to investigate how such distributions interact with one another and with the tree prior. Often, this is referred to as "running the analysis under the prior". In RevBayes, however, calibration distributions are treated as part of the data rather than of the prior. Accordingly, when setting `underPrior=TRUE` in `mcmcmc.burnin()` or `mcmcmc.run()`, both character data and calibration distributions are disregarded. To disregard the character data but not the calibrations, use `suppressCharacterData=TRUE` instead.

## authors
Michael Landis
Sebastian Hoehna
## see_also
mcmc
## example
	# Create a simple model (unclamped)
	a ~ exponential(1)
	mymodel = model(a)
	
	# Create a move vector and a monitor vector
	moves[1] = mvScale(a, lambda=1.0, weight=1.0)
	monitors[1] = mnFile(a,"output/out.log")
	
	# Create an mcmcmc object
	myMcmcmcObject = mcmcmc( mymodel, monitors, moves, nchains=4, deltaHeat=5)
	
	# Run a short analysis
	myMcmcmcObject.burnin( generations = 400, tuningInterval = 100)
	myMcmcmcObject.run( generations = 400)
	
	# print the summary of the operators (now tuned)
	myMcmcmcObject.operatorSummary()
	
## references
- citation: "Geyer,C.J. (1991) Markov chain Monte Carlo maximum likelihood. In Keramidas\
    \ (ed.), Computing Science and Statistics: Proceedings of the 23rd Symposium on\
    \ the Interface. Interface Foundation, Fairfax Station, pp. 156\u2013163."
  doi: null
  url: null
- citation: "Gilks,W.R. and Roberts,G.O. (1996) Strategies for improving MCMC. In\
    \ Gilks,W.R., Richardson,S. and Spiegelhalter (eds) Markov chain Monte Carlo in\
    \ Practice. Chapman&Hall, London, 89\u2013114."
  doi: null
  url: null
- citation: Altekar, G.; Dwarkadas, S.; Huelsenbeck, J. P. & Ronquist, F. Parallel
    metropolis coupled Markov chain Monte Carlo for Bayesian phylogenetic inference
    Bioinformatics, Oxford Univ Press, 2004, 20, 407-415.
  doi: null
  url: null

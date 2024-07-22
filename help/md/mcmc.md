## name
mcmc
## title
MCMC analysis object
## description
The MCMC analysis object keeps a model and the associated moves and monitors. The object is used to run Markov chain Monte Carlo (MCMC) simulation on the model, using the provided moves, to obtain a sample of the posterior probability distribution. During the analysis, the monitors are responsible for sampling model parameters of interest.
## details
The MCMC analysis object produced by a call to this function keeps copies of the model and the associated moves and monitors. The MCMC analysis object is used to run Markov chain Monte Carlo (MCMC) simulation on the model, using the provided moves, to obtain a sample of the posterior probability distribution. During the analysis, the monitors are responsible for sampling model parameters of interest.

The `mcmc.run()` method begins or continues an MCMC analysis. 

During each iteration of an analysis, moves are selected from those listed in the `moves` parameter.  With the default `moveschedule = "random"`, or `moveschedule = "sequential"`, moves will be attempted, on average, `weight` times per iteration.  If `moveschedule = "single"`, RevBayes will attempt exactly one move per iteration, corresponding to the behavior of software like BEAST or MrBayes.

The run will continue for `generations` iterations, or until a stopping rule is triggered: perhaps once the run has attained convergence, or after a certain amount of time has passed.  The run will be terminated once *all* convergence rules ([`srGelmanRubin()`], [`srGeweke()`], [`srMinESS()`], [`srStationarity()`]) in its `StoppingRule[]` argument have been fulfilled; or once *any* threshold rules ([`srMaxTime()`], [`srMaxIteration()`]) are met.

The parameters `checkpointFile` and `checkpointInterval` generate snapshots of the current state of the MCMC run from which the run can be continued if interrupted using the `mcmc.initializeFromCheckpoint()` method.

The `mcmc.initializeFromCheckpoint()` method allows an analysis to be continued from a checkpoint file.
New generations will be appended to existing monitor files.

## authors
Sebastian Hoehna
## see_also
mcmcmc
## example
    # Create a simple model (unclamped)
    a ~ dnExponential(1)
    mymodel = model(a)
    
    # Create a move vector and a monitor vector
    moves[1] = mvScale( a, lambda = 1.0, weight = 1.0 )
    monitors[1] = mnFile( a, filename = "output/out.log" )
    
    # Create an mcmc object
    mymcmcObject = mcmc( mymodel, monitors, moves )
    
    # Run a short analysis
    mymcmcObject.burnin( generations = 400, tuningInterval = 100 )
    mymcmcObject.run( generations = 400, checkpointFile = "output/out.ckp", checkpointInterval = 100 )
    
    # print the summary of the operators (now tuned)
    mymcmcObject.operatorSummary()

    # Resume analysis from the checkpoint file
    mymcmcObject.initializeFromCheckpoint( "output/out.ckp" )

    # Conduct an additional 400 generations
    mymcmcObject.run( generations = 400 )

    # Stopping rules are defined on the total number of generations
    # This command will have no effect, as 400 generations have already been performed.
    mymcmcObject.run( rules = [ srMaxIteration(400) ] )
	
## references
- citation: Metropolis N, AW Rosenbluth, MN Rosenbluth, AH Teller, E Teller (1953).
    Equation of state calculations by fast computing machines. Journal of Chemical
    Physics, 21:1087-1092.
  doi: 10.1063/1.1699114
  url: null
- citation: Hastings WK (1970) Monte Carlo sampling methods using Markov chains and
    their applications. Biometrika, 57:97-109.
  doi: 10.2307/2334940
  url: null
- citation: Höhna S, Landis MJ, Heath TA (2017). 
  Phylogenetic inference using `RevBayes`.
  Current Protocols in Bioinformatics.
  doi: 10.1002/cpbi.22
  url: null

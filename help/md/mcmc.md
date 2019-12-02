## name
mcmc
## title
MCMC analysis object
## description
The MCMC analysis object keeps a model and the associated moves and monitors. The object is used to run Markov chain Monte Carlo (MCMC) simulation on the model, using the provided moves, to obtain a sample of the posterior probability distribution. During the analysis, the monitors are responsible for sampling model parameters of interest.
## details
The MCMC analysis object produced by a call to this function keeps copies of the model and the associated moves and monitors. The MCMC analysis object is used to run Markov chain Monte Carlo (MCMC) simulation on the model, using the provided moves, to obtain a sample of the posterior probability distribution. During the analysis, the monitors are responsible for sampling model parameters of interest.
## authors
Sebastian Hoehna
## see_also
mcmcmc
## example
	# Create a simple model (unclamped)
	a ~ exponential(1)
	mymodel = model(a)
	
	# Create a move vector and a monitor vector
	moves[1] = mvScale(a, lambda=1.0, weight=1.0)
	monitors[1] = mnFile(a,"output/out.log")
	
	# Create an mcmc object
	mymcmcObject = mcmc( mymodel, monitors, moves)
	
	# Run a short analysis
	mymcmcObject.burnin( generations = 400, tuningInterval = 100)
	mymcmcObject.run( generations = 400)
	
	# print the summary of the operators (now tuned)
	mymcmcObject.operatorSummary()
	
## references
- citation: Metropolis N, AW Rosenbluth, MN Rosenbluth, AH Teller, E Teller (1953).
    Equation of state calculations by fast computing machines. Journal of Chemical
    Physics, 21:1087-1092.
  doi: null
  url: null
- citation: Hastings WK (1970) Monte Carlo sampling methods using Markov chains and
    their applications. Biometrika, 57:97-109.
  doi: null
  url: null

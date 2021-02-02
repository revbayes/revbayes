## name
dnPoisson
## title
Poisson Distribution
## description
A Poisson distribution defines probabilities for natural numbers. It is defined as the number of exponentially distributed events in a given interval.
## details
## authors
Sebastian Hoehna
## see_also
dnGeom
## example
	l ~ dnUniform(0.0,100.0)
	x ~ dnPoisson(l)
	x.clamp(10)
	moves[1] = mvSlide(l, delta=0.1, weight=1.0)
	monitors[1] = mnScreen(printgen=1000, separator = "	", l)
	mymodel = model(l)
	mymcmc = mcmc(mymodel, monitors, moves)
	mymcmc.burnin(generations=20000,tuningInterval=100)
	mymcmc.run(generations=200000)
	
## references

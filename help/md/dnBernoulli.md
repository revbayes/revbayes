## name
dnBernoulli
## title
Bernoulli Distribution
## description
A Bernoulli-distributed random variable takes the value 1 with probability p and the value 0 with probability 1-p.
## details
## authors
John Huelsenbeck
## see_also
dnBinomial
## example
	p ~ dnBeta(1.0,1.0)
	x ~ dnBernoulli(p)
	x.clamp(1)
	moves[1] = mvSlide(p, delta=0.1, weight=1.0)
	monitors[1] = screenmonitor(printgen=1000, separator = "	", x)
	mymodel = model(p)
	mymcmc = mcmc(mymodel, monitors, moves)
	mymcmc.burnin(generations=20000,tuningInterval=100)
	mymcmc.run(generations=200000)
	
## references

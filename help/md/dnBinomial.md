## name
dnBinomial
## title
Binomial Distribution
## description
Binomial probability distribution of x successes in size trials.
## details
The binomial probability distribution defines the number of success in size trials, where each trial has the same success probability p. The probability is given by (size choose x) p^(x) * (1-p)^(size-p)
## authors
Sebastian Hoehna
## see_also
dnBernoulli
## example
	p ~ dnBeta(1.0,1.0)
	x ~ dnBinomial(size=10,p)
	x.clamp(8)
	moves[1] = mvSlide(p, delta=0.1, weight=1.0)
	monitors[1] = screenmonitor(printgen=1000, separator = "	", x)
	mymodel = model(p)
	mymcmc = mcmc(mymodel, monitors, moves)
	mymcmc.burnin(generations=20000,tuningInterval=100)
	mymcmc.run(generations=200000)
	
## references

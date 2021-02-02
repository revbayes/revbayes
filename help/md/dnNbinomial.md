## name
dnNbinomial
## title
Negative Binomial Distribution
## description
Negative binomial probability distribution of x successes before r failures.
## details
The negative binomial probability distribution describes the number of successes before r failures, where the success probability is p. The probability is given by (x + r - 1 choose x) p^(x) * (1-p)^r
## authors
Walker Pett
## see_also
dnBinomial
## example
	p ~ dnBeta(1.0,1.0)
	x ~ dnNegativeBinomial(r=10,p)
	x.clamp(8)
	moves[1] = mvSlide(p, delta=0.1, weight=1.0)
	monitors[1] = screenmonitor(printgen=1000, separator = "	", x)
	mymodel = model(p)
	mymcmc = mcmc(mymodel, monitors, moves)
	mymcmc.burnin(generations=20000,tuningInterval=100)
	mymcmc.run(generations=200000)
	
## references

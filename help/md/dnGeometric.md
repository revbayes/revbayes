## name
dnGeometric
## title
Geometric Distribution
## description
A geometric distribution defines a random variable on natural numbers. The geometric distribution describes the number of success until the first failure, with success probability p.
## details
## authors
Sebastian Hoehna
## see_also
dnPoisson
mvRandomIntegerWalk
## example
	p ~ dnBeta(1.0,1.0)
	x ~ dnGeom(p)
	x.clamp(10)
	moves[1] = mvSlide(p, delta=0.1, weight=1.0)
	monitors[1] = screenmonitor(printgen=1000, separator = "	", p)
	mymodel = model(p)
	mymcmc = mcmc(mymodel, monitors, moves)
	mymcmc.burnin(generations=20000,tuningInterval=100)
	mymcmc.run(generations=200000)
	
## references

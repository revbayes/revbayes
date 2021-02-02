## name
dnBeta
## title
Beta Distribution
## description
The Beta probability distribution.
## details
The Beta distribution is a probability distribution defined on 0 to 1.
## authors
Sebastian Hoehna
## see_also
dnDirichlet
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

## name
dnBimodalLognormal
## title
Bimodal Lognormal Distribution
## description
A bimodal lognormal distribution, that is, with probability p a value is distributed according to the first lognormal distribution and with probability 1-p from the second lognormal distribution.
## details
## authors
Sebastian Hoehna
## see_also
dnBimodalNormal
dnLognormal
## example
	p ~ dnBeta(1.0,1.0)
	x ~ dnBimodalLognormal(mean1=-1,mean2=1,sd1=0.1,sd2=0.1,p=p)
	x.clamp( exp(1) )
	moves[1] = mvSlide(p, delta=0.1, weight=1.0)
	monitors[1] = screenmonitor(printgen=1000, separator = "	", x)
	mymodel = model(p)
	mymcmc = mcmc(mymodel, monitors, moves)
	mymcmc.burnin(generations=20000,tuningInterval=100)
	mymcmc.run(generations=200000)
	
## references

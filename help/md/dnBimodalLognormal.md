## name
dnBimodalLognormal
## title
Bimodal Lognormal Distribution
## description
The Bimodal Lognormal distribution unites two separate lognormal distributions.
## details
The bimodal lognormal distribution takes on five parameters: mean1 and mean2 (the means of two lognormal distributions), sd1 and sd2 (the standard deviations of two lognormal distributions), and p (between 0 and 1). The value will be distributed according to the first lognormal distribution with probability p, and distributed according to the second lognormal distribution with probability (1 - p).
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

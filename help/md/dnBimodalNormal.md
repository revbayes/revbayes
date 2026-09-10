## name
dnBimodalNormal
## title
Bimodal Normal dsitribution
## description
The Bimodal Normal distribution unites two separate normal distributions.
## details
The bimodal normal distribution takes five parameters: mean1 and mean2 (the means of two normal distributions), sd1 and sd2 (the standard deviations of two normal distributions), and p (between 0 and 1). The value will be distributed according to the first normal distribution with probability p, and distributed according to the second normal distribution with probability (1 - p).
## authors
Sebastian Hoehna
## see_also
dnBimodalLognormal
dnNormal
## example
	p ~ dnBeta(1.0,1.0)
	x ~ dnBimodalNormal(mean1=-1,mean2=1,sd1=0.1,sd2=0.1,p=p)
	x.clamp( 1 )
	moves[1] = mvSlide(p, delta=0.1, weight=1.0)
	monitors[1] = screenmonitor(printgen=1000, separator = "	", x)
	mymodel = model(p)
	mymcmc = mcmc(mymodel, monitors, moves)
	mymcmc.burnin(generations=20000,tuningInterval=100)
	mymcmc.run(generations=200000)
	
## references

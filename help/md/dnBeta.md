## name
dnBeta
## title
Beta Distribution
## description
The Beta probability distribution is a flexible distribution that returns a number between 0 and 1, so it is often used as a distribution for probabilities themselves.
## details
The Beta distribution takes on two parameters, alpha and beta. It is equivalent to the uniform when alpha = beta = 1. 

The probability density function is f(x) = x^(alpha - 1) * (1 - x)^(beta - 1) * Gamma(alpha + beta) / (Gamma(alpha) * Gamma(beta)), where Gamma is the gamma function.
## authors
Sebastian Hoehna
## see_also
dnDirichlet
gamma
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

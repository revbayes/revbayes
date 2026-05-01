## name
dnBinomial
## title
Binomial Distribution
## description
The Binomial probability distribution describes the probability of a number of successes for an experiment with a certain number of trials and probability of success per trial.
## details
The binomial distribution takes two parameters, p and size. It defines the number of successes in size trials, where each trial has the same success probability p. 

The probability density function is f(x) =  choose(size,x) * p^(x) * (1-p)^(size-p).
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
